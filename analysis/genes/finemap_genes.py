#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Identify, cluster, and fine-map all significant gene associations per HPO
"""


from os import path
import gzip
import csv
import pybedtools as pbt
import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod import families
from functools import reduce
from scipy.stats import norm
import argparse
from sys import stdout


def format_stat(x, is_phred=False, na_val=0):
    """
    Helper function to format & convert p-values as needed
    """

    if x == 'NA':
        return float(na_val)
    else:
        if is_phred:
            return 10 ** -float(x)
        else:
            return float(x)


def is_gene_sig(primary_p, secondary_p, n_nominal, primary_p_cutoff, 
                secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                secondary_or_nominal=True):
    """
    Checks if a gene should be considered exome-wide significant
    """

    if primary_p >= primary_p_cutoff:
        return False
    else:
        secondary = (secondary_p < secondary_p_cutoff)
        n_nom = (n_nominal >= n_nominal_cutoff)
        if secondary_or_nominal:
            if not any([secondary, n_nom]):
                return False
        elif not all([secondary, n_nom]):
            return False

    return True


def ci2se(ci):
    """
    Converts a tuple of (lower, upper) confidence interval bounds to standard error
    """

    ci = sorted(ci)

    return (ci[1] - ci[0]) / (2 * 1.96)


def finemap(gene_priors, gene_info, null_variance=0.42 ** 2):
    """
    Fine-maps a set of genes given their prior probability and effect size estimates
    Inputs:
        gene_priors : dict of genes with prior probabilities
        gene_info : dict of gene association stats as processed by process_hpo()
        null_se : float, variance of OR estimates under the null. By default,
                  this value is set to a 5% chance that ORs are > 2, as suggested
                  by Wakefield 2009
    """

    genes = list(gene_priors.keys())

    if len(genes) == 1:
        finemap_res = {genes[0] : {'ABF' : None, 'PIP' : 1}}
    else:
        finemap_res = {}

        # Compute ABF per gene
        for gene in genes:
            # From Wakefield, 2009
            theta = gene_info[gene]['lnOR']
            se = ci2se((gene_info[gene]['lnOR_upper'], gene_info[gene]['lnOR_lower']))
            V = se ** 2
            if V > 0:
                zsq = (theta ** 2) / V
                W = null_variance
                ABF = np.sqrt((V+W) / V) * np.exp((-zsq / 2) * (W / (V+W)))

                # Wakefield 2009 formulates BF relative to H0. We need to invert to 
                # obtain evidence & posterior for H1 (i.e., true non-zero effect)
                # However, we also are only testing for a _positive_ effect in cases,
                # so we will only invert ABF if theta >= 0. This is necessary to
                # prevent sites with enrichments in controls being scored with high ABF
                if theta >= 0:
                    ABF = 1 / ABF

            else:
                ABF = 0

            finemap_res[gene] = {'ABF' : ABF}

        # Compute PIP per gene as fraction of total BF, adjusted for prior probs
        # As per Mahajan 2018 (T2D fine-mapping paper, Nat. Genet.)
        posteriors = {}
        for gene in genes:
            posteriors[gene] = finemap_res[gene]['ABF'] * gene_priors[gene]
        posterior_sum = np.sum(list(posteriors.values()))
        for gene in genes:
            finemap_res[gene]['PIP'] = posteriors[gene] / posterior_sum

    return finemap_res


def parse_stats(stats_in, primary_p_cutoff, p_is_phred=True, 
                secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                secondary_or_nominal=True, sig_only=False, keep_genes=None,
                finemap_secondary=False):
    """
    Input: csv.reader of meta-analysis association stats
    Output: dict of gene stats (either sig_only or all genes in reader)
    """

    if path.splitext(stats_in)[1] in '.gz .bgz .bgzip'.split():
        csvin = gzip.open(stats_in, 'rt')
    else:
        csvin = open(stats_in)
    reader = csv.reader(csvin, delimiter='\t')

    stats_dict = {}

    for chrom, start, end, gene, n_nominal, top_cohort, case_freq, control_freq, \
        lnOR, lnOR_lower, lnOR_upper, zscore, primary_p, secondary_lnOR, \
        secondary_lnOR_lower, secondary_lnOR_upper, secondary_zscore, secondary_p \
        in reader:

        # Skip header line
        if chrom.startswith('#'):
            continue

        # If optioned, restrict on gene name
        if keep_genes is not None:
            if gene not in keep_genes:
                continue

        # Clean up gene data
        primary_p = format_stat(primary_p, p_is_phred, 1)
        secondary_p = format_stat(secondary_p, p_is_phred, 1)
        n_nominal = int(n_nominal)
        case_freq = format_stat(case_freq)
        control_freq = format_stat(control_freq)
        if finemap_secondary:
            use_lnOR = format_stat(secondary_lnOR)
            use_lnOR_lower = format_stat(secondary_lnOR_lower)
            use_lnOR_upper = format_stat(secondary_lnOR_upper)
            use_zscore = format_stat(secondary_zscore)
        else:
            use_lnOR = format_stat(lnOR)
            use_lnOR_lower = format_stat(lnOR_lower)
            use_lnOR_upper = format_stat(lnOR_upper)
            use_zscore = format_stat(zscore)

        # Store gene association stats
        if sig_only:
            if is_gene_sig(primary_p, secondary_p, n_nominal, primary_p_cutoff,
                           secondary_p_cutoff, n_nominal_cutoff, secondary_or_nominal):
                gene_bt = pbt.BedTool('\t'.join([chrom, start, end, gene]), 
                                      from_string=True)
                gene_stats = {'case_freq' : case_freq, 'control_freq' : control_freq,
                              'lnOR' : use_lnOR, 
                              'lnOR_lower' : use_lnOR_lower,
                              'lnOR_upper' : use_lnOR_upper, 
                              'zscore' : use_zscore,
                              'primary_p' : primary_p, 'secondary_p' : secondary_p,
                              'n_nominal' : n_nominal, 'gene_bt' : gene_bt}
                stats_dict[gene] = gene_stats
        else:
            gene_bt = pbt.BedTool('\t'.join([chrom, start, end, gene]), 
                                  from_string=True)
            gene_stats = {'case_freq' : case_freq, 'control_freq' : control_freq,
                          'lnOR' : use_lnOR, 
                          'lnOR_lower' : use_lnOR_lower,
                          'lnOR_upper' : use_lnOR_upper, 
                          'zscore' : use_zscore,
                          'primary_p' : primary_p, 'secondary_p' : secondary_p,
                          'n_nominal' : n_nominal, 'gene_bt' : gene_bt}
            stats_dict[gene] = gene_stats

    csvin.close()

    return stats_dict


def process_hpo(hpo, stats_in, primary_p_cutoff, p_is_phred=True, 
                secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                secondary_or_nominal=True, block_merge_dist=500000, 
                block_prefix='gene_block', null_variance=0.42 ** 2,
                finemap_secondary=False):
    """
    Loads & processes all necessary data for a single phenotype
    Returns a dict with the following entries:
        sig_genes : dict of sig genes, with each entry corresponding to stats
                    for a single significant gene
        all_genes : dict of sig genes + all genes within block_merge_dist with
                    stats per gene
        blocks : pbt.BedTool of all clustered genes to be fine-mapped
    """

    hpo_info = {'blocks' : {}}
    se_by_chrom = {}

    # First pass: parse data for significant genes only
    hpo_info['sig_genes'] = parse_stats(stats_in, primary_p_cutoff, p_is_phred, 
                                        secondary_p_cutoff, n_nominal_cutoff, 
                                        secondary_or_nominal, sig_only=True,
                                        finemap_secondary=finemap_secondary)

    # Second pass: parse data for all genes within block_merge_dist of sig_genes
    if len(hpo_info['sig_genes']) > 0:
        # Make bt of significant genes
        sig_gene_bts = [g['gene_bt'] for g in hpo_info['sig_genes'].values()]
        if len(sig_gene_bts) > 1:
            sig_genes_bt = sig_gene_bts[0].cat(*sig_gene_bts[1:], postmerge=False).sort()
        else:
            sig_genes_bt = sig_gene_bts[0]

        # Intersect sig genes with all genes
        all_genes_bt = pbt.BedTool(stats_in).cut(range(4)).sort()
        nearby_genes = all_genes_bt.closest(sig_genes_bt.sort(), d=True).\
                           filter(lambda x: int(x[8]) > -1 and \
                                            int(x[8]) <= block_merge_dist).\
                           saveas().to_dataframe().loc[:, 'name'].values.tolist()

        # Gather gene stats
        hpo_info['all_genes'] = parse_stats(stats_in, primary_p_cutoff, p_is_phred, 
                                            secondary_p_cutoff, n_nominal_cutoff, 
                                            secondary_or_nominal, sig_only=False,
                                            keep_genes=nearby_genes, 
                                            finemap_secondary=finemap_secondary)

        # Cluster significant genes into blocks to be fine-mapped
        gene_bts = [g['gene_bt'] for g in hpo_info['all_genes'].values()]
        if len(gene_bts) > 1:
            genes_bt = gene_bts[0].cat(*gene_bts[1:], postmerge=False).sort()
        else:
            genes_bt = gene_bts[0]
        blocks =  genes_bt.merge(d=block_merge_dist, c=4, o='distinct')

        # Perform initial genetic fine-mapping of each block with flat prior
        k = 0
        for block in blocks:
            k += 1
            block_id = '_'.join([hpo, block_prefix, str(k)])
            genes = block[3].split(',')
            gene_priors = {gene : 1 / len(genes) for gene in genes}
            finemap_res = finemap(gene_priors, hpo_info['all_genes'], null_variance)
            hpo_info['blocks'][block_id] = {'coords' : block,
                                            'finemap_res' : finemap_res}

    # If no genes are significant, add empty placeholder dict for all genes
    else:
        hpo_info['all_genes'] = {}

    return hpo_info


def load_all_hpos(statslist, secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                  secondary_or_nominal=True, block_merge_dist=500000,
                  block_prefix='gene_block', finemap_secondary=False):
    """
    Wrapper function to process each HPO with process_hpo()
    Returns a dict with one entry per HPO
    """

    hpo_data = {}

    with open(statslist) as infile:
        reader = csv.reader(infile, delimiter='\t')
        for hpo, stats_in, pval, in reader:
            primary_p_cutoff = float(pval)
            hpo_data[hpo] = process_hpo(hpo, stats_in, primary_p_cutoff, 
                                        p_is_phred=True, 
                                        secondary_p_cutoff=secondary_p_cutoff, 
                                        n_nominal_cutoff=n_nominal_cutoff, 
                                        secondary_or_nominal=secondary_or_nominal,
                                        block_merge_dist=block_merge_dist,
                                        block_prefix=block_prefix,
                                        finemap_secondary=finemap_secondary)

    return hpo_data


def estimate_null_variance_basic(hpo_data):
    """
    Estimates null variance per phenotype from average of all significant genes
    """

    vardict_all = {hpo : {} for hpo in hpo_data.keys()}
    vardict_sig = {hpo : [] for hpo in hpo_data.keys()}
    vardict_best = {hpo : [] for hpo in hpo_data.keys()}

    for hpo, dat in hpo_data.items():
        
        for gene, gdat in dat['all_genes'].items():
            # Estimate null variance from effect size per Wakefield, AJHG, 2007
            var = (float(gdat['lnOR']) / 1.96) ** 2
            vardict_all[hpo][gene] = var
            if gene in dat['sig_genes'].keys():
                vardict_sig[hpo].append(var)
        
        for bdat in dat['blocks'].values():
            bpvals = [(gene, dat['all_genes'][gene]['primary_p']) for gene \
                          in bdat['finemap_res'].keys()]
            best_gene = sorted(bpvals, key=lambda x: x[1])[0][0]
            vardict_best[hpo].append(vardict_all[hpo][best_gene])

    # Compute 2 null variance estimates for fine-mapping:
    # 1. Mean of all significant genes
    # 2. Mean of all top genes (one per block)

    v1 = np.nanmean([x for l in vardict_sig.values() for x in l if x > 0])
    v2 = np.nanmean([x for l in vardict_best.values() for x in l if x > 0])

    return [v1, v2]


def estimate_null_variance_gs(gs_lists, statslist, finemap_secondary=False):
    """
    Estimates null variance for all cases from the average of a list of known causal genes
    """

    # Read effect sizes per gene from highest level phenotype
    if finemap_secondary:
        cols = 'gene meta_lnOR_secondary'.split()
    else:
        cols = 'gene meta_lnOR'.split()
    stats = pd.read_csv(open(statslist).readline().rstrip().split('\t')[1], 
                        delimiter='\t').loc[:, cols]
    stats.columns = 'gene lnOR'.split()

    # Iterate over lists of known causal genes and compute mean variance
    var = []
    for gslist in gs_lists:
        gs_genes = open(gslist).read().splitlines()
        gs_vars = (stats.lnOR[stats.gene.isin(gs_genes)].astype(float) / 1.96) ** 2
        var.append(float(np.nanmean(gs_vars)))

    return var


def update_finemap(hpo_data, W):
    """
    Update initial finemapping results (flat prior) with a new null variance
    """

    for hpo, hdat in hpo_data.items():
        for block_id, bdat in hdat['blocks'].items():
            genes = list(bdat['finemap_res'].keys())
            gene_priors = {gene : 1 / len(genes) for gene in genes}
            finemap_res = finemap(gene_priors, hpo_data[hpo]['all_genes'], W)
            hpo_data[hpo]['blocks'][block_id]['finemap_res'] = finemap_res

    return hpo_data


def make_cs(finemap_res, cs_val=0.95):
    """
    Aggregate genes into credible set based on sum of ranked PIPs
    """

    cs_genes = []

    vals = pd.DataFrame.from_dict(finemap_res, orient='index')
    vals = vals.sort_values(by='PIP', ascending=False)

    cs_sum = 0
    for i in range(len(vals)):
        cs_genes.append(vals.index[i])
        cs_sum += vals.PIP[i]

        if cs_sum >= cs_val:
            break
    
    return sorted(list(set(cs_genes)))


def make_sig_genes_df(hpo_data, naive=False, sig_only=False, cs_val=0.95):
    """
    Makes pd.DataFrame of all significant genes, CSs, HPOs, ABFs, and PIPs
    """

    sig_df = pd.DataFrame(columns='HPO gene ABF PIP credible_set'.split())

    for hpo in hpo_data.keys():
        for block_id, block in hpo_data[hpo]['blocks'].items():
            cs_genes = make_cs(block['finemap_res'], cs_val)
            for gene, stats in block['finemap_res'].items():
                if not naive:
                    ABF = stats['ABF']
                    PIP = stats['PIP']
                else:
                    ABF = np.nan
                    PIP = 1 / len(block['finemap_res'].keys())
                if gene in cs_genes:
                    cs_assign = block_id
                else:
                    cs_assign = None
                if sig_only:
                    if gene in hpo_data[hpo]['sig_genes'].keys():
                        sig_df = sig_df.append(pd.Series([hpo, gene, ABF, PIP, 
                                                          cs_assign],
                                                         index=sig_df.columns), 
                                                         ignore_index=True)
                else:
                    sig_df = sig_df.append(pd.Series([hpo, gene, ABF, PIP,
                                                      cs_assign],
                                                     index=sig_df.columns), 
                                                     ignore_index=True)

    return sig_df.sort_values('PIP', ascending=False)


def rmse(pairs):
    """
    Compute root mean-squared error for a list of tuples
    """

    mse = [(a - b) ** 2 for a, b in pairs]
    
    return np.sqrt(np.sum(mse) / len(pairs))


def functional_finemap(hpo_data, gene_features_in, l1_l2_mix, logit_alpha, 
                       cs_val=0.95, null_variance=0.42 ** 2, converge_rmse=10e-8, 
                       quiet=False):
    """
    Conduct E-M optimized functional fine-mapping for all gene blocks & HPOs
    Two returns:
        - pd.DataFrame of all genes with their HPOs, ABFs, and PIPs
        - statsmodels.iolib.table.SimpleTable of final logit coefficients
    """

    if not quiet:
        msg = '\nStarting fine-mapping with null variance (W) = {:.5}\n' + \
              '  Iter.\tPIP RMSE\tCoeff. RMSE'
        print(msg.format(null_variance))

    # Re-fine-map each block with specified null variance (for BMA)
    for hpo in hpo_data.keys():
        for block_id, block_data in hpo_data[hpo]['blocks'].items():
            genes = list(block_data['finemap_res'].keys())
            gene_priors = {gene : 1 / len(genes) for gene in genes}
            updated_finemap_res = finemap(gene_priors, hpo_data[hpo]['all_genes'], 
                                          null_variance)
            hpo_data[hpo]['blocks'][block_id]['finemap_res'] = updated_finemap_res

    # Load gene features, standardize, and subset to genes present in any gene block
    features = pd.read_csv(gene_features_in, delimiter='\t')
    for colname in '#chr chr start end'.split():
        if colname in features.columns:
            features.drop(labels=colname, axis=1, inplace=True)
    for col in features.columns.tolist()[1:]:
        missing = (features[col] == '.')
        if missing.sum() > 0:
            fmean = np.nanmean(features.loc[~missing, col].astype(float))
            features.loc[missing, col] = fmean
    sig_df = make_sig_genes_df(hpo_data, cs_val=cs_val)
    features.iloc[:, 1:] = scale(features.iloc[:, 1:])
    features = features.loc[features.gene.isin(sig_df.gene), :]

    # Iterate until convergence
    coeffs = [0 for x in features.columns.tolist()[1:]]
    k = 0
    rmse_PIP, rmse_coeffs = 100, 100
    while rmse_PIP >= converge_rmse or rmse_coeffs >= converge_rmse:
        k += 1
        # Join sig_df with features for logistic regression
        logit_df = sig_df.loc[:, 'gene PIP'.split()].merge(features, how='left', 
                                                           on='gene')

        # When fitting regression, take mean of PIPs for genes appearing multiple times
        # This can happen due to multiple HPO associations with the same gene
        logit_df = logit_df.groupby('gene').mean()

        # Fit logit GLM & predict new priors
        glm = GLM(logit_df.PIP, logit_df.drop(labels='PIP', axis=1),
                  family=families.Binomial())
        if logit_alpha is None:
            logit = glm.fit()
        else:
            logit = glm.fit_regularized(L1_wt = 1 - l1_l2_mix, alpha=logit_alpha)
        pred_priors = logit.predict(features.drop(labels='gene', axis=1))
        new_priors = features.gene.to_frame().join(pred_priors.to_frame(name='prior'))
        
        # Re-finemap all blocks with new priors
        for hpo in hpo_data.keys():
            for block_id, block_data in hpo_data[hpo]['blocks'].items():
                genes = list(block_data['finemap_res'].keys())
                gene_priors = {gene : new_priors[new_priors.gene == gene].iloc[0]['prior'] \
                               for gene in genes}
                # Normalize new priors within each block such that sum(priors) = 1
                sum_priors = np.sum(list(gene_priors.values()))
                gene_priors = {g : p / sum_priors for g, p in gene_priors.items()}
                new_finemap_res = finemap(gene_priors, hpo_data[hpo]['all_genes'],
                                          null_variance)
                hpo_data[hpo]['blocks'][block_id]['finemap_res'] = new_finemap_res

        # Compute RMSE for PIPs and coefficients
        new_sig_df = make_sig_genes_df(hpo_data, cs_val=cs_val)
        PIPs_oldnew = sig_df.merge(new_sig_df, on='HPO gene'.split(), 
                                   suffixes=('_old', '_new'))\
                            .loc[:, 'PIP_old PIP_new'.split()]\
                            .to_records(index=False).tolist()
        rmse_PIP = rmse(PIPs_oldnew)
        new_coeffs = logit.params.to_list()
        coeffs_oldnew = tuple(zip(coeffs, new_coeffs))
        rmse_coeffs = rmse(coeffs_oldnew)

        # Update coeffs and sig_df
        coeffs = new_coeffs
        sig_df = new_sig_df

        # Print iteration info
        if not quiet:
            print('  {:,}\t{:.3E}\t{:.3E}'.format(k, rmse_PIP, rmse_coeffs))

    # Report completion
    print('Converged after {:,} iterations'.format(k))

    # Save table
    if logit_alpha is None:
        tab_out = logit.summary().tables[1]
    else:
        tab_out = logit.params

    return sig_df.sort_values('PIP', ascending=False), tab_out


def bmavg(sig_dfs, hpo_data, outfile, sig_only=False, cs_val=0.95):
    """
    Average ABFs and PIPs for all genes across models (list of sig_dfs)
    """

    if sig_only:
        sig_genes = []
        for hpo, info in hpo_data.items():
            for gene in info['sig_genes'].keys():
                sig_genes.append((hpo, gene))
    else:
        sig_genes = sig_dfs[0].loc[:, 'HPO gene'.split()].to_records(index=False).tolist()

    sig_df = pd.DataFrame(columns=sig_dfs[0].columns)

    for hpo, gene in sig_genes:
        ABFs = [df.loc[(df.HPO == hpo) & (df.gene == gene), :].iloc[0]['ABF'] \
                for df in sig_dfs]
        if all(a is None for a in ABFs):
            avg_ABF = None
        else:
            avg_ABF = np.nanmean(ABFs)
        PIPs = [df.loc[(df.HPO == hpo) & (df.gene == gene), :].iloc[0]['PIP'] \
                for df in sig_dfs]
        avg_PIP = np.nanmean(PIPs)
        sig_df = sig_df.append(pd.Series([hpo, gene, avg_ABF, avg_PIP, None], 
                                         index=sig_df.columns),
                               ignore_index=True)

    sig_df = sig_df.sort_values('PIP', ascending=False)

    for hpo in hpo_data:
        for bid, binfo in hpo_data[hpo]['blocks'].items():
            bgenes = list(binfo['finemap_res'].keys())
            block_df = sig_df.loc[(sig_df.HPO == hpo) & (sig_df.gene.isin(bgenes)), :]
            block_pips = [tuple(x) for x in block_df.loc[:, 'gene PIP'.split()].to_numpy()]
            cs_sum = 0
            for gene, pip in block_pips:
                sig_df.loc[(sig_df.HPO == hpo) & (sig_df.gene == gene), 'credible_set'] \
                    = bid
                cs_sum += pip
                if cs_sum >= cs_val:
                    break

    sig_df.rename(columns={'HPO' : '#HPO'}).\
           to_csv(outfile, sep='\t', index=False, na_rep='NA')
    outfile.close()

    return sig_df


def coeff_avg(coeff_tables, outfile, logit_alpha):
    """
    Average logit coefficients across all models & write as tsv
    Input: list of statsmodels.iolib.table.SimpleTable
    """

    if logit_alpha is None:

        cols = 'feature coeff stderr zscore pvalue lower upper'.split()

        def ct2df(ct, k):
            df = pd.read_html(ct.as_html(), header=0)[0]
            fmtcols = ' '.join([cols[0]] + [x + '_{0}' for x in cols[1:]])
            df.columns = fmtcols.format(k).split()
            return df

        ct_dfs = [ct2df(ct, k+1) for k, ct in enumerate(coeff_tables)]

    else:

        cols = 'feature coeff'.split()

        def ct2df(ct, k):
            features = ct.index.tolist()
            values = ct.values.tolist()
            df = pd.DataFrame(tuple(zip(features, values)), columns = cols)
            return df

        ct_dfs = [ct2df(ct, k+1) for k, ct in enumerate(coeff_tables)]

    merged_df = reduce(lambda left, right: pd.merge(left, right, on='feature'), 
                       ct_dfs)
    avg_df = pd.DataFrame(merged_df.feature, columns=['feature'])
    for col in cols[1:]:
        subdf = merged_df.loc[:, [x for x in merged_df.columns if col in x]]
        avg_df[col] = subdf.mean(axis=1)
    
    # Recompute P-values according to averaged Z-score
    if logit_alpha is None:
        avg_df['pvalue'] = [2 * norm.sf(abs(z)) for z in avg_df.zscore]

    avg_df.rename(columns={'feature' : '#feature'}).\
           to_csv(outfile, sep='\t', index=False, na_rep='NA')
    outfile.close()


def output_assoc_bed(hpo_data, sig_df, outfile, conf_pip=0.1, vconf_pip=0.9, cnv='NS'):
    """
    Format final list of significant gene-phenotype pairs with summary statistics
    """

    cols = 'chr start end gene cnv hpo control_freq case_freq ln_or ln_or_ci_lower ' + \
            'ln_or_ci_upper pvalue pip credible_set_id'
    outfile.write('#' + '\t'.join(cols.split()) + '\n')

    # Iterate over each gene with PIP >= conf_pip
    for gdict in sig_df.loc[sig_df.PIP >= conf_pip, :].to_dict(orient='index').values():
        # Get gene info
        hpo = gdict['HPO']
        gene = gdict['gene']
        pip = gdict['PIP']
        cred = gdict['credible_set']
        chrom = str(hpo_data[hpo]['all_genes'][gene]['gene_bt'][0].chrom)
        start = str(hpo_data[hpo]['all_genes'][gene]['gene_bt'][0].start)
        end = str(hpo_data[hpo]['all_genes'][gene]['gene_bt'][0].end)
        case_freq = hpo_data[hpo]['all_genes'][gene]['case_freq']
        control_freq = hpo_data[hpo]['all_genes'][gene]['control_freq']
        lnor = hpo_data[hpo]['all_genes'][gene]['lnOR']
        lnor_lower = hpo_data[hpo]['all_genes'][gene]['lnOR_lower']
        lnor_upper = hpo_data[hpo]['all_genes'][gene]['lnOR_upper']
        pval = hpo_data[hpo]['all_genes'][gene]['primary_p']
        # If p-value was rounded to zero, modify by taking the unadjusted normal P-value
        if pval == 0:
            pval = norm.sf(hpo_data[hpo]['all_genes'][gene]['zscore'])

        # Write gene-phenotype pair to file
        outline = '\t'.join([chrom, start, end, gene, cnv, hpo])
        outnums_fmt = '\t{:.3E}\t{:.3E}\t{:.3}\t{:.3}\t{:.3}\t{:.3E}\t{:.3}'
        outline += outnums_fmt.format(control_freq, case_freq, lnor, lnor_lower, 
                                      lnor_upper, pval, pip)
        outline += '\t' + cred + '\n'
        outfile.write(outline)

    outfile.close()


def iv_mean(values, variances, conf=0.95):
    """
    Returns inverse-variance weighted mean of values and conf% confidence interval
    """

    weights = [1 / v for v in variances]

    wsum = np.nansum(weights)

    numerator = np.nansum([x / v for x, v in zip(values, variances)])

    ivm = numerator / wsum

    pooled_se = np.sqrt(1 / wsum)

    ci_dist = norm.ppf(conf) * pooled_se

    return ivm, (ivm - ci_dist, ivm + ci_dist)


def output_genes_bed(hpo_data, sig_df, outfile, ncase_dict,
                     conf_pip=0.1, vconf_pip=0.9, cnv='NS'):
    """
    Format final list of significant genes and compute pooled summary statistics
    """

    cols = 'chr start end gene cnv pooled_control_freq pooled_case_freq ' + \
           'pooled_ln_or pooled_ln_or_ci_lower pooled_ln_or_ci_upper ' + \
           'best_pip n_hpos all_hpos vconf_hpos conf_hpos credible_set_ids'
    outfile.write('#' + '\t'.join(cols.split()) + '\n')

    for gene in set(sig_df.loc[sig_df.PIP >= conf_pip, ].gene.tolist()):
        # Get gene info
        genedf = sig_df.loc[(sig_df.gene == gene) & (sig_df.PIP >= conf_pip), :]
        all_hpos = sorted(genedf.HPO.tolist())
        n_hpos = len(all_hpos)
        vconf_hpos = sorted(genedf.HPO[sig_df.PIP >= vconf_pip].tolist())
        conf_hpos = sorted(genedf.HPO[sig_df.PIP < vconf_pip].tolist())
        chrom = str(hpo_data[all_hpos[0]]['all_genes'][gene]['gene_bt'][0].chrom)
        start = str(hpo_data[all_hpos[0]]['all_genes'][gene]['gene_bt'][0].start)
        end = str(hpo_data[all_hpos[0]]['all_genes'][gene]['gene_bt'][0].end)
        control_freq = hpo_data[all_hpos[0]]['all_genes'][gene]['control_freq']
        best_pip = np.nanmax(genedf.PIP)
        cred_ids = sorted(genedf.credible_set.tolist())

        if n_hpos > 1:
            # Compute pooled effect size as inverse-variance weighted average
            lnors = []
            variances = []
            for hpo in all_hpos:
                lnors.append(hpo_data[hpo]['all_genes'][gene]['lnOR'])
                lower = hpo_data[hpo]['all_genes'][gene]['lnOR_lower']
                upper = hpo_data[hpo]['all_genes'][gene]['lnOR_upper']
                variances.append(ci2se((lower, upper)) ** 2)
            lnor, lnor_ci = iv_mean(lnors, variances)
            lnor_lower, lnor_upper = sorted(lnor_ci)

            # Compute pooled case frequency as mean weighted by np.sqrt(N_cases)
            cfreqs = [hpo_data[hpo]['all_genes'][gene]['case_freq'] for hpo in all_hpos]
            cfreq_weights = [np.sqrt(ncase_dict[hpo]) for hpo in all_hpos]
            case_freq = np.average(cfreqs, weights=cfreq_weights)

        else:
            lnor = hpo_data[all_hpos[0]]['all_genes'][gene]['lnOR']
            lnor_lower = hpo_data[all_hpos[0]]['all_genes'][gene]['lnOR_lower']
            lnor_upper = hpo_data[all_hpos[0]]['all_genes'][gene]['lnOR_upper']
            case_freq = hpo_data[all_hpos[0]]['all_genes'][gene]['case_freq']

        # Reassign empty variables for writing to outfile, if necessary
        if len(vconf_hpos) == 0:
            vconf_hpos = ['NA']
        if len(conf_hpos) == 0:
            conf_hpos = ['NA']

        # Write gene stats to file
        outline = '\t'.join([chrom, start, end, gene, cnv])
        outnums_fmt = '\t{:.3E}\t{:.3E}\t{:.3}\t{:.3}\t{:.3}\t{:.3}'
        outline += outnums_fmt.format(control_freq, case_freq, lnor, lnor_lower, 
                                      lnor_upper, best_pip)
        outline += '\t' + '\t'.join([str(n_hpos), ';'.join(all_hpos), ';'.join(vconf_hpos),
                                     ';'.join(conf_hpos), ';'.join(cred_ids)]) + \
                   '\n'
        outfile.write(outline)

    outfile.close()


def output_credsets_bed(hpo_data, sig_df, outfile, conf_pip=0.1, vconf_pip=0.9, cnv='NS'):
    """
    Format final list of credible sets with summary statistics
    """
    
    cols = 'chr start end credible_set_id cnv hpo mean_control_freq mean_case_freq ' + \
           'pooled_ln_or pooled_ln_or_ci_lower pooled_ln_or_ci_upper ' + \
           'best_pvalue n_genes all_genes top_gene vconf_genes conf_genes'
    outfile.write('#' + '\t'.join(cols.split()) + '\n')

    for cred in [x for x in set(sig_df.credible_set.tolist()) if x is not None]:
        # Get basic credible set info
        cred_df = sig_df.loc[sig_df.credible_set == cred, :]
        hpo = cred_df.HPO.tolist()[0]
        genes = sorted(list(set(cred_df.gene.tolist())))
        vconf_genes = sorted(list(set(cred_df.gene[cred_df.PIP >= vconf_pip].tolist())))
        conf_genes = sorted(list(set(cred_df.gene[(cred_df.PIP < vconf_pip) & \
                                                  (cred_df.PIP >= conf_pip)].tolist())))
        top_pip = np.nanmax(cred_df.PIP)
        top_gene = sorted(list(set(cred_df.gene[cred_df.PIP == top_pip].tolist())))
        top_gene = [g for g in top_gene if g in vconf_genes + conf_genes]
        n_genes = len(genes)
        gdat = hpo_data[hpo]['all_genes']
        chrom = str(gdat[genes[0]]['gene_bt'][0].chrom)
        start = str(np.nanmin([gdat[g]['gene_bt'][0].start for g in genes]))
        end = str(np.nanmax([gdat[g]['gene_bt'][0].end for g in genes]))
        best_p = np.nanmin([gdat[g]['primary_p'] for g in genes])
        if best_p == 0:
            best_z = np.nanmax([gdat[g]['zscore'] for g in genes])
            best_p = norm.sf(best_z)
        control_freq = np.nanmean([gdat[g]['control_freq'] for g in genes])
        case_freq = np.nanmean([gdat[g]['case_freq'] for g in genes])

        # Get pooled effect size as inverse-variance weighted mean of all genes
        if n_genes > 1:
            gors = [gdat[g]['lnOR'] for g in genes]
            gvars = [ci2se((gdat[g]['lnOR_lower'], gdat[g]['lnOR_upper'])) ** 2 for g in genes]
            lnor, lnor_ci = iv_mean(gors, gvars)
            lnor_lower, lnor_upper = sorted(lnor_ci)

        else:
            lnor = gdat[genes[0]]['lnOR']
            lnor_lower = gdat[genes[0]]['lnOR_lower']
            lnor_upper = gdat[genes[0]]['lnOR_upper']

        # Reassign empty variables for writing to outfile, if necessary
        if len(vconf_genes) == 0:
            vconf_genes = ['NA']
        if len(conf_genes) == 0:
            conf_genes = ['NA']
        if len(top_gene) == 0:
            top_gene = ['NA']

        # Write gene stats to file
        outline = '\t'.join([chrom, start, end, cred, cnv, hpo])
        outnums_fmt = '\t{:.3E}\t{:.3E}\t{:.3}\t{:.3}\t{:.3}\t{:.3E}'
        outline += outnums_fmt.format(control_freq, case_freq, lnor, lnor_lower, 
                                      lnor_upper, best_p)
        outline += '\t' + '\t'.join([str(n_genes), ';'.join(genes), ';'.join(top_gene),
                                     ';'.join(vconf_genes), ';'.join(conf_genes)]) + \
                   '\n'
        outfile.write(outline)

    outfile.close()


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('statslist', help='tsv of metadata per phenotype. Three ' +
                        'required columns: HPO, path to meta-analysis stats, ' +
                        'and primary P-value cutoff.')
    parser.add_argument('gene_features', help='tsv of gene features to use ' +
                        'for functional fine-mapping. First column = gene name. ' +
                        'All other columns must be numeric features. Optionally, ' + 
                        'first three columns can be BED-like, and will be dropped.')
    parser.add_argument('hpos_by_cohort', help='tsv of sample sizes per HPO.')
    parser.add_argument('--secondary-p-cutoff', help='Maximum secondary P-value to ' + 
                        'consider as significant. [default: 1]', default=1, type=float)
    parser.add_argument('--cnv', help='Indicate CNV type. [default: NS]', default='NS')
    parser.add_argument('--min-nominal', help='Minimum number of individual cohorts ' + 
                        'required to be nominally significant to consider ' +
                        'significant. [default: 1]', default=1, type=int)
    parser.add_argument('--secondary-or-nominal', dest='secondary_or_nom', 
                        help='Allow genes to meet either --secondary-p-cutoff ' +
                        'or --min-nominal, but do not require both. ' +
                        '[default: require both]', default=False, action='store_true')
    parser.add_argument('--credible-sets', dest='cs_val', type=float, default=0.95,
                        help='Credible set value. [default: 0.95]')
    parser.add_argument('--regularization-alpha', dest='logit_alpha', type=float,
                        help='Regularization penalty weight for logit glm. Must ' +
                        'be in ~ [0, 1]. [default: no regularization]', default=0.2)
    parser.add_argument('--regularization-l1-l2-mix', dest='l1_l2_mix', type=float,
                        help='Regularization parameter (elastic net alpha) for ' +
                        'logit glm. 0 = L1, 1 = L2, (0, 1) = elastic net. ' +
                        '[default: L2 regularization]', default=1)
    parser.add_argument('--distance', help='Distance to pad each significant gene ' +
                        'prior to fine-mapping. [default: 1Mb]', default=1000000, 
                        type=int)
    parser.add_argument('--confident-pip', default=0.1, type=float, help='Minimum ' +
                        'PIP for classifying a gene as confidently fine-mapped. ' +
                        '[default: 0.1]')
    parser.add_argument('--very-confident-pip', default=0.9, type=float, help='Minimum ' +
                        'PIP for classifying a gene as very confidently fine-mapped. ' +
                        '[default: 0.9]')
    parser.add_argument('--known-causal-gene-lists', help='.tsv list of paths to ' +
                        '.txt lists of known causal genes. Used for estimating null ' +
                        'variance. Can be specified multiple times. [default: ' +
                        'no known causal genes]', dest='gs_list')
    parser.add_argument('--finemap-secondary', action='store_true', default=False,
                        help='Use secondary P-values for fine-mapping priors ' + 
                        '[default: use primary P-values]')
    parser.add_argument('-o', '--outfile', default='stdout', help='Output tsv of ' +
                        'final fine-mapping results for significant genes and ' + 
                        'phenotypes.')
    parser.add_argument('--sig-genes-bed', help='Output BED of significant genes ' +
                        'and their overall association statistics.')
    parser.add_argument('--sig-assoc-bed', help='Output BED of significant ' +
                        'gene-phenotype pairs and their corresponding association ' +
                        'statistics.')
    parser.add_argument('--sig-credsets-bed', help='Output BED of credible sets ' +
                        ' comprised of significant genes, and their corresponding ' +
                        'association statistics.')
    parser.add_argument('--all-genes-outfile', help='Output tsv of final ' +
                        'fine-mapping results for all genes and phenotypes.')
    parser.add_argument('--naive-outfile', help='Output tsv of naive results ' +
                        'before fine-mapping for all genes and phenotypes.')
    parser.add_argument('--genetic-outfile', help='Output tsv of genetics-only ' +
                        'fine-mapping results for all genes and phenotypes.')
    parser.add_argument('--coeffs-out', help='Output tsv of logit coefficients ' +
                        'from fine-mapping model.')
    args = parser.parse_args()

    # Open connections to output files
    if args.outfile in 'stdout - /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Process data per hpo
    hpo_data = load_all_hpos(args.statslist, args.secondary_p_cutoff, 
                             args.min_nominal, args.secondary_or_nom, 
                             args.distance, '{}_gene_block'.format(args.cnv), 
                             args.finemap_secondary)

    # Estimate null variance based on:
    #   1. most significant gene from each block
    #   2. all significant genes
    #   3. known causal genes (optional; can be multiple lists)
    Wsq = estimate_null_variance_basic(hpo_data)
    if args.gs_list is not None:
        with open(args.gs_list) as gsf:
            Wsq += estimate_null_variance_gs(gsf.read().splitlines(), args.statslist)
    Wsq = sorted(Wsq)
    print('Null variance estimates: ' + ', '.join([str(round(x, 3)) for x in Wsq]))

    # Update original finemapping results with re-estimated mean null variance
    hpo_data = update_finemap(hpo_data, np.nanmean(Wsq))

    # Write naive and/or genetics-only fine-mapping results (for ROC comparisons)
    if args.naive_outfile is not None:
        naive_outfile = open(args.naive_outfile, 'w')
        # Note: cs_val set to 1.0 fixed value here such that all genes will be assigned to their original credible set
        # rather than randomly throwing out ~5% of genes
        make_sig_genes_df(hpo_data, naive=True, sig_only=True, cs_val=1.0l).\
            rename(columns={'HPO' : '#HPO'}).\
            to_csv(naive_outfile, sep='\t', index=False, na_rep='NA')
        naive_outfile.close()

    if args.genetic_outfile is not None:
        genetic_outfile = open(args.genetic_outfile, 'w')
        make_sig_genes_df(hpo_data, sig_only=True, cs_val=args.cs_val).\
            rename(columns={'HPO' : '#HPO'}).\
            to_csv(genetic_outfile, sep='\t', index=False, na_rep='NA')
        genetic_outfile.close()

    # Perform functional fine-mapping with Bayesian model averaging
    finemap_res = [functional_finemap(hpo_data, args.gene_features, args.l1_l2_mix, 
                                      args.logit_alpha, args.cs_val, w) for w in Wsq]
    
    # Average models across Wsq priors and write to --outfile
    bms = [x[0] for x in finemap_res]
    bmavg(bms, hpo_data, outfile, sig_only=True)
    if args.all_genes_outfile is not None:
        all_genes_outfile = open(args.all_genes_outfile, 'w')
        final_sig_df = bmavg(bms, hpo_data, all_genes_outfile, cs_val=args.cs_val)

    # If optioned, average logit coefficients across models and write to --coeffs-out
    if args.coeffs_out is not None:
        coeffs_out = open(args.coeffs_out, 'w')
        coeff_tables = [x[1] for x in finemap_res]
        coeff_avg(coeff_tables, coeffs_out, args.logit_alpha)

    # Format & write final table of significant associations
    if args.sig_assoc_bed is not None:
        sig_assoc_bed = open(args.sig_assoc_bed, 'w')
        output_assoc_bed(hpo_data, final_sig_df, sig_assoc_bed,
                         args.confident_pip, args.very_confident_pip, args.cnv)

    # Read dict of N_case per HPO
    ncase_df = pd.read_csv(args.hpos_by_cohort, sep='\t').loc[:, '#HPO Total'.split()]
    ncase_df.index = ncase_df.iloc[:, 0]
    ncase_dict = ncase_df.drop(columns='#HPO').transpose().to_dict(orient='records')[0]

    # Format & write final table of significant genes
    if args.sig_genes_bed is not None:
        sig_genes_bed = open(args.sig_genes_bed, 'w')
        output_genes_bed(hpo_data, final_sig_df, sig_genes_bed, ncase_dict,
                         args.confident_pip, args.very_confident_pip, args.cnv)

    # Format & write final table of credible sets (including sig. genes only)
    if args.sig_credsets_bed is not None:
        sig_credsets_bed = open(args.sig_credsets_bed, 'w')
        output_credsets_bed(hpo_data, final_sig_df, sig_credsets_bed,
                            args.confident_pip, args.very_confident_pip, args.cnv)


if __name__ == '__main__':
    main()
