#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
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
from copy import deepcopy
from sklearn.preprocessing import scale
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod import families
from functools import reduce
from scipy.stats import norm
import networkx as nx
import random
import argparse
from sys import stdout


def format_stat(x, is_neg_log10=False, na_val=0):
    """
    Helper function to format & convert p-values as needed
    """

    if x == 'NA':
        return float(na_val)
    else:
        if is_neg_log10:
            return 10 ** -float(x)
        else:
            return float(x)


def get_sig_label(primary_p, secondary_p, n_nominal, primary_q, primary_p_cutoff, 
                  secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                  secondary_or_nominal=True, fdr_q_cutoff=0.05,
                  secondary_for_fdr=False):
    """
    Checks if a gene should be considered exome-wide or FDR significant
    """

    # Run all comparisons
    primary_p_is_sig = (primary_p < primary_p_cutoff)
    secondary_p_is_sig = (secondary_p < secondary_p_cutoff)
    n_nominal_is_sig = (n_nominal >= n_nominal_cutoff)
    fdr_q_is_sig = (primary_q < fdr_q_cutoff)

    # Determine secondary criteria
    if secondary_p_is_sig and n_nominal_is_sig:
        secondary_is_sig = True
    elif secondary_or_nominal and (secondary_p_is_sig or n_nominal_is_sig):
        secondary_is_sig = True
    else:
        secondary_is_sig = False

    # First consider exome-wide significance
    if primary_p_is_sig and secondary_is_sig:
        return 'EWS'
    
    # Second consider FDR significance
    elif fdr_q_is_sig and secondary_for_fdr and secondary_is_sig:
        return 'FDR'
    elif fdr_q_is_sig and not secondary_for_fdr:
        return 'FDR'

    # Otherwise, non-significant
    else:
        return 'NS'

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
        null_se : float or list of floats, variance of OR estimates under the null. 
                  By default, this value is set to a 5% chance that ORs are > 2, 
                  as suggested by Wakefield 2009. If multiple values are provided, 
                  outputs will be averaged across models.
    """

    genes = list(gene_priors.keys())

    if isinstance(null_variance, float):
        null_variance = [null_variance]

    finemap_res = {W : {} for W in null_variance}

    for W in null_variance:

        # Compute ABF per gene
        for gene in genes:
            # From Wakefield, 2009
            theta = gene_info[gene]['lnOR']
            se = ci2se((gene_info[gene]['lnOR_upper'], gene_info[gene]['lnOR_lower']))
            V = se ** 2
            if V > 0:
                zsq = (theta ** 2) / V
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

            finemap_res[W][gene] = {'ABF' : ABF}

        # Compute PIP per gene as fraction of total BF, adjusted for prior probs
        # As per Mahajan 2018 (T2D fine-mapping paper, Nat. Genet.)
        posteriors = {}
        for gene in genes:
            posteriors[gene] = finemap_res[W][gene]['ABF'] * gene_priors[gene]
        posterior_sum = np.sum(list(posteriors.values()))
        for gene in genes:
            finemap_res[W][gene]['PIP'] = posteriors[gene] / posterior_sum

    # Average across null variance estimates, if necessary
    if len(null_variance) == 1:
        finemap_res_bma = list(finemap_res.values())[0]
    else:
        finemap_res_bma = average_finemap_results(finemap_res)

    return finemap_res_bma


def average_finemap_results(finemap_res, ncase_df=None):
    """
    Perform Bayesian model averaging across multiple sets of finemapping results

    By default, computes arithmetic mean of PIPs

    For averaging across HPOs, computes weighted mean of ABFs and PIPs while 
    weighting by sqrt(N_cases)
    """

    genes = set()
    for subres in finemap_res.values():
        genes.update(subres.keys())

    finemap_res_bma = {g : {'ABF' : [], 'PIP': [], 'weights': []} for g in genes}
    
    for key, subres in finemap_res.items():

        # Collect information for all genes, and fill with zeros if not included in credset
        for gene in genes:
            for field in 'ABF PIP'.split():
                finemap_res_bma[gene][field].append(subres.get(gene, {field : 0})[field])

            # Average PIPs by sqrt(N_case) if ncase_df supplied
            if ncase_df is not None:
                # Assumes HPO for each individual result is the first underscore-delimited element of the key
                hpo = key.split('_')[0]
                weight = np.sqrt(ncase_df.Total[ncase_df['#HPO'] == hpo].values[0])

            # Otherwise, use uniform weights
            else:
                weight = 1

            finemap_res_bma[gene]['weights'].append(weight)

    # Compute weighted PIP per gene
    for gene in genes:
        for field in 'ABF PIP'.split():
            newval = np.average(finemap_res_bma[gene][field], 
                                weights=finemap_res_bma[gene]['weights'])
            finemap_res_bma[gene][field] = newval
        finemap_res_bma[gene].pop('weights')

    # Re-scale PIPs after weighting
    pip_sum = np.nansum([x['PIP'] for x in finemap_res_bma.values()])
    for gene in genes:
        finemap_res_bma[gene]['PIP'] = finemap_res_bma[gene]['PIP'] / pip_sum

    return finemap_res_bma


def parse_stats(stats_in, primary_p_cutoff, p_is_neg_log10=True, 
                secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                secondary_or_nominal=True, fdr_q_cutoff=0.05, 
                secondary_for_fdr=False, sig_only=False, keep_genes=None, 
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

    for chrom, start, end, gene, n_nominal, top_cohort, excluded_cohorts, \
        case_freq, control_freq, lnOR, lnOR_lower, lnOR_upper, zscore, \
        primary_p, primary_q, secondary_lnOR, secondary_lnOR_lower, \
        secondary_lnOR_upper, secondary_zscore, secondary_p, secondary_q \
        in reader:

        # Skip header line
        if chrom.startswith('#'):
            continue

        # If optioned, restrict on gene name
        if keep_genes is not None:
            if gene not in keep_genes:
                continue

        # Clean up gene data
        primary_p = format_stat(primary_p, p_is_neg_log10, 1)
        primary_q = format_stat(primary_q, p_is_neg_log10, 1)
        secondary_p = format_stat(secondary_p, p_is_neg_log10, 1)
        n_nominal = int(n_nominal)
        sig_label = get_sig_label(primary_p, secondary_p, n_nominal, primary_q, 
                                  primary_p_cutoff, secondary_p_cutoff, 
                                  n_nominal_cutoff, secondary_or_nominal, 
                                  fdr_q_cutoff, secondary_for_fdr)
        ew_sig = False
        fdr_sig = False
        if sig_label == 'EWS':
            ew_sig = True
            fdr_sig = True
        elif sig_label == 'FDR':
            fdr_sig = True
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
        
        gene_stats = {'case_freq' : case_freq, 'control_freq' : control_freq,
                      'lnOR' : use_lnOR, 
                      'lnOR_lower' : use_lnOR_lower,
                      'lnOR_upper' : use_lnOR_upper, 
                      'zscore' : use_zscore,
                      'primary_p' : primary_p, 'primary_q' : primary_q, 
                      'secondary_p' : secondary_p, 'n_nominal' : n_nominal, 
                      'fdr_q' : primary_q, 'ew_sig' : ew_sig,
                      'fdr_sig' : fdr_sig}

        # Store gene association stats
        if sig_only:
            if sig_label in 'EWS FDR'.split():
                gene_bt = pbt.BedTool('\t'.join([chrom, start, end, gene]), 
                                      from_string=True)
                gene_stats['gene_bt'] = gene_bt
                stats_dict[gene] = gene_stats
        else:
            gene_bt = pbt.BedTool('\t'.join([chrom, start, end, gene]), 
                                  from_string=True)
            gene_stats['gene_bt'] = gene_bt
            stats_dict[gene] = gene_stats

    csvin.close()

    return stats_dict


def _get_genes_in_ld(baits_list, cov_df, min_covariance=0.2):
    """
    Extract a list of all genes with covariance >= min_covariance
    vs. at least one gene in baits_list
    """

    hits = (cov_df.jaccard >= min_covariance) & \
           ((cov_df.geneA.isin(baits_list)) | (cov_df.geneB.isin(baits_list)))
    ld_friends = set(list(cov_df.geneA[hits]) + list(cov_df.geneB[hits]))
    ld_friends.update(baits_list)

    return sorted(list(ld_friends))


def cluster_genes_by_cov(hpo_info, all_genes_bt, cov_df, min_covariance=0.2, 
                         block_merge_dist=1000000):
    """
    Make clusters of genes for fine-mapping based on covariance and distance
    """

    sig_genes = list(hpo_info['sig_genes'].keys())

    # Make graph of all significant genes
    gg = nx.Graph()
    gg.add_nodes_from(sig_genes)

    # Add edges to all genes "in LD" with at least one significant gene
    for sg in sig_genes:
        for og in _get_genes_in_ld([sg], cov_df, min_covariance):
            if og != sg:
                if og not in gg.nodes:
                    gg.add_node(og)
                gg.add_edge(sg, og)

    # Add edges to all genes within block_merge_dist of a gene already in gg
    genes_in_gg_bt = all_genes_bt.filter(lambda x: x[3] in gg.nodes).saveas()
    other_genes_bt = all_genes_bt.filter(lambda x: x[3] not in gg.nodes).saveas()
    other_near_gg = other_genes_bt.closest(genes_in_gg_bt, d=True, 
                                           t='first', k=1, N=True).\
                                   filter(lambda x: float(x[-1]) <= block_merge_dist and \
                                                    float(x[-1]) >= 0).\
                                   cut(range(4))
    for hit in other_near_gg.closest(genes_in_gg_bt, d=True, k=10e10, N=True):
        dist = float(hit[-1])
        og = hit[3]
        sg = hit[7]
        if dist <= block_merge_dist and dist >= 0:
            if og not in gg.nodes:
                gg.add_node(og)
            gg.add_edge(og, sg)

    # Convert clustered genes into a pbt.BedTool of blocks
    blocks_bt_str = ''
    for cluster in nx.connected_components(gg):
        # Only process clusters that have at least one significant gene
        if len(cluster.intersection(set(sig_genes))) == 0:
            continue

        # Get chromosome and min/max coordinates
        gdf = all_genes_bt.filter(lambda x: x[3] in cluster).saveas().\
                           to_dataframe(names='chrom start end gene'.split())
        chrom = str(gdf.chrom[0])
        start = str(gdf.start.min())
        end = str(gdf.end.max())
        genes_str = ','.join(cluster)
        blocks_bt_str += '\t'.join([chrom, start, end, genes_str]) + '\n'

    return pbt.BedTool(blocks_bt_str, from_string=True).sort()


def process_hpo(hpo, stats_in, primary_p_cutoff, p_is_neg_log10=True, 
                secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                secondary_or_nominal=True, fdr_q_cutoff=0.05, 
                secondary_for_fdr=False, block_merge_dist=1000000, 
                nonsig_distance=1000000, cov_df=None, min_covariance=0.2,
                block_prefix='gene_block', null_variance=0.42 ** 2, 
                finemap_secondary=False, include_non_sig=True,
                force_load_genes=[]):
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

    # First pass: parse data for exome-wide or FDR significant genes only
    hpo_info['sig_genes'] = parse_stats(stats_in, primary_p_cutoff, p_is_neg_log10, 
                                        secondary_p_cutoff, n_nominal_cutoff, 
                                        secondary_or_nominal, fdr_q_cutoff, 
                                        secondary_for_fdr, sig_only=True, 
                                        finemap_secondary=finemap_secondary)

    # Second pass: parse data for all genes within nonsig_distance of sig_genes
    # (Or, if specified, will treat all genes with at least min_covariance of a 
    #  significant gene as significant before loading more gene stats)
    if len(hpo_info['sig_genes']) > 0:

        # Make pbt.BedTool of gene coordinates, for reference
        all_genes_bt = pbt.BedTool(stats_in).cut(range(4)).sort()

        if include_non_sig:
            # Make bt of significant genes
            sig_gene_bts = [g['gene_bt'] for g in hpo_info['sig_genes'].values()]
            if len(sig_gene_bts) > 1:
                sig_genes_bt = sig_gene_bts[0].cat(*sig_gene_bts[1:], postmerge=False).sort()
            else:
                sig_genes_bt = sig_gene_bts[0]

            # Intersect sig genes with all genes
            if cov_df is not None and len(sig_gene_bts) > 0:
                # If CNV covariance is provided, supplement sig_genes_bt with all
                # other genes >= min_covariance with at least one sig gene
                sig_genes = list(hpo_info['sig_genes'].keys())
                all_cov_genes = _get_genes_in_ld(sig_genes, cov_df, min_covariance)
                sig_genes_bt = all_genes_bt.filter(lambda x: x[3] in all_cov_genes).saveas()

            nearby_genes = all_genes_bt.closest(sig_genes_bt.sort(), d=True).\
                               filter(lambda x: int(x[8]) > -1 and \
                                                int(x[8]) <= nonsig_distance).\
                               saveas().to_dataframe().loc[:, 'name'].values.tolist()
            nearby_genes = list(set(nearby_genes + list(force_load_genes))) # Deduplicates

            # Gather gene stats
            hpo_info['all_genes'] = parse_stats(stats_in, primary_p_cutoff, p_is_neg_log10, 
                                                secondary_p_cutoff, n_nominal_cutoff, 
                                                secondary_or_nominal, fdr_q_cutoff, 
                                                secondary_for_fdr, sig_only=False, 
                                                keep_genes=nearby_genes, 
                                                finemap_secondary=finemap_secondary)

        # If --no-non-significant is passed, do not include NS genes when fine-mapping
        else:
            hpo_info['all_genes'] = hpo_info['sig_genes'].copy()

        # Cluster genes into blocks to be fine-mapped
        gene_bts = [g['gene_bt'] for g in hpo_info['all_genes'].values()]
        if len(gene_bts) > 1:
            genes_bt = gene_bts[0].cat(*gene_bts[1:], postmerge=False).sort()
        else:
            genes_bt = gene_bts[0]
        # Cluster by distance or covariance if cov_df provided
        if cov_df is None:
            blocks = genes_bt.merge(d=block_merge_dist, c=4, o='distinct')
        else:
            blocks = cluster_genes_by_cov(hpo_info, all_genes_bt, cov_df, min_covariance,
                                          block_merge_dist)

        # Ensure each block has at least one significant gene
        sig_set = set(hpo_info['sig_genes'].keys())
        blocks = blocks.filter(lambda x: len(set(x[3].split(',')).intersection(sig_set)) > 0)

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
                  secondary_or_nominal=True, fdr_q_cutoff=0.05, 
                  secondary_for_fdr=False, block_merge_dist=1000000, 
                  nonsig_distance=1000000, cov_df=None, min_covariance=0.2,
                  block_prefix='gene_block', finemap_secondary=False, 
                  include_non_sig=True, force_load_genes=[], verbose=True):
    """
    Wrapper function to process each HPO with process_hpo()
    Returns a dict with one entry per HPO
    """

    hpo_data = {}

    with open(statslist) as infile:
        reader = csv.reader(infile, delimiter='\t')
        for hpo, stats_in, pval, in reader:
            if verbose:
                print('Loading association statistics from {}'.format(hpo))
            primary_p_cutoff = float(pval)
            hpo_data[hpo] = process_hpo(hpo, stats_in, primary_p_cutoff, 
                                        p_is_neg_log10=True, 
                                        secondary_p_cutoff=secondary_p_cutoff, 
                                        n_nominal_cutoff=n_nominal_cutoff, 
                                        secondary_or_nominal=secondary_or_nominal,
                                        fdr_q_cutoff=fdr_q_cutoff,
                                        secondary_for_fdr=secondary_for_fdr,
                                        block_merge_dist=block_merge_dist,
                                        nonsig_distance=nonsig_distance,
                                        cov_df=cov_df,
                                        min_covariance=min_covariance,
                                        block_prefix=block_prefix,
                                        finemap_secondary=finemap_secondary,
                                        include_non_sig=include_non_sig,
                                        force_load_genes=force_load_genes)

    return hpo_data


def estimate_null_variance_basic(hpo_data):
    """
    Estimates null variance per phenotype from average of all significant genes (exome-wide + FDR)
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
    for gs_genes in gs_lists.values():
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


def functional_finemap(hpo_data_orig, gene_features_in, l1_l2_mix, logit_alpha, 
                       cs_val=1.0, null_variance=0.42 ** 2, exclude_genes=None,
                       use_max_pip=False, converge_rmse=10e-8, max_iter=100, 
                       quiet=False):
    """
    Conduct E-M optimized functional fine-mapping for all gene blocks & HPOs
    Two returns:
        - pd.DataFrame of HPO data updated with new finemapping results
        - statsmodels.iolib.table.SimpleTable of final logit coefficients
    """

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
    all_genes = make_sig_genes_df(hpo_data_orig, cs_val=1.0).gene.tolist()
    features.iloc[:, 1:] = scale(features.iloc[:, 1:])
    features = features.loc[features.gene.isin(all_genes), :]

    # Load list of genes to be excluded from training, if optioned
    if exclude_genes is not None:
        with open(exclude_genes) as xin:
            xgenes = list(set([g.rstrip() for g in xin.readlines()]))
    else:
        xgenes = []

    # Iterate over all null variance values provided and finemap each separately
    if isinstance(null_variance, float):
        null_variance = [null_variance]

    finemap_res = {}

    for W in null_variance:

        # Make a deep copy of HPO data prior to finemapping
        hpo_data = deepcopy(hpo_data_orig)

        if not quiet:
            msg = '\nStarting fine-mapping with null variance (W) = {:.5}\n' + \
                  '  Iter.\tPIP RMSE\tCoeff. RMSE'
            print(msg.format(W))

        # Re-finemap each block with specified null variance (for BMA)
        for hpo in hpo_data.keys():
            for block_id, block_data in hpo_data[hpo]['blocks'].items():
                genes = list(block_data['finemap_res'].keys())
                gene_priors = {gene : 1 / len(genes) for gene in genes}
                updated_finemap_res = finemap(gene_priors, hpo_data[hpo]['all_genes'], W)
                hpo_data[hpo]['blocks'][block_id]['finemap_res'] = updated_finemap_res

        # Iterate until convergence
        coeffs = [0 for x in features.columns.tolist()[1:]]
        k = 0
        rmse_PIP, rmse_coeffs, prev_rmse_PIP, prev_rmse_coeffs = 100, 100, 100, 100
        while rmse_PIP >= converge_rmse or rmse_coeffs >= converge_rmse:
            k += 1
            # Join sig_df with features for logistic regression
            sig_df = make_sig_genes_df(hpo_data, cs_val=1.0)
            logit_df = sig_df.loc[:, 'gene PIP'.split()].\
                              merge(features, how='left', on='gene')
            logit_df['PIP'] = logit_df['PIP'].astype('float64')

            # Drop excluded genes from training data
            logit_df = logit_df.loc[~logit_df['gene'].isin(xgenes), :]

            # When fitting regression, take mean (or max) of PIPs for genes appearing multiple times
            # This can happen due to multiple HPO associations with the same gene
            if use_max_pip:
                logit_df = logit_df.groupby('gene').max()
            else:
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
                    new_finemap_res = finemap(gene_priors, hpo_data[hpo]['all_genes'], W)
                    hpo_data[hpo]['blocks'][block_id]['finemap_res'] = new_finemap_res

            # Compute RMSE for PIPs and coefficients
            new_sig_df = make_sig_genes_df(hpo_data, cs_val=1.0)
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

            # Check for convergence above RMSE target by comparing RMSE from previous iteration to this iteration
            d_rmse_PIP = prev_rmse_PIP - rmse_PIP
            d_rmse_coeffs = prev_rmse_coeffs - rmse_coeffs

            if d_rmse_PIP <= converge_rmse and d_rmse_coeffs < converge_rmse:
                print('RMSE stablized below target after {:,} iterations'.format(k))
                break

            else:
                prev_rmse_PIP = rmse_PIP
                prev_rmse_coeffs = rmse_coeffs

            # Force exit from E-M after max_iter iterations
            if k > max_iter:
                print('Failed to converge after {:,} iterations'.format(k))
                break

        # Report completion & store results prior to BMA
        print('Finished after {:,} iterations'.format(k))
        if logit_alpha is None:
            tab_out = logit.summary().tables[1]
        else:
            tab_out = logit.params
        finemap_res[W] = {'hpo_data' : hpo_data, 'coeffs' : tab_out}

    # Update original copy of hpo_data by averaging results for each block 
    # across null variance estimates
    for hpo, hdat in hpo_data_orig.items():
        for bid in hdat['blocks'].keys():
            finemap_res_dict = {W : finemap_res[W]['hpo_data'][hpo]['blocks'][bid]['finemap_res'] \
                                for W in null_variance}
            finemap_res_avg = average_finemap_results(finemap_res_dict)
            hpo_data_orig[hpo]['blocks'][bid]['finemap_res'] = finemap_res_avg

    # Average logit coefficients
    coeff_tables = [finemap_res[W]['coeffs'] for W in null_variance]
    logit_coeffs = coeff_avg(coeff_tables, logit_alpha)

    return hpo_data_orig, logit_coeffs


def coeff_avg(coeff_tables, logit_alpha):
    """
    Average logit coefficients across all models
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

    return avg_df


def output_assoc_bed(hpo_data, sig_df, prejoint_sig_df, outfile, 
                     conf_pip=0.1, vconf_pip=0.9, cnv='NS'):
    """
    Format final list of significant gene-phenotype pairs with summary statistics
    """

    cols = 'chr start end gene cnv hpo sig_level control_freq case_freq ln_or ' + \
           'ln_or_ci_lower ln_or_ci_upper pvalue pip_final pip_HPO_specific ' + \
           'top_gene credible_set_id'
    outfile.write('#' + '\t'.join(cols.split()) + '\n')

    outlines_presort = []

    # Iterate over each significant gene with PIP >= conf_pip
    for gdict in sig_df.loc[sig_df.PIP >= conf_pip, :].to_dict(orient='index').values():
        # Skip if gene wasn't at least FDR significant in original analysis
        hpo = gdict['HPO']
        gene = gdict['gene']
        if not hpo_data[hpo]['all_genes'][gene]['fdr_sig']:
            continue

        # Get gene info
        pip = gdict['PIP']
        prejoint_pip = prejoint_sig_df.PIP[(prejoint_sig_df.gene == gene) \
                                            & (prejoint_sig_df.HPO == hpo)].\
                                       values[0]
        cred = gdict['credible_set']
        finemap_res = hpo_data[hpo]['blocks'][cred]['finemap_res']
        best_pip_in_cs = np.nanmax([x['PIP'] for x in finemap_res.values()])
        top_gene = int(pip == best_pip_in_cs)
        chrom = str(hpo_data[hpo]['all_genes'][gene]['gene_bt'][0].chrom)
        start = str(hpo_data[hpo]['all_genes'][gene]['gene_bt'][0].start)
        end = str(hpo_data[hpo]['all_genes'][gene]['gene_bt'][0].end)
        if hpo_data[hpo]['all_genes'][gene]['ew_sig']:
            sig_level = 'exome_wide'
        elif hpo_data[hpo]['all_genes'][gene]['fdr_sig']:
            sig_level = 'FDR'
        else:
            sig_level = 'not_significant'
        case_freq = hpo_data[hpo]['all_genes'][gene]['case_freq']
        control_freq = hpo_data[hpo]['all_genes'][gene]['control_freq']
        lnor = hpo_data[hpo]['all_genes'][gene]['lnOR']
        lnor_lower = hpo_data[hpo]['all_genes'][gene]['lnOR_lower']
        lnor_upper = hpo_data[hpo]['all_genes'][gene]['lnOR_upper']
        pval = hpo_data[hpo]['all_genes'][gene]['primary_p']
        # If p-value was rounded to zero, modify by taking the unadjusted normal P-value
        if pval == 0:
            pval = norm.sf(hpo_data[hpo]['all_genes'][gene]['zscore'])

        # Prepare line to be written out and add keys to be sorted
        outline = '\t'.join([chrom, start, end, gene, cnv, hpo, sig_level])
        outnums_fmt = '\t{:.3E}\t{:.3E}\t{:.3}\t{:.3}\t{:.3}\t{:.3E}\t{:.3}\t{:.3}\t{}'
        outline += outnums_fmt.format(control_freq, case_freq, lnor, lnor_lower, 
                                      lnor_upper, pval, pip, prejoint_pip, top_gene)
        outline += '\t' + cred + '\n'
        outlines_presort.append([int(chrom), int(start), hpo, outline])
        
    # Sort lines by coordinate and write to outfile
    for vals in sorted(outlines_presort, key=lambda x: x[:3]):
        outfile.write(vals[3])
    outfile.close()


def iv_mean(values, variances, conf=0.95):
    """
    Returns inverse-variance weighted mean of values and conf% confidence interval
    """

    # Drop values where variance == 0
    keep_idx = [i for i, v in enumerate(variances) if v > 0]
    values = [v for i, v in enumerate(values) if i in keep_idx]
    variances = [v for i, v in enumerate(variances) if i in keep_idx]

    weights = [1 / v for v in variances]

    wsum = np.nansum(weights)

    numerator = np.nansum([x / v for x, v in zip(values, variances)])

    ivm = numerator / wsum

    pooled_se = np.sqrt(1 / wsum)

    ci_dist = norm.ppf(conf) * pooled_se

    return ivm, (ivm - ci_dist, ivm + ci_dist)


def output_credsets_bed(hpo_data, sig_df, outfile, conf_pip=0.1, vconf_pip=0.9, cnv='NS'):
    """
    Format final list of credible sets with summary statistics
    """
    
    cols = 'chr start end credible_set_id cnv hpo best_sig_level mean_control_freq ' + \
           'mean_case_freq pooled_ln_or pooled_ln_or_ci_lower pooled_ln_or_ci_upper ' + \
           'best_pvalue n_genes all_genes n_sig_genes sig_genes top_gene ' +\
           'vconf_genes conf_genes'
    outfile.write('#' + '\t'.join(cols.split()) + '\n')

    outlines_presort = []

    for cred in [x for x in set(sig_df.credible_set.tolist()) if x is not None]:
        # Get basic credible set info
        cred_df = sig_df.loc[sig_df.credible_set == cred, :]
        hpo = cred_df.HPO.tolist()[0]
        hpo_res_df = sig_df[sig_df.HPO == hpo]
        genes = sorted(list(set(cred_df.gene.tolist())))

        # Gather lists of sig genes *IN* credible set and list of all sig genes
        # (including those who were fine-mapped out of the credible set)
        sig_genes = [g for g in genes if g in hpo_data[hpo]['sig_genes'].keys()]
        all_sig_genes = [g for g in hpo_data[hpo]['blocks'][cred]['finemap_res'].keys() \
                         if g in hpo_data[hpo]['sig_genes'].keys()]
        if any([hpo_data[hpo]['all_genes'][g]['ew_sig'] for g in genes]):
            best_sig_level = 'exome_wide'
        elif any([hpo_data[hpo]['all_genes'][g]['fdr_sig'] for g in genes]):
            best_sig_level = 'FDR'
        else:
            best_sig_level = 'not_significant'

        # Candidate genes must be significantly associated with HPO in question
        vconf_genes = sorted(list(set(cred_df.gene[(cred_df.gene.isin(sig_genes)) & \
                                                   (cred_df.PIP >= vconf_pip)].tolist())))
        conf_genes = sorted(list(set(cred_df.gene[(cred_df.gene.isin(sig_genes)) & \
                                                  (cred_df.PIP < vconf_pip) & \
                                                  (cred_df.PIP >= conf_pip)].tolist())))
        top_pip = np.nanmax(hpo_res_df.PIP[(hpo_res_df.gene.isin(all_sig_genes))])
        top_gene = sorted(list(set(hpo_res_df.gene[hpo_res_df.PIP == top_pip].tolist())))
        n_genes = len(genes)
        n_sig_genes = len(sig_genes)
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

        # Prepare lines to be sorted and written to outfile
        outline = '\t'.join([chrom, start, end, cred, cnv, hpo, best_sig_level])
        outnums_fmt = '\t{:.3E}\t{:.3E}\t{:.3}\t{:.3}\t{:.3}\t{:.3E}'
        outline += outnums_fmt.format(control_freq, case_freq, lnor, lnor_lower, 
                                      lnor_upper, best_p)
        outline += '\t' + '\t'.join([str(n_genes), ';'.join(genes), str(n_sig_genes),
                                     ';'.join(sig_genes), ';'.join(top_gene),
                                     ';'.join(vconf_genes), ';'.join(conf_genes)]) + '\n'
        outlines_presort.append([int(chrom), int(start), hpo, outline])
        
    # Sort lines by coordinate and write to outfile
    for vals in sorted(outlines_presort, key=lambda x: x[:3]):
        outfile.write(vals[3])
    outfile.close()


def cluster_credsets(hpo_data, sig_only=False, seed=2021):
    """
    Link credible sets into clusters based on overlapping genes
    Returns a list of clusters
    """

    # Build dictionary of genes per credset
    sig_df = make_sig_genes_df(hpo_data, naive=False, sig_only=sig_only, cs_val=1.0)
    cs_ids = set(sig_df.credible_set[~sig_df.credible_set.isna()].to_list())
    cs_dict = {cs : set(sig_df.gene[sig_df.credible_set == cs].to_list()) for cs in cs_ids}

    # Add each credible set as a node to nx.Graph()
    # Add edge between any two nodes that share at least one gene in common
    csg = nx.Graph()
    for csA, genesA in cs_dict.items():
        csg.add_node(csA)
        for csB in csg.nodes():
            if csA == csB:
                continue
            if len(cs_dict[csA].intersection(cs_dict[csB])) > 0:
                csg.add_edge(csA, csB)

    random.seed(seed)

    return sorted([sorted(x) for x in nx.connected_components(csg)])


def joint_finemap(hpo_data, cred_clusters, ncase_df, seed=2021):
    """
    Average finemapping results across all HPOs per cluster of credsets
    """

    random.seed(seed)

    itervals = [(i, cluster) for i, cluster in enumerate(cred_clusters)]
    removed_cs = {}

    for i, cluster in itervals:

        # Get basic info for each credset in cluster
        cs_info = {cs : {'hpo' : cs.split('_')[0]} for cs in cluster}
        for cs in cluster:
            fres = hpo_data[cs_info[cs]['hpo']]['blocks'][cs]['finemap_res']
            cs_info[cs]['finemap_res'] = fres

        # Re-finemap each cluster based on average ABFs across all credible sets
        all_finemap_res = {cs : info['finemap_res'] for cs, info in cs_info.items()}
        avg_finemap_res = average_finemap_results(all_finemap_res, ncase_df)

        # Update each credible set in hpo_data with new finemapping results
        used_hpos = set()
        for cs in cluster:
            hpo = cs_info[cs]['hpo']
            if hpo in used_hpos:
                if i in removed_cs.keys():
                    removed_cs[i].add(cs)
                else:
                    removed_cs[i] = set([cs])
            else:
                hpo_data[hpo]['blocks'][cs]['finemap_res'] = avg_finemap_res
                used_hpos.add(hpo)

    # Remove redundant credible sets
    for i, cs_list in removed_cs.items():
        for cs in cs_list:
            hpo = cs.split('_')[0]
            hpo_data[hpo]['blocks'].pop(cs)
            cred_clusters[i].remove(cs)

    return hpo_data, cred_clusters


def output_joint_credsets_bed(hpo_data, sig_df, cs_clusters, genes_bt, outfile, 
                              ncase_dict, conf_pip=0.1, vconf_pip=0.9, 
                              cnv='NS', prefix=None):
    """
    Format final list of credible sets collapsed across HPOs with summary statistics
    """
    
    cols = 'chr start end credible_set_id cnv best_sig_level mean_control_freq ' + \
           'mean_case_freq pooled_ln_or pooled_ln_or_ci_lower pooled_ln_or_ci_upper ' + \
           'best_pvalue n_genes all_genes n_sig_genes sig_genes top_gene ' + \
           'vconf_genes conf_genes n_hpos hpos'
    cols = cols.split()

    if prefix is None:
        prefix = cnv + '_gene_block'

    out_df = pd.DataFrame(columns=cols)

    for cluster in cs_clusters:
        
        # Get basic credible set info
        cs_hpo_dict = {cs: cs.split('_')[0] for cs in cluster}
        hpos = set(cs_hpo_dict.values())
        n_hpos = len(hpos)

        # Assign a single member of the cluster as a "proxy" (since all finemap results have already been averaged)
        all_cs_data = {cs : hpo_data[hpo]['blocks'].get(cs, None) for cs, hpo in cs_hpo_dict.items()}
        cs_proxy, cs_proxy_data = [(cs, dat) for cs, dat in all_cs_data.items() \
                                   if any(sig_df.credible_set == cs)][0]

        # Slice sig_df to supserset of all genes from any HPO
        # This is necessary since some genes will be significant in certain HPOs but not others
        cred_df = sig_df[sig_df.credible_set.isin(cluster)]['gene PIP'.split()].drop_duplicates()

        # Use proxy to get more credible set info
        genes = sorted(list(set(cred_df.gene.tolist())))
        n_genes = len(genes)
        sig_genes = set()
        for hpo in hpos:
            sig_genes.update(set(hpo_data[hpo]['sig_genes'].keys()).intersection(set(genes)))
        n_sig_genes = len(sig_genes)
        all_sig_genes = set()
        for hpo in hpos:
            all_sig_genes.update([g for g in hpo_data[hpo]['sig_genes'].keys() \
                                  if g in cs_proxy_data['finemap_res'].keys()])

        # Candidate genes must be significant in at least one HPO
        vconf_genes = sorted(list(set(cred_df.gene[(cred_df.gene.isin(sig_genes)) & \
                                                   (cred_df.PIP >= vconf_pip)].tolist())))
        conf_genes = sorted(list(set(cred_df.gene[(cred_df.gene.isin(sig_genes)) & \
                                                  (cred_df.PIP < vconf_pip) & \
                                                  (cred_df.PIP >= conf_pip)].tolist())))
        top_pip = np.nanmax(sig_df.PIP[(sig_df.gene.isin(all_sig_genes)) & \
                                        (sig_df.HPO.isin(hpos))])
        top_gene = sorted(list(set(sig_df.gene[(sig_df.PIP == top_pip) & \
                                               (sig_df.gene.isin(all_sig_genes) & \
                                               (sig_df.HPO.isin(hpos)))].tolist())))
        cred_genes_bt = genes_bt.filter(lambda x: x.name in genes).saveas()
        chrom = str(cred_genes_bt[0].chrom)
        start = str(np.nanmin([x.start for x in cred_genes_bt]))
        end = str(np.nanmax([x.end for x in cred_genes_bt]))

        # Build pd.DataFrame of all gene stats across all HPOs
        gdat = None
        for hpo in hpos:
            gdat_h = {}
            for gene in genes:
                gdat_h_i = hpo_data[hpo]['all_genes'].get(gene, None)
                if gdat_h_i is not None:
                    gdat_h[gene] = {k : v for k, v in gdat_h_i.items() if k != 'gene_bt'}
            gdf_h = pd.DataFrame.from_dict(gdat_h, orient='index')
            gdf_h['gene'] = gdf_h.index
            gdf_h['HPO'] = hpo
            if gdat is None:
                gdat = gdf_h.copy()
            else:
                gdat = gdat.append(gdf_h, ignore_index=True)

        # Get summary information across all HPOs in cluster
        if any(gdat.ew_sig):
            best_sig_level = 'exome_wide'
        elif any(gdat.fdr_sig):
            best_sig_level = 'FDR'
        else:
            best_sig_level = 'not_significant'
        best_p = np.nanmin(gdat.primary_p)
        if best_p == 0:
            best_z = np.nanmax(gdat.zscore)
            best_p = norm.sf(best_z)
        control_freq = np.nanmean(gdat.control_freq)

        # Compute pooled case frequency as mean weighted by np.sqrt(N_cases)
        cfreqs = gdat.case_freq
        cfreq_weights = np.sqrt(gdat.HPO.map(ncase_dict))
        case_freq = np.average(cfreqs, weights=cfreq_weights)

        # Get pooled effect size as inverse-variance weighted mean of all genes across all HPOs
        gors = gdat.lnOR.tolist()
        gvars = [ci2se(ci) ** 2 for ci in \
                 gdat.loc[:, 'lnOR_lower lnOR_upper'.split()].itertuples(index=False)]
        lnor, lnor_ci = iv_mean(gors, gvars)
        lnor_lower, lnor_upper = sorted(lnor_ci)

        # Reassign empty variables for writing to outfile, if necessary
        if len(vconf_genes) == 0:
            vconf_genes = ['NA']
        if len(conf_genes) == 0:
            conf_genes = ['NA']
        if len(top_gene) == 0:
            top_gene = ['NA']

        # Add stats to out_df
        outvals = [chrom, start, end, 'tmp', cnv, best_sig_level, control_freq, 
                   case_freq, lnor, lnor_lower, lnor_upper, best_p, n_genes, 
                   ';'.join(sorted(genes)), n_sig_genes, ';'.join(sorted(sig_genes)), 
                   ';'.join(sorted(top_gene)), ';'.join(sorted(vconf_genes)), 
                   ';'.join(sorted(conf_genes)), n_hpos, ';'.join(sorted(hpos))]
        out_df = out_df.append(pd.DataFrame([outvals], columns=cols), ignore_index=True)

    # Sort entries in out_df and write to file
    for col in 'chr start end'.split():
        out_df[col] = pd.to_numeric(out_df[col])
    out_df.sort_values(by=['chr', 'start'], inplace=True)
    out_df['credible_set_id'] = ['{}_{}'.format(prefix, i + 1) for i in range(len(out_df))]
    out_df.rename(columns={'chr' : '#chr'}).\
           to_csv(outfile, index=False, na_rep='.', sep='\t')
    outfile.close()

    # Return dict mapping genes to merged credsets for writing out final gene file
    gene_cs_map = {}
    gcmap_vals = out_df['credible_set_id all_genes top_gene'.split()].\
                       itertuples(index=False, name=None)
    for cs_id, gstr, topstr in gcmap_vals:
        top_genes = topstr.split(';')
        for gene in gstr.split(';'):
            gene_cs_map[gene] = {'credset' : cs_id, 'top' : gene in top_genes}

    all_genes = set()
    for hdat in hpo_data.values():
        for bdat in hdat['blocks'].values():
            all_genes.update(set(bdat['finemap_res'].keys()))

    return gene_cs_map


def output_genes_bed(hpo_data, sig_df, prejoint_sig_df, joint_gene_cs_map, 
                     outfile, ncase_dict, conf_pip=0.1, vconf_pip=0.9, cnv='NS'):
    """
    Format final list of significant genes and compute pooled summary statistics
    """

    cols = 'chr start end gene cnv best_sig_level pooled_control_freq pooled_case_freq ' + \
           'pooled_ln_or pooled_ln_or_ci_lower pooled_ln_or_ci_upper ' + \
           'pip_final top_gene mean_pip_HPO_specific n_hpos all_hpos vconf_hpos ' + \
           'conf_hpos credible_set_ids'
    cols = cols.split()

    out_df = pd.DataFrame(columns=cols)

    for gene in set(sig_df.loc[sig_df.PIP >= conf_pip, ].gene.tolist()):
        # Get gene info
        genedf = sig_df.loc[(sig_df.gene == gene) & (sig_df.PIP >= conf_pip), :]
        all_hpos = sorted(genedf.HPO.tolist())
        if any([hpo_data[h]['all_genes'][gene]['ew_sig'] for h in all_hpos]):
            best_sig_level = 'exome_wide'
        elif any([hpo_data[h]['all_genes'][gene]['fdr_sig'] for h in all_hpos]):
            best_sig_level = 'FDR'
        else:
            best_sig_level = 'not_significant'
        n_hpos = len(all_hpos)
        prejoint_genedf = prejoint_sig_df.loc[(prejoint_sig_df.gene == gene) \
                                              & (prejoint_sig_df.PIP >= conf_pip), :]
        vconf_hpos = sorted(prejoint_genedf.HPO[prejoint_genedf.PIP >= vconf_pip].tolist())
        conf_hpos = sorted(prejoint_genedf.HPO[prejoint_genedf.PIP < vconf_pip].tolist())
        chrom = str(hpo_data[all_hpos[0]]['all_genes'][gene]['gene_bt'][0].chrom)
        start = str(hpo_data[all_hpos[0]]['all_genes'][gene]['gene_bt'][0].start)
        end = str(hpo_data[all_hpos[0]]['all_genes'][gene]['gene_bt'][0].end)
        control_freq = hpo_data[all_hpos[0]]['all_genes'][gene]['control_freq']
        pip = np.nanmax(genedf.PIP)
        top_gene = int(joint_gene_cs_map[gene]['top'])
        mean_prejoint_pip = np.nanmean(prejoint_genedf.PIP)
        cred_ids = [joint_gene_cs_map[gene]['credset']] + sorted(genedf.credible_set.tolist())

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

        # Add gene values to out_df
        outvals = [chrom, start, end, gene, cnv, best_sig_level, control_freq, 
                   case_freq, lnor, lnor_lower, lnor_upper, pip, top_gene, 
                   mean_prejoint_pip, str(n_hpos), ';'.join(all_hpos), 
                   ';'.join(vconf_hpos), ';'.join(conf_hpos), ';'.join(cred_ids)]
        out_df = out_df.append(pd.DataFrame([outvals], columns=cols), ignore_index=True)

    # Sort entries in out_df and write to file
    for col in 'chr start end'.split():
        out_df[col] = pd.to_numeric(out_df[col])
    out_df.sort_values(by=['chr', 'start'], inplace=True)
    out_df.rename(columns={'chr' : '#chr'}).\
           to_csv(outfile, index=False, na_rep='.', sep='\t')
    outfile.close()


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('statslist', help='.tsv of metadata per phenotype. Three ' +
                        'required columns: HPO, path to meta-analysis stats, ' +
                        'and primary P-value cutoff.')
    parser.add_argument('gene_features', help='.tsv of gene features to use ' +
                        'for functional fine-mapping. First column = gene name. ' +
                        'All other columns must be numeric features. Optionally, ' + 
                        'first three columns can be BED-like, and will be dropped.')
    parser.add_argument('hpos_by_cohort', help='.tsv of sample sizes per HPO.')
    parser.add_argument('--cnv', help='Indicate CNV type. [default: NS]', default='NS')
    parser.add_argument('--secondary-p-cutoff', help='Maximum secondary P-value to ' + 
                        'consider as exome-wide significant. [default: 1]', 
                        default=1, type=float)
    parser.add_argument('--min-nominal', help='Minimum number of individual cohorts ' + 
                        'required to be nominally significant to consider ' +
                        'exome-wide significant. [default: 1]', default=1, type=int)
    parser.add_argument('--secondary-or-nominal', dest='secondary_or_nom', 
                        help='Allow genes to meet either --secondary-p-cutoff ' +
                        'or --min-nominal, but do not require both for exome-wide ' +
                        'significance. [default: require both]', default=False, 
                        action='store_true')
    parser.add_argument('--fdr-q-cutoff', help='Maximum FDR Q-value to ' + 
                        'consider as FDR significant. [default: 0.05]', 
                        default=0.05, type=float)
    parser.add_argument('--secondary-for-fdr', action='store_true', 
                        default=False, help='Apply sample secondary and/or min. ' +
                        'nominal criteria when evaluating FDR significance. ' +
                        '[default: apply no secondary criteria to FDR genes]')
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
    parser.add_argument('--nonsig-distance', help='Distance to pad each significant ' +
                        'gene when searching for non-sig genes to include in ' +
                        'fine-mapping. [defaults to value of --distance]', type=int)
    parser.add_argument('--covariance-tsv', help='.tsv of gene-gene CNV covariance ' +
                        'generated by count_cnvs_per_gene.py. If specified, will ' +
                        'be used in addition to --distance when clustering gene blocks.')
    parser.add_argument('--min-covariance', type=float, default=0.2, 
                        help='Minimum gene-gene covariance required when ' +
                        'clustering gene blocks.')
    parser.add_argument('--no-non-significant', action='store_true', default=False,
                        help='Do not include non-significant genes when fine-mapping. ' +
                        '[default: include nearby non-sig genes]')
    parser.add_argument('-x', '--training-exclusion', help='List of genes to exclude ' +
                        'in the training step of the functional fine-mapping model. ' +
                        '[default: include all genes]')
    parser.add_argument('--use-max-pip-per-gene', action='store_true', default=False,
                        help='Use maximum PIP when training model for genes with ' + 
                        'multiple associations [default: use mean PIP]')
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
    parser.add_argument('-o', '--outfile', help='Output tsv of final fine-mapping ' +
                        'results for significant genes and phenotypes in any ' +
                        'credible set.')
    parser.add_argument('--all-genes-outfile', help='Output tsv of final ' +
                        'fine-mapping results for all genes and phenotypes.')
    parser.add_argument('--naive-outfile', help='Output tsv of naive results ' +
                        'before fine-mapping for all genes and phenotypes.')
    parser.add_argument('--genetic-outfile', help='Output tsv of genetics-only ' +
                        'fine-mapping results for all genes and phenotypes.')
    parser.add_argument('--coeffs-out', help='Output tsv of logit coefficients ' +
                        'from fine-mapping model.')
    parser.add_argument('--sig-assoc-bed', help='Output BED of significant ' +
                        'gene-phenotype pairs and their corresponding association ' +
                        'statistics.')
    parser.add_argument('--sig-genes-bed', help='Output BED of all prioritized genes ' +
                        'significant  in for at least one HPO, as well as their ' +
                        'overall association statistics.')
    parser.add_argument('--hpo-credsets-bed', help='Output BED of credible sets ' +
                        'comprised of significant genes per HPO after fine-mapping, ' +
                        'and their corresponding association statistics.')
    parser.add_argument('--joint-credsets-bed', help='Output BED of credible sets ' +
                        'comprised of significant genes after joint fine-mapping, ' +
                        'and their corresponding association statistics.')
    parser.add_argument('--prefix', help='Prefix for naming loci & credible sets.')
    args = parser.parse_args()

    # Sets value of nonsig distance if not specified
    if args.nonsig_distance is None:
        args.nonsig_distance = args.distance

    # Set block prefix
    if args.prefix is None:
        block_prefix = '{}_gene_block'.format(args.cnv)
    else:
        block_prefix = '{}_{}_gene_block'.format(args.cnv, args.prefix)

    # Read dict of N_case per HPO
    ncase_df = pd.read_csv(args.hpos_by_cohort, sep='\t').loc[:, '#HPO Total'.split()]
    ncase_df.index = ncase_df.iloc[:, 0]
    ncase_dict = ncase_df.drop(columns='#HPO').transpose().to_dict(orient='records')[0]

    # Load CNV covariance, if optioned
    if args.covariance_tsv is not None:
        cov_df = pd.read_csv(args.covariance_tsv, sep='\t').\
                    rename(columns={"#geneA" : 'geneA'})
    else:
        cov_df = None

    # Load gold-standard causal gene lists to ensure their summary stats are loaded
    if args.gs_list is not None:
        gs_lists = {}
        all_gs_genes = set()
        with open(args.gs_list) as gsf:
            for gslist in gsf.read().splitlines():
                gs_genes = open(gslist).read().splitlines()
                gs_lists[gslist] = gs_genes
                all_gs_genes.update(gs_genes)

    # Process data per hpo
    hpo_data = load_all_hpos(args.statslist, args.secondary_p_cutoff, 
                             args.min_nominal, args.secondary_or_nom, 
                             args.fdr_q_cutoff, args.secondary_for_fdr,
                             args.distance, args.nonsig_distance,
                             cov_df, args.min_covariance,
                             block_prefix, args.finemap_secondary, 
                             not args.no_non_significant, all_gs_genes)

    # Estimate null variance based on:
    #   1. all significant genes
    #   2. most significant gene from each block
    #   3. known causal genes (optional; can be multiple lists)
    Wsq = estimate_null_variance_basic(hpo_data)
    if args.gs_list is not None:
        Wsq += estimate_null_variance_gs(gs_lists, args.statslist)
    Wsq = sorted(Wsq)
    print('Null variance estimates: ' + ', '.join([str(round(x, 3)) for x in Wsq]))

    # Update original finemapping results with BMA across re-estimated null variances
    hpo_data = update_finemap(hpo_data, Wsq)

    # Write naive and/or genetics-only fine-mapping results (for ROC comparisons)
    if args.naive_outfile is not None:
        # Note: cs_val set to 1.0 fixed value here such that all genes will be assigned to their original credible set
        # rather than randomly throwing out ~5% of genes
        make_sig_genes_df(hpo_data, naive=True, sig_only=True, cs_val=1.0).\
            sort_values(by='PIP gene HPO'.split(), ascending=(False, True, True, )).\
            rename(columns={'HPO' : '#HPO'}).\
            to_csv(args.naive_outfile, sep='\t', index=False, na_rep='NA')

    if args.genetic_outfile is not None:
        make_sig_genes_df(hpo_data, sig_only=True, cs_val=args.cs_val).\
            sort_values(by='PIP gene HPO'.split(), ascending=(False, True, True, )).\
            rename(columns={'HPO' : '#HPO'}).\
            to_csv(args.genetic_outfile, sep='\t', index=False, na_rep='NA')

    # Perform functional fine-mapping
    hpo_df, logit_coeffs = functional_finemap(hpo_data, args.gene_features, 
                                              args.l1_l2_mix, args.logit_alpha, 
                                              1.0, Wsq, args.training_exclusion,
                                              args.use_max_pip_per_gene)

    # If optioned, write average of logit coefficients to --coeffs-out
    if args.coeffs_out is not None:
        logit_coeffs.rename(columns={'feature' : '#feature'}).\
                     to_csv(args.coeffs_out, sep='\t', index=False, na_rep='NA')

    # Prepare data for writing out HPO-specific association information prior to joint finemapping
    prejoint_sig_df = make_sig_genes_df(hpo_data, sig_only=True, cs_val=args.cs_val)
    prejoint_allgenes_df = make_sig_genes_df(hpo_data, sig_only=False, cs_val=args.cs_val)
    prejoint_hpo_data = deepcopy(hpo_data)

    # Write out initial credible sets per HPO prior to joint finemapping, if optioned
    if args.hpo_credsets_bed is not None:
        hpo_credsets_bed = open(args.hpo_credsets_bed, 'w')
        output_credsets_bed(hpo_data, prejoint_allgenes_df, hpo_credsets_bed,
                            args.confident_pip, args.very_confident_pip, args.cnv)

    # Cluster credible sets across HPOs 
    cred_clusters = cluster_credsets(hpo_data, sig_only=True)

    # Joint refinement of overlapping credible sets across HPOs
    hpo_data, cred_clusters = joint_finemap(hpo_data, cred_clusters, ncase_df)

    # Write out final finemapping results for significant genes with both prejoint and joint pips
    final_sig_df = make_sig_genes_df(hpo_data, naive=False, sig_only=True, 
                                     cs_val=args.cs_val)
    pd.merge(final_sig_df, prejoint_sig_df.drop(columns='credible_set'), 
             on='HPO gene'.split(), suffixes='_final _HPO_specific'.split()).\
       sort_values(by='PIP_final gene PIP_HPO_specific'.split(), ascending=False).\
       rename(columns={'HPO' : '#HPO'}).\
       to_csv(args.outfile, sep='\t', index=False, na_rep='NA')

    # Write out final finemapping results for all genes with both prejoint and joint pips
    final_allgenes_df = make_sig_genes_df(hpo_data, naive=False, sig_only=False, 
                                          cs_val=args.cs_val)
    if args.all_genes_outfile is not None:
        pd.merge(final_allgenes_df, prejoint_sig_df.drop(columns='credible_set'), 
                 on='HPO gene'.split(), suffixes='_final _HPO_specific'.split()).\
           sort_values(by='PIP_final gene PIP_HPO_specific'.split(), ascending=False).\
           rename(columns={'HPO' : '#HPO'}).\
           to_csv(args.all_genes_outfile, sep='\t', index=False, na_rep='NA')

    # Write out significant gene-HPO associations, if optioned
    if args.sig_assoc_bed is not None:
        sig_assoc_bed = open(args.sig_assoc_bed, 'w')
        output_assoc_bed(hpo_data, final_sig_df, prejoint_sig_df, sig_assoc_bed,
                         args.confident_pip, args.very_confident_pip, args.cnv)

    # Format & write final table of credible sets (including sig. genes only)
    if args.joint_credsets_bed is not None:
        # Read pbt.BedTool of gene coordinates for convenience
        with open(args.statslist) as fin:
            genes_bt = pbt.BedTool(fin.readline().rstrip().split('\t')[1]).cut(range(4)).saveas()

        # Format & output final joint credsets
        joint_credsets_bed = open(args.joint_credsets_bed, 'w')
        joint_gene_cs_map = \
            output_joint_credsets_bed(hpo_data, final_allgenes_df, cred_clusters, genes_bt, 
                                      joint_credsets_bed, ncase_dict, args.confident_pip, 
                                      args.very_confident_pip, args.cnv, block_prefix)

    # # Format & write final table of significant genes
    if args.sig_genes_bed is not None:
        sig_genes_bed = open(args.sig_genes_bed, 'w')
        output_genes_bed(hpo_data, final_sig_df, prejoint_sig_df, joint_gene_cs_map,
                         sig_genes_bed, ncase_dict, args.confident_pip, 
                         args.very_confident_pip, args.cnv)


if __name__ == '__main__':
    main()
