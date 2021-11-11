#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Compile BED of final segments, other known genomic disorders, and other NAHR segments
Outputs summary table for rCNV2 paper analyses
"""


cnvtypes = 'DEL DUP'.split()


import pandas as pd
import pybedtools as pbt
from itertools import combinations
from numpy import nansum, log10, NaN, nanmax
from scipy.stats import hmean
import re
import argparse
from os.path import splitext
from sys import stdout
import csv
import subprocess


def read_bed_as_df(bed, output_hpos=False, output_final_assocs=False):
    """
    Reads a BED file as a pd.DataFrame and reformats data as needed
    """

    df = pd.read_csv(bed, sep='\t')
    c1h = list(df.columns)[0]
    df.rename(columns={c1h : c1h.replace('#', '')}, inplace=True)

    # Extract hpos as dictionary, if optioned
    if output_hpos:
        hpo_df = df['hpos']
        hpo_df.index = df['region_id']
        hpo_dict = {r : h.split(';') for r, h in hpo_df.to_dict().items()}

    # Extract constituent associations, if optioned
    if output_final_assocs:
        assoc_df = df['constituent_assocs']
        assoc_df.index = df['region_id']
        assoc_dict = {r : a.split(';') for r, a in assoc_df.to_dict().items()}

    # Add dummy column for significance level if not already present
    if 'best_sig_level' not in df.columns:
        df['best_sig_level'] = 'not_significant'

    # Subset columns of interest
    keep_cols = 'best_sig_level cred_interval_coords cred_intervals_size n_genes genes'.split()
    df = pd.concat([df.iloc[:, :5], df.loc[:, df.columns.isin(keep_cols)]], axis=1)

    # Rename columns
    old_cnames = list(df.columns)
    df.rename(columns={old_cnames[0] : 'chr',
                       old_cnames[1] : 'start',
                       old_cnames[2] : 'end',
                       old_cnames[3] : 'region_id',
                       old_cnames[4] : 'cnv',
                       'best_sig_level' : 'sig_level',
                       'cred_interval_coords' : 'coords',
                       'cred_intervals_size' : 'size'}, 
              inplace=True)

    # Add interval coordinates and size (if not optioned)
    if 'coords' not in df.columns:
        df['coords'] = df.apply(lambda x: '{}:{}-{}'.format(*x[['chr', 'start', 'end']]), axis=1)
    if 'size' not in df.columns:
        df['size'] = df.end - df.start

    cols_ordered = 'chr start end region_id cnv sig_level coords size n_genes genes'.split()
    if output_hpos and output_final_assocs:
        return df.loc[:, cols_ordered], hpo_dict, assoc_dict
    elif output_hpos:
        return df.loc[:, cols_ordered], hpo_dict
    elif output_final_assocs:
        return df.loc[:, cols_ordered], assoc_dict
    else:
        return df.loc[:, cols_ordered]


def find_overlap(bta, btb, r=0.000001):
    """
    Overlaps two pbt.BedTool objects based on reciprocal overlap (r)
    Matches on CNV type (assumes the fifth column = CNV type)
    Returns: two lists of interval IDs of the overlapping set
    """

    ibt = bta.cut(range(5)).intersect(btb.cut(range(5)), f=r, r=True, wa=True, wb=True).\
              filter(lambda x: x[4] == x[9]).saveas()

    ids_a = [x[3] for x in ibt]
    ids_b = [x[8] for x in ibt]

    return ids_a, ids_b


def classify_nahr(all_df, nahr_bt, cnv_bed, r=0.5, frac_nahr_cnvs=0.2):
    """
    Classify segments in all_df as NAHR or not based on their underlying case CNVs
    Returns a list of IDs from all_bt meeting NAHR criteria
    """

    # Create BedTool of segment credible intervals
    ci_bt_str = ''
    for seg_id, cnv, cstr in all_df['region_id cnv coords'.split()].values.tolist():
        for substr in cstr.split(';'):
            ci_bt_str += '\t'.join(re.split(':|-', substr) + [seg_id, cnv]) + '\n'
    ci_bt = pbt.BedTool(ci_bt_str, from_string=True)

    # Pre-filter CNVs to only retain cases
    cnv_bt = pbt.BedTool(cnv_bed).filter(lambda x: 'HEALTHY_CONTROL' not in x[5]).saveas()

    # Extract list of CNV IDs meeting NAHR overlap criteria
    nahr_cnv_ids=[cnv[3] for cnv in cnv_bt.intersect(nahr_bt, wa=True, f=r, r=True)]

    # Overlap CNVs with each region and make dict of CNVs per region while matching on CNV type
    seg_cnv_dict = {}
    for hit in cnv_bt.intersect(ci_bt, wa=True, wb=True, f=r, r=True):

        cnv_id, cnv_type, seg_id, seg_type = map(hit.fields.__getitem__, [3, 4, 9, 10])
        
        # Only keep CNV type matches
        if cnv_type != seg_type:
            continue

        if seg_id in seg_cnv_dict.keys():
            seg_cnv_dict[seg_id].add(cnv_id)
        else:
            seg_cnv_dict[seg_id] = set([cnv_id])

    # Compute NAHR fraction per segment
    seg_nahr_frac = {seg_id : len([x for x in cnvs if x in nahr_cnv_ids]) / len(cnvs) 
                     for seg_id, cnvs in seg_cnv_dict.items()}

    # Return list of segments meeting minimum NAHR CNV fraction
    return [seg_id for seg_id, frac in seg_nahr_frac.items() if frac >= frac_nahr_cnvs]


def make_terminal_bt(genomefile, dist=1000000):
    """
    Create a pbt.BedTool of terminal search windows from a BEDTools-styled .genome file
    """

    term_str = ''
    line_fmt = '{}\t{}\t{}\t{}\t{}\t\n'

    with open(genomefile) as gfile:
        for chrom, clen in csv.reader(gfile, delimiter='\t'):
            for cnv in 'DEL DUP'.split():
                term_str += line_fmt.format(chrom, 0, dist, chrom + '_p', cnv)
                term_str += line_fmt.format(chrom, max([0, int(clen)-dist]), 
                                            clen, chrom + '_q', cnv)

    return pbt.BedTool(term_str, from_string=True)


def overlap_common_cnvs(all_bt, common_cnvs, cnv, min_cov):
    """
    Overlap regions with a BED of common CNVs, and return list of IDs ≥ min_cov
    """

    if common_cnvs is not None:
        bt_cov = all_bt.coverage(common_cnvs)
        ids = [x[3] for x in bt_cov if x[4] == cnv and float(x[-1]) >= min_cov]

    else:
        ids = []

    return ids


def score_pleiotropy(df, loci_hpos, jaccard_matrix_in, min_jaccard_sum=2.0):
    """
    Score genome-wide significant loci for pleiotropic effects based on the
    sum of (1 - Jaccard index) for all pairs of associated HPOs

    Note that 1 - Jaccard is necessary because we are interested in the fraction 
    of samples that _don't_ overlap (whereas Jaccard index measures fraction overlap)

    Returns: a list of dummy indicator variables
    """

    jac = pd.read_csv(jaccard_matrix_in, sep='\t')
    if 'HPO' in jac.columns:
        jac.index = jac['HPO']
        jac.drop('HPO', axis=1, inplace=True)

    def _calc_pleio(rid):
        pleio_sum = 0
        if rid in loci_hpos.keys():
            for hpoA, hpoB in combinations(loci_hpos[rid], 2):
                pleio_sum += (1 - jac.loc[hpoA, hpoB])
            if pleio_sum >= min_jaccard_sum:
                return 1
            else:
                return 0
        else:
            return 0

    return [_calc_pleio(rid) for rid in df['region_id'].tolist()]


def _subset_genes(gstr, genes):
    """
    Helper function to subset a semicolon-delimited string of genes on a second set of genes
    """

    hits = [g for g in str(gstr).split(';') if g in genes]
    
    if len(hits) == 0 or hits[0] == '':
        return 'NA'
    
    else:
        return ';'.join(sorted(hits))



def annotate_genelist(df, genes, glname):
    """
    Annotates overlap of regions with a specified genelist
    Returns: input df, with two extra columns:
               1. n_glname_genes = number of genes in region also present in genelist
               2. glname_genes = gene symbols of genes in region also present in genelist
    """

    df[glname + '_genes'] = df.genes.apply(_subset_genes, genes=genes)

    df['n_' + glname + '_genes'] = \
        df[glname + '_genes'].apply(lambda gstr: len([g for g in str(gstr).split(';') if g != 'NA']))

    return df


def load_hpo_genelists(tsv_in):
    """
    Load HPO-specific gene lists
    Returns: dict of HPO : [genes] pairings
    """

    hpo_genes = {}

    with open(tsv_in) as fin:
        reader = csv.reader(fin, delimiter='\t')
        for hpo, path in reader:
            hpo_genes[hpo] = [g.rstrip() for g in open(path).readlines()]

    return hpo_genes


def annotate_hpo_genes(df, loci_hpos, hpo_genes):
    """
    Annotates overlap of regions with HPO-specific genelists
    Returns: input df, with two extra columns:
               1. n_HPOmatched_genes = number of genes in region also present in HPO-matched genelist(s)
               2. HPOmatched_genes = gene symbols of genes in region also present in HPO-matched genelist(s)
    """

    # Make list of eligible OMIM genes per region
    region_targets = {}
    for rid, hpos in loci_hpos.items():
        region_targets[rid] = \
            sorted(list(set([g for s in [hpo_genes[h] for h in hpos] for g in s])))

    # Find list of hits for each region
    region_hits = {}
    hit_counts = {}
    for rid, targets in region_targets.items():
        region_hits[rid] = \
            _subset_genes(df.genes[df['region_id'] == rid].tolist()[0], targets)
        hit_counts[rid] = \
            len([g for g in str(region_hits[rid]).split(';') if g != 'NA'])

    df['HPOmatched_genes'] = df.region_id.map(region_hits)
    df['n_HPOmatched_genes'] = df.region_id.map(hit_counts)

    return df


def annotate_constraint(df, constraint_path):
    """
    Annotates each segment with mutational constraint info from gnomAD
    """

    if constraint_path.endswith('.bgz'):
        constr_df = pd.read_csv(constraint_path, sep='\t', compression='gzip')
    else:
        constr_df = pd.read_csv(constraint_path, sep='\t')

    # Map minimum constraint metric among all genes in region
    def _min_metric(genes_str, constr_df, col):
        genes = str(genes_str[0]).split(';')
        return constr_df.loc[constr_df.gene.isin(genes), col].min()

    df['min_LOEUF'] \
        = pd.DataFrame(df.genes).apply(_min_metric, axis=1, raw=True, 
                                       constr_df=constr_df, col='oe_lof_upper')

    df['min_MisOEUF'] \
        = pd.DataFrame(df.genes).apply(_min_metric, axis=1, raw=True, 
                                       constr_df=constr_df, col='oe_mis_upper')

    # Compute total LOF and missense O/E from all genes in region
    def _sum_metric(genes_str, constr_df, csq):
        genes = str(genes_str[0]).split(';')
        obs = constr_df.loc[constr_df.gene.isin(genes), 'obs_' + csq].sum()
        exp = constr_df.loc[constr_df.gene.isin(genes), 'exp_' + csq].sum()
        if exp > 0:
            return obs / exp
        else:
            return NaN

    df['total_LoF_OE'] \
        = pd.DataFrame(df.genes).apply(_sum_metric, axis=1, raw=True, 
                                       constr_df=constr_df, csq='lof')

    df['total_mis_OE'] \
        = pd.DataFrame(df.genes).apply(_sum_metric, axis=1, raw=True, 
                                       constr_df=constr_df, csq='mis')

    return df


def count_dnms(df, dnm_path, col_prefix, mu_df=None, xgenes=None):
    """
    Annotate segments based on sum of observed de novo coding mutations
    """

    # Load DNMs
    dnms = pd.read_csv(dnm_path, sep='\t').\
              rename(columns={'#gene' : 'gene'})

    # Drop excluded genes, if optioned
    if xgenes is not None:
        dnms = dnms.loc[~dnms.gene.isin(xgenes), :]

    def _sum_dnms(genes_str, dnms, key, include=None):
        if include is not None:
            genes = [x for x in str(genes_str[0]).split(';') if x in include]
        else:
            genes = str(genes_str[0]).split(';')
        if len(genes) == 0 or 'NaN' in genes:
            return 0
        else:
            return dnms.loc[dnms.gene.isin(genes), key].sum()
            
    for csq in 'lof mis syn'.split():
        df['_'.join([col_prefix, 'dnm', csq])] \
            = pd.DataFrame(df.genes).apply(_sum_dnms, axis=1, raw=True, dnms=dnms, key=csq)
        if mu_df is not None:
            mu_df_x = mu_df.copy(deep=True)
            if xgenes is not None:
                mu_df_x = mu_df_x.loc[~mu_df_x.gene.isin(xgenes), :]
            n_dnms = df['_'.join([col_prefix, 'dnm', csq])].sum()
            mu_df_x['mu_' + csq] = mu_df_x['mu_' + csq] * n_dnms
            genes_with_mus = mu_df_x.dropna()['gene'].tolist()
            obs_with_mus = pd.DataFrame(df.genes).\
                                    apply(_sum_dnms, axis=1, raw=True, dnms=dnms, 
                                          key=csq, include=genes_with_mus)
            expected = pd.DataFrame(df.genes).\
                          apply(_sum_dnms, axis=1, raw=True, dnms=mu_df_x, 
                                key='mu_' + csq, include=genes_with_mus)
            oe_ratio = obs_with_mus / expected
            df['_'.join([col_prefix, 'dnm', csq, 'obs_wMu'])] \
                = (obs_with_mus)
            df['_'.join([col_prefix, 'dnm', csq, 'exp_wMu'])] \
                = (expected).round(decimals=6)

    return df


def count_bcas(all_df, bca_tsv):
    """
    Count number of BCA carriers with a breakpoint in each region
    """

    # Convert CNV loci to pbt.BedTool and intersect with bca_tsv
    bca_hits = {i : set() for i in all_df.region_id.tolist()}
    all_bt = pbt.BedTool.from_dataframe(all_df.iloc[:, :4])
    for hit in all_bt.intersect(bca_tsv, wa=True, wb=True):
        bca_hits[hit.name].add(hit[-1])

    # Count number of unique BCA carriers per region
    bca_k = {i : len(bca_hits[i]) for i in bca_hits.keys()}

    # Map counts onto all_df
    all_df['Redin_BCAs'] = all_df['region_id'].map(bca_k)

    return all_df


def annotate_expression(all_df, gtex_in, min_expression=1):
    """
    Annotate segments with various expression-based metrics
    Adds the following columns to all_df:
        1. gene_expression_harmonic_mean: harmonic mean of expression levels across tissues for all genes per segment
        2. n_ubiquitously_expressed_genes: count of genes expressed >= min_expression for all tissues
    """

    # Load GTEx stats
    gtex = pd.read_csv(gtex_in, sep='\t').rename(columns = {'#gene' : 'gene'})
    gtex.index = gtex['gene']
    gtex.drop(columns='gene', inplace=True)

    # Compute harmonic mean of expression levels for all genes per segment
    def _gtex_geneset_hmean(genes_str, gtex_means):
        elig_genes = gtex_means.index.tolist()
        genes = [x for x in str(genes_str[0]).split(';') if x in elig_genes]
        if len(genes) > 0:
            return hmean(gtex_means[gtex_means.index.isin(genes)])
        else:
            return 0

    gtex_means = gtex.mean(axis=1, skipna=True)
    all_df['gene_expression_harmonic_mean'] \
        = pd.DataFrame(all_df.genes).apply(_gtex_geneset_hmean, axis=1, raw=True, 
                                           gtex_means=gtex_means)

    # Compute number of tissues with expression > min_expression per gene
    # Note that GTEx matrix has been log10(x+1) transformed
    def _count_ubi(genes_str, gtex_ubi):
        elig_genes = gtex_means.index.tolist()
        genes = [x for x in str(genes_str[0]).split(';') if x in elig_genes]
        if len(genes) > 0:
            return gtex_ubi[gtex_ubi.index.isin(genes)].sum()
        else:
            return 0

    trans_min = log10(min_expression + 1)
    gtex_ubi = gtex.apply(lambda vals: all(vals.dropna() >= trans_min), axis=1)
    all_df['n_ubiquitously_expressed_genes'] \
        = pd.DataFrame(all_df.genes).apply(_count_ubi, axis=1, raw=True, gtex_ubi=gtex_ubi)

    return all_df


def annotate_meta_sumstats(all_df, sumstats_in, loci_assocs, assocs_in=None, 
                           add_nom_neuro=False, neuro_hpos=None, 
                           add_nom_dev=False, dev_hpos=None):
    """
    Annotate segments with sumstats for best phenotype association while matching on CNV type
    Adds the following columns to all_df:
        1. meta_best_p: top P-value for any phenotype
            Note that if assocs_in is provided, P-values will be extracted from there
            where possible to ensure consistency with individual association reporting
        2. nom_sig: dummy indicator if the region is nominally significant in at least one phenotype
        3. nom_neuro [optional]: dummy indicator if the region is nominally significant in at least one neurological phenotype
        4. nom_dev [optional]: dummy indicator if the region is nominally significant in at least one developmental phenotype
        5. meta_best_lnor: top lnOR for any phenotype for any window in the credible interval(s) 
    """

    # Load sumstats
    ss = pd.read_csv(sumstats_in, sep='\t')
    numeric_cols = 'lnor lnor_lower lnor_upper pvalue pvalue_secondary'.split()
    ss[numeric_cols] = ss[numeric_cols].apply(pd.to_numeric)

    # Load final associations, if optioned
    if assocs_in is not None:
        assoc_df = pd.read_csv(assocs_in, sep='\t')

    # Get best P-value per segment
    best_ps = []
    nomsig = []
    nomneuro = []
    nomdev = []
    best_lnORs = []
    for rid, cnv in all_df.loc[:, 'region_id cnv'.split()].itertuples(index=False, name=None):
        if rid in loci_assocs.keys() and assocs_in is not None:
            pvals = -log10(assoc_df.best_pvalue[assoc_df.credible_set_id.isin(loci_assocs.get(rid))])
        else:
            pvals = ss.loc[(ss.region_id == rid) & (ss.cnv == cnv), 'pvalue']
        if len(pvals) > 0:
            top_p = float(nanmax(pvals.astype(float)))
            best_ps.append(top_p)
            if top_p >= -log10(0.05):
                nomsig.append(1)
            else:
                nomsig.append(0)
            if add_nom_neuro:
                nom_hpos = ss.loc[(ss.region_id == rid) & (ss.cnv == cnv) & (ss.pvalue >= -log10(0.05)), 'hpo']
                nom_neuro_hpos = [h for h in nom_hpos.tolist() if h in neuro_hpos]
                if len(nom_neuro_hpos) > 0:
                    nomneuro.append(1)
                else:
                    nomneuro.append(0)
            if add_nom_dev:
                nom_hpos = ss.loc[(ss.region_id == rid) & (ss.cnv == cnv) & (ss.pvalue >= -log10(0.05)), 'hpo']
                nom_dev_hpos = [h for h in nom_hpos.tolist() if h in dev_hpos]
                if len(nom_dev_hpos) > 0:
                    nomdev.append(1)
                else:
                    nomdev.append(0)
        else:
            best_ps.append(NaN)
            nomsig.append(0)
            nomneuro.append(0)
            nomdev.append(0)
        lnors = ss.loc[(ss.region_id == rid) & (ss.cnv == cnv), 'lnor']
        if len(lnors) > 0:
            top_lnOR = float(nanmax(lnors.astype(float)))
            best_lnORs.append(top_lnOR)
        else:
            best_lnORs.append(NaN)

    all_df['meta_best_p'] = best_ps
    all_df['nom_sig'] = nomsig
    if add_nom_neuro:
        all_df['nom_neuro'] = nomneuro
    if add_nom_dev:
        all_df['nom_dev'] = nomdev
    all_df['meta_best_lnor'] = best_lnORs

    return all_df


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--final-loci', required=True, help='BED of final large ' +
                        'segment loci from association testing. Required.')
    parser.add_argument('--hc-gds', required=True, help='BED of high-confidence GDs. ' +
                        'Required.')
    parser.add_argument('--mc-gds', required=True, help='BED of medium-confidence GDs. ' +
                        'Required.')
    parser.add_argument('--lc-gds', required=True, help='BED of low-confidence GDs. ' +
                        'Required.')
    parser.add_argument('--nahr-cnvs', required=True, help='BED of predicted NAHR-' +
                        'mediated CNVs. Required.')
    parser.add_argument('--final-associations', help='BED of final associations ' +
                        'corresponding to --final-loci. If provided, will be used ' +
                        'to annotate best P-value to ensure consistency with ' +
                        'original associations.')
    parser.add_argument('--loose-nahr-mask', help='BED of more lenient possible NAHR ' +
                        'events. Will only be used for labeling significant loci and/' +
                        'or GDs as NAHR; will not be included as rows in final output. ' +
                        'If not provided, will default to using --nahr-cnvs for labeling.')
    parser.add_argument('--cnv-bed', help='BED file of all CNVs. Required.')
    parser.add_argument('-g', '--genome-file', help='BEDTools-style .genome file.' +
                        'If provided, will be used to annotate segments as terminal.')
    parser.add_argument('--terminal-window', type=int, default=1000000, 
                        help='Distance to search from telomeres when labeling ' +
                        'segments as terminal.')
    parser.add_argument('--common-dels', help='BED of common deletions.')
    parser.add_argument('--common-dups', help='BED of common duplications.')
    parser.add_argument('--common-cnv-cov', type=float, default=0.3, help='Coverage ' +
                        'of common CNVs required to label as "benign". [default: 0.3]')
    parser.add_argument('--hpo-jaccard-matrix', help='Tsv of Jaccard indexes ' +
                        'between HPOs. Will be used to label pleiotropic loci.')
    parser.add_argument('--genelists', help='Tsv of genelists to annotate. Two ' +
                        'columns expected: genelist name, and path to genelist.')
    parser.add_argument('--hpo-genelists', help='Tsv of HPO-specific genelists to ' +
                        'annotate. Will match on segment HPO. Two columns ' +
                        'expected: HPO, and path to genelist.')
    parser.add_argument('--gnomad-constraint-tsv', help='Tsv of gene-level ' +
                        'constraint metadata. If provided, will annotate each ' +
                        'segment with constraint metrics from gnomAD.')
    parser.add_argument('--dnm-tsvs', help='Tsv of de novo mutation counts ' +
                        ' to annotate. Three columns expected: study prefix, ' +
                        'path to tsv with dnm counts, and path to list of exome-' +
                        'wide significant genes from that study.')
    parser.add_argument('--snv-mus', help='Tsv of snv mutation rates per gene. ' +
                        'Four columns expected: gene, and relative mutation rates ' +
                        'for lof, missense, and synonymous mutations.')
    parser.add_argument('--bca-tsv', help='Tsv of BCA breakpoint coordinates ' +
                        'from Redin et al. Four columns expected: chrom, start, ' +
                        'end, and sample ID.')
    parser.add_argument('--gtex-matrix', help='Tsv gene X tissue expression levels ' +
                        'from GTEx. Will be used for various expression-based ' +
                        'gene annotations.')
    parser.add_argument('--meta-sumstats', help='Tsv of meta-analysis summary statistics ' +
                        'per phenotype per region. Computed with calc_all_seg_stats.py')
    parser.add_argument('--neuro-hpos', help='List of neurological HPOs.')
    parser.add_argument('--dev-hpos', help='List of developmental HPOs.')
    parser.add_argument('--min-expression', default=1, help='Minimum expression ' +
                        'level (in unscaled TPM) to consider a gene as "expressed". ' +
                        '[default: 1]')
    parser.add_argument('--gd-recip', type=float, default=0.2, help='Reciprocal ' +
                        'overlap required for GD match. [default: 0.2]')
    parser.add_argument('--nahr-recip', type=float, default=0.5, help='Reciprocal ' +
                        'overlap required for NAHR CNV match. [default: 0.5]')
    parser.add_argument('--min-jaccard-sum', type=float, default=2.0, help='Minimum ' +
                        'sum of Jaccard indexes to consider a locus as pleiotropic. ' +
                        '[default: 2.0]')
    parser.add_argument('--keep-all-gds', default=False, action='store_true',
                        help='Keep all HC/MC/LC GDs. Default behavior is to remove ' +
                        'lower-confidence GDs that overlap with any higher-confidence ' +
                        'segment. Not a recommended option for most analyses.')
    parser.add_argument('-o', '--outfile', help='Path to output BED file. [default: ' +
                        'stdout]', default='stdout')
    parser.add_argument('-z', '--bgzip', action='store_true', help='Compress ' + 
                        'output BED files with bgzip. [Default: do not compress]')
    args = parser.parse_args()

    # Open connection to outfile
    if args.outfile in '- stdout'.split():
        outfile = stdout
    else:
        gz_suf = 'gz gzip bgzip bgz bz'.split()
        if splitext(args.outfile)[-1].replace('.' ,'') in gz_suf:
            outfile_path = splitext(args.outfile)[0]
        else:
            outfile_path = args.outfile
        outfile = open(outfile_path, 'w')

    # Load final significant loci
    loci_bt = pbt.BedTool(args.final_loci).cut(range(5))
    loci_df, loci_hpos, loci_assocs = \
        read_bed_as_df(args.final_loci, output_hpos=True, output_final_assocs=True)
    loci_ids = loci_df.iloc[:, 3].tolist()
    all_df = loci_df.copy(deep=True)
    all_bt = loci_bt.saveas()

    # Load high-confidence GDs and merge into master df
    hc_gd_bt = pbt.BedTool(args.hc_gds).cut(range(5))
    hc_gd_df = read_bed_as_df(args.hc_gds)
    hc_gd_ids = hc_gd_df.iloc[:, 3].tolist()
    all_ids_in_hc_gd, hc_gd_ids_in_all = find_overlap(all_bt, hc_gd_bt, r=args.gd_recip)
    if args.keep_all_gds:
        hc_gd_ids_in_all = []
    all_df = pd.concat([all_df, hc_gd_df.loc[~hc_gd_df.region_id.isin(hc_gd_ids_in_all), :]], axis=0)
    all_bt = pbt.BedTool().from_dataframe(all_df)

    # Load medium-confidence GDs and merge into master df
    mc_gd_bt = pbt.BedTool(args.mc_gds).cut(range(5))
    mc_gd_df = read_bed_as_df(args.mc_gds)
    mc_gd_ids = mc_gd_df.iloc[:, 3].tolist()
    all_ids_in_mc_gd, mc_gd_ids_in_all = find_overlap(all_bt, mc_gd_bt, r=args.gd_recip)
    if args.keep_all_gds:
        mc_gd_ids_in_all = []
    all_df = pd.concat([all_df, mc_gd_df.loc[~mc_gd_df.region_id.isin(mc_gd_ids_in_all), :]], axis=0)
    all_bt = pbt.BedTool().from_dataframe(all_df)

    # Load low-confidence GDs and merge into master df
    lc_gd_bt = pbt.BedTool(args.lc_gds).cut(range(5))
    lc_gd_df = read_bed_as_df(args.lc_gds)
    lc_gd_ids = lc_gd_df.iloc[:, 3].tolist()
    all_ids_in_lc_gd, lc_gd_ids_in_all = find_overlap(all_bt, lc_gd_bt, r=args.gd_recip)
    if args.keep_all_gds:
        lc_gd_ids_in_all = []
    all_df = pd.concat([all_df, lc_gd_df.loc[~lc_gd_df.region_id.isin(lc_gd_ids_in_all), :]], axis=0)
    all_bt = pbt.BedTool().from_dataframe(all_df)

    # Load predicted NAHR-mediated CNVs
    nahr_bt = pbt.BedTool(args.nahr_cnvs).cut(range(5))
    nahr_df = read_bed_as_df(args.nahr_cnvs)
    nahr_ids = nahr_df.iloc[:, 3].tolist()
    all_ids_in_nahr, nahr_ids_in_all = find_overlap(all_bt, nahr_bt, r=args.nahr_recip)
    if args.cnv_bed is not None:
        # If optioned, reclassify NAHR mechanism for significant and/or GD segments based on underlying CNV evidence
        if args.loose_nahr_mask is not None:
            loose_nahr_bt = pbt.BedTool(args.loose_nahr_mask)
        else:
            loose_nahr_bt = nahr_bt
        all_ids_in_nahr = classify_nahr(all_df, loose_nahr_bt, args.cnv_bed, r=args.nahr_recip)
    all_df = pd.concat([all_df, nahr_df.loc[~nahr_df.region_id.isin(nahr_ids_in_all), :]], axis=0)
    all_bt = pbt.BedTool().from_dataframe(all_df)

    # Annotate terminal CNV segments
    if args.genome_file is not None:
        terminal_bt = make_terminal_bt(args.genome_file, args.terminal_window)
        terminal_ids, not_terminal_ids = find_overlap(all_bt, terminal_bt)

    # Find regions overlapping common CNVs, if optioned
    benign_del_ids = overlap_common_cnvs(all_bt, args.common_dels, "DEL", args.common_cnv_cov)
    benign_dup_ids = overlap_common_cnvs(all_bt, args.common_dups, "DUP", args.common_cnv_cov)

    # Annotate merged df
    all_df['any_sig'] = all_df.sig_level.isin(['genome_wide', 'FDR']).astype(int)
    all_df['gw_sig'] = all_df.sig_level.isin(['genome_wide']).astype(int)
    all_df['fdr_sig'] = all_df.sig_level.isin(['FDR']).astype(int)
    all_df.drop(columns='sig_level', inplace=True)
    all_df['hc_gd'] = pd.get_dummies(all_df.region_id.isin(set(hc_gd_ids + all_ids_in_hc_gd)), 
                                     drop_first=True)
    all_df['mc_gd'] = pd.get_dummies(all_df.region_id.isin(set(mc_gd_ids + all_ids_in_mc_gd)), 
                                     drop_first=True)
    all_df['lc_gd'] = pd.get_dummies(all_df.region_id.isin(set(lc_gd_ids + all_ids_in_lc_gd)), 
                                     drop_first=True)
    all_df['any_gd'] = pd.get_dummies(all_df.region_id.isin(set(hc_gd_ids + all_ids_in_hc_gd + \
                                                                mc_gd_ids + all_ids_in_mc_gd + \
                                                                lc_gd_ids + all_ids_in_lc_gd)), 
                                      drop_first=True)
    benign_ids = benign_del_ids + benign_dup_ids
    cand_path_ids = loci_ids + hc_gd_ids + mc_gd_ids + lc_gd_ids
    path_ids = list(set([i for i in cand_path_ids if i not in benign_ids]))
    if all_df.region_id.isin(path_ids).all():
        all_df['pathogenic'] = 1
    else:
        all_df['pathogenic'] = pd.get_dummies(all_df.region_id.isin(path_ids), drop_first=True)
    if args.common_dels is not None and args.common_dups is not None:
        all_df['benign'] = pd.get_dummies(all_df.region_id.isin(benign_ids), drop_first=True)
    all_df['nahr'] = pd.get_dummies(all_df.region_id.isin(set(nahr_ids + all_ids_in_nahr)), 
                                    drop_first=True)
    if args.genome_file is not None:
        all_df['terminal'] = pd.get_dummies(all_df.region_id.isin(terminal_ids), drop_first=True)

    # Label pleiotropic loci based on HPO overlap, if optioned
    if args.hpo_jaccard_matrix is not None:
        all_df['pleiotropic'] = score_pleiotropy(all_df, loci_hpos, 
                                                 args.hpo_jaccard_matrix,
                                                 args.min_jaccard_sum)

    # Annotate genelists
    if args.genelists is not None:
        with open(args.genelists) as genetsv:
            reader = csv.reader(genetsv, delimiter='\t')
            for glname, glpath in reader:
                genes = [g.rstrip() for g in open(glpath).readlines()]
                all_df = annotate_genelist(all_df, genes, glname)

    # Annotate with HPO-specific gene lists, if optioned
    if args.hpo_genelists is not None:
        hpo_genes = load_hpo_genelists(args.hpo_genelists)
        all_df = annotate_hpo_genes(all_df, loci_hpos, hpo_genes)

    # Annotate with gene-level constraint metadata, if optioned
    if args.gnomad_constraint_tsv is not None:
        all_df = annotate_constraint(all_df, args.gnomad_constraint_tsv)

    # Annotate with de novo mutations, if optioned
    if args.dnm_tsvs is not None:
        # Also, load snv mutation rates, if optioned
        if args.snv_mus is not None:
            mu_df = pd.read_csv(args.snv_mus, sep='\t').\
                       rename(columns={'#gene' : 'gene'}).\
                       dropna()
        else:
            mu_df = None
        with open(args.dnm_tsvs) as dnm_ins:
            for study, dnm_path, xgene_list in csv.reader(dnm_ins, delimiter='\t'):
                all_df = count_dnms(all_df, dnm_path, study, mu_df)
                xgenes = [g.rstrip() for g in open(xgene_list).readlines()]
                all_df = count_dnms(all_df, dnm_path, study + '_noSig', mu_df, xgenes)

    # Annotate with breakpoints from Redin et al., if optioned
    if args.bca_tsv is not None:
        all_df = count_bcas(all_df, args.bca_tsv)

    # Annotate with various expression-based metrics
    if args.gtex_matrix is not None:
        all_df = annotate_expression(all_df, args.gtex_matrix, args.min_expression)

    # Annotate with best P-value from meta-analysis summary statistics, if optioned
    if args.meta_sumstats is not None:
        if args.neuro_hpos is not None:
            add_nom_neuro = True
            neuro_hpos = [x.rstrip() for x in open(args.neuro_hpos).readlines()]
        else:
            add_nom_neuro = False
            neuro_hpos = None
        if args.dev_hpos is not None:
            add_nom_dev = True
            dev_hpos = [x.rstrip() for x in open(args.dev_hpos).readlines()]
        else:
            add_nom_dev = False
            dev_hpos = None
        all_df = annotate_meta_sumstats(all_df, args.meta_sumstats, loci_assocs,
                                        args.final_associations, add_nom_neuro, 
                                        neuro_hpos, add_nom_dev, dev_hpos)
<<<<<<< HEAD

    # Annotate GDs passing Bonferroni cutoff for total number of lit GDs tested
    gd_bonf_cutoff = -log10(0.05 / len(hc_gd_ids + mc_gd_ids + lc_gd_ids))
    all_df['bonf_sig_gd'] = pd.get_dummies((all_df.any_gd == 1) & \
                                           (all_df.meta_best_p >= gd_bonf_cutoff), 
                                           drop_first=True)
=======
>>>>>>> df5ba8e245e8774b1c1a11c8de746305e02d49db

    # Sort & write out merged BED
    all_df.sort_values(by='chr start end cnv region_id'.split(), inplace=True)
    all_df.rename(columns={list(all_df.columns)[0] : '#chr'}).\
           to_csv(outfile, sep='\t', na_rep='NA', index=False)
    outfile.close()
    if args.outfile not in '- stdout'.split():
        if args.bgzip:
            subprocess.run(['bgzip', '-f', outfile_path])


if __name__ == '__main__':
    main()
