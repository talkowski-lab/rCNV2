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
from numpy import nansum
import argparse
from os.path import splitext
from sys import stdout
import csv
from athena.utils import bgzip


def read_bed_as_df(bed, output_hpos=False):
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

    # Subset columns of interest
    keep_cols = 'cred_interval_coords cred_intervals_size n_genes genes'.split()
    df = pd.concat([df.iloc[:, :5], df.loc[:, df.columns.isin(keep_cols)]], axis=1)

    # Rename columns
    old_cnames = list(df.columns)
    df.rename(columns={old_cnames[0] : 'chr',
                      old_cnames[1] : 'start',
                      old_cnames[2] : 'end',
                      old_cnames[3] : 'region_id',
                      old_cnames[4] : 'cnv',
                      'cred_interval_coords' : 'coords',
                      'cred_intervals_size' : 'size'}, 
              inplace=True)

    # Add interval coordinates and size (if not optioned)
    if 'coords' not in df.columns:
        df['coords'] = df.apply(lambda x: '{}:{}-{}'.format(*x[['chr', 'start', 'end']]), axis=1)
    if 'size' not in df.columns:
        df['size'] = df.end - df.start

    cols_ordered = 'chr start end region_id cnv coords size n_genes genes'.split()
    if output_hpos:
        return df.loc[:, cols_ordered], hpo_dict
    else:
        return df.loc[:, cols_ordered]


def find_overlap(bta, btb, r=0.01):
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


def count_dnms(df, dnm_path, col_prefix):
    """
    Annotate segments based on sum of observed de novo coding mutations
    """

    # Load DNMs
    dnms = pd.read_csv(dnm_path, sep='\t').\
              rename(columns={'#gene' : 'gene'})

    def _sum_dnms(genes_str, dnms, key):
        genes = str(genes_str[0]).split(';')
        if len(genes) == 0 or 'NaN' in genes:
            return 0
        else:
            return dnms.loc[dnms.gene.isin(genes), key].sum()
            
    for csq in 'lof mis syn'.split():
        df['_'.join([col_prefix, 'dnm', csq])] \
            = pd.DataFrame(df.genes).apply(_sum_dnms, axis=1, raw=True, dnms=dnms, key=csq)

    return df


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
    parser.add_argument('-o', '--outfile', help='Path to output BED file. [default: ' +
                        'stdout]', default='stdout')
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
    parser.add_argument('--dnm-tsvs', help='Tsv of de novo mutation counts ' +
                        ' to annotate. Two columns expected: study prefix, and ' +
                        'path to tsv with dnm counts.')
    parser.add_argument('--gd-recip', type=float, default=0.2, help='Reciprocal ' +
                        'overlap required for GD match. [default: 0.2]')
    parser.add_argument('--nahr-recip', type=float, default=0.5, help='Reciprocal ' +
                        'overlap required for NAHR CNV match. [default: 0.5]')
    parser.add_argument('--min-jaccard-sum', type=float, default=2.0, help='Minimum ' +
                        'sum of Jaccard indexes to consider a locus as pleiotropic. ' +
                        '[default: 2.0]')
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
    loci_df, loci_hpos = read_bed_as_df(args.final_loci, output_hpos=True)
    loci_ids = loci_df.iloc[:, 3].tolist()
    all_df = loci_df.copy(deep=True)
    all_bt = loci_bt.saveas()

    # Load high-confidence GDs and merge into master df
    hc_gd_bt = pbt.BedTool(args.hc_gds).cut(range(5))
    hc_gd_df = read_bed_as_df(args.hc_gds)
    hc_gd_ids = hc_gd_df.iloc[:, 3].tolist()
    all_ids_in_hc_gd, hc_gd_ids_in_all = find_overlap(all_bt, hc_gd_bt, r=args.gd_recip)
    all_df = pd.concat([all_df, hc_gd_df.loc[~hc_gd_df.region_id.isin(hc_gd_ids_in_all), :]], axis=0)
    all_bt = pbt.BedTool().from_dataframe(all_df)

    # Load medium-confidence GDs and merge into master df
    mc_gd_bt = pbt.BedTool(args.mc_gds).cut(range(5))
    mc_gd_df = read_bed_as_df(args.mc_gds)
    mc_gd_ids = mc_gd_df.iloc[:, 3].tolist()
    all_ids_in_mc_gd, mc_gd_ids_in_all = find_overlap(all_bt, mc_gd_bt, r=args.gd_recip)
    all_df = pd.concat([all_df, mc_gd_df.loc[~mc_gd_df.region_id.isin(mc_gd_ids_in_all), :]], axis=0)
    all_bt = pbt.BedTool().from_dataframe(all_df)

    # Load low-confidence GDs and merge into master df
    lc_gd_bt = pbt.BedTool(args.lc_gds).cut(range(5))
    lc_gd_df = read_bed_as_df(args.lc_gds)
    lc_gd_ids = lc_gd_df.iloc[:, 3].tolist()
    all_ids_in_lc_gd, lc_gd_ids_in_all = find_overlap(all_bt, lc_gd_bt, r=args.gd_recip)
    all_df = pd.concat([all_df, lc_gd_df.loc[~lc_gd_df.region_id.isin(lc_gd_ids_in_all), :]], axis=0)
    all_bt = pbt.BedTool().from_dataframe(all_df)

    # Load predicted NAHR-mediated CNVs
    nahr_bt = pbt.BedTool(args.nahr_cnvs).cut(range(5))
    nahr_df = read_bed_as_df(args.nahr_cnvs)
    nahr_ids = nahr_df.iloc[:, 3].tolist()
    all_ids_in_nahr, nahr_ids_in_all = find_overlap(all_bt, nahr_bt, r=args.nahr_recip)
    all_df = pd.concat([all_df, nahr_df.loc[~nahr_df.region_id.isin(nahr_ids_in_all), :]], axis=0)
    all_bt = pbt.BedTool().from_dataframe(all_df)

    # Find regions overlapping common CNVs, if optioned
    benign_del_ids = overlap_common_cnvs(all_bt, args.common_dels, "DEL", args.common_cnv_cov)
    benign_dup_ids = overlap_common_cnvs(all_bt, args.common_dups, "DUP", args.common_cnv_cov)

    # Annotate merged df
    all_df['gw_sig'] = pd.get_dummies(all_df.region_id.isin(loci_ids), drop_first=True)
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
    path_ids = [i for i in cand_path_ids if i not in benign_ids]
    all_df['pathogenic'] = pd.get_dummies(all_df.region_id.isin(path_ids), drop_first=True)
    if args.common_dels is not None and args.common_dups is not None:
        all_df['benign'] = pd.get_dummies(all_df.region_id.isin(benign_ids), drop_first=True)
    all_df['nahr'] = pd.get_dummies(all_df.region_id.isin(set(nahr_ids + all_ids_in_nahr)), 
                                    drop_first=True)

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

    # Annotate with de novo mutations, if optioned
    if args.dnm_tsvs is not None:
        with open(args.dnm_tsvs) as dnm_ins:
            for study, dnm_path in csv.reader(dnm_ins, delimiter='\t'):
                all_df = count_dnms(all_df, dnm_path, study)

    
    # Sort & write out merged BED
    all_df.sort_values(by='chr start end cnv region_id'.split(), inplace=True)
    all_df.rename(columns={list(all_df.columns)[0] : '#chr'}).\
           to_csv(outfile, sep='\t', na_rep='NA', index=False)
    outfile.close()
    if args.outfile not in '- stdout'.split():
        if args.bgzip:
            bgzip(outfile_path)


if __name__ == '__main__':
    main()
