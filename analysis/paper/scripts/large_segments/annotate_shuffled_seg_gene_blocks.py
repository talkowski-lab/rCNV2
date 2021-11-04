#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Annotates a tsv of gene blocks (as produced by shuffle_gene_blocks.py)
Designed for large segment permutation test matching on # of genes
"""


import csv
import argparse
from os.path import splitext
import gzip
from sys import stdout
import pandas as pd
from numpy import log10
from scipy.stats import hmean
import subprocess


dnm_csqs = 'lof mis syn'.split()


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


def preprocess_gtex(gtex_in, min_expression=1):
    """
    Preprocess GTEx expression matrix for segment annotation
    """

    # Load GTEx stats
    gtex = pd.read_csv(gtex_in, sep='\t').rename(columns = {'#gene' : 'gene'})
    gtex.index = gtex['gene']
    gtex.drop(columns='gene', inplace=True)

    # Compute arithmetic mean per gene
    gtex_means = gtex.mean(axis=1, skipna=True)

    # Collect ubiquitous expressor labels
    trans_min = log10(min_expression + 1)
    gtex_ubi = gtex.apply(lambda vals: all(vals.dropna() >= trans_min), axis=1)

    return gtex_means, gtex_ubi


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('blocks', help='tsv of gene blocks to be annotated. ' +
    					'Expects at least three columns: segment id, cnv type, and ' +
    					'semicolon-delimited list of genes to be annotated.')
    parser.add_argument('--gene-sets', help='tsv of gene sets to annotate. ' +
                        'Requires two columns: gene set name, and path to gene list')
    parser.add_argument('--hpo-genelists', help='tsv of HPO-specific genelists to ' +
                        'annotate. Will match on segment HPO. Two columns ' +
                        'expected: HPO, and path to genelist. Also requires ' +
                        'providing --segment-hpos.')
    parser.add_argument('--segment-hpos', help='tsv of segment ids and semicolon-' +
                        'delimited list of associated HPOs.')
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
    parser.add_argument('--gtex-matrix', help='Tsv gene X tissue expression levels ' +
                        'from GTEx. Will be used for various expression-based ' +
                        'gene annotations.')
    parser.add_argument('--min-expression', default=1, help='Minimum expression ' +
                        'level (in unscaled TPM) to consider a gene as "expressed". ' +
                        '[default: 1]')
    parser.add_argument('-o', '--outfile', help='Path to output tsv file. [default: ' +
                        'stdout]', default='stdout')
    parser.add_argument('-z', '--gzip', action='store_true', help='Compress ' + 
                        'output tsv file with gzip. [Default: do not compress]')
    args = parser.parse_args()

    gz_suf = 'gz gzip bgzip bgz bz'.split()

    # Open connection to infile
    if splitext(args.blocks)[-1].replace('.' ,'') in gz_suf:
        blocks_in = gzip.open(args.blocks, 'rt')
    else:
        blocks_in = open(args.blocks)

    # Open connection to outfile
    if args.outfile in '- stdout'.split():
        outfile = stdout
    else:
        
        if splitext(args.outfile)[-1].replace('.' ,'') in gz_suf:
            outfile_path = splitext(args.outfile)[0]
        else:
            outfile_path = args.outfile
        outfile = open(outfile_path, 'w')

    # Prep header vector
    header_cols = '#region_id cnv perm_idx n_genes'.split()

    # Read gene sets for annotation, if optioned
    genesets = {}
    if args.gene_sets is not None:
        with open(args.gene_sets) as gslin:
            reader = csv.reader(gslin, delimiter='\t')
            for gsname, gspath in reader:
                genesets[gsname] = [x.rstrip() for x in open(gspath).readlines()]
        header_cols += list(['_'.join(['n', x, 'genes']) for x in genesets.keys()])

    # Read info for HPO-matched gene lists, if optioned
    if args.hpo_genelists is not None and args.segment_hpos is not None:
        hpo_genes = load_hpo_genelists(args.hpo_genelists)
        with open(args.segment_hpos) as fin:
            segment_hpos = {bid : hpos.split(';') for bid, hpos \
                            in csv.reader(fin, delimiter='\t')}
        header_cols += ['n_HPOmatched_genes']

    # Load gene-level constraint metadata, if optioned
    if args.gnomad_constraint_tsv is not None:
        if args.gnomad_constraint_tsv.endswith('.bgz'):
            constr_df = pd.read_csv(args.gnomad_constraint_tsv, sep='\t', 
                                    compression='gzip')
        else:
            constr_df = pd.read_csv(args.gnomad_constraint_tsv, sep='\t')
        header_cols += ['min_LOEUF', 'min_MisOEUF', 'total_LoF_OE', 'total_mis_OE']

    # Load snv mutation rates, if optioned
    if args.snv_mus is not None:
        mu_df = pd.read_csv(args.snv_mus, sep='\t').\
                   rename(columns={'#gene' : 'gene'}).\
                   dropna()
    else:
        mu_df = None

    # Read DNM counts, if optioned
    if args.dnm_tsvs is not None:
        dnms = {}
        mus = {}
        xgenes = {}
        with open(args.dnm_tsvs) as dnm_in:
            reader = csv.reader(dnm_in, delimiter='\t')
            for study, dnm_path, xgene_list in reader:
                dnm_df = pd.read_csv(dnm_path, sep='\t').\
                            rename(columns={'#gene' : 'gene'})
                dnms[study] = dnm_df
                xgenes[study] = list(set([g.rstrip() for g in open(xgene_list).readlines()]))
                if mu_df is not None:
                    mus[study] = mu_df.copy(deep=True)
                for csq in dnm_csqs:
                    header_cols.append('_'.join([study, 'dnm', csq]))
                    if mu_df is not None:
                        n_dnms = dnms[study][csq].sum()
                        mus[study]['mu_' + csq] = mus[study]['mu_' + csq] * n_dnms
                        header_cols += ['_'.join([study, 'dnm', csq, 'obs_wMu']),
                                        '_'.join([study, 'dnm', csq, 'exp_wMu'])]
                for csq in dnm_csqs:
                    header_cols.append('_'.join([study, 'noSig', 'dnm', csq]))
                    if mu_df is not None:
                        header_cols += ['_'.join([study, 'noSig', 'dnm', csq, 'obs_wMu']),
                                        '_'.join([study, 'noSig', 'dnm', csq, 'exp_wMu'])]


    # Preprocess GTEx expression matrix, if optioned
    if args.gtex_matrix is not None:
        gtex_means, gtex_ubi = preprocess_gtex(args.gtex_matrix, min_expression=args.min_expression)
        header_cols += 'gene_expression_harmonic_mean n_ubiquitously_expressed_genes'.split()

    # Write header to outfile
    outfile.write('\t'.join(header_cols) + '\n')

    # Parse each row from input tsv
    reader = csv.reader(blocks_in, delimiter='\t')
    for bid, cnv, genes_str in reader:

        if '#' in bid:
            continue

        # Get basic region info
        orig_bid = '_'.join(bid.split('_')[1:])
        perm_idx = bid.split('_')[0].replace('perm', '')
        genes = genes_str.split(';')
        ngenes = len([g for g in genes if g != 'NA' and g != ''])
        outvals = [orig_bid, cnv, perm_idx, ngenes]

        # Annotate with --gene-sets
        for gsname, gsmems in genesets.items():
            outvals.append(len([g for g in genes if g in gsmems]))

        # Annotate with HPO-specific gene lists, if optioned
        if args.hpo_genelists is not None and args.segment_hpos is not None:
            hgenes = [g for l in [hpo_genes[h] for h in segment_hpos[orig_bid]] for g in l]
            outvals.append(len([g for g in genes if g in hgenes]))

        # Annotate with gene-level constraint metadata, if optioned
        if args.gnomad_constraint_tsv is not None:
            min_loeuf = constr_df.loc[constr_df.gene.isin(genes), 'oe_lof_upper'].min()
            min_misoeuf = constr_df.loc[constr_df.gene.isin(genes), 'oe_mis_upper'].min()
            outvals += [min_loeuf, min_misoeuf]
            for csq in 'lof mis'.split():
                obs_snv = constr_df.loc[constr_df.gene.isin(genes), 'obs_' + csq].sum()
                exp_snv = constr_df.loc[constr_df.gene.isin(genes), 'exp_' + csq].sum()
                if exp_snv > 0:
                    snv_oe = round(obs_snv / exp_snv, 6)
                else:
                    snv_oe = 'NA'
                outvals.append(snv_oe)

        # Annotate with count of de novo mutations per study, if optioned
        if args.dnm_tsvs is not None:
            for study in dnms.keys():
                dnm_df = dnms[study]
                # First annotate without excluding exome-wide significant genes from each study
                for csq in dnm_csqs:
                    outvals.append(dnm_df.loc[dnm_df.gene.isin(genes), csq].sum())
                    if mu_df is not None:
                        mu_df_x = mus[study]
                        genes_with_mus = mu_df.dropna()['gene'].tolist()
                        obs = dnm_df.loc[dnm_df.gene.isin(genes) & dnm_df.gene.isin(genes_with_mus), csq].sum()
                        exp = mu_df_x.loc[mu_df_x.gene.isin(genes) & mu_df_x.gene.isin(genes_with_mus), 'mu_' + csq].sum()
                        outvals += [obs, round(exp, 6)]
                # Annotate again after holding out exome-wide significant genes from each study
                for csq in dnm_csqs:
                    genes_x = [g for g in genes if g not in xgenes[study]]
                    outvals.append(dnm_df.loc[dnm_df.gene.isin(genes_x), csq].sum())
                    if mu_df is not None:
                        mu_df_x2 = mu_df_x.loc[~mu_df_x.gene.isin(xgenes[study]), :]
                        obs = dnm_df.loc[dnm_df.gene.isin(genes_x) & dnm_df.gene.isin(genes_with_mus), csq].sum()
                        exp = mu_df_x2.loc[mu_df_x2.gene.isin(genes_x) & mu_df_x2.gene.isin(genes_with_mus), 'mu_' + csq].sum()
                        outvals += [obs, round(exp, 6)]

        # Annotate with expression-based features, if optioned
        if args.gtex_matrix is not None:
            elig_genes = gtex_means.index.tolist()
            gtex_genes = [g for g in genes if g in elig_genes]
            if len(gtex_genes) > 0:
                mean_expr = hmean(gtex_means[gtex_means.index.isin(gtex_genes)])
                n_ubi = gtex_ubi[gtex_ubi.index.isin(gtex_genes)].sum()
                outvals += [round(mean_expr, 6), n_ubi]
            else:
                outvals += [0, 0]

        outfile.write('\t'.join([str(x) for x in outvals]) + '\n')

    outfile.close()

    # Gzip, if optioned
    if args.gzip:
        subprocess.run(['gzip', '-f', outfile_path])


if __name__ == '__main__':
    main()

