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
    parser.add_argument('--dnm-tsvs', help='Tsv of de novo mutation counts ' +
                        ' to annotate. Two columns expected: study prefix, and ' +
                        'path to tsv with dnm counts.')
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

    # Read DNM counts, if optioned
    if args.dnm_tsvs is not None:
        dnms = {}
        with open(args.dnm_tsvs) as dnm_in:
            reader = csv.reader(dnm_in, delimiter='\t')
            for study, dnm_path in reader:
                dnm_df = pd.read_csv(dnm_path, sep='\t').\
                            rename(columns={'#gene' : 'gene'})
                dnms[study] = dnm_df
                header_cols += ['_'.join([study, 'dnm', csq]) for csq in dnm_csqs]

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

        # Annotate with count of de novo mutations per study, if optioned
        if args.dnm_tsvs is not None:
            for dnm_df in dnms.values():
                for csq in dnm_csqs:
                    outvals.append(dnm_df.loc[dnm_df.gene.isin(genes), csq].sum())

        outfile.write('\t'.join([str(x) for x in outvals]) + '\n')

    outfile.close()

    # Gzip, if optioned
    if args.gzip:
        subprocess.run(['gzip', '-f', outfile_path])


if __name__ == '__main__':
    main()

