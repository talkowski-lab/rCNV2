#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Convert optimization data to summary statistics for gene scoring parameterization
"""


import pandas as pd
from os import path
import gzip
import csv
import argparse
import subprocess


def parse_cnvs(tsv_in, hpos, counts_out, cnvs_out, gpos=[], gneg=[], gexcl=[], 
               control_hpo='HEALTHY_CONTROL', max_genes_for_summary=20000):
    """
    Parse tsv of genes per CNV
    """

    # Open CNV file
    if path.splitext(tsv_in)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
        fin = gzip.open(tsv_in, 'rt')
    else:
        fin = open(tsv_in)

    # Prepare counter for cnv hits per HPO
    hpo_counts = {h : {'pos' : {'hit' : 0, 'no_hit' : 0},
                       'neg' : {'hit' : 0, 'no_hit' : 0}} for h in hpos}

    # Iterate over CNVs and process one at a time
    for cnvid, phenos, ngenes, gstr in csv.reader(fin, delimiter='\t'):
        if cnvid in 'cnvid #cnvid'.split():
            continue

        if control_hpo in phenos:
            cc = 'control'
        else:
            cc = 'case'

        genes = [g for g in gstr.split(';') if g not in gexcl]
        ngenes = len(genes)

        hits_pos = any(g in gpos for g in genes)
        npos = len([g for g in genes if g in gpos])
        hits_neg = any(g in gneg for g in genes)
        nneg = len([g for g in genes if g in gneg])
        noth = ngenes - (npos + nneg)

        # Write CNV stats to cnvs_out
        cnv_line = '\t'.join([str(x) for x in [cnvid, cc, ngenes, npos, nneg, noth]])
        cnvs_out.write(cnv_line + '\n')

        # Add counts to hpos_counts
        if ngenes <= max_genes_for_summary:
            for hpo in phenos.split(';'):
                if hpo not in hpos:
                    continue
                if hits_pos:
                    hpo_counts[hpo]['pos']['hit'] += 1
                else:
                    hpo_counts[hpo]['pos']['no_hit'] += 1
                if hits_neg:
                    hpo_counts[hpo]['neg']['hit'] += 1
                else:
                    hpo_counts[hpo]['neg']['no_hit'] += 1

    # Write counts to counts_out
    for hpo in hpos:
        for gset in 'pos neg'.split():
            hit, nohit = hpo_counts[hpo][gset].values()
            cline = '\t'.join([str(x) for x in [hpo, gset, hit, nohit]])
            counts_out.write(cline + '\n')

    # Close output files
    cnvs_out.close()
    counts_out.close()


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--genes-per-cnv', required=True, 
                        help='tsv of cnv id, phenos, ngenes, genes.')
    parser.add_argument('--hpos', required=True, help='list of hpos.')
    parser.add_argument('--positive-truth-genes', required=True, 
                        help='list of positive truth genes.')
    parser.add_argument('--negative-truth-genes', required=True, 
                        help='list of negative truth genes.')
    parser.add_argument('--exclude-genes',  help='list of genes to ignore.')
    parser.add_argument('--control-hpo', default='HEALTHY_CONTROL', help='HPO code ' +
                        'to use for control CNV counts. [default: HEALTHY_CONTROL]')
    parser.add_argument('--max-genes-for-summary', default=20000, type=int, 
                        help='Maximum number of genes per CNV to include in ' +
                        '--summary-counts.')
    parser.add_argument('--summary-counts', required=True, help='Output .tsv of ' +
                        'CNV counts per HPO.')
    parser.add_argument('--cnv-stats', required=True, help='Output .tsv of ' +
                        'CNV stats.')
    parser.add_argument('-z', '--gzip', dest='gzip', action='store_true',
                        help='Compress output tsvs with gzip.')
    args = parser.parse_args()

    # Open connection to output files and write headers
    if path.splitext(args.summary_counts)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
        counts_out_path = path.splitext(args.summary_counts)[0]
    else:
        counts_out_path = args.summary_counts
    counts_out = open(counts_out_path, 'w')
    counts_out.write('\t'.join('#hpo gset hit no_hit'.split()) + '\n')
    if path.splitext(args.cnv_stats)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
        cnvs_out_path = path.splitext(args.cnv_stats)[0]
    else:
        cnvs_out_path = args.cnv_stats
    cnvs_out = open(cnvs_out_path, 'w')
    cnvs_out.write('\t'.join('#cnvid pheno n_genes n_pos n_neg n_oth'.split()) + '\n')

    # Load positive, negative, and exclusion gene sets
    gpos = [x.rstrip() for x in open(args.positive_truth_genes).readlines()]
    gneg = [x.rstrip() for x in open(args.negative_truth_genes).readlines()]
    if args.exclude_genes is not None:
        gexcl = [x.rstrip() for x in open(args.exclude_genes).readlines()]
    else:
        gexcl = []

    # Load list of HPOs
    hpos = [x.rstrip() for x in open(args.hpos).readlines()]
    if args.control_hpo not in hpos:
        hpos.append(args.control_hpo)

    # Process CNV file
    parse_cnvs(args.genes_per_cnv, hpos, counts_out, cnvs_out, gpos, gneg, gexcl, 
               args.control_hpo, args.max_genes_for_summary)

    # Gzip, if optioned
    if args.gzip:
        subprocess.run(['gzip', '-f', counts_out_path])
        subprocess.run(['gzip', '-f', cnvs_out_path])


if __name__ == '__main__':
    main()

