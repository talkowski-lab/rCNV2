#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Compute counts of CNV carriers vs non-carriers per phenotype based on overlap
with a prespecified gene list
"""


import pandas as pd
from os import path
import gzip
import csv
import argparse
import subprocess


def parse_cnvs(tsv_in, cohort, hpos, counts_out, gpass=[], 
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
    hpo_counts = {h : {'all' : 0,
                       'not_unconstrained' : 0,
                       'unconstrained_only' : 0} for h in hpos}

    # Iterate over CNVs and process one at a time
    for cnvid, phenos, ngenes, gstr in csv.reader(fin, delimiter='\t'):
        if cnvid in 'cnvid #cnvid'.split():
            continue

        if gstr == '':
            genes = []
        else:
            genes = gstr.split(';')
        ngenes = len(genes)

        if ngenes > 0:
            nonpass_hit = any(g not in gpass for g in genes)
        else:
            nonpass_hit = False

        # Add counts to hpos_counts
        if ngenes <= max_genes_for_summary:
            for hpo in phenos.split(';'):
                if hpo not in hpos:
                    continue

                hpo_counts[hpo]['all'] += 1
                if nonpass_hit:
                    hpo_counts[hpo]['not_unconstrained'] += 1
                else:
                    hpo_counts[hpo]['unconstrained_only'] += 1

    # Write counts to counts_out
    for hpo in hpos:
        for gset in 'all not_unconstrained unconstrained_only'.split():
            hit = hpo_counts[hpo][gset]
            cline = '\t'.join([str(x) for x in [cohort, hpo, gset, hit]])
            counts_out.write(cline + '\n')


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--genes-per-cnv', required=True, 
                        help='list of tsvs of cnv id, phenos, ngenes, genes. One ' +
                        'per cohort. First column is cohort id.')
    parser.add_argument('--hpos', required=True, help='list of hpos.')
    parser.add_argument('--unconstrained-genes', required=True, 
                        help='list of unconstrained genes.')
    parser.add_argument('--control-hpo', default='HEALTHY_CONTROL', help='HPO code ' +
                        'to use for control CNV counts. [default: HEALTHY_CONTROL]')
    parser.add_argument('--max-genes-for-summary', default=20000, type=int, 
                        help='Maximum number of genes per CNV to include in ' +
                        '--summary-counts.')
    parser.add_argument('--summary-counts', required=True, help='Output .tsv of ' +
                        'CNV counts per HPO.')
    parser.add_argument('-z', '--gzip', dest='gzip', action='store_true',
                        help='Compress output tsvs with gzip.')
    args = parser.parse_args()

    # Open connection to output files and write headers
    if path.splitext(args.summary_counts)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
        counts_out_path = path.splitext(args.summary_counts)[0]
    else:
        counts_out_path = args.summary_counts
    counts_out = open(counts_out_path, 'w')
    counts_out.write('\t'.join('#cohort hpo gset cnvs'.split()) + '\n')

    # Load unconstrained genes
    gpass = [x.rstrip() for x in open(args.unconstrained_genes).readlines()]

    # Load list of HPOs
    hpos = [x.rstrip() for x in open(args.hpos).readlines()]
    if args.control_hpo not in hpos:
        hpos.append(args.control_hpo)

    # Process CNV files for each cohort
    with open(args.genes_per_cnv) as tsvin:
        for cohort, cnvin in csv.reader(tsvin, delimiter='\t'):
            parse_cnvs(cnvin, cohort, hpos, counts_out, gpass, 
                       args.control_hpo, args.max_genes_for_summary)
    counts_out.close()

    # Gzip, if optioned
    if args.gzip:
        subprocess.run(['gzip', '-f', counts_out_path])


if __name__ == '__main__':
    main()

