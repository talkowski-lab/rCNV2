#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Match two BED files of CNV segments bsed on gene overlap
"""


import argparse
from sys import stdout
from os import path
import pandas as pd


def intersect_genes(genes_a, genes_b, recip=False, frac=0.5):
    """
    Assess overlap between a set of genes from -a and all sets of genes in -b
    """

    ga = set([x for x in genes_a.split(';') if x != 'nan'])
    na = len(genes_a)

    if na < 1:
        return False

    for gb in genes_b:
        inter = len(ga.intersection(gb))
        if recip:
            denom = len(ga.union(gb))
        else:
            denom = na
        if inter / denom >= frac:
            return True

    return False


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-a', required=True, help='First BED file to be intersected')
    parser.add_argument('-b', required=True, help='Second BED file to be intersected')
    parser.add_argument('-v', '--invert', action='store_true', default=False, 
                        help='Return BED of intervals from -a that _don\'t_ ' +
                        'match any element in -b. Analogous to grep -v.')
    parser.add_argument('-f', '--frac', default=0.5, type=float, help='Fraction ' +
                        'of genes that must overlap between -a and -b to be ' +
                        'considered a match')
    parser.add_argument('-r', '--reciprocal', action='store_true', default=False,
                        help='Require reciprocal overlap between -a and -b. ' +
                        'Analogous to bedtools intersect -r. Default: only ' +
                        'consider overlap vs. genes in -a.')
    parser.add_argument('-o', '--outbed', default='stdout', help='Output .bed of ' +
                        'curated cnvs. [default: stdout]')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output .bed with bgzip.')
    args = parser.parse_args()

    # Open connections to output files
    if args.outbed in 'stdout - /dev/stdout'.split():
        outbed = stdout
    else:
        if path.splitext(args.outbed)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outbed_path = path.splitext(args.outbed)[0]
            bgzip = True
        else:
            outbed_path = args.outbed
            bgzip = args.bgzip
        outbed = open(outbed_path, 'w')

    # Load -a and -b as pd.DataFrames
    dfa = pd.read_csv(args.a, sep='\t')
    dfb = pd.read_csv(args.b, sep='\t')

    # Extract genes from dfb for comparison
    if 'genes' in dfb.columns:
        dfb_gcol = dfb.genes
    else:
        dfb_gcol = dfb.iloc[:, -1]
    genes_b = [set(x.split(';')) for x in dfb_gcol.astype(str).tolist()]

    # Compute inclusion/exclusion mask for dfa
    if 'genes' in dfa.columns:
        dfa_gcol = dfa.genes
    else:
        dfa_gcol = dfa.iloc[:, -1]
    include = dfa_gcol.astype(str).\
                       apply(intersect_genes, genes_b=genes_b, 
                             recip=args.reciprocal, frac=args.frac)
    if args.invert:
        include = ~include

    # Write subset to outfile
    dfa[include].to_csv(outbed, sep='\t', na_rep='NA', index=False)
    if outbed != stdout:
        outbed.close()
        if bgzip:
            subprocess.run(['bgzip', '-f', outbed_path])

’’’
if __name__ == '__main__':
    main()

