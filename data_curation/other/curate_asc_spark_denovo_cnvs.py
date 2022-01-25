#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Reformat ASC/SPARK de novo CNVs for rCNV analyses
"""


import csv
import pandas as pd
import argparse
from sys import stdout
from os import path
import subprocess


def load_cnvs(cnvs_in, max_size=20000000):
    """
    Load and filter de novo CNVs
    """

    elig_chroms = ['chr{}'.format(k) for k in range(1, 23)]

    cnvs = pd.read_csv(cnvs_in, sep='\t')
    
    keep_rows = ((cnvs.chr.isin(elig_chroms)) & (cnvs.end - cnvs.start <= max_size))
    keep_cols = cnvs.columns.isin('chr start end sample call Affected_Status'.split())

    cnvs = cnvs.loc[keep_rows, keep_cols]

    cnvs['pheno'] = cnvs.Affected_Status.map({1 : 'Control', 2 : 'ASD'})
    cnvs = cnvs['chr start end sample call pheno'.split()]
    return cnvs.rename(columns={'chr' : '#chr', 'call' : 'cnv'})


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cnvs', help='.tsv of de novo cnvs.')
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

    # Load and filter CNVs
    cnvs = load_cnvs(args.cnvs)

    # Write CNVs to outfile
    cnvs.to_csv(outbed, sep='\t', index=False)
    if outbed != stdout:
        outbed.close()
        if bgzip:
            subprocess.run(['bgzip', '-f', outbed_path])


if __name__ == '__main__':
    main()


