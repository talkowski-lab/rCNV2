#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Reformat ASC/SPARK rare CNVs for rCNV replication analysis
"""


import pandas as pd
import argparse
from sys import stdout
from os import path
import subprocess


def load_cnvs(cnvs_in):
    """
    Load and filter CNVs
    """

    cnvdf = pd.read_csv(cnvs_in, sep='\t')

    # Filter on:
    # - either unaffected siblings (1) or ASD probands (2)
    # - parental site frequency < 1%
    # - QS >= QS threshold
    # - >= 3 coding exons
    keepers = (cnvdf.asd.isin([1, 2]) & 
               (cnvdf.sf <= 0.01) &
               (cnvdf.QS >= cnvdf.QS_thresh) &
               (cnvdf.NEC >= 3))
    cnvdf_filt = cnvdf[keepers].loc[:, 'chr start end call sample asd NEC'.split()]

    # Rename columns
    cnvdf_filt.rename(columns={'chr' : '#chrom', 'call' : 'CNV', 
                               'asd' : 'phenotype'}, inplace=True)
    cnvdf_filt.phenotype = cnvdf_filt.phenotype.map({2 : 'ASD', 1 : 'control'})
    
    return cnvdf_filt


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cnvs', help='.bed of rare cnvs.')
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

