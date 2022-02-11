#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Collect CNV counts for replication attempt
"""


import argparse
from sys import stdout
from os import path
import subprocess


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('segments', help='.bed of final discovery segments')
    parser.add_argument('cnvs', help='.bed of cnvs to use for replication.')
    parser.add_argument('--case-phenotype', default='ASD', help='phenotype to ' +
                        'consider a case for purposes of replication')
    parser.add_argument('--control-phenotype', default='control', help='phenotype ' +
                        'to consider a case for purposes of replication')
    parser.add_argument('-o', '--outbed', default='stdout', help='Output .bed of ' +
                        'segments with formatted counts. [default: stdout]')
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

    # # Load and filter CNVs
    # cnvs = load_cnvs(args.cnvs)

    # # Write CNVs to outfile
    # cnvs.to_csv(outbed, sep='\t', index=False)
    # if outbed != stdout:
    #     outbed.close()
    #     if bgzip:
    #         subprocess.run(['bgzip', '-f', outbed_path])


if __name__ == '__main__':
    main()


