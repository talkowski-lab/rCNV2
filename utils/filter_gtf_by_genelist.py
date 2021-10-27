#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021- Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Filter features in GTF by gene list
"""


import pybedtools as pbt
import argparse
from sys import stdout
from os import path
import subprocess
        

def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gtf', help='GTF of genes to consider.')
    parser.add_argument('genelist', help='List of genes to be included.')
    parser.add_argument('-v', '--invert', action='store_true', help='Invert filter ' +
                        'such that all genes from genelist will be excluded ' +
                        '(analogous to grep -v).')
    parser.add_argument('--field', default='gene_name', 
                        help='Field to evaluate when filtering [default: "gene_name"].')
    parser.add_argument('-o', '--outgtf', help='Path to output GTF file. ' +
                        '[default: stdout]')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output GTF with bgzip.')
    
    args = parser.parse_args()

    # Load genelist
    genes = [x.rstrip() for x in open(args.genelist).readlines()]

    # Open connection to output file
    if args.outgtf is None \
    or args.outgtf in 'stdout -'.split():
        outgtf = stdout
    else:
        if path.splitext(args.outgtf)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outgtf_path = path.splitext(args.outgtf)[0]
        else:
            outgtf_path = args.outgtf

    # Filter GTF
    gtfbt = pbt.BedTool(args.gtf)
    if args.invert:
        filtered_gtf = gtfbt.filter(lambda f: f.attrs.get(args.field, None) not in genes)
    else:
        filtered_gtf = gtfbt.filter(lambda f: f.attrs.get(args.field, None) in genes)

    # Save filtered GTF and bgzip (if optioned)
    filtered_gtf.saveas(outgtf_path)
    if args.outgtf is not None \
    and args.outgtf not in 'stdout -'.split() \
    and args.bgzip:
        subprocess.run(['bgzip', '-f', outgtf_path])


if __name__ == '__main__':
    main()

