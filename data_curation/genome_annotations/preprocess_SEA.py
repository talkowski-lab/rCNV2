#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Parse simple SEA super-enhancer BED by cell types
"""


import argparse
import csv
import subprocess


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bed', help='Path to BED4 of super enhancers')
    parser.add_argument('outdir', help='Output directory')
    args = parser.parse_args()

    outfiles = {}

    with open(args.bed) as fin:
        for chrom, start, end, source in csv.reader(fin, delimiter='\t'):
            source = source.replace(' ', '_').replace('(', '').replace(')', '')
            if source not in outfiles.keys():
                outfiles[source] = open('{}/SEA.{}.bed'.format(args.outdir, source), 'w')
            outfiles[source].write('\t'.join([chrom, start, end]) + '\n')

    for outfile in outfiles.values():
        outpath = outfile.name
        outfile.close()
        subprocess.run(['bgzip', '-f', outpath])


if __name__ == '__main__':
    main()

