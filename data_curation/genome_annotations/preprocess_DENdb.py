#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Parse DENdb enhancer track
"""


import argparse
import subprocess


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('csv', help='Path to DENdb enhancers.csv')
    parser.add_argument('outdir', help='Output directory')
    args = parser.parse_args()

    outfiles = {}

    for line in open(args.csv):
        eid, chrom, start, end, clid = line.rstrip().split(',')[0:5]
        if clid not in outfiles.keys():
            outfiles[clid] = open('{}/DENdb.{}.bed'.format(args.outdir, clid), 'w')
        outfiles[clid].write('\t'.join([chrom, start, end]) + '\n')

    for outfile in outfiles.values():
        outpath = outfile.name
        outfile.close()
        subprocess.run(['bgzip', '-f', outpath])


if __name__ == '__main__':
    main()

