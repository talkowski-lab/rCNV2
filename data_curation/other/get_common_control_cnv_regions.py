#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Gather common CNV regions in controls from sliding window analysis summary stats
"""


import argparse
from os.path import splitext
from sys import stdout
import pandas as pd
import pybedtools as pbt
from athena.utils import bgzip


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('stats', help='BED of sliding window summary stats')
    parser.add_argument('-f', '--min-freq', type=float, default=0.01,
                        help='Minimum control CNV frequency')
    parser.add_argument('-n', '--n-controls', type=int, help='Total number of ' +
                        'control samples. Only used if input BED does not have ' +
                        'control_freq column.')
    parser.add_argument('-o', '--outfile', help='Path to output BED file. [default: ' +
                        'stdout]', default='stdout')
    parser.add_argument('-g', '--genome', help='Genome file (for sorting output).')
    parser.add_argument('-z', '--bgzip', action='store_true', help='Compress ' + 
                        'output BED files with bgzip. [Default: do not compress]')
    args = parser.parse_args()

    # Open connection to outfile
    if args.outfile in '- stdout'.split():
        outfile = stdout
    else:
        gz_suf = 'gz gzip bgzip bgz bz'.split()
        if splitext(args.outfile)[-1].replace('.' ,'') in gz_suf:
            outfile_path = splitext(args.outfile)[0]
        else:
            outfile_path = args.outfile
        outfile = open(outfile_path, 'w')

    # Load stats & calculate control frequency (if not provided)
    df = pd.read_csv(args.stats, delimiter='\t')
    if 'control_freq' not in df.columns:
        if args.n_controls is None or 'cnvs' not in df.columns:
            print('Must provide stats BED with control_freq column, or provide ' + \
                  'BED with cnvs column and specify --n-controls\n')
            exit()
        else:
            df['control_freq'] = df.cnvs / args.n_controls

    # Filter to sites with common control CNVs
    bt = pbt.BedTool().from_dataframe(df.loc[df.control_freq >= args.min_freq, :].iloc[:, :3])
    if args.genome is not None:
        bt = bt.sort(g=args.genome).merge()
    else:
        bt = bt.sort().merge()

    # Write to outfile
    for x in bt:
        outfile.write('\t'.join([x[0], str(x[1]), str(x[2])]) + '\n')
    outfile.close()
    if args.bgzip:
        bgzip(outfile_path)


if __name__ == '__main__':
    main()
