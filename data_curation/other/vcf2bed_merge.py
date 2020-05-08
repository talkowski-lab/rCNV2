#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Merge VCFs of common CNVs and convert to BED format
"""


import argparse
from sys import stdout
import pysam
from athena.utils import vcf2bed
import pybedtools as pbt


def main():
    """
    Main block
    """

    # Read input VCFs and output BED (if optioned)
    parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcfs', nargs='*', help='Input VCFs [can accept unlimited ' +
                        'VCFs]')
    parser.add_argument('-g', '--genome', help='BEDTools-style genome file ' +
                        '(for sorting).')
    parser.add_argument('-o', '--outfile', help='Path to output BED file [default: ' +
                        'stdout]')
    args = parser.parse_args()

    # Open connection to output BED
    if args.outfile in '- stdout':
        outfile = stdout
    else:
        outfile = args.outfile

    # Read input VCFs and convert to BED
    bts = [vcf2bed(pysam.VariantFile(x)) for x in args.vcfs]

    # Collapse BEDs and write to outfile
    if len(bts) > 1:
        pooled = bts[0].cat(*bts[0:], postmerge=False).sort(g=args.genome)
    else:
        pooled = bts[0].sort(g=args.genome)
    pooled.cut(range(3)).merge().saveas(outfile, trackline='#chr\tstart\tend')


if __name__ == '__main__':
    main()
