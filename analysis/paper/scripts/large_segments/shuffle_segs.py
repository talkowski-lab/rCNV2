#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Randomly shuffles segments in the genome
"""


import argparse


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('segs', help='BED of segments to be ' +
                        'shuffled. Requires at least six columns: chrom, start, ' +
                        'end, id, cnv type, and coords listing the exact ' +
                        'intervals to be shuffled in chr:start-end format.')
    parser.add_argument('-w', '--whitelist', help='BED of regions to allow ' +
                        'during shuffling.')
    args = parser.parse_args()


    # Load segments and restrict to first six columns
    segs = pbt.BedTool(args.segs).cut(range(6)).saveas()

    import pdb; pdb.set_trace()




