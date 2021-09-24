#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Determine effective number of samples matching at least one of a list of HPOs
"""


import argparse
import csv


def main():
    """
    Main block
    """

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('samples', help='tsv of sample ID & HPO terms.')
    parser.add_argument('hpos', help='List of HPO terms to keep.')
    args = parser.parse_args()

    keepers = set([x.rstrip() for x in open(args.hpos).readlines()])

    k = 0
    with open(args.samples) as fin:
        for sample, hpos in csv.reader(fin, delimiter='\t'):
            if len(set(hpos.split(';')).intersection(set(keepers))) > 0:
                k += 1

    print(str(k))


if __name__ == '__main__':
    main()

