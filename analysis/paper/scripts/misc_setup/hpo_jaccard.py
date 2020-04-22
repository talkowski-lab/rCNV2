#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Compute Jaccard similarity index for all pairs of HPOs
"""


import numpy as np
import pandas as pd
import csv
from collections import Counter
from itertools import chain, product
import argparse
from sys import stdout


def count_samples(hpos, samples_tsv):
    """
    Count raw sample overlap for all pairs of HPOs
    """

    # Read all sample phenotypes
    phenos = []
    with open(samples_tsv) as fin:
        for sid, shpos in csv.reader(fin, delimiter='\t'):
            shpos_cleaned = set([h for h in shpos.split(';') if h in hpos])
            if len(shpos_cleaned) > 0:
                phenos.append(shpos_cleaned)

    # Compute counts of all pairs of hpos
    counts_dict = Counter(chain.from_iterable(product(x, repeat=2) for x in phenos))

    # Convert pairwise counts to pd.DataFrame of HPOs x HPOs
    zero_counts = np.zeros(shape=(len(hpos), len(hpos)))
    counts = pd.DataFrame(zero_counts, columns=hpos, index=hpos)
    counts.index.name = 'HPO'
    for hpoA in hpos:
        for hpoB in hpos:
            counts.loc[hpoA, hpoB] = counts_dict[(hpoA, hpoB)]

    return counts


def compute_jaccards(counts, hpos, asymmetric=False):
    """
    Normalize matrix of sample counts 
    """

    jaccards = counts.copy(deep=True)

    for hpoA in hpos:
        for hpoB in hpos:
            if hpoA == hpoB:
                jaccards.loc[hpoA, hpoB] = 1.0
            else:
                nA = counts.loc[hpoA, hpoA]
                nB = counts.loc[hpoB, hpoB]
                nAB = counts.loc[hpoA, hpoB]
                if asymmetric:
                    j = nAB / nA
                else:
                    j = nAB / (nA + nB - nAB)
                jaccards.loc[hpoA, hpoB] = j
    
    return jaccards


def main():
    """
    Main command-line block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('hpos', help='tsv of HPOs (assumes HPOs in last column if ' +
                        'multiple columns provided).')
    parser.add_argument('samples', help='Two-column tsv of sample ID & HPOs.')
    parser.add_argument('--jaccardfile', default='stdout', help='Output tsv of ' +
                        'reordered phenotypes. [default: stdout]')
    parser.add_argument('--countsfile', help='Optional tsv output of raw ' +
                        'sample counts matrix.')
    parser.add_argument('--asymfile', help='Optional tsv output of asymmetric ' + 
                        'sample overlap matrix.')
    args = parser.parse_args()
    
    # Open connections to output files
    if args.jaccardfile in 'stdout - /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.jaccardfile, 'w')
    if args.countsfile is not None:
        counts_out = open(args.countsfile, 'w')
    if args.asymfile is not None:
        asym_out = open(args.asymfile, 'w')

    # Read HPOs
    hpos = [x.rstrip().split('\t')[-1] for x in open(args.hpos).readlines()]

    # Read sample phenotypes and count raw overlaps 
    counts = count_samples(hpos, args.samples)

    # Normalize counts to compute Jaccard indexes and asymmetric sample overlaps
    jaccards = compute_jaccards(counts, hpos)
    asym = compute_jaccards(counts, hpos, asymmetric=True)

    # Write to outfile
    jaccards.to_csv(outfile, sep='\t', na_rep='NA')
    if args.countsfile is not None:
        counts.to_csv(counts_out, sep='\t', na_rep='NA')
    if args.asymfile is not None:
        asym.to_csv(asym_out, sep='\t', na_rep='NA')


if __name__ == '__main__':
    main()

