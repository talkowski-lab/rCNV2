#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Compute Bayes factors and false discovery probabilties given assumed effect sizes
"""


from os import path
import numpy as np
import pandas as pd
from scipy.stats import norm
import argparse
from sys import stdout


def load_sumstats(stats_in):
    """
    Load & format association summary statistics
    """

    # Read data
    keep_cols = '#chr gene meta_lnOR meta_lnOR_lower meta_lnOR_upper'.split()
    ss = pd.read_csv(stats_in, delimiter='\t').loc[:, keep_cols]
    ss.set_axis(ss.gene, axis=0, inplace=True)
    ss.drop(labels='gene', axis=1, inplace=True)

    # Convert odds ratios to numeric
    ss = ss.apply(pd.to_numeric, errors='coerce')

    ss.rename(columns={'#chr' : 'chrom',
                       'meta_lnOR' : 'lnOR', 
                       'meta_lnOR_lower' : 'lnOR_lower',
                       'meta_lnOR_upper' : 'lnOR_upper'},
              inplace=True)

    return ss


def ci2se(ci):
    """
    Converts a tuple of (lower, upper) confidence interval bounds to standard error
    """

    ci = sorted(ci)

    return (ci[1] - ci[0]) / (2 * 1.96)


def calc_bf(lnOR, lnOR_lower, lnOR_upper, theta0, theta1, var0):
    """
    Calculate Bayes Factor for a single gene
    """

    # H0 : true effect size ≤ null effect size
    pr0 = 1 - norm.cdf(lnOR, theta0, np.sqrt(var0))

    # H1 : true effect size ≥ alternative effect size
    se = ci2se((lnOR_lower, lnOR_upper))
    pr1 = norm.cdf(lnOR, theta1, se)

    BF = pr0 / pr1

    return BF


def calc_bfdp(bf, prior):
    """
    Calculate BFDP for a single gene per Wakefield, AJHG, 2007
    """

    # Wakefield phrases prior as probability of no effect, so we need to invert
    prior = 1 - prior

    po = prior / (1 - prior)
    bftpo = bf * po

    bfdp = bftpo / (1 + bftpo)

    return bfdp


def calc_all_bfs(sumstats, theta0, theta1, var0, prior):
    """
    Calculate BFs and BFDPs for all genes
    """

    bfs = {}
    bfdps = {}
    
    for gene, stats in sumstats.transpose().to_dict().items():
        bf = calc_bf(stats['lnOR'], stats['lnOR_lower'], stats['lnOR_upper'], 
                     theta0, theta1, var0)
        bfdp = calc_bfdp(bf, prior)
        bfs[gene] = bf
        bfdps[gene] = bfdp

    sumstats['bf'] = sumstats.index.map(bfs)
    sumstats['bfdp'] = sumstats.index.map(bfdps)

    return sumstats


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('stats', help='bed.gz of meta-analysis association stats.')
    parser.add_argument('-o', '--outfile', default='stdout', help='Output .tsv of ' +
                        'gene BF/BFDPs per gene. [default: stdout]')
    parser.add_argument('--theta0', help='Null effect size (ln-scaled). ' +
                        '[default: 1.167]', default=1.167, type=float)
    parser.add_argument('--theta1', help='Alternative effect size (ln-scaled). ' +
                        '[default: 2.373]', default=2.373, type=float)
    parser.add_argument('--var0', help='Variance of null effect size (ln-scaled). ' +
                        '[default: 1.467]', default=1.467, type=float)
    parser.add_argument('-p', '--prior', help='Prior probability of intolerance ' +
                        'per gene. [default: 0.115]', default=0.115, type=float)
    args = parser.parse_args()

    # Open connections to output files
    if args.outfile in 'stdout - /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Import association summary stats
    sumstats = load_sumstats(args.stats)

    # Calculate BFs and BFDPs
    res = calc_all_bfs(sumstats, args.theta0, args.theta1, args.var0, args.prior)

    res.to_csv(outfile, sep='\t', index=True, na_rep='NA')


if __name__ == '__main__':
    main()
