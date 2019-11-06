#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Compute quantiles of input data
"""


import argparse
import sys
import numpy as np


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='input values [default: stdin]', 
                        default='stdin', nargs='?')
    parser.add_argument('output', help='path to output [default: stdout].',
                        default='stdout', nargs='?')
    parser.add_argument('-q', '--quantiles', help='quantiles to compute ' + 
                        '(comma-delimited) [default: quartiles]', 
                        default='0,0.25,0.5,0.75,1.0')
    parser.add_argument('--no-header', action='store_true', default=False, 
                        help='suppress header [default: print header]')
    args = parser.parse_args()

    if args.input in 'stdin /dev/stdin -'.split():
        input = sys.stdin
    else:
        input = open(args.input)

    quantiles = np.array([float(x) for x in args.quantiles.split(',')])

    vals = np.array([float(x.rstrip()) for x in input.readlines()])

    quant_vals = np.quantile(a=vals, q=quantiles)

    if args.output in 'stdout /dev/stdout -'.split():
        output = sys.stdout
    else:
        output = open(args.output, 'w')

    if not args.no_header:
        output.write('#quantile\tvalue\n')

    for q, v in tuple(zip(quantiles, quant_vals)):
        output.write('\t'.join([str(q), str(int(np.round(v)))]) + '\n')


if __name__ == '__main__':
    main()

