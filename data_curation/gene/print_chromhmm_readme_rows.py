#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Format README rows for chromatin-based gene features
"""

import argparse
from sys import stdout
import csv

def main():
    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('manifest', help='tsv of state manifest')
    parser.add_argument('-o', '--outfile', help='path to outfile [default: stdout]',
    					default='stdout')
    args = parser.parse_args()

    # Open connection to outfile
    if args.outfile in 'stdout /dev/stdout -'.split():
        fout = stdout
    else:
        fout = open(args.outfile, 'w')

    # Print rows for states, one at a time
    with open(args.manifest) as fin:
        reader = csv.reader(fin, delimiter='\t'):
        for i, state, descrip, color, hexcode in reader:
            



if __name__ == '__main__':
	main()
