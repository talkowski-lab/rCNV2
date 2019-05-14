#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Filter HPO terms for a list of sample-HPO pairs
"""


import csv
import argparse
from sys import stdout


def filter_phenotypes(phenos, pass_terms, outfile):
    """
    Read sample-HPO pairings from file, and filter HPO terms
    """

    with open(phenos) as infile:

        reader = csv.reader(infile, delimiter='\t')

        for sample, terms in reader:
            keep_terms = [t for t in terms.split(';') if t in pass_terms]
            outfile.write('{0}\t{1}\n'.format(sample, ';'.join(keep_terms)))


def main():
    """
    Main block
    """

    # Parse arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('phenos', help='Test file of phenotypes. One line per ' +
                        'patient. Two tab-delimited columns: patient ID and ' +
                        'string of semicolon-delimited HPO terms.', 
                        metavar='file')
    parser.add_argument('eligible_terms', help='File containing eligible HPO ' + \
                        'terms to be retained. Header lines with hashes will ' + \
                        'be ignored, and HPO terms must be in first column ' + \
                        'if multiple are provided.', metavar='file')
    parser.add_argument('-o', '--outfile', help='Path to outfile. ' +
                        '[default: stdout]', metavar='file')
    args = parser.parse_args()

    # Read list of eligible terms
    pass_terms = []
    for line in open(args.eligible_terms):
        li = line.strip()
        if not li.startswith('#'):
            pass_terms.append(li.split('\t')[0])

    # Open connection to outfile
    if args.outfile is None:
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Filter hpo terms per sample and write to outfile
    filter_phenotypes(args.phenos, pass_terms, outfile)


if __name__ == '__main__':
    main()

