#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Subset an HPO co-occurrence table
"""


import argparse
import pandas as pd


def main():
    
    # Parse arguments & options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pairs_in', help='Input .tsv of HPO pairs.')
    parser.add_argument('pass_list', help='Input list of eligible HPO terms.')
    parser.add_argument('pairs_out', help='Output .tsv of filtered HPO pairs.')
    args = parser.parse_args()

    # Read input data
    pairs = pd.read_csv(args.pairs_in, sep='\t')
    elig_terms = [x.rstrip() for x in open(args.pass_list).readlines()]

    # Filter to pairs that are wholly within pass_list
    keepers = (pairs.iloc[:, 0].isin(elig_terms) & pairs.iloc[:, 1].isin(elig_terms))
    pairs = pairs.loc[keepers, :]

    # Write to outfile
    pairs.to_csv(args.pairs_out, sep='\t', index=False)


if __name__ == '__main__':
    main()

