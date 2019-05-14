#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Reads an HPO .obo file and converts it to a two-column tsv matching all HPO
terms to their plain-text descriptions
"""


import networkx
import obonet
import argparse
from sys import stdout


def obo2tsv(hpo_g, outfile):
    """
    Iterate over all nodes in an obonet graph and write two-column tsv out
    """

    template = '{0}\t{1}\n'

    for term, data in hpo_g.nodes(data=True):

        newline = template.format(term, data.get('name', 'NA'))

        outfile.write(newline)


def main():
    """
    Main block
    """

    # Parse arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--obo', help='Path to HPO .obo file. [default: ' + 
                        'http://purl.obolibrary.org/obo/hp.obo]',
                        default='http://purl.obolibrary.org/obo/hp.obo',
                        metavar='(file|url)')
    parser.add_argument('-o', '--outfile', help='Path to outfile. ' +
                        '[default: stdout]', metavar='file')
    args = parser.parse_args()

    # Open connection to obo
    hpo_g = obonet.read_obo(args.obo)

    # Open connection to outfile
    if args.outfile is None:
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Add header to outfile
    outfile.write('#HPO\tdescription\n')

    # Convert obo to tsv
    obo2tsv(hpo_g, outfile)
    

if __name__ == '__main__':
    main()

