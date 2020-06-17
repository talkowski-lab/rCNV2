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

rep = 'Roadmap Epigenomics Project'

def format_line(i, state, descrip, prefix, suffix):
    
    name = descrip.lower().replace('tss', 'TSS')
    
    l1 = '{} {} coverage'.format(prefix, name)
    l2 = '_'.join(['chromhmm', str(i), state, suffix])
    l3 = '{} gene coverage by {} from ChromHMM across 98 tissues from the {}'.format(prefix, name, rep)
    l4 = rep + ' [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248)'

    line = '| ' + ' | '.join([l1, l2, l3, l4]) + ' |  \n'

    return line

    # | Mean active TSS coverage | `chromhmm_1_TssA_mean` | Mean gene coverage by active TSS regions from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
    # | Standard deviation of active TSS coverage | `chromhmm_1_TssA_sd` | Standard deviation of gene coverage by active TSS regions from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  


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
        reader = csv.reader(fin, delimiter='\t')
        for i, state, descrip, color, hexcode in reader:
            if '#' in i:
                continue
            fout.write(format_line(i, state, descrip, "Mean", "mean"))
            fout.write(format_line(i, state, descrip, "Standard deviation of", "sd"))
    fout.close()


if __name__ == '__main__':
	main()
