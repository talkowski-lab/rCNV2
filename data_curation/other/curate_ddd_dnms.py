#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Curates de novo mutations from DDD supplemental table
Expects: https://www.biorxiv.org/content/biorxiv/early/2020/04/01/797787/DC4/embed/media-4.txt
"""


import argparse
from os.path import splitext
import csv
import gzip
import subprocess

lof_csqs = 'frameshift_variant stop_gained splice_donor_variant splice_acceptor_variant'.split()


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dnm-tsv', help='tsv of de novo mutations from Kaplanis ' +
                        'et al., bioRxiv, 2020. Required.', required=True)
    parser.add_argument('--genes', help='list of gene symbols to consider. Required.',
                        required=True)
    parser.add_argument('-o', '--outfile', help='Path to output tsv file. [default: ' +
                        'stdout]', default='stdout')
    parser.add_argument('-z', '--gzip', action='store_true', help='Compress ' + 
                        'output tsv file with gzip. [Default: do not compress]')
    args = parser.parse_args()

    gz_suf = 'gz gzip'.split()

    # Open connection to infile
    if splitext(args.dnm_tsv)[-1].replace('.' ,'') in gz_suf:
        dnms_in = csv.reader(gzip.open(args.dnm_tsv, 'rt'), delimiter='\t')
    else:
        dnms_in = csv.reader(open(args.dnm_tsv), delimiter='\t')

    # Open connection to outfile
    if args.outfile in '- stdout'.split():
        outfile = stdout
    else:
        if splitext(args.outfile)[-1].replace('.' ,'') in gz_suf:
            outfile_path = splitext(args.outfile)[0]
        else:
            outfile_path = args.outfile
        outfile = open(outfile_path, 'w')

    # Build dict of genes & counts
    zeros = {'lof' : 0, 'mis' : 0, 'syn' : 0}
    gdict = {g.rstrip() : zeros.copy() for g in open(args.genes).readlines()}

    # Parse each line of DNM file
    for var, chrom, pos, ref, alt, csq, gene, study, prop, hgnc in dnms_in:
        
        if gene not in gdict.keys():
            continue

        if csq in lof_csqs:
            gdict[gene]['lof'] += 1
        elif 'missense' in csq:
            gdict[gene]['mis'] += 1
        elif 'synonymous' in csq:
            gdict[gene]['syn'] += 1

    # Write to outfile
    header = '\t'.join('#gene lof mis syn'.split())
    outfile.write(header + '\n')
    for gene, gdat in gdict.items():
         gline = gene + '\t' + '\t'.join([str(x) for x in gdat.values()])
         outfile.write(gline + '\n')
    outfile.close()

    # Gzip output, if optioned
    if args.gzip:
        subprocess.run(['gzip', '-f', outfile_path])


if __name__ == '__main__':
    main()

