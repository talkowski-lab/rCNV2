#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Create whitelist for CRB clustering
"""


import csv
import pybedtools as pbt
import argparse
from os import path
import subprocess


def load_genome(gfile):
    """
    Load BEDTools-style genome file (gfile) as pbt.BedTool of whole-chromosome intervals
    """

    gstr = ''

    with open(gfile) as fin:
        for chrom, end in csv.reader(fin, delimiter='\t'):
            gstr += '\t'.join([chrom, '0', end]) + '\n'

    return pbt.BedTool(gstr, from_string=True)


def process_gtf(gtf_in, whitelist=[]):
    """
    Read gtf and extract exons as pbt.BedTool
    Returns: pbt.BedTool of all exons for genes _not_ in whitelist
    """

    gtfbt = pbt.BedTool(gtf_in)

    # Build lists of eligible gene names and transcript IDs
    genes, transcripts = [], []

    for f in gtfbt:
        if f.fields[2] == 'transcript':
            gname = f.attrs['gene_name']
            tname = f.attrs['transcript_id']
            if gname not in genes:
                genes.append(gname)
            if tname not in transcripts:
                transcripts.append(tname)

    def _filter_gtf(feature):
        """
        Restrict GTF features to desired elements
        """
        if feature.fields[2] in 'exon'.split() \
        and feature.attrs['gene_name'] in genes \
        and feature.attrs['transcript_id'] in transcripts:
            return True
        else:
            return False

    def _simplify(feature):
        """
        Simplify feature to BED4
        """
        return '\t'.join([str(x) for x in [feature.chrom, feature.start, feature.end, 
                                           feature.attrs['gene_name']]])

    exonbt = pbt.BedTool('\n'.join([_simplify(x) for x in gtfbt.filter(_filter_gtf)]),
                         from_string=True).\
                 filter(lambda x: x[3] not in whitelist).\
                 saveas()

    return exonbt


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gtf', help='GTF of genes to consider.')
    parser.add_argument('whitelist', help='List of genes to allow in whitelist.')
    parser.add_argument('genome', help='BEDTools-style genome file.')
    parser.add_argument('outbed', help='Path to output BED.')
    parser.add_argument('-m', '--min-size', default=100000, type=int,
                        help='Minimum interval size to allow in final whitelist.')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output BED(s) with bgzip.')
    args = parser.parse_args()

    # Opens connections to output files
    if args.outbed is not None:
        if path.splitext(args.outbed)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            out_path = path.splitext(args.outbed)[0]
        else:
            out_path = args.outbed

    # Loads gene whitelist (if optioned)
    if args.whitelist is not None:
        whitelist = [g.rstrip() for g in open(args.whitelist).readlines()]
    else:
        whitelist = []

    # Read genome as pbt.BedTool
    gbt = load_genome(args.genome)

    # Extract pbt.BedTool of exons from input GTF
    exonbt = process_gtf(args.gtf, whitelist)
    
    # Filter regions of the genome not meeting criteria
    whitebt = gbt.subtract(exonbt).filter(lambda x: len(x) >= args.min_size)

    # Write filtered whitelist to outfile
    if args.outbed is not None:
        whitebt.saveas(out_path)
        if args.bgzip:
            subprocess.run(['bgzip', '-f', out_path])


if __name__ == '__main__':
    main()

