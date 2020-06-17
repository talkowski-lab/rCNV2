#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Count & weight CNVs overlapping genes
"""


import pybedtools as pbt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
import argparse
from sys import stdout
from os import path
import subprocess
import gzip


def process_gtf(gtf_in, blacklist=[]):
    """
    Read gtf and extract exons as pbt.BedTool
    Returns:
        1. pbt.BedTool of all exons
        2. pbt.BedTool of all exons after excluding genes in blacklist
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
                         from_string=True)
    sub_exonbt = exonbt.filter(lambda x: x[3] in blacklist).saveas()

    return exonbt, sub_exonbt


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cnvs', help='CNV BED file to compare vs GTF.')
    parser.add_argument('gtf', help='GTF of genes to consider.')
    parser.add_argument('--whitelist', help='List of genes to allow for more ' +
                        'lenient subset of noncoding CNVs.')
    parser.add_argument('--strict-outbed', help='Path to output BED for strictly ' +
                        'noncoding CNVs.')
    parser.add_argument('--loose-outbed', help='Path to output BED for more ' +
                        'lenient noncoding CNVs.')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output BED(s) with bgzip.')
    args = parser.parse_args()

    # Opens connections to output files
    if args.strict_outbed is not None:
        if path.splitext(args.strict_outbed)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            strict_out_path = path.splitext(args.strict_outbed)[0]
        else:
            strict_out_path = args.strict_outbed
    if args.loose_outbed is not None:
        if path.splitext(args.loose_outbed)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            loose_out_path = path.splitext(args.loose_outbed)[0]
        else:
            loose_out_path = args.loose_outbed

    # Loads CNV bed
    cnvbt = pbt.BedTool(args.cnvs)

    # Loads gene whitelist (if optioned)
    if args.whitelist is not None:
        whitelist = [g.rstrip() for g in open(args.whitelist).readlines()]
    else:
        whitelist = []

    # Extract relevant data from input GTF
    strict_exonbt, loose_exonbt \
        = process_gtf(args.gtf, whitelist)

    # Filter CNVs based on exon overlap
    strict_cnvbt = cnvbt.intersect(strict_exonbt, v=True, wa=True, header=True)
    loose_cnvbt = cnvbt.intersect(loose_exonbt, v=True, wa=True, header=True)

    # Write filtered CNVs to outfile
    if args.strict_outbed is not None:
        strict_cnvbt.saveas(strict_out_path)
        if args.bgzip:
            subprocess.run(['bgzip', '-f', strict_out_path])
    if args.loose_outbed is not None:
        loose_cnvbt.saveas(loose_out_path)
        if args.bgzip:
            subprocess.run(['bgzip', '-f', loose_out_path])


if __name__ == '__main__':
    main()
