#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Collect lists of genes overlapped by a BED file of regions
"""


import pybedtools as pbt
import argparse
from sys import stdout
from os import path
import gzip
import csv
import subprocess
import gzip


def process_gtf(gtf_in):
    """
    Read & process gtf
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

    # Filter & clean records in gtf
    def _filter_gtf(feature):
        """
        Restrict GTF features to desired elements
        """
        if feature.fields[2] == 'transcript' \
        and feature.attrs['gene_name'] in genes \
        and feature.attrs['transcript_id'] in transcripts:
            return True
        else:
            return False

    attrs_to_drop = 'gene_id gene_type gene_status transcript_type ' + \
                    'transcript_status transcript_name protein_id ' + \
                    'tag ccdsid havana_gene havana_transcript'
    attrs_to_drop = attrs_to_drop.split()

    def _clean_feature(feature):
        """
        Clean unnecessary fields & info from GTF features
        """
        for key in attrs_to_drop:
            if key in feature.attrs.keys():
                feature.attrs.pop(key)
        return feature

    gtfbt = gtfbt.filter(_filter_gtf).filter(_clean_feature).saveas()

    return gtfbt, genes, transcripts


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('regions', help='CNV BED file to compare vs windows.')
    parser.add_argument('gtf', help='GTF of genes to consider.')
    parser.add_argument('-o', '--outbed', help='Path to output file. ' +
                        '[default: stdout]')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output BED with bgzip.')
    
    args = parser.parse_args()

    # Open connection to output file
    if args.outbed is None \
    or args.outbed in 'stdout -'.split():
        outbed = stdout
    else:
        if path.splitext(args.outbed)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outbed_path = path.splitext(args.outbed)[0]
        else:
            outbed_path = args.outbed
        outbed = open(outbed_path, 'w')

    # Extract canonical transcripts from input GTF
    gtfbt, genes, transcripts = process_gtf(args.gtf)

    # Intersect transcripts with regions
    regbt = pbt.BedTool(args.regions).each(lambda x: x[0:4]).saveas()
    intbt = regbt.intersect(gtfbt, wa=True, wb=True).saveas()
    genedict = {}
    for reg in regbt:
        genes = []
        for info in [x.fields[-1] for x in intbt if x.name == reg.name]:
            genes.append(str([x.split()[1].replace('"', '') for x in info.split(';') \
                         if x.startswith('gene_name')][0]))
        genes = sorted(list(set(genes)))
        genedict[reg.name] = genes

    # Write original bed out to file with extra column for genes
    if path.splitext(args.regions)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
        csvin = gzip.open(args.regions, 'rt')
    else:
        csvin = open(args.regions)
    reader = csv.reader(csvin, delimiter='\t')
    for row in reader:
        if row[0].startswith('#'):
            newline = '\t'.join(row + ['genes']) + '\n'
        else:
            newline = '\t'.join(row + [';'.join(genedict.get(row[3], []))]) + '\n'
        outbed.write(newline)

    # Format output table and write to outfile
    if args.outbed is not None \
    and args.outbed not in 'stdout -'.split() \
    and args.bgzip:
        subprocess.run(['bgzip', '-f', outbed_path])


if __name__ == '__main__':
    main()

