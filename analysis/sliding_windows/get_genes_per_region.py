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
import pandas as pd
import re
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


def annotate_regions(regions_path, gtfbt):
    """
    Load regions and annotate with genes
    """

    regions = pd.read_csv(regions_path, sep='\t')
    # regions.rename(columns={regions.columns.tolist()[0] : \
    #                         regions.columns.tolist()[0].replace('#', '')},
    #                inplace=True)

    intervals_dict = {rid : i.split(';') for rid, i in regions.iloc[:, [3, -2]].values}

    intervals_str = ''
    for rid, ints in intervals_dict.items():
        for i in ints:
            intervals_str += '\t'.join([re.sub(':|-', '\t', i), rid]) + '\n'
    intervals_bt = pbt.BedTool(intervals_str, from_string=True)

    genes_dict = {x : [] for x in regions.iloc[:, 3].values}
    for x in intervals_bt.intersect(gtfbt, wa=True, wb=True):
        gene = str([f.split()[1].replace('"', '') for f in x[-1].split(';') \
                    if f.startswith('gene_name')][0])
        genes_dict[x.name].append(gene)
    
    genes_dict = {rid : sorted(list(set(g))) for rid, g in genes_dict.items()}
    ngenes_dict = {rid : len(g) for rid, g in genes_dict.items()}
    genes_str_dict = {rid : ';'.join(g) for rid, g in genes_dict.items()}

    regions['n_genes'] = regions.iloc[:, 3].map(ngenes_dict)
    regions['genes'] = regions.iloc[:, 3].map(genes_str_dict)

    return regions.sort_values(by=regions.columns[0:3].tolist(), axis=0)


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('regions', help='BED file of final regions after refinement. ' +
                        'Fourth column must be unique row name and wecond-to-last ' +
                        'column must be a semicolon-delimited list of chr:start-end ' +
                        'intervals to compare vs gtf.')
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

    # Load regions & annotate with genes
    regions = annotate_regions(args.regions, gtfbt)

    # Sort & write original bed out to file with extra columns for n_genes and genes
    regions.to_csv(outbed, sep='\t', na_rep='NA', header=True, index=False)
    outbed.close()

    # Bgzip output, if optioned
    if args.outbed is not None \
    and args.outbed not in 'stdout -'.split() \
    and args.bgzip:
        subprocess.run(['bgzip', '-f', outbed_path])


if __name__ == '__main__':
    main()

