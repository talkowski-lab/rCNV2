#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Collect features per gene for an input gtf
"""


import pybedtools as pbt
import argparse
from sys import stdout
from os import path
from pysam import TabixFile
from pandas import to_numeric
from numpy import nanmean, isnan
import subprocess
import gzip


def load_gtf(gtf_in):
    """
    Read gtf & filter to minimal info required
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
        if feature.attrs['gene_name'] in genes \
        and feature.attrs['transcript_id'] in transcripts:
            return True
        else:
            return False

    gtfbt = gtfbt.filter(_filter_gtf).saveas()

    return gtfbt, genes, transcripts


def pext_filter_feature(feature, pext, min_pext):
    """
    Subroutine to apply pext filter to single features from a BedTool
    """

    geneid = feature.attrs['gene_name']

    if feature.fields[2] == 'exon':
        pext_data = [x.split('\t')[-1] for x in \
                     pext.fetch(feature.chrom, feature.start, feature.end) \
                     if x.split('\t')[3] == geneid]
        pext_data = to_numeric(pext_data, errors='coerce')
        pext_data = pext_data[~isnan(pext_data)]
        return len(pext_data) == 0 or nanmean(pext_data) > min_pext
    else:
        return True


def pext_filter(gtfbt, pext, genes, min_pext):
    """
    Filter exons with mean pext below min_pext
    """

    # Filter exons
    gtfbt_filt = gtfbt.filter(pext_filter_feature, pext=pext, 
                              min_pext=min_pext).saveas()
    exons_before = len(gtfbt.filter(lambda x: x[2] == 'exon'))
    exons_after = len(gtfbt_filt.filter(lambda x: x[2] == 'exon'))
    n_exons_lost =  exons_before - exons_after

    # Count remaining exons per gene
    exons_per_gene = {gene : 0 for gene in genes}
    for exon in gtfbt_filt.filter(lambda x: x[2] == 'exon'):
        exons_per_gene[exon.attrs['gene_name']] += 1

    # Remove genes with no remaining exons
    genes_to_drop = [gene for gene, k in exons_per_gene.items() if k == 0]
    n_genes_lost = len(genes_to_drop)
    gtfbt_filt = gtfbt_filt.filter(lambda x: x.attrs['gene_name'] not in genes_to_drop).saveas()

    # Compile filtering stats
    filter_stats = {'n_genes_lost' : n_genes_lost,
                    'n_exons_lost' : n_exons_lost,
                    'genes_lost' : genes_to_drop}

    return gtfbt_filt, filter_stats
        

def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gtf', help='GTF of genes to consider.')
    parser.add_argument('pext', help='BED of pext scores. Must be tabixed.')
    parser.add_argument('--min-pext', default=0.1, type=float, 
                        help='Minimum mean pext score to retain exon. ' +
                        '[default: 0.1]')
    parser.add_argument('-o', '--outgtf', help='Path to output GTF file. ' +
                        '[default: stdout]')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output GTF with bgzip.')
    parser.add_argument('--lost-genes', help='Path to output file listing genes ' +
                        'lost due to pext filtering.')
    
    args = parser.parse_args()

    # Open connection to output file
    if args.outgtf is None \
    or args.outgtf in 'stdout -'.split():
        outgtf = stdout
    else:
        if path.splitext(args.outgtf)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outgtf_path = path.splitext(args.outgtf)[0]
        else:
            outgtf_path = args.outgtf

    # Load GTF & pext data
    gtfbt, genes, transcripts = load_gtf(args.gtf)
    pext = TabixFile(args.pext)

    # Apply pext filter
    gtfbt, filter_stats = pext_filter(gtfbt, pext, genes, args.min_pext)
    gtfbt.saveas(outgtf_path)
    filt_msg = 'Finished. Removed {:,} exons, resulting in the loss of {:,} genes.'
    print(filt_msg.format(filter_stats['n_exons_lost'], filter_stats['n_genes_lost']) + '\n')

    # Bgzip output GTF, if optioned
    if args.outgtf is not None \
    and args.outgtf not in 'stdout -'.split() \
    and args.bgzip:
        subprocess.run(['bgzip', '-f', outgtf_path])

    # Write out list of lost genes, if optioned
    if args.lost_genes is not None:
        with open(args.lost_genes, 'w') as lost_out:
            for gene in filter_stats['genes_lost']:
                lost_out.write(gene + '\n')


if __name__ == '__main__':
    main()

