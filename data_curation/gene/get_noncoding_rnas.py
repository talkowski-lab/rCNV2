#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Extract noncoding RNAs as BEDs from Gencode GTF
"""


import pybedtools as pbt
import argparse
from os import path, makedirs
import gzip
import subprocess


def extract_noncoding_bts(gtf_in, pc_txs):
    """
    Iterate over a GTF and extract noncoding RNAs by class
    """

    nc_bts = {}

    for x in pbt.BedTool(gtf_in):

        gtype = x.attrs['gene_type']

        # Skip protein-coding transcripts
        if gtype == 'protein_coding' or x.attrs['transcript_id'] in pc_txs:
            continue

        # Get basic gene features
        chrom, start, end = tuple([str(x) for x in [x.chrom, x.start, x.end]])
        gname = x.attrs['gene_name']
        feature = '\t'.join([chrom, start, end, gname]) + '\n'

        # Add element to pbt.BedTool corresponding to gtype
        if gtype in nc_bts.keys():
            nc_bts[gtype] += feature
        else:
            nc_bts[gtype] = feature

    for gtype, bt_str in nc_bts.items():
        nc_bts[gtype] = pbt.BedTool(bt_str, from_string=True)

    return nc_bts


# def process_gtf(gtf_in):
#     """
#     Read gtf & filter to minimal info required
#     """

#     gtfbt = pbt.BedTool(gtf_in)

#     # Build lists of eligible gene names and ensembl IDs
#     genes, ensg_ids, transcripts = [], [], []
#     ensg_to_gene, gene_to_ensg = {}, {}

#     for f in gtfbt:
#         if f.fields[2] == 'transcript':
#             gname = f.attrs['gene_name']
#             ensg_id = f.attrs['gene_id']
#             tname = f.attrs['transcript_id']
#             if gname not in genes:
#                 genes.append(gname)
#             if ensg_id not in ensg_ids:
#                 ensg_ids.append(ensg_id)
#             if tname not in transcripts:
#                 transcripts.append(tname)
#             if ensg_id not in ensg_to_gene.keys():
#                 ensg_to_gene[ensg_id] = gname
#             if gname not in gene_to_ensg.keys():
#                 gene_to_ensg[gname] = ensg_id

#     # Filter & clean records in gtf
#     def _filter_gtf(feature):
#         """
#         Restrict GTF features to desired elements
#         """
#         if feature.fields[2] in 'exon transcript'.split() \
#         and feature.attrs['gene_name'] in genes \
#         and feature.attrs['transcript_id'] in transcripts:
#             return True
#         else:
#             return False

#     attrs_to_drop = 'gene_id gene_type gene_status transcript_type ' + \
#                     'transcript_status transcript_name protein_id ' + \
#                     'tag ccdsid havana_gene havana_transcript'
#     attrs_to_drop = attrs_to_drop.split()

#     def _clean_feature(feature):
#         """
#         Clean unnecessary fields & info from GTF features
#         """
#         for key in attrs_to_drop:
#             if key in feature.attrs.keys():
#                 feature.attrs.pop(key)
#         return feature

#     gtfbt = gtfbt.filter(_filter_gtf).filter(_clean_feature).saveas()

#     # Make separate BedTools for exons and transcripts
#     txbt = gtfbt.filter(lambda x: x.fields[2] == 'transcript').saveas()
#     exonbt = gtfbt.filter(lambda x: x.fields[2] == 'exon').saveas()

#     return gtfbt, txbt, exonbt, genes, ensg_ids, transcripts, ensg_to_gene, gene_to_ensg


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gtf', help='Gencode gene annotations.')
    parser.add_argument('pc_transcripts', help='Gencode-style fasta of protein-' +
                        'coding transcripts. These will be excluded.')
    parser.add_argument('--outdir', default='./', help='Directory to which output ' +
                        'BED files will be written.')
    parser.add_argument('-p', '--prefix', default='noncoding_rnas', help='Prefix ' +
                        'to append to output BED files.')
    args = parser.parse_args()

    # Load coding transcripts
    if path.splitext(args.pc_transcripts)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
        pc_in = gzip.open(args.pc_transcripts, 'rt')
    else:
        pc_in = open(args.pc_transcripts)
    pc_txs = [x.rstrip().split('|')[0].replace('>', '') for x in pc_in]

    # Prep output directory, if needed
    if not path.exists(args.outdir):
        makedirs(args.outdir)

    # Process GTF
    nc_bts = extract_noncoding_bts(args.gtf, pc_txs)

    # Write out one BED per noncoding transcript class
    bheader = '#chrom\tstart\tend\tgenes\n'
    for gtype, ncbt in nc_bts.items():
        outpath = '{}/{}.{}.bed'.format(args.outdir, args.prefix, gtype)
        ncbt.sort().\
             merge(c=4, o='distinct').\
             saveas(outpath, trackline=bheader)
        subprocess.run(['bgzip', '-f', outpath])


if __name__ == '__main__':
    main()

