#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Reformats de novo translocations & inversions from Redin 2017 supplementary table
"""


import pybedtools as pbt
import pandas as pd
import re
import argparse
from os.path import splitext
import csv
import gzip
import subprocess


def process_gtf(gtf_in):
    """
    Parse & filter GTF
    Returns: pbt.BedTool of gene bodies to consider
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
        Restrict GTF features to canonical transcripts
        """
        if feature.fields[2] == 'transcript' \
        and feature.attrs['gene_name'] in genes \
        and feature.attrs['transcript_id'] in transcripts:
            return True
        else:
            return False

    # Cut to simple BED of gene bodies
    tx_str = ''
    for tx in gtfbt.filter(_filter_gtf):
        tx_str += '\t'.join([tx.chrom, str(tx.start), str(tx.end), tx.attrs['gene_name']]) + '\n'
    txbt = pbt.BedTool(tx_str, from_string=True)

    return txbt, genes, transcripts


def load_bcas(redin_tsv, genome=None):
    """
    Preprocess tsv of Redin BCAs
    """

    elig_clusters = '+/+ +/- -/+ -/-'.split()

    bcadf = pd.read_csv(redin_tsv, sep='\t')

    # Subset to resolved BCA breakpoints that are either confirmed de novo or
    # segregate with disease in the family
    keepers = (bcadf['cluster'].isin(elig_clusters)) & \
              ((bcadf['inheritance'] == 'de novo') | \
               (bcadf['inheritance'].str.contains('Segregates')))
    bcadf = bcadf[keepers]

    # Clean sample IDs
    bcadf['sample'] = bcadf['sample'].map(lambda x: x.strip().replace('*', '').replace(' ', '_'))

    # Convert all breakpoints to pbt.BedTool for intersection
    bca_bt_str = ''
    for bca in bcadf.values.tolist():
        sid, cyto1, pos1, cyto2, pos2 = tuple(bca[i] for i in [0, 3, 4, 5, 6])
        chrom1, chrom2 = tuple(re.split('p|q', str(c))[0] for c in [cyto1, cyto2])
        bp1 = '\t'.join([chrom1, str(int(pos1)), str(int(pos1) + 1), str(sid)])
        bp2 = '\t'.join([chrom2, str(int(pos2)), str(int(pos2) + 1), str(sid)])
        bca_bt_str += bp1 + '\n' + bp2 + '\n'
    bca_bt = pbt.BedTool(bca_bt_str, from_string=True)
    if genome is None:
        bca_bt = bca_bt.sort().saveas()
    else:
        bca_bt = bca_bt.sort(g=genome).saveas()

    # Map BCA classifications to samples
    samples = list(set(bcadf['sample'].values.tolist()))
    samp_dict = {}
    for sid in samples:
        n_breaks = sum(bcadf['sample'] == sid)
        if n_breaks > 2:
            samp_dict[sid] = 'cpx'
        elif n_breaks == 2:
            bkpts = sorted(bcadf[bcadf['sample'] == sid].cluster.values.tolist())
            if bkpts == ['+/-', '-/+']:
                samp_dict[sid] = 'inv'
            else:
                samp_dict[sid] = 'tloc'
        elif n_breaks == 1:
            bkpt = bcadf[bcadf['sample'] == sid].cluster.values.tolist()[0]
            if bkpt in ['+/-', '-/+']:
                samp_dict[sid] = 'inv'
            else:
                samp_dict[sid] = 'tloc'

    return bca_bt, samp_dict


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--redin-tsv', help='tsv of BCA breakpoints from Redin et al., ' +
                        'Nat. Genet., 2017. Required.', required=True)
    parser.add_argument('--gtf', help='GTF of genes to consider. Required.',
                        required=True)
    parser.add_argument('--genes', help='list of gene symbols to consider. Required.',
                        required=True)
    parser.add_argument('-o', '--outfile', help='Path to output tsv file of gene ' +
                        'counts. [default: stdout]', default='stdout')
    parser.add_argument('-b', '--outbed', help='Path to output BED file of breakpoints.')
    parser.add_argument('-g', '--genome', help='BEDTools-style genome file to use ' +
                        'for sorting --outbed.')
    parser.add_argument('-z', '--gzip', action='store_true', help='Compress ' + 
                        'output tsv file with gzip. [Default: do not compress]')
    args = parser.parse_args()

    gz_suf = 'gz gzip'.split()

    # Open connection to outfiles
    if args.outfile in '- stdout'.split():
        outfile = stdout
    else:
        if splitext(args.outfile)[-1].replace('.' ,'') in gz_suf:
            outfile_path = splitext(args.outfile)[0]
        else:
            outfile_path = args.outfile
        outfile = open(outfile_path, 'w')
    if args.outbed is None:
        outbed_path = None
    else:
        if splitext(args.outbed)[-1].replace('.', '') in gz_suf:
            outbed_path = splitext(args.outbed)[0]
        else:
            outbed_path = args.outbed

    # Load GTF as pbt.BedTool of gene boundaries
    txbt, genes, transcripts = process_gtf(args.gtf)

    # Preprocess Redin BCAs
    bca_bt, samp_dict = load_bcas(args.redin_tsv, args.genome)

    # Build dicts of genes & BCA breakpoint sample IDs
    genes = [g.rstrip() for g in open(args.genes).readlines()]
    bca_classes = 'tloc inv cpx any_bca'.split()
    counts = {g : {c : 0 for c in bca_classes} for g in genes}
    samples = {g : {c : set() for c in bca_classes} for g in genes}

    # Intersect breakpoints with genes
    for hit in txbt.intersect(bca_bt, wa=True, wb=True):
        gene = hit[3]
        sample = hit[-1]
        bca_type = samp_dict[sample]
        samples[gene][bca_type].add(sample)
        samples[gene]['any_bca'].add(sample)

    # Count unique samples per gene
    counts = {g : {c : len(samples[g][c]) for c in bca_classes} for g in genes}

    # Write to outfiles
    header = '\t'.join(['#gene'] + bca_classes)
    outfile.write(header + '\n')
    for gene, gdat in counts.items():
         gline = gene + '\t' + '\t'.join([str(x) for x in gdat.values()])
         outfile.write(gline + '\n')
    outfile.close()
    bca_bt.saveas(outbed_path, trackline='\t'.join('#chrom start end sample'.split()))

    # Gzip output, if optioned
    if args.gzip:
        subprocess.run(['gzip', '-f', outfile_path])
        subprocess.run(['gzip', '-f', outbed_path])


if __name__ == '__main__':
    main()

