#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Create single BED file of pext scores per base per gene
"""


import json
from pandas import to_numeric
from numpy import nanmax
import argparse
from sys import stdout
from os import path
from pysam import TabixFile
import subprocess
import gzip


def process_pext_line(line, gene_field, pext_field, pan_tissue):
    """
    Extract single line from data in pext tsv and reformat as BED
    """

    dat = line.split('\t')

    chrom = dat[0]
    start = int(dat[1])
    end = start + 1
    annos = json.loads(dat[-1])[0]
    gene = annos[gene_field]

    if pan_tissue:
        ignore_fields = 'ensg csq symbol lof lof_flag'.split()
        vals = [v for k, v in annos.items() if k not in ignore_fields]
        val = str(nanmax(to_numeric(vals, errors='coerce')))
    else:
        val = str(annos[pext_field])
    if val in 'NaN nan'.split():
        val = '.'

    return chrom, start, end, gene, val


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pext', help='TSV of all pext scores')
    parser.add_argument('--gene-field', help='Field name in pext tsv to use as ' +
                        'gene identifier [default: "symbol"]', default='symbol')
    parser.add_argument('--pext-field', help='Field name in pext tsv to report as ' +
                        'pext value [default: "mean_proportion"]', 
                        default='mean_proportion')
    parser.add_argument('--pan-tissue', action='store_true', help='Report max pext ' +
                        'across all tissues. Ignores --pext-field. [Default: false]')
    parser.add_argument('--contig', help='Specify a single contig to process. ' +
                        '[default: process all autosomes]')
    parser.add_argument('-o', '--outbed', help='Path to output BED file. ' +
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
    outbed_header = '\t'.join('#chr start end gene'.split())

    # Iterate over chromosomes with tabix and process pext data
    with TabixFile(args.pext) as pextfile:
        if args.contig is None:
            contigs = [i + 1 for i in range(22)]
        else:
            contigs = [args.contig]
        for contig in contigs:
            processed = {}
            for line in pextfile.fetch(contig):
                chrom, start, end, gene, val \
                    = process_pext_line(line, args.gene_field, args.pext_field,
                                        args.pan_tissue)
                if gene not in processed.keys():
                    processed[gene] = [start]
                else:
                    if start not in processed[gene]:
                        pext_data = [chrom, start, end, gene, val]
                        outbed.write('\t'.join([str(x) for x in pext_data]) + '\n')
                        processed[gene].append(start)

    # Bgzip output, if optioned
    if args.outbed is not None \
    and args.outbed not in 'stdout -'.split() \
    and args.bgzip:
        subprocess.run(['bgzip', '-f', outbed_path])


if __name__ == '__main__':
    main()

