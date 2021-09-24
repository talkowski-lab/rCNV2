#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Replace coordinates of SVs in a VCF from a precomputed BED file
"""


import csv
import argparse
from sys import stdin, stdout
import pysam
from numpy import array
from os import path


def read_new_coords(inbed):
    """
    Read new coordinates to be replaced in VCF
    """

    coords = {}

    with open(inbed) as fin:
        for chrom, start, end, svid, cnv in csv.reader(fin, delimiter='\t'):
            coords[svid] = {'chrom' : chrom, 'start' : int(start), 'end' : int(end)}

    return coords


def main():
    # Command-line args
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('invcf', help='Input VCF (supports "stdin")')
    parser.add_argument('inbed', help='Input BED of new coords')
    parser.add_argument('newheader', help='New header')
    parser.add_argument('outvcf', help='Output VCF (supports "stdout")')
    args = parser.parse_args()

    # Load new coordinates as dict
    coords = read_new_coords(args.inbed)

    # Open connection to input VCF
    if args.invcf in '- stdin /dev/stdin'.split():
        vcf = pysam.VariantFile(stdin)
    else:
        vcf = pysam.VariantFile(args.invcf)

    # Load correct header
    header = pysam.VariantFile(args.newheader).header

    # Open connection to output VCF
    if args.outvcf in '- stdout /dev/stdout':
        outfile = pysam.VariantFile(stdout, 'w', header=header)
    else:
        outfile = pysam.VariantFile(args.outvcf, 'w', header=header)

    # Process each record in input VCF
    for record in vcf.fetch():
        if 'AF' in record.info.keys():
            if record.info['AF'] == 0:
                continue
        rid = record.id
        if rid in coords.keys():
            chrom = coords[rid]['chrom']
            if chrom not in header.contigs.keys():
                continue
            start = coords[rid]['start']
            end = coords[rid]['end']
            newrec = header.new_record(contig=chrom, start=start, stop=end, 
                                       alleles=record.alleles, id=rid, 
                                       qual=record.qual, filter=record.filter,
                                       info=record.info)
            newrec.stop = end
            outfile.write(newrec)

    # Close output
    outfile.close()


if __name__ == '__main__':
    main()
