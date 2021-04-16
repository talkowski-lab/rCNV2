#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Convert an HGDP SV VCF (from Almarri, Cell, 2020) to conventional format 
able to be processed by GATK-SV, svtk, and related tools
"""


import argparse
from sys import stdin, stdout
import pysam
from os import path


def process_manta_record(record, cnv, k):
    """
    Reformat a single HGDP manta record
    """

    chrom = record.chrom
    start = record.pos
    altinfo = record.alleles[1].split(':')
    svlen = int(altinfo[1].split('=')[1])
    end = start + svlen

    newrec = record.copy()

    newrec.id = 'HGDP_{}_{}_{}'.format(cnv, chrom, k)
    newrec.info['SVTYPE'] = cnv
    newrec.stop = end
    newrec.info['SVLEN'] = svlen
    newrec.alleles = ('N', '<{}>'.format(cnv))
    newrec.qual = None

    return newrec


def main():
    # Command-line args
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('invcf', help='Input VCF (supports "stdin")')
    parser.add_argument('outvcf', help='Output VCF (supports "stdout")')
    parser.add_argument('-a', '--algorithm', default='manta', 
                        help='Specify algorithm [default: "manta"]')
    parser.add_argument('-c', '--cnv', help='Specify CNV type (manta only)')
    args = parser.parse_args()

    # Check if CNV type has been specified for Manta
    if args.cnv is None and args.algorithm == 'manta':
        from sys import exit
        exit('Must provide --cnv when using --algorithm manta')

    # Open connection to input VCF
    if args.invcf in '- stdin /dev/stdin'.split():
        vcf = pysam.VariantFile(stdin)
    else:
        vcf = pysam.VariantFile(args.invcf)
    header = vcf.header

    # Open connection to output VCF
    if args.outvcf in '- stdout /dev/stdout':
        outfile = pysam.VariantFile(stdout, 'w', header=header)
    else:
        outfile = pysam.VariantFile(args.outvcf, 'w', header=header)

    # Process each record in input VCF
    k = 0
    for record in vcf.fetch():
        if args.algorithm == 'manta':
            if 'PASS' in record.filter.keys():
                k += 1
                newrec = process_manta_record(record, args.cnv, k)
                outfile.write(newrec)


    # Close output and bgzip (if optioned)
    outfile.close()


if __name__ == '__main__':
    main()

