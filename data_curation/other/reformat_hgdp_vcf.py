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
from numpy import array
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

    newrec.id = 'HGDP_manta_{}_{}_{}'.format(cnv, chrom, k)
    newrec.info['SVTYPE'] = cnv
    newrec.stop = end
    newrec.info['SVLEN'] = svlen
    newrec.alleles = ('N', '<{}>'.format(cnv))
    newrec.qual = None

    return newrec


def process_gs_record(record, k):
    """
    Reformat a single GenomeSTRiP record
    """

    chrom = record.chrom
    start = record.pos
    altinfo = record.id.split('_')
    end = int(altinfo[3])
    svlen = end - start

    newrec = record.copy()

    nonref_cns = array([int(a.split('CN')[1].replace('>', '')) for a in record.alleles[1:] if a != '<CN2>'])
    loss, gain = False, False
    if any(nonref_cns < 2):
        loss = True
    if any(nonref_cns > 2):
        gain = True
    if loss and gain:
        cnv = 'CNV'
        filt = 'MULTIALLELIC'
    elif loss:
        cnv = 'DEL'
        filt = 'PASS'
    elif gain:
        cnv = 'DUP'
        filt = 'PASS'

    newrec.id = 'HGDP_GS_{}_{}_{}'.format(cnv, chrom, k)
    newrec.stop = end
    newrec.info['SVLEN'] = svlen
    newrec.info['SVTYPE'] = cnv
    newrec.alleles = ('N', '<{}>'.format(cnv))
    newrec.qual = None
    newrec.filter.clear()
    newrec.filter.add(filt)

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
    header.add_line('##FILTER=<ID=MULTIALLELIC,Description="mCNV site">')

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
        if args.algorithm in 'GS strip GenomeSTRiP'.split():
            if 'PASS' in record.filter.keys() \
            and len([a for a in record.alleles[1:] if a != '<CN2>']) > 0:
                k += 1
                newrec = process_gs_record(record, k)
                outfile.write(newrec)

    # Close output
    outfile.close()


if __name__ == '__main__':
    main()

