#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Filter an input CNV BED file
"""


import argparse
from sys import stdin, stdout
from os import path
import pybedtools
import pysam
from math import ceil
import subprocess


# Generic BED vs BED frequency filtering function
def freq_filter(cnvsA, cnvsB, nsamp, maxFreq=0.01, ro=0.5, dist=50000):
    """
    Compare two sets of CNVs, and retain those in set A that have fewer than
    a maximum specified number of hits in set B
    """

    # Intersect
    xbed = cnvsA.intersect(cnvsB, wa=True, wb=True, f=ro, r=True)
    
    # Throw out self-hits
    xbed = xbed.filter(lambda x: str(x[3]) != str(x[9]))

    # Throw out CNV type mismatches
    def _cnv_match(feature):
        if str(feature[4]) == str(feature[10]) or str(feature[10]) == 'MCNV':
            return feature
    xbed = xbed.filter(_cnv_match)

    # Apply breakpoint distance requirement
    def _bp_dist(feature, dist):
        if abs(int(feature[1]) - int(feature[7])) <= dist \
        and abs(int(feature[2]) - int(feature[8])) <= dist:
            return feature[0:6]
    xbed = xbed.filter(_bp_dist, dist)

    # Get dictionary of hits per ID
    hits = {}
    for feature in xbed:
        fname = str(feature[3])
        if fname not in hits.keys():
            hits[fname] = 0
        hits[fname] += 1

    # Get list of CNVs to keep based on freq in xbed
    cutoff = ceil(maxFreq * nsamp)
    def _remove_fails(feature, hits, cutoff):
        fname = str(feature[3])
        if fname not in hits.keys() \
        or hits[fname] <= cutoff:
            return feature

    return cnvsA.filter(_remove_fails, hits, cutoff).saveas()


# Read gnomAD sites VCF for frequency filtering
def read_gnomad(vcfin, maxfreq):
    """
    Reads a gnomAD-SV sites vcf and converts it to a BedTool
    """
    vcf = pysam.VariantFile(vcfin)
    header = vcf.header

    def _filter_vcf(vcf, maxfreq):
        vcf_str = ''

        for record in vcf:

            if record.info['SVTYPE'] not in 'DEL DUP MCNV'.split():
                continue

            af = record.info['AF']
            if len(af) > 2:
                af = sum(af)
            else:
                af = af[0]
            if af <= maxfreq:
                continue

            flist = [str(record.chrom), str(record.pos), str(record.stop), 
                     str(record.id), str(record.info['SVTYPE']), 'CTRL']

            fstr = ' '.join(flist) + '\n'

            vcf_str = vcf_str + fstr

        return pybedtools.BedTool(vcf_str, from_string=True)

    gbed = _filter_vcf(vcf, maxfreq)

    return gbed


# Main block
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inbed', help='Input BED (supports "stdin").')
    parser.add_argument('outbed', help='Output BED (supports "stdout").')
    parser.add_argument('--minsize', type=int, help='Minimum CNV size. ' +
                        '[50kb]', default=50000)
    parser.add_argument('--maxsize', type=int, help='Maximum CNV size. ' +
                        '[10Mb]', default=10000000)
    parser.add_argument('-N', '--nsamp', type=int, help='Number of samples ' +
                        'represented in input BED.')
    parser.add_argument('--maxfreq', type=float, help='Maximum CNV ' + 
                        'frequency. [0.01]', default=0.01)
    parser.add_argument('--recipoverlap', type=float, help='Reciprocal overlap ' + 
                        'cutoff after which to consider two CNVs matching. [0.5]',
                        default=0.5)
    parser.add_argument('--dist', type=int, help='Maximum distance permitted ' + 
                        'between breakpoints to consider two CNVs matching. [50kb]',
                        default=50000)
    parser.add_argument('-x', '--blacklist', action='append', help='BED file ' +
                        'containing regions to blacklist based on CNV coverage. ' +
                        'May be specified multiple times.')
    parser.add_argument('--xcov', type=float, help='Maximum coverage ' + 
                        'by any blacklist before excluding a CNV. [0.3]',
                        default=0.3)
    parser.add_argument('--allcohorts', help='BED file of all unfiltered ' + 
                        'CNVs across all cohorts, to be used for freq filtering.')
    parser.add_argument('--allcohorts_nsamp', type=int, help='Number of samples ' +
                        'combined across all cohorts. Only used if --allcohorts ' +
                        'is specified.')
    parser.add_argument('-g', '--gnomad', help='gnomAD-SV VCF to use for ' +
                        'frequency filtering.')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output BED with bgzip.')

    args = parser.parse_args()

    # Read input bed
    if args.inbed in 'stdin -'.split():
        cnvs = pybedtools.BedTool(stdin)
    else:
        cnvs = pybedtools.BedTool(args.inbed)

    # Restrict to autosomes
    autosomes = [i for subl in ['{0} chr{1}'.format(c, c).split() for c in range(1, 23)] for i in subl]
    cnvs = cnvs.filter(lambda x: x.chrom in autosomes).saveas()

    # Restrict on self-intersection
    if args.nsamp is not None:
        cnvs = freq_filter(cnvs, cnvs, args.nsamp, maxFreq=args.maxfreq, 
                           ro=args.recipoverlap, dist=args.dist)

    # Restrict on blacklist coverage
    if args.blacklist is not None:
        def _bl_filter(feature, xcov):
            if float(feature[-1]) <= xcov:
                return feature
        for blpath in args.blacklist:
            bl = pybedtools.BedTool(blpath)
            cnvs = cnvs.coverage(bl).filter(_bl_filter, args.xcov).cut(range(0, 6))

    # Restrict on size
    cnvs = cnvs.filter(lambda x: len(x) >= args.minsize and len(x) <= args.maxsize).saveas()

    # Restrict on gnomAD 
    if args.gnomad is not None:
        gbed = read_gnomad(args.gnomad, args.maxfreq)
        cnvs = freq_filter(cnvs, gbed, 0, maxFreq=args.maxfreq, 
                           ro=args.recipoverlap, dist=args.dist)
    
    # Restrict on global intersection
    if args.allcohorts is not None \
    and args.allcohorts_nsamp is not None:
        cnvs = freq_filter(cnvs, args.allcohorts, args.allcohorts_nsamp,
                           maxFreq=args.maxfreq, ro=args.recipoverlap, 
                           dist=args.dist)

    # Write filtered CNVs out to file
    outheader = '#chr\tstart\tend\tname\tcnv\tpheno'
    if args.outbed in 'stdout -'.split():
        cnvs.saveas(stdout, trackline=outheader)
    else:
        outbed = args.outbed
        if '.gz' in outbed:
            outbed = path.splitext(outbed)[0]
        cnvs.saveas(outbed, trackline=outheader)
        if args.bgzip:
            subprocess.run(['bgzip', '-f', outbed])

if __name__ == '__main__':
    main()
