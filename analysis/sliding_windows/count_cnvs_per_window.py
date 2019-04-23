#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Count CNVs overlapping a set of bins
"""


import pybedtools
import argparse
from sys import stdout
from os import path
import subprocess
import gzip


# Parse list of input BED files
def parse_bed_list(bedlist, bednames=None):

    paths = [line.rstrip('\n') for line in open(bedlist)]

    if bednames is None:
        names = [path.split('/')[-1].split('.bed')[0] for path in paths]
    else:
        names = [line.rstrip('\n') for line in open(bednames)]

    if len(names) != len(paths):
        from sys import exit
        exit('ERROR: number of BED files does not match number of names')

    bed_table = list(zip(paths, names))

    return bed_table


# Get header from BED file
def get_bed_header(bedpath):

    if path.splitext(bedpath)[1] in '.gz .bgz .bgzip'.split():
        header = gzip.GzipFile(bedpath).readline().decode('utf-8').rstrip('\n')
    else:
        header = open(bedpath).readline().rstrip('\n')

    return header


# Main block
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bedlist', help='List of CNV BED files to compare vs windows.')
    parser.add_argument('windows', help='BED file of windows.')
    parser.add_argument('-f', '--fraction', help='Minimum fraction of window ' +
                        'covered by a CNV to count as overlapping. [default: 0.75]',
                        type=float, default=0.75)
    parser.add_argument('-c', '--cnv', help='Type of CNV to include (DEL/DUP). ' +
                        '[default: all]')
    parser.add_argument('-n', '--names', help='File listing names for each set ' +
                        'of CNVs to be appended to the header of the output file. ' +
                        'Order must match input bedlist. Otherwise, will be ' +
                        'inferred from the input filenames.')
    parser.add_argument('-o', '--outbed', help='Path to output file. ' +
                        '[default: stdout]')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output BED with bgzip.')
    
    args = parser.parse_args()

    cnv_bed_table = parse_bed_list(args.bedlist, args.names)

    bins = pybedtools.BedTool(args.windows)

    header = get_bed_header(args.windows)

    for cnvs, setname in cnv_bed_table:
        cnvbt = pybedtools.BedTool(cnvs)
        if args.cnv is not None:
            cnvbt = cnvbt.filter(lambda x: args.cnv in x)
        bins = bins.intersect(cnvbt, f=args.fraction, wa=True, c=True)
        header = '\t'.join([header, setname])

    if args.outbed is None \
    or args.outbed in 'stdout -'.split():
        bins.saveas(stdout, trackline=header)
    else:
        outbed = args.outbed
        if path.splitext(outbed)[1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outbed = path.splitext(outbed)[0]
        bins.saveas(outbed, trackline=header)
        if args.bgzip:
            subprocess.run(['bgzip', '-f', outbed])


if __name__ == '__main__':
    main()
