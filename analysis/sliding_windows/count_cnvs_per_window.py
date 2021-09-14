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


def get_bed_header(bedpath):
    """
    Get header from BED file
    """

    if path.splitext(bedpath)[1] in '.gz .bgz .bgzip'.split():
        header = gzip.GzipFile(bedpath).readline().decode('utf-8').rstrip('\n')
    else:
        header = open(bedpath).readline().rstrip('\n')

    return header


def process_cnvs(bedpath, pad_controls, control_hpo):
    """
    Read CNVs & extend control CNV breakpoints by a specified distance
    """

    cnvbt = pybedtools.BedTool(bedpath)

    def _pad_control_cnv(feature, pad_controls, control_hpo):
        """
        Add a fixed distance to control CNV breakpoints
        """

        if feature[5] == control_hpo:
            feature.start = max([0, feature.start-pad_controls])
            feature.stop = feature.stop+pad_controls
        
        return feature

    cnvbt = cnvbt.each(_pad_control_cnv, pad_controls, control_hpo).saveas()

    return cnvbt


def main():
    """
    Main block
    """

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cnvs', help='CNV BED file to compare vs windows.')
    parser.add_argument('windows', help='BED file of windows.')
    parser.add_argument('-f', '--fraction', help='Minimum fraction of window ' +
                        'covered by a CNV to count as overlapping. [default: 0.75]',
                        type=float, default=0.75)
    parser.add_argument('--pad-controls', help='Distance to be added to control ' +
                        'breakpoints. [default: 0]',
                        type=float, default=0)
    parser.add_argument('-t', '--type', help='Type of CNV to include (DEL/DUP). ' +
                        '[default: all]')
    parser.add_argument('--hpo', help='HPO term to consider for case samples. ' +
                        'If no --hpo is supplied, will count all CNVs ' +
                        'irrespective of phenotype.')
    parser.add_argument('--control-hpo', default='HEALTHY_CONTROL', help='HPO code ' +
                        'to use for control CNV counts. [default: HEALTHY_CONTROL]',
                        dest='control_hpo')
    parser.add_argument('-o', '--outbed', help='Path to output file. ' +
                        '[default: stdout]')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output BED with bgzip.')
    
    args = parser.parse_args()

    bins = pybedtools.BedTool(args.windows)

    header = get_bed_header(args.windows)

    cnvbt = process_cnvs(args.cnvs, args.pad_controls, args.control_hpo)

    if args.type is not None:
        if args.type != 'CNV':
            cnvbt = cnvbt.filter(lambda x: args.type in x.fields).saveas()

    if args.hpo is None:
        bins = bins.intersect(cnvbt, f=args.fraction, wa=True, c=True)
        header = '\t'.join([header, 'cnvs'])
    else:
        casebt = cnvbt.filter(lambda x: args.hpo in x.fields[5].split(';'))
        ctrlbt = cnvbt.filter(lambda x: args.control_hpo in x.fields[5].split(';'))
        bins = bins.intersect(casebt, f=args.fraction, wa=True, c=True)
        bins = bins.intersect(ctrlbt, f=args.fraction, wa=True, c=True)
        header = '\t'.join([header, 'case_cnvs', 'control_cnvs'])

    if args.outbed is None \
    or args.outbed in 'stdout -'.split():
        stdout.write(header + '\n')
        for rec in bins:
            stdout.write('\t'.join([str(x) for x in rec.fields]) + '\n')
    else:
        outbed = args.outbed
        if path.splitext(outbed)[1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outbed = path.splitext(outbed)[0]
        bins.saveas(outbed, trackline=header)
        if args.bgzip:
            subprocess.run(['bgzip', '-f', outbed])


if __name__ == '__main__':
    main()
