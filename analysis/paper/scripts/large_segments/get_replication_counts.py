#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Collect CNV counts for replication attempt
"""


import pybedtools as pbt
import pandas as pd
import argparse
from sys import stdout
from os import path
import subprocess


def get_counts(segments_in, cnvs_in, excl_samples=[], prefix='replication', 
               case_phenotype='ASD', control_phenotype='control', ro=0.25):
    """
    Collect counts of case & control CNVs for each segment
    """

    cnvs_per_seg = {}
    rids_all = []

    # Load CNVs
    cnvs_bt = pbt.BedTool(cnvs_in).filter(lambda x: x[5] not in excl_samples).saveas()

    # Load credible intervals from segments
    segs_bt_l = []
    for sdat in pd.read_csv(segments_in, sep='\t').\
                   loc[:, 'region_id cnv cred_interval_coords'.split()].\
                   itertuples(index=False):
        rids_all.append(sdat.region_id)
        for cstr in sdat.cred_interval_coords.split(';'):
            chrom = cstr.split(':')[0]
            start = cstr.split(':')[1].split('-')[0]
            stop = cstr.split(':')[1].split('-')[1]
            segs_bt_l.append('\t'.join([chrom, start, stop, sdat.region_id, sdat.cnv]))
    segs_bt = pbt.BedTool('\n'.join(segs_bt_l), from_string=True).saveas()

    # Find all qualifying samples with a CNV per segment
    for hit in segs_bt.intersect(cnvs_bt, f=ro, r=True, wa=True, wb=True).\
                       filter(lambda x: x[4] == x[9]):
        rid = hit[3]
        sample_id = hit[10]
        pheno = hit[11]
        if rid not in cnvs_per_seg.keys():
            cnvs_per_seg[rid] = {case_phenotype : set(), control_phenotype : set()}
        cnvs_per_seg[rid][pheno].add(sample_id)

    # Count samples with CNVs
    case_counts = {}
    control_counts = {}
    for rid, sdicts in cnvs_per_seg.items():
        case_counts[rid] = len(sdicts[case_phenotype])
        control_counts[rid] = len(sdicts[control_phenotype])

    # Format as simple table
    rids = pd.Series(rids_all)
    case_series = rids.map(case_counts).fillna(0).astype(int)
    control_series = rids.map(control_counts).fillna(0).astype(int)
    out_df = pd.concat([rids, case_series, control_series], axis=1)
    out_df.columns = ['#region_id', 
                      '_'.join([prefix, case_phenotype]), 
                      '_'.join([prefix, control_phenotype])]
    
    return out_df


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('segments', help='.bed of final discovery segments')
    parser.add_argument('cnvs', help='.bed of cnvs to use for replication.')
    parser.add_argument('-x', '--exclude-samples', help='list of sample IDs to exclude')
    parser.add_argument('-p', '--prefix', default='replication', help='prefix ' +
                        'to prepend to column names of replication counts in --outfile')
    parser.add_argument('--case-phenotype', default='ASD', help='phenotype to ' +
                        'consider a case for purposes of replication')
    parser.add_argument('--control-phenotype', default='control', help='phenotype ' +
                        'to consider a case for purposes of replication')
    parser.add_argument('-r', '--recip', default=0.25, type=float, 
                        help='reciprocal overlap to require')
    parser.add_argument('-o', '--outfile', default='stdout', help='Output .bed of ' +
                        'segments with formatted counts. [default: stdout]')
    parser.add_argument('-z', '--gzip', dest='gzip', action='store_true',
                        help='Compress output .bed with gzip.')
    args = parser.parse_args()

    # Open connections to output files
    if args.outfile in 'stdout - /dev/stdout'.split():
        outfile = stdout
    else:
        if path.splitext(args.outfile)[-1] in '.gz .bz .bgz .gzip .gzip'.split():
            outfile_path = path.splitext(args.outfile)[0]
            gzip = True
        else:
            outfile_path = args.outfile
            gzip = args.gzip
        outfile = open(outfile_path, 'w')

    # Load list of excuded samples, if optioned
    if args.exclude_samples is not None:
        with open(args.exclude_samples) as fin:
            excl_samples = [s.rstrip() for s in fin.readlines()]
    else:
        excl_samples = []

    # Compile table of counts per segment
    count_df = get_counts(args.segments, args.cnvs, excl_samples, args.prefix,
                          args.case_phenotype, args.control_phenotype, args.recip)

    # Write counts table to outfile
    count_df.to_csv(outfile, sep='\t', index=False)
    if outfile != stdout:
        outfile.close()
        if gzip:
            subprocess.run(['gzip', '-f', outfile_path])


if __name__ == '__main__':
    main()

