#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Compile BED of final segments, other known genomic disorders, and other NAHR segments
Outputs summary table for rCNV2 paper analyses
"""


cnvtypes = 'DEL DUP'.split()


import pandas as pd
import pybedtools as pbt
import argparse
from os.path import splitext
from sys import stdout
from athena.utils import bgzip


def read_bed_as_df(bed):
    """
    Reads a BED file as a pd.DataFrame and reformats data as needed
    """

    df = pd.read_csv(bed, sep='\t')
    c1h = list(df.columns)[0]
    df.rename(columns={c1h : c1h.replace('#', '')}, inplace=True)

    # Subset columns of interest
    keep_cols = 'cred_interval_coords cred_intervals_size n_genes genes'.split()
    df = pd.concat([df.iloc[:, :5], df.loc[:, df.columns.isin(keep_cols)]], axis=1)

    # Rename columns
    old_cnames = list(df.columns)
    df.rename(columns={old_cnames[0] : 'chr',
                      old_cnames[1] : 'start',
                      old_cnames[2] : 'end',
                      old_cnames[3] : 'region_id',
                      old_cnames[4] : 'cnv',
                      'cred_interval_coords' : 'coords',
                      'cred_intervals_size' : 'size'}, 
              inplace=True)

    # Add interval coordinates and size (if not optioned)
    if 'coords' not in df.columns:
        df['coords'] = df.apply(lambda x: '{}:{}-{}'.format(*x[['chr', 'start', 'end']]), axis=1)
    if 'size' not in df.columns:
        df['size'] = df.end - df.start

    cols_ordered = 'chr start end region_id cnv coords size n_genes genes'.split()    
    return df.loc[:, cols_ordered]


def find_overlap(bta, btb, r=0.01):
    """
    Overlaps two pbt.BedTool objects based on reciprocal overlap (r)
    Matches on CNV type (assumes the fifth column = CNV type)
    Returns: two lists of interval IDs of the overlapping set
    """

    ibt = bta.cut(range(5)).intersect(btb.cut(range(5)), f=r, r=True, wa=True, wb=True).\
              filter(lambda x: x[4] == x[9]).saveas()

    ids_a = [x[3] for x in ibt]
    ids_b = [x[8] for x in ibt]

    return ids_a, ids_b


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--final-loci', required=True, help='BED of final large ' +
                        'segment loci from association testing. Required.')
    parser.add_argument('--hc-gds', required=True, help='BED of high-confidence GDs. ' +
                        'Required.')
    parser.add_argument('--lc-gds', required=True, help='BED of low-confidence GDs. ' +
                        'Required.')
    parser.add_argument('--nahr-cnvs', required=True, help='BED of predicted NAHR-' +
                        'mediated CNVs. Required.')
    parser.add_argument('-o', '--outfile', help='Path to output BED file. [default: ' +
                        'stdout]', default='stdout')
    parser.add_argument('--gd-recip', type=float, default=0.2, help='Reciprocal ' +
                        'overlap required for GD match. [default: 0.2]')
    parser.add_argument('--nahr-recip', type=float, default=0.5, help='Reciprocal ' +
                        'overlap required for NAHR CNV match. [default: 0.5]')
    parser.add_argument('-z', '--bgzip', action='store_true', help='Compress ' + 
                        'output BED files with bgzip. [Default: do not compress]')
    args = parser.parse_args()

    # Open connection to outfile
    if args.outfile in '- stdout'.split():
        outfile = stdout
    else:
        gz_suf = 'gz gzip bgzip bgz bz'.split()
        if splitext(args.outfile)[-1].replace('.' ,'') in gz_suf:
            outfile_path = splitext(args.outfile)[0]
        else:
            outfile_path = args.outfile
        outfile = open(outfile_path, 'w')

    # Load final significant loci
    loci_bt = pbt.BedTool(args.final_loci).cut(range(5))
    loci_df = read_bed_as_df(args.final_loci)
    loci_ids = loci_df.iloc[:, 3].tolist()
    all_df = loci_df.copy(deep=True)
    all_bt = loci_bt.saveas()

    # Load high-confidence GDs and merge into master df
    hc_gd_bt = pbt.BedTool(args.hc_gds).cut(range(5))
    hc_gd_df = read_bed_as_df(args.hc_gds)
    hc_gd_ids = hc_gd_df.iloc[:, 3].tolist()
    all_ids_in_hc_gd, hc_gd_ids_in_all = find_overlap(all_bt, hc_gd_bt, r=args.gd_recip)
    all_df = pd.concat([all_df, hc_gd_df.loc[~hc_gd_df.region_id.isin(hc_gd_ids_in_all), :]], axis=0)
    all_bt = pbt.BedTool().from_dataframe(all_df)

    # Load low-confidence GDs and merge into master df
    lc_gd_bt = pbt.BedTool(args.lc_gds).cut(range(5))
    lc_gd_df = read_bed_as_df(args.lc_gds)
    lc_gd_ids = lc_gd_df.iloc[:, 3].tolist()
    all_ids_in_lc_gd, lc_gd_ids_in_all = find_overlap(all_bt, lc_gd_bt, r=args.gd_recip)
    all_df = pd.concat([all_df, lc_gd_df.loc[~lc_gd_df.region_id.isin(lc_gd_ids_in_all), :]], axis=0)
    all_bt = pbt.BedTool().from_dataframe(all_df)

    # Load predicted NAHR-mediated CNVs
    nahr_bt = pbt.BedTool(args.nahr_cnvs).cut(range(5))
    nahr_df = read_bed_as_df(args.nahr_cnvs)
    nahr_ids = nahr_df.iloc[:, 3].tolist()
    all_ids_in_nahr, nahr_ids_in_all = find_overlap(all_bt, nahr_bt, r=args.nahr_recip)
    all_df = pd.concat([all_df, nahr_df.loc[~nahr_df.region_id.isin(nahr_ids_in_all), :]], axis=0)

    # Annotate merged df
    all_df['gw_sig'] = pd.get_dummies(all_df.region_id.isin(loci_ids), drop_first=True)
    all_df['hc_gd'] = pd.get_dummies(all_df.region_id.isin(set(hc_gd_ids + all_ids_in_hc_gd)), 
                                     drop_first=True)
    all_df['lc_gd'] = pd.get_dummies(all_df.region_id.isin(set(lc_gd_ids + all_ids_in_lc_gd)), 
                                     drop_first=True)
    all_df['any_gd'] = pd.get_dummies(all_df.region_id.isin(set(hc_gd_ids + all_ids_in_hc_gd + \
                                                                lc_gd_ids + all_ids_in_lc_gd)), 
                                      drop_first=True)
    all_df['nahr'] = pd.get_dummies(all_df.region_id.isin(set(nahr_ids + all_ids_in_nahr)), 
                                    drop_first=True)
    
    # Sort & write out merged BED
    all_df.sort_values(by='chr start end cnv region_id'.split(), inplace=True)
    all_df.rename(columns={list(all_df.columns)[0] : '#chr'}).\
           to_csv(outfile, sep='\t', na_rep='NA', index=False)
    outfile.close()
    if args.outfile not in '- stdout'.split():
        if args.bgzip:
            bgzip(outfile_path)


if __name__ == '__main__':
    main()
