#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Determine the minimal set of most specific HPO terms to retain for analysis 
based on a list of patients and corresponding HPO terms
"""


import pybedtools as pbt
import numpy as np
from athena.utils.misc import determine_filetype, bgzip
from athena.utils.dfutils import float_cleanup
import gzip
import csv
import pandas as pd
import argparse
from athena.mutrate import annotate_bins
from os import path
from sys import stdout


def expand_interval(interval, min_size):
    """
    Uniformly expand intervals smaller than min_size
    """

    if len(interval) < min_size:
        midpoint = (interval.end - interval.start) / 2
        pad = np.ceil(min_size / 2)
        new_start = max([0, interval.start - pad])
        d_left = midpoint - new_start
        new_end = np.ceil(midpoint + min_size - d_left)
        interval.start = new_start
        interval.end = new_end

    return interval


def load_intervals(bed_in, min_size, cols_to_keep):
    """
    Load and expand intervals
    """

    if 'compressed' in determine_filetype(bed_in):
        fin = gzip.open(bed_in, 'rt')
    else:
        fin = open(bed_in)
    header = fin.readline().rstrip().split('\t')[:cols_to_keep]
    fin.close()

    bt = pbt.BedTool(bed_in).each(expand_interval, min_size=min_size).\
                             cut(range(cols_to_keep)).\
                             saveas('resized_intervals.bed', 
                                    trackline='\t'.join(header))

    return bt, header


def parse_probesets(probes_tsv):
    """
    Parse probes_tsv input
    """

    tracks, tnames = [], []

    with open(probes_tsv) as fin:
        for bedpath, tname in csv.reader(fin, delimiter='\t'):
            tracks.append(bedpath)
            tnames.append(tname)

    return tracks, tnames


def replace_coords(intervals, orig_bed, cols_to_keep):
    """
    Replace expanded coordinates for all intervals with their original coordinates
    """

    old_df = pd.read_csv(orig_bed, sep='\t')
    new_df = intervals.to_dataframe(names=range(intervals.field_count()))

    fixed_df = pd.concat([old_df.iloc[:, range(cols_to_keep)],
                          new_df.iloc[:, cols_to_keep:(new_df.shape[1])]],
                         ignore_index=True, axis=1)

    return pbt.BedTool().from_dataframe(fixed_df).saveas()


def label_array_fails(intervals, min_probes, bed_header, tnames):
    """
    Generate pd.DataFrame of pass/fail labels per interval per array
    """

    bed_df = intervals.to_dataframe(names=bed_header + tnames, comment='#')

    for track in tnames:
        bed_df[track] = (bed_df[track] >= min_probes)
    
    return bed_df


def load_array_counts(samples_tsv, cohorts, elig_arrays):
    """
    Load dict of number of samples per array per cohort
    """

    samples_df = pd.read_csv(samples_tsv, sep='\t')
    samples_df.index = samples_df['array']
    samples_df.drop('array', axis=1, inplace=True)

    counts = {}

    with open(cohorts) as fin:
        for meta, cohorts in csv.reader(fin, delimiter='\t'):
            counts[meta] = {}

            for cohort in cohorts.split(';'):
                if cohort not in samples_df.columns:
                    continue

                for array in samples_df.index[samples_df[cohort] > 0]:
                    if array not in elig_arrays:
                        continue

                    n = samples_df.loc[array, cohort]
                    if array not in counts[meta]:
                        counts[meta][array] = n
                    else:
                        counts[meta][array] += n

    return counts


def get_passing_fracs(array_counts, array_labels_df, keep_n_columns):
    """
    Compute fraction of passing samples per cohort for all intervals
    """

    cohorts = array_counts.keys()
    cohort_totals = {c : sum(v.values()) for c, v in array_counts.items()}

    fracs_df = array_labels_df.iloc[:, :keep_n_columns]

    for cohort in cohorts:
        npass = pd.Series(np.zeros(len(fracs_df)), index=fracs_df.index)
        for array in array_counts[cohort].keys():
            n_samps = array_counts[cohort][array]
            npass[array_labels_df[array]] += n_samps
        fracs_df[cohort] = npass / cohort_totals[cohort]

    return float_cleanup(fracs_df, 6, keep_n_columns)


def label_cohort_fails(cohort_fracs_df, min_frac, keep_n_columns):
    """
    Label all intervals with failing cohorts, if any
    """

    cohorts = list(cohort_fracs_df.columns[keep_n_columns:])

    fails_df = cohort_fracs_df.iloc[:, :keep_n_columns].copy()
    fails_df['exclude_cohorts'] = ''

    def _append_cohort(x, cohort):
        if x == '':
            return cohort
        else:
            return ';'.join([x, cohort])

    for cohort in cohorts:
        fail_labs = (cohort_fracs_df[cohort] < min_frac)
        fails_df.loc[fail_labs, 'exclude_cohorts'] = \
            fails_df.loc[fail_labs, 'exclude_cohorts'].\
                apply(_append_cohort, cohort=cohort)

    return fails_df


def main():
    """
    Main block
    """

    # Parse arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bed', help='BED file of intervals to be evaluated')
    parser.add_argument('probes_tsv', help='.tsv of path to probeset BED and ' +
                        'array name')
    parser.add_argument('samples_tsv', help='.tsv matrix of sample counts for ' +
                        'arrays (rows) X cohorts (columns)')
    parser.add_argument('cohorts', help='.tsv of metacohort assignments')
    parser.add_argument('-o', '--outfile', help='Output BED file annotated with ' +
                        'exclusion criteria [default: stdout]', default='stdout')
    parser.add_argument('--probecounts-outfile', help='Output BED file annotated with ' +
                        'number of probes per interval [optional]')
    parser.add_argument('--frac-pass-outfile', help='Output BED file annotated with ' +
                        'fraction of passing samples per cohort per per interval [optional]')
    parser.add_argument('--min-interval-size', dest='min_size', type=int, 
                        default=100000, help='Uniformly expand small intervals to ' +
                        'be at least this size [default: 100,000]')
    parser.add_argument('--min-probes', type=int, default=10, help='Minimum number ' +
                        'of probes required per interval to pass [default: 10]')
    parser.add_argument('--min-frac-samples', type=float, default=0.9, 
                        help='Minimum fraction of samples required per interval ' +
                        'to pass [default: 0.9]')
    parser.add_argument('-k', '--keep-n-columns', type=int, default=3, 
                        help='Number of columns from input BED to keep in ' +
                        '--outfile [default: 3]')
    parser.add_argument('-z', '--bgzip', action='store_true', default=False,
                        help='compress --outfile and --probecounts-outfile with ' +
                        'bgzip')
    args = parser.parse_args()


    # Step 1. Load intervals, and expand (as necessary)
    intervals, header = load_intervals(args.bed, args.min_size, args.keep_n_columns)

    # Step 2. Annotate intervals with probe counts using athena
    tracks, tnames = parse_probesets(args.probes_tsv)
    intervals = annotate_bins(bins=intervals.fn, chroms=None, ranges=None, 
                              tracks=tracks, ucsc_tracks=[], ucsc_ref=None, 
                              actions=['count' for i in range(len(tracks))], 
                              fasta=None, snv_mus=None, maxfloat=8, 
                              ucsc_chromsplit=False, quiet=False)
    intervals = replace_coords(intervals, args.bed, args.keep_n_columns)
    counts_outfile = args.probecounts_outfile
    if counts_outfile is not None:
        if 'compressed' in determine_filetype(counts_outfile):
            counts_outfile = path.splitext(counts_outfile)[0]
        intervals.saveas(counts_outfile, trackline='\t'.join(header + tnames))
        if args.bgzip:
            bgzip(counts_outfile)

    # Step 3. Determine pass/fail labels per interval per array
    array_labels_df = label_array_fails(intervals, args.min_probes, header, tnames)

    # Step 4. Compute fraction of passing samples per interval per cohort
    array_counts = load_array_counts(args.samples_tsv, args.cohorts, tnames)
    cohort_fracs_df = get_passing_fracs(array_counts, array_labels_df, 
                                        args.keep_n_columns)
    fracs_outfile = args.frac_pass_outfile
    if fracs_outfile is not None:
        if 'compressed' in determine_filetype(fracs_outfile):
            fracs_outfile = path.splitext(fracs_outfile)[0]
        cohort_fracs_df.to_csv(fracs_outfile, sep='\t', index=False)
        if args.bgzip:
            bgzip(fracs_outfile)

    # Step 5. Label each interval with cohorts to be excluded
    cohort_labels_df = label_cohort_fails(cohort_fracs_df, args.min_frac_samples,
                                          args.keep_n_columns)

    # Step 6. Format output file and write out
    if args.outfile in '- stdout /dev/stdout'.split():
        cohort_labels_df.to_csv(stdout, sep='\t', index=False)
    else:
        outfile = args.outfile
        if 'compressed' in determine_filetype(outfile):
            outfile = path.splitext(outfile)[0]
        cohort_labels_df.to_csv(outfile, sep='\t', index=False)
        if args.bgzip:
            bgzip(outfile)


if __name__ == '__main__':
    main()

