#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Randomly shuffles segments in the genome
"""


import pandas as pd
import pybedtools as pbt
import argparse
from os.path import splitext
from sys import stdout
from athena.utils import bgzip


def load_segs(segs_in):
    """
    Load and reformat input segments to prepare for permutation
    """

    segs_df = pd.read_csv(segs_in,  delimiter='\t').iloc[:, :6]
    coord_colname = segs_df.columns.tolist()[5]

    def _normalize_intervals(row):
        start = int(row[1])
        end = int(row[2])
        intervals_str = row[-1]
        norm_intervals = []
        for x in intervals_str.split(';'):
            icoords = [int(n) for n in x.split(':')[1].split('-')]
            norm_intervals.append('-'.join([str(icoords[0] - start), str(icoords[1] - start)]))
        return ';'.join(norm_intervals)

    segs_df[coord_colname] = segs_df.apply(_normalize_intervals, axis=1)

    return pbt.BedTool().from_dataframe(segs_df), coord_colname


def custom_shuffle(segs, seed, genome, whitelist=None, blacklist=None, coords_colname='coords'):
    """
    Customized implementation of bedtools shuffle
    """

    def _seed_prefix(x, seed):
        x.name = '_'.join(['perm' + str(seed), x.name])
        return x

    # Shuffle DEL and DUP separately using +seed for DEL and -seed for DUP
    # This ensures non-overlapping intervals within CNV types (DEL vs DEL)
    # But allows for overlapping CNV intervalls between CNV types (DEL vs DUP)
    if whitelist is not None:
        shuf_del = segs.filter(lambda x: x[4] == 'DEL').\
                        shuffle(g=genome, incl=whitelist, excl=blacklist, seed=seed, f=1e-9, 
                                noOverlapping=True, allowBeyondChromEnd=False)
        shuf_dup = segs.filter(lambda x: x[4] == 'DUP').\
                        shuffle(g=genome, incl=whitelist, excl=blacklist, seed=-seed, f=1e-9, 
                                noOverlapping=True, allowBeyondChromEnd=False)
    else:
        shuf_del = segs.filter(lambda x: x[4] == 'DEL').\
                        shuffle(g=genome, seed=seed, f=1.0, noOverlapping=True, 
                                allowBeyondChromEnd=False)
        shuf_dup = segs.filter(lambda x: x[4] == 'DUP').\
                        shuffle(g=genome, seed=-seed, f=1.0, noOverlapping=True,
                                allowBeyondChromEnd=False)
    shuf = shuf_del.cat(shuf_dup, postmerge=False).sort(g=genome).\
                    each(_seed_prefix, seed=seed).saveas()
                        
    def _unnormalize_intervals(row):
        chrom = row.chrom
        start = row.start
        end = row.end
        intervals_str = row[-1]
        unnorm_intervals = []
        for x in intervals_str.split(';'):
            icoords = [int(n) for n in x.split('-')]
            newstart = str(start + icoords[0])
            newend = str(start + icoords[1])
            unnorm_intervals.append('{}:{}-{}'.format(chrom, newstart, newend))
        return ';'.join(unnorm_intervals)

    unnorm_intervals = [_unnormalize_intervals(x) for x in shuf]

    bt_colnames = 'chr start end region_id cnv'.split() + [coords_colname]
    shuf_bt = shuf.to_dataframe(names=bt_colnames)
    shuf_bt[coords_colname] = unnorm_intervals
    shuf_bt['perm'] = seed

    return shuf_bt


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('segs', help='BED of segments to be ' +
                        'shuffled. Requires at least six columns: chrom, start, ' +
                        'end, id, cnv type, and coords listing the exact ' +
                        'intervals to be shuffled in chr:start-end format.')
    parser.add_argument('-g', '--genome', required=True, help='BEDTools-style ' +
                        'genome file.')
    parser.add_argument('-w', '--whitelist', help='BED of regions to allow ' +
                        'during shuffling.')
    parser.add_argument('-b', '--blacklist', help='BED of regions to exclude ' +
                        'during shuffling.')
    parser.add_argument('-n', '--n-perms', type=int, default=1, help='Number of ' +
                        'permutations to perform.')
    parser.add_argument('-s', '--first-seed', type=int, default=1, help='Integer ' +
                        'seed passed to BEDTools shuffle for first permutation. ' +
                        'Will be incremented successively for each subsequent ' +
                        'permutation.')
    parser.add_argument('-o', '--outfile', help='Path to output BED file. [default: ' +
                        'stdout]', default='stdout')
    parser.add_argument('-z', '--bgzip', action='store_true', help='Compress ' + 
                        'output BED files with bgzip. [Default: do not compress]')
    parser.add_argument('-q', '--quiet', action='store_true', help='Suppress ' +
                        'verbose output.')
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

    # Load and reformat segments to prepare for permutation
    segs, coords_colname = load_segs(args.segs)

    # Shuffle intervals for each of --n-perms permutations
    for seed in range(args.first_seed, args.first_seed + args.n_perms):
        if not args.quiet:
            print('Starting permutation {}'.format(seed))
        shuffled_bt = custom_shuffle(segs, seed, args.genome, args.whitelist, 
                                     args.blacklist, coords_colname=coords_colname)
        if seed == args.first_seed:
            shuffled_bt.to_csv(outfile, sep='\t', na_rep='NA', header=True, 
                               index=False, mode='w')
        else:
            shuffled_bt.to_csv(outfile, sep='\t', na_rep='NA', header=False, 
                               index=False, mode='a')
    outfile.close()

    # Bgzip, if optioned
    if args.bgzip:
        bgzip(outfile_path)


if __name__ == '__main__':
    main()
