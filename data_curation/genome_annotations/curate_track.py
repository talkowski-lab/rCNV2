#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Download and curate a single genome annotation track
"""


import pybedtools as pbt
import pandas as pd
import numpy as np
import argparse
from os import path
import subprocess


def download_track(trackpath):
    """
    Download track (or locate its local copy, if already downloaded)
    """

    if trackpath.startswith('gs://'):
        localpath = path.basename(trackpath)
        if not any([path.exists(localpath), path.isfile(localpath)]):
            subprocess.run(['gsutil', '-m', 'cp', trackpath, './'], check=True)

    elif any(trackpath.startswith(prefix) for prefix in 'http:// https:// ftp://'.split()):
        localpath = path.basename(trackpath)
        if not any([path.exists(localpath), path.isfile(localpath)]):
            subprocess.run(['wget', '--no-check-certificate', trackpath], check=True)

    else:
        localpath = trackpath

    return localpath


def get_track(trackpath, n_download_retries=5, prefix=None):
    """
    Download (if necessary) and load annotation track
    """

    # Download track, if necessary, trying a total of n_download_retries times
    for tries in range(1, n_download_retries+1):
        try:
            localpath = download_track(trackpath)
            print('\nSuccessfully localized {} after {} attempts\n'.format(trackpath, tries))
            break
        except:
            print('\nAttepmt {} to localize {} was unsuccessful. Retrying...\n'.format(tries, trackpath))
            continue

    # Set track name equal to base filename without bed suffix
    if '.bed' in localpath:
        trackname = path.basename(localpath).split('.bed')[0]
    elif '.txt' in localpath:
        trackname = path.basename(localpath).split('.txt')[0]
    elif '.tsv' in localpath:
        trackname = path.basename(localpath).split('.tsv')[0]

    if prefix is not None:
        if prefix != '':
            trackname = '.'.join([prefix, trackname])

    return pbt.BedTool(localpath), trackname


def build_blacklist(blacklists):
    """
    Build universal blacklist from one or more inputs
    """

    if blacklists is None:
        xbt = pbt.BedTool('', from_string=True)

    else:
        xlist = [x for s in blacklists for x in s]
        n_bl = len(xlist)
        if n_bl == 0:
            xbt = pbt.BedTool()
        elif n_bl == 1:
            xbt = pbt.BedTool(xlist[0])
        else:
            xbt = pbt.BedTool(xlist[0]).\
                      cat(*[pbt.BedTool(bl) for bl in xlist[1:]])

        xbt = xbt.sort().merge().saveas()

    return xbt


def clean_track(track, trackname, xbt, xcov=0.5, genome=None, 
                min_size=10, max_size=100000):
    """
    Clean & curate pbt.BedTool of elements
    """

    def _clean_contig(feature):
        newchrom = str(feature.chrom).replace('chr', '')
        feature.chrom = newchrom
        return feature

    def _is_autosomal(feature):
        return any(str(feature.chrom) == str(x) for x in range(1, 23))

    track = track.cut(range(3)).\
                  each(_clean_contig).\
                  filter(_is_autosomal).\
                  coverage(xbt).\
                  filter(lambda x: float(x[-1]) < xcov).\
                  cut(range(3)).\
                  saveas()

    def _size_filter(feature, min_size, max_size):
        return all([len(feature) >= min_size, len(feature) <= max_size])

    if genome is None:
        track_df = track.sort().merge().\
                         filter(_size_filter, min_size, max_size).\
                         saveas().to_dataframe()
    else:
        track_df = track.sort(g=genome).merge().\
                         filter(_size_filter, min_size, max_size).\
                         saveas().to_dataframe()

    track_df['name'] = trackname

    return pbt.BedTool.from_dataframe(track_df)


def get_track_stats(track, trackname, orig_trackpath, outfile):
    """
    Compute basic statistics for elements in track
    """

    # Add header to outfile
    hcols = '#trackname n_elements min_size median_size mean_size max_size total_bp original_path'
    outfile.write('\t'.join(hcols.split()) + '\n')

    # Compute stats
    sizes = np.array([len(x) for x in track], dtype=float)
    n_ele = len(track)
    min_size = np.nanmin(sizes)
    med_size = np.nanmedian(sizes)
    mean_size = np.nanmean(sizes)
    max_size = np.nanmax(sizes)
    total_bp = np.nansum(sizes)

    # Write stats to outfile
    svals = [n_ele, min_size, med_size, mean_size, max_size, total_bp]
    outfile.write('\t'.join([trackname] + [str(round(x, 1)) for x in svals] + [orig_trackpath]) + '\n')
    outfile.close()


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('path', help='Path to annotation track. Can be local, ' +
                        'gs://, http://, https://, or ftp://')
    parser.add_argument('-g', '--genome', help='BEDTools-style genome file to ' +
                        'be used when sorting.')
    parser.add_argument('-x', '--blacklist', nargs='*', action='append',
                        help='Blacklist BED files to exclude elements. May be ' +
                        'specified multiple times.')
    parser.add_argument('--blacklist-cov', default=0.5, type=float, 
                        help='Minimum fraction of element that must be covered ' +
                        'by any blacklist before being excluded.')
    parser.add_argument('--min-size', default=10, type=int, 
                        help='Minimum size of element to retain after merging.')
    parser.add_argument('--max-size', default=100000, type=int, 
                        help='Minimum size of element to retain after merging.')
    parser.add_argument('-o', '--outbed', help='Path to output BED file. ' +
                        '[default: trackname + .curated.bed.gz]')
    parser.add_argument('-s', '--stats', action='store_true', help='Compute ' +
                        'track statistics.')
    parser.add_argument('--n-download-retries', default=5, type=int, 
                        help='Number of times to retry downloading each track.')
    parser.add_argument('-p', '--prefix',  help='Prefix for track and element names.')
    parser.add_argument('-z', '--bgzip', action='store_true', help='Compress ' +
                        'output BED file with bgzip.')
    args = parser.parse_args()

    # Open connections to output files (unless not specified)
    if args.outbed is not None:
        if path.splitext(args.outbed)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outbed_path = path.splitext(args.outbed)[0]
        else:
            outbed_path = args.outbed

    # Load annotation track
    track, trackname = get_track(args.path, args.n_download_retries, args.prefix)

    # Open connections to unspecified output files
    if args.outbed is None:
        outbed_path = trackname + '.curated.bed'
        bgzip = True
    else:
        bgzip = args.bgzip
    if args.stats:
        stats_out = open(trackname + '.curated.stats.tsv', 'w')

    # Load blacklists
    xbt = build_blacklist(args.blacklist)

    # Clean annotation track
    track = clean_track(track, trackname, xbt, args.blacklist_cov, args.genome,
                        args.min_size, args.max_size)

    # Compute track stats, if optioned
    get_track_stats(track, trackname, args.path, stats_out)
    
    # Write cleaned track to file
    track.saveas(outbed_path)
    if bgzip:
        subprocess.run(['bgzip', '-f', outbed_path])


if __name__ == '__main__':
    main()

