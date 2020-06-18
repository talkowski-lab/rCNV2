#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Clusters elements from input annotation tracks into cis-regulatory blocks (CRBs)
"""


import pandas as pd
import pybedtools as pbt
import argparse


def load_tracks(tracklist, genome=None):
    """
    Loads all tracks as pbt.BedTool and suffixes each with unique element ID
    """

    def _load_and_suffix_bt(btpath):
        tdf = pd.read_csv(btpath, sep='\t', names='chrom start end name'.split())
        tdf['name'] = ['{}_e{}'.format(n, i) for n, i in zip(tdf.name, tdf.index + 1)]
        return pbt.BedTool.from_dataframe(tdf)

    tracks = [_load_and_suffix_bt(t) for t in tracklist]

    if len(tracks) == 1:
        all_bt = tracks[0]
    else:
        all_bt = tracks[0].cat(*tracks[1:], postmerge=False)

    if genome is None:
        return all_bt.sort().saveas()
    else:
        return all_bt.sort(g=genome).saveas()


def build_dist_matrix(ebt):
    """
    Builds an all X all distance matrix from a pbt.BedTool of elements
    """

    for feature in ebt:
        fstr = '\t'.join([str(x) for x in [feature.chrom, feature.start, feature.end]]) + '\n'
        # pbt.BedTool(fstr, from_string=True).


def cluster_chrom(elements, chrom):
    """
    Cluster all elements for a single chromosome
    """

    # Filter to elements on chromosome
    ebt = elements.filter(lambda x: x.chrom == chrom)

    # Build distance matrix of all elements
    dmat = 


    import pdb; pdb.set_trace()


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tracks', nargs='+', help='Annotation ' +
                        'tracks to be included in CRB clustering. Must specify ' +
                        'at least one track')
    parser.add_argument('-g', '--genome', help='BEDTools-style genome file to ' +
                        'be used when sorting.')
    args = parser.parse_args()

    # Load all tracks
    elements = load_tracks(args.tracks, args.genome)

    # Build list of contigs to consider
    if args.genome is None:
        contigs = list(set([x.chrom for x in elements]))
    else:
        contigs = [x.split('\t')[0] for x in open(args.genome).readlines()]

    # Cluster each chromosome
    crbs_split = [cluster_chrom(elements, chrom) for chrom in contigs]


    import pdb; pdb.set_trace()


if __name__ == '__main__':
    main()
