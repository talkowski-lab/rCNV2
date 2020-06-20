#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Clusters elements from input annotation tracks into cis-regulatory blocks (CRBs)
"""


from pysam import TabixFile
import pandas as pd
from io import StringIO
import pybedtools as pbt
import numpy as np
from sklearn.cluster import DBSCAN
import argparse
from os import path
import subprocess


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
            xbt = pbt.BedTool('', from_string=True)
        elif n_bl == 1:
            xbt = pbt.BedTool(xlist[0])
        else:
            xbt = pbt.BedTool(xlist[0]).\
                      cat(*[pbt.BedTool(bl) for bl in xlist[1:]])

        xbt = xbt.sort().merge().saveas()

    return xbt


def load_tracks(tracklist, chrom, genome=None):
    """
    Loads all tracks for a single chromsome as pbt.BedTool
    Suffixes each element with unique ID
    """

    def _load_and_suffix_bt(btpath, chrom):
        invals = StringIO('\n'.join([x for x in TabixFile(btpath).fetch(chrom)]))
        tdf = pd.read_csv(invals, sep='\t', names='chrom start end name'.split())
        tdf['name'] = ['{}.{}_e{}'.format(n, c, i) for n, c, i \
                       in zip(tdf.name, tdf.chrom, tdf.index + 1)]
        return pbt.BedTool.from_dataframe(tdf)

    tracks = []
    for t in tracklist:
        try:
            tracks.append(_load_and_suffix_bt(t, chrom))
        except:
            msg = 'WARNING: zero elements loaded for chromosome {} from track {}.\n' + \
                  'This may indicate an error with this script, or with data preprocessing.'
            print(msg.format(chrom, t))

    if len(tracks) == 1:
        all_bt = tracks[0]
    else:
        all_bt = tracks[0].cat(*tracks[1:], postmerge=False)

    if genome is None:
        return all_bt.sort().saveas()
    else:
        return all_bt.sort(g=genome).saveas()


def make_clusters(ebt, pos_df, min_elements=1, neighborhood_dist=10000, genome=None):
    """
    Clusters elements based on 1D position
    Returns:
        1. pbt.BedTool of maximal coordinates per cluster
        2. dict of {cluster_id : [element_ids]}
    """

    edf = ebt.to_dataframe()

    # Clusters elements with DBSCAN
    clustering = DBSCAN(eps=neighborhood_dist, min_samples=min_elements, 
                        metric='euclidean').fit(pos_df)
    pos_df['cluster'] = clustering.labels_
    cluster_members = {k : pos_df.index[pos_df.cluster == k].tolist() \
                       for k in list(set(clustering.labels_.tolist())) \
                       if k != -1}

    # Gets maximal coordinates for each cluster
    clust_bt_str = ''
    for k, eids in cluster_members.items():
        kchrom = edf.chrom[edf.name.isin(eids)].iloc[0]
        kstart = np.nanmin(edf.start[edf.name.isin(eids)])
        kend = np.nanmax(edf.end[edf.name.isin(eids)])
        clust_bt_str += '\t'.join(str(x) for x in [kchrom, kstart, kend, k]) + '\n'
    clust_bt = pbt.BedTool(clust_bt_str, from_string=True).sort(g=genome)

    return clust_bt, cluster_members


def refine_clusters(clust_bt, clust_members, ebt, blacklist, xcov=0.5, genome=None, 
                    min_crb_separation=10000, prefix='CRB'):
    """
    Refine & reformat final clusters & their constituent elements
    """

    crb_bt_str = ''
    crb_ele_bt_str = ''

    clust_df = clust_bt.to_dataframe()
    edf = ebt.to_dataframe()

    # Aggregate clusters within min_crb_separation
    clust_groups = clust_bt.merge(d=min_crb_separation, c=4, o='distinct')

    # Blacklist final CRB clusters
    clust_groups.coverage(blacklist).\
                 filter(lambda x: float(x[-1]) < xcov).\
                 cut(range(4))

    # Iterate over cluster groups and reformat CRB & elements from each
    k = 0
    for g in clust_groups:
        k += 1

        # Gather final CRB info
        crb_id = '_'.join([prefix, g.chrom, str(k)])
        cidxs = [int(x) for x in g[3].split(',')]
        cstart = str(np.nanmin(clust_df.start[clust_df.name.isin(cidxs)]))
        cend = str(np.nanmax(clust_df.end[clust_df.name.isin(cidxs)]))

        # Gather elements belonging to CRB and format them
        elists = [clust_members[x] for x in cidxs]
        crb_eles = [i for s in elists for i in s]
        n_ele = str(len(crb_eles))
        eles_df = edf[edf.name.isin(crb_eles)]
        for edat in eles_df.values.tolist():
            ename = '.'.join(edat[3].split('.')[:-1])
            crb_ele_bt_str += '\t'.join([str(x) for x in edat[:-1] + [ename, crb_id]]) + '\n'

        # Format final CRB info & write to file
        crb_bt_str += '\t'.join([g.chrom, cstart, cend, crb_id, n_ele]) + '\n'


    crb_bt = pbt.BedTool(crb_bt_str, from_string=True).sort(g=genome)
    crb_ele_bt = pbt.BedTool(crb_ele_bt_str, from_string=True).sort(g=genome)

    return crb_bt, crb_ele_bt


def cluster_chrom(tracklist, chrom, genome, blacklist, xcov=0.5, min_elements=None,
                  n_ele_prop=0.1, neighborhood_dist=10000, min_crb_separation=10000, 
                  prefix='CRB'):
    """
    Load & cluster all elements for a single chromosome
    Returns:
        1. pbt.BedTool of final clusters
        2. pbt.BedTool of all elements belonging to final clusters
    """

    # Count number of input tracks & determine proportional min_elements
    n_tracks = len(tracklist)
    if min_elements is None:
        min_elements = np.nanmax([np.floor(n_tracks * n_ele_prop), 1])

    # Load elements for chromosome
    ebt = load_tracks(tracklist, chrom, genome)

    # Build pd.DataFrame of positions
    pos_df = pd.DataFrame([(x.start + x.end)/2 for x in ebt],
                          columns=['pos'], 
                          index=[x.name for x in ebt])

    # Initial clustering of CRBs
    clust_bt, clust_members = make_clusters(ebt, pos_df, min_elements, 
                                            neighborhood_dist, genome)

    # Refine & annotate clusters
    crb_bt, crb_ele_bt = refine_clusters(clust_bt, clust_members, ebt, blacklist, 
                                         xcov, genome, min_crb_separation, prefix)

    return crb_bt, crb_ele_bt


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
                        'at least one track.')
    parser.add_argument('-g', '--genome', required=True, help='BEDTools-style ' +
                        'genome file to be used when sorting.')
    parser.add_argument('-x', '--blacklist', nargs='*', action='append',
                        help='Blacklist BED files to exclude CRBs. May be ' +
                        'specified multiple times.')
    parser.add_argument('--blacklist-cov', default=0.5, type=float, 
                        help='Minimum fraction of CRB that must be covered ' +
                        'by any blacklist before being excluded.')
    parser.add_argument('--min-elements', default=None, type=int, help='Minimum ' +
                        'number of elements in CRB.')
    parser.add_argument('--prop-min-elements', default=0.1, type=float, help='If ' +
                        '--min-elements is not supplied, will automatically set ' +
                        'minimum number of elements as a proportion of the total ' +
                        'number of tracks.')
    parser.add_argument('--neighborhood-dist', default=10000, type=int, help='Maximum ' +
                        'distance between two elements to allow.')
    parser.add_argument('--min-crb-separation', default=10000, type=int, help='Minimum ' +
                        'distance between two CRBs before merging them.')
    parser.add_argument('-p', '--crb-prefix', default='CRB', help='Prefix for ' +
                        'CRB names.')
    parser.add_argument('--crb-outbed', required=True, help='Output BED for final CRBs.')
    parser.add_argument('--element-outbed', required=True, help='Output BED for ' +
                        'elements belonging to final CRBs.')
    parser.add_argument('-z', '--bgzip', action='store_true', help='Compress ' +
                        'output BEDs with bgzip.')
    args = parser.parse_args()

    # Clean suffixes for outfiles
    outpaths = {'crbs' : args.crb_outbed, 'elements' : args.element_outbed}
    for key, outpath in outpaths.items():
        if path.splitext(outpath)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outpaths[key] = path.splitext(outpath)[0]
            bgzip_out = True
        else:
            bgzip_out = args.gzip

    # Build list of contigs to consider
    contigs = [x.split('\t')[0] for x in open(args.genome).readlines()]
    eligible_contigs = [str(x) for x in range(1, 23)]
    contigs = [x for x in contigs if x in eligible_contigs]

    # Load blacklists
    blacklist = build_blacklist(args.blacklist)

    # Process elements from each chromosome in serial
    final_elements = pbt.BedTool('', from_string=True)
    final_crbs = pbt.BedTool('', from_string=True)
    for chrom in contigs:
        print('Clustering CRBs on chromosome {}...'.format(chrom))
        new_crbs, new_elements = \
            cluster_chrom(args.tracks, chrom, args.genome, blacklist, 
                          args.blacklist_cov, args.min_elements, 
                          args.prop_min_elements, args.neighborhood_dist, 
                          args.min_crb_separation, args.crb_prefix)
        final_crbs = final_crbs.cat(new_crbs, postmerge=False)
        final_elements = final_elements.cat(new_elements, postmerge=False)

    # Write final tracks to file
    final_crbs.saveas(outpaths['crbs'], 
                      trackline='\t'.join('#chrom start end crb_id n_elements'.split()))
    final_elements.saveas(outpaths['elements'],
                          trackline='\t'.join('#chrom start end element_id crb_id'.split()))

    # Bgzip, if optioned
    if bgzip_out:
        for outpath in outpaths.values():
            subprocess.run(['bgzip', '-f', outpath])


if __name__ == '__main__':
    main()
