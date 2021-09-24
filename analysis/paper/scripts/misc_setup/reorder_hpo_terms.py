#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Reorder HPO terms by hierarchy for rCNV plots
"""


import pandas as pd
import scipy.cluster.hierarchy as shc
import numpy as np
import argparse
from sys import stdout


def load_hpos(hpo_table):
    """
    Load & clean HPO table. Returns: pd.DataFrame
    """

    hpo_df = pd.read_csv(hpo_table, sep='\t').rename(columns={'#HPO_term' : 'HPO'})
    hpo_df.index = hpo_df.HPO
    
    # Drop controls, sort by sample size, and move UNKNOWN to last
    hpo_df = hpo_df.drop('HEALTHY_CONTROL', axis=0).drop('HPO', axis=1)
    hpo_df.sort_values(by='samples', axis=0, ascending=False, inplace=True)
    hpo_df = hpo_df.loc[hpo_df.index != 'UNKNOWN', :].\
                 append(hpo_df.loc[hpo_df.index == 'UNKNOWN', :])

    return hpo_df


def load_sample_overlaps(overlap_matrix):
    """
    Load & clean sample overlap matrix
    """

    ovr_df = pd.read_csv(overlap_matrix, sep='\t')
    ovr_df.index = ovr_df.HPO

    return ovr_df.drop('HPO', axis=1)


def reorder_subset(hpos, hpo_df, ovr_df):
    """
    Reorder a list of hpos based on:
    1. Hierarchical clustering of ovr_df
    2. Sample size from hpo_df
    """

    hpo_df_sub = hpo_df.loc[hpo_df.index.isin(hpos), :]
    ovr_df_sub = ovr_df.loc[ovr_df.index.isin(hpos), :]

    # Hierarchical clustering of hpos by ovr_df
    Z = shc.linkage(ovr_df_sub, method='average')

    # Iterate over hierarchical clustering results
    clusters = {}
    nhpos = len(hpos)
    k = len(hpos) - 1
    for idxA, idxB in Z[:, :2].tolist():
        k += 1
        # If both items in the next pair are original HPOs, sort them by sample
        # size and add them as a new cluster to clusters dict
        if idxA < nhpos and idxB < nhpos:
            pair_hpos = ovr_df_sub.index[[int(idxA), int(idxB)]].to_list()
            ordered_hpos = hpo_df_sub.index[hpo_df_sub.index.isin(pair_hpos)].tolist()

        # If one (or both) items in the next pair are clustered pairs, create a new 
        # cluster by concatenating them, with the items first being sorted by max 
        # sample size
        else:
            if idxA < nhpos:
                hposA = [ovr_df_sub.index[int(idxA)]]
            else:
                hposA = clusters[int(idxA)]['hpos']
            if idxB < nhpos:
                hposB = [ovr_df_sub.index[int(idxB)]]
            else:
                hposB = clusters[int(idxB)]['hpos']
            # Determine which item has larger max sample size
            if len([x for x in hposA if x != 'UNKNOWN']) > 0:
                nmaxA = np.nanmax(hpo_df_sub.samples[(hpo_df_sub.index.isin(hposA)) \
                                                     & (hpo_df_sub.index != 'UNKNOWN')])
            else:
                nmaxA = 0
            if len([x for x in hposA if x != 'UNKNOWN']) > 0:
                nmaxB = np.nanmax(hpo_df_sub.samples[(hpo_df_sub.index.isin(hposB)) \
                                                     & (hpo_df_sub.index != 'UNKNOWN')])
            else:
                nmaxB = 0
            if nmaxA >= nmaxB:
                ordered_hpos = hposA + hposB
            else:
                ordered_hpos = hposB + hposA

        if 'UNKNOWN' in ordered_hpos:
            ordered_hpos = [h for h in ordered_hpos if h != 'UNKNOWN'] + ['UNKNOWN']

        max_n = np.nanmax(hpo_df_sub.samples[(hpo_df_sub.index.isin(ordered_hpos)) \
                                             & (hpo_df_sub.index != 'UNKNOWN')])
        clusters[k] = {'hpos' : ordered_hpos, 'max_n' : max_n}

    # The last cluster contains the full set of reordered hpos
    return clusters[k]['hpos']


def get_children(hpo, hpo_df):
    """
    Returns a list of direct child HPOs for the input HPO
    """

    parent_tier = hpo_df.HPO_tier[hpo_df.index == hpo][0]
    children = hpo_df.index[(hpo_df.parent_terms.str.contains(hpo, na=False).values) \
                            & (hpo_df.HPO_tier == parent_tier + 1)]

    return list(children)


def get_parent(hpo, hpo_df):
    """
    Returns the direct parent HPO for the input HPO
    """

    child_tier = hpo_df.HPO_tier[hpo_df.index == hpo][0]
    parents_df = hpo_df.loc[(hpo_df.child_terms.str.contains(hpo, na=False).values), :]
    parents_df = parents_df.sort_values(by='HPO_tier', ascending=False)
    if len(parents_df) > 0:
        parent = parents_df.index.to_list()[0]
    else:
        parent = hpo_df.index[hpo_df.HPO_tier == 1].to_list()[0]

    return parent


def reorder_hpos(hpo_df, ovr_df):
    """
    Reorder phenotypes by hierarchical HPO relationship, sample size, and sample overlap
    """

    # Descend HPO tiers, in order, while reordering each tier using reorder_subset()
    tiers = {t : [] for t in range(np.nanmin(hpo_df.HPO_tier), np.nanmax(hpo_df.HPO_tier) + 1)}
    for t in tiers.keys():
        tier_hpos = list(hpo_df.index[hpo_df.HPO_tier == t])
        if len(tier_hpos) > 1:
            ordered_hpos = reorder_subset(tier_hpos, hpo_df, ovr_df)
        else:
            ordered_hpos = tier_hpos
        tiers[t] = ordered_hpos

    # Flatten tier dict
    flat_tiers = [h for s in [t for t in tiers.values()] for h in s]

    # Ascend tiers in reverse (from lowest to highest), rolling lower HPO terms up
    # into their parent terms at each tier while preserving order
    all_terms = {h : [h] for h in flat_tiers}
    for t in reversed(range(2, len(tiers) + 1)):
        tier_parent_map = {h : get_parent(h, hpo_df) for h in tiers[t]}
        for hpo, parent in tier_parent_map.items():
            all_terms[parent] += all_terms[hpo]

    # The top-tier term will contain the reordered list of all terms
    return all_terms[tiers[1][0]]


def main():
    """
    Main command-line block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('hpo_table', help='tsv of HPOs with parent/child mappings.')
    parser.add_argument('overlap_matrix', help='tsv of asymmetrical sample overlap ' +
                        'between HPOs.')
    parser.add_argument('-o', '--outfile', default='stdout', help='Output tsv of ' +
                        'reordered phenotypes. [default: stdout]')
    args = parser.parse_args()
    
    # Open connections to output files
    if args.outfile in 'stdout - /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Load HPOs, sorted by sample size
    hpo_df = load_hpos(args.hpo_table)

    # Load sample overlap matrix
    ovr_df = load_sample_overlaps(args.overlap_matrix)

    # Reorder HPOs
    reordered_hpos = reorder_hpos(hpo_df, ovr_df)

    # Write to outfile
    outfile.writelines([h + '\n' for h in reordered_hpos])
    outfile.close()


if __name__ == '__main__':
    main()

