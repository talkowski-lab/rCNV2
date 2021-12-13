#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Cluster final segments prior to plotting (for clarity)
"""


import pandas as pd
import re
import pybedtools as pbt
import networkx as nx
import argparse
from sys import stdout


def load_credints(loci_in, sig=None, pad=100000):
    """
    Load segment credible intervals as pbt.BedTool
    """

    # Load loci table as pd.DataFrame
    loci_df = pd.read_csv(loci_in, sep='\t')
    loci_df.index = loci_df.region_id

    # Subset based on significance, if optioned
    if sig is not None:
        loci_df = loci_df[loci_df.best_sig_level == sig]

    # Parse each row to build pbt.BedTool of all intervals
    ci_bt_str = ''
    for rid, coords in loci_df['cred_interval_coords'].to_dict().items():
        for cstr in coords.split(';'):
            parts = re.split('[:|-]', cstr)
            start = max([0, int(parts[1]) - pad])
            end = int(parts[2]) + pad
            ci_bt_str += '\t'.join([parts[0], str(start), str(end), rid]) + '\n'

    # Return pbt.BedTool of credible intervals
    return pbt.BedTool(ci_bt_str, from_string=True)


def cluster_credints(cred_bt):
    """
    Cluster credible intervals based on simple overlap
    Annotate each edge with the number of base pairs overlapping
    """

    G = nx.Graph()

    for hit in cred_bt.intersect(cred_bt, wo=True):

        rid1, rid2, ovr = hit[3], hit[7], hit[-1]

        for rid in [rid1, rid2]:
            if rid not in G.nodes():
                G.add_node(rid)
        
        if rid1 == rid2:
            continue

        G.add_edge(rid1, rid2, ovr=int(ovr))

    return G


def prune_clusters(G):
    """
    Prune a cluster graph such that each segment has no more than one edge for each CNV type
    Determine best edge based on greatest overlap
    """

    g_dat = nx.to_dict_of_dicts(G)

    for node in G.nodes():

        for cnv in 'DEL DUP'.split():
        
            mates = [m for m in g_dat[node] if cnv in m]

            if len(mates) <= 1:
                continue

            m_dat = tuple({m : g_dat[node][m]['ovr'] for m in mates}.items())
            best_mate = sorted(m_dat, key=lambda x: x[1], reverse=True)[0][0]

            for mate in mates:
                # Need to check for existing edges that may have been removed earlier in loop
                if mate != best_mate and G.has_edge(node, mate):
                    G.remove_edge(node, mate)

    return G


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('loci', help='BED file of final segments.')
    parser.add_argument('-s', '--significance', help='Restrict to specified ' +
                        'significance level.')
    parser.add_argument('-o', '--outfile', default='stdout', help='Output .txt ' +
                        'file of all locus clusters.')
    args = parser.parse_args()

    # Open connections to output files
    if args.outfile in 'stdout - /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Load credible intervals
    cred_bt = load_credints(args.loci, args.significance)

    # Cluster credible intervals into graph
    clust_G = cluster_credints(cred_bt)

    # Prune cluster graph
    pruned_G = prune_clusters(clust_G)

    # Write clusters to outfile
    for subg in nx.connected_components(pruned_G):
        outfile.write(','.join(subg) + '\n')
    outfile.close()


if __name__ == '__main__':
    main()
