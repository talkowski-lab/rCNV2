#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Collapse segdup pairs based on pairwise distance
"""


import csv
import numpy as np
import networkx as nx
import pybedtools as pbt
import argparse
from sys import stdout


def load_bedpe(pairs, distance, recip):
    """
    Load bedpe of segdup pairs and self-intersect
    """
    
    hits = pbt.BedTool(pairs).pairtopair(pairs, type='both', slop=distance)

    def _check_overlap(x, recip):
        bfmt = '{}\t{}\t{}\n'
        b1 = pbt.BedTool(bfmt.format(x[0], x[2], x[4]), from_string=True)
        b2 = pbt.BedTool(bfmt.format(x[7], x[9], x[11]), from_string=True)
        return len(b1.intersect(b2, r=True, f=recip)) > 0

    if recip > 0:
        hits = hits.filter(_check_overlap, recip=recip).saveas()

    return hits


def load_coords(fin):
    """
    Load inner min/max coordinates of segdup pairs from input BEDPE
    """

    coords = {}

    with open(fin) as bedpe:
        reader = csv.reader(bedpe, delimiter='\t')
        for chrA, startA, endA, chrB, startB, endB, name in reader:
            imin = np.nanmin([int(endA), int(endB)])
            imax = np.nanmax([int(startA), int(startB)])
            coords[name] = {'chrom' : chrA,
                            'min' : imin,
                            'max' : imax,
                            'size' : imax - imin}

    return coords


def make_graph(hits):
    """
    Create nx.Graph of intersected segdup pairs
    """

    G = nx.Graph()

    for pair in hits:
        id1 = pair[6]
        id2 = pair[13]
        if id1 not in G.nodes():
            G.add_node(id1)
        if id2 != id1:
            if id2 not in G.nodes():
                G.add_node(id2)
            G.add_edge(id1, id2)

    return G


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pairs', help='BEDPE of segdup pairs to consider.')
    parser.add_argument('-d', '--distance', type=int, default=0,
                        help='Slop distance for bedtools pairtopair.')
    parser.add_argument('-r', '--recip', type=float, default=0,
                        help='Reciprocal overlap for bedtools intersect.')
    parser.add_argument('-o', '--outfile', default='stdout',
                        help='Outfile [default: stdout].')
    args = parser.parse_args()

    # Open connection to outfile
    if args.outfile in '- stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Load bedpe and self-intersect
    hits = load_bedpe(args.pairs, args.distance, args.recip)

    # Load dictionary of coordinates
    coords = load_coords(args.pairs)

    # Make graph of intersected pairs
    hits_g = make_graph(hits)

    # Write one line to output BED file for each subgraph of clustered segdup pairs
    for subgraph in nx.connected_components(hits_g):
        # Retain segdup pair with smallest intervening distance among linked pairs
        smallest = sorted([(n, coords[n]['size']) for n in subgraph], key=lambda x: x[1])[0][0]
        chrom = coords[smallest]['chrom']
        newmin = coords[smallest]['min']
        newmax = coords[smallest]['max']
        # newmin = np.nanmax([coords[x]['min'] for x in subgraph])
        # newmax = np.nanmin([coords[x]['max'] for x in subgraph])
        outfile.write('\t'.join([chrom, str(newmin), str(newmax)]) + '\n')

    outfile.close()


if __name__ == '__main__':
    main()

