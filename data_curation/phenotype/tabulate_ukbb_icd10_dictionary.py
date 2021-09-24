#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Prepare master UKBB ICD-10 manifest for manual review
"""


import csv
import networkx as nx
import argparse
from sys import stdout


def make_icd_graph(icd_10_map):
    """
    Convert ICD-10 dictionary into directed graph
    """

    # Make directed graph
    G = nx.DiGraph()
    icd_to_node_dict = {}

    # Iterate over all terms in ICD-10 dictionary
    with open(icd_10_map) as infile:

        reader = csv.reader(infile, delimiter='\t')

        for icd, descrip, node_id, parent_id, selectable in reader:

            # Skip header line
            if icd.startswith('coding'):
                continue

            # Add node for term & parent
            G.add_edge(parent_id, node_id)

            # Add node info
            G.nodes[node_id]['icd'] = icd
            G.nodes[node_id]['descrip'] = descrip
            if selectable == 'Y':
                G.nodes[node_id]['selectable'] = True
            else:
                G.nodes[node_id]['selectable'] = False
            G.nodes[node_id]['samples'] = 0

            # Add icd to node_id conversion into helper dict
            icd_to_node_dict[icd] = node_id

    return G, icd_to_node_dict


def get_children(G, node_id):
    """
    Return all generations of children for a node
    """

    children = list(set([node for node in nx.dfs_preorder_nodes(G, node_id)]))

    return children


def get_parents(G, node_id):
    """
    Return all generations of parents for a node
    """

    parents = list(set(nx.ancestors(G, node_id)))

    return parents


def annotate_sample_counts(samples, G, icd_to_node_dict):
    """
    Annotate ICD-10 graph with count of samples per node
    """

    # Iterate over samples listed in args.samples
    with open(samples) as infile:

        reader = csv.reader(infile, delimiter='\t')

        for samp_id, icd10s in reader:
            if icd10s == '':
                if G.nodes['0'].get('samples', None) is None:
                    G.nodes['0']['samples'] = 1
                else:
                    G.nodes['0']['samples'] += 1
            else:
                icds = icd10s.split(';')
                nodes = [icd_to_node_dict.get(i, None) for i in icds]
                parents = [get_parents(G, n) for n in nodes]
                parents = [item for sublist in parents for item in sublist]

                all_nodes = list(set(nodes + parents))

                for node in all_nodes:
                    G.nodes[node]['samples'] += 1

    return G


def drop_small_nodes(G, minsamp):
    """
    Prune nodes from graph that are below minimum number of samples
    """

    for node in list(nx.dfs_postorder_nodes(G, '0')):

        if G.nodes[node].get('samples', 0) < minsamp:
            G.remove_node(node)

    return G


def write_final_table(G, outfile):
    """
    Write out formatted table of curated results
    """

    header = '#samples\tdescription\tICD\tnode_id\tparent_id\n'
    outfile.write(header)

    for node in list(nx.dfs_preorder_nodes(G, '0')):

        samples = G.nodes[node].get('samples', 0)
        descrip = G.nodes[node].get('descrip', 'NA')
        ICD = G.nodes[node].get('icd', 'NA')
        node_id = node
        parent_id = ';'.join(list(G.predecessors(node)))
        if parent_id == '':
            parent_id = 'NA'

        outline = f'{samples}\t{descrip}\t{ICD}\t{node_id}\t{parent_id}\n'

        outfile.write(outline)


def main():
    """
    Main block
    """

    # Parse arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('icd_10_map', help='Mapping of ICD-10 codes to ' + 
                        'descriptions. Follows format provided by UKBB.', 
                        metavar='file')
    parser.add_argument('samples', help='Test file of samples & ICD10 codes. ' + 
                        'One line per sample. Two tab-delimited columns: ' + 
                        'sample ID and string of semicolon-delimited ICD-10s.', 
                        metavar='file')
    parser.add_argument('-m', '--minsamp', default=100, type=int, help='Minimum count ' +
                        'of samples required per term to be included in output. ' +
                        '[default: 100]', metavar='int')
    parser.add_argument('-o', '--outfile', help='Path to outfile. ' +
                        '[default: stdout]', metavar='file')
    args = parser.parse_args()

    # Open connection to outfile
    if args.outfile is None:
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Build tree of ICD-10 codes & conversion dict
    icd_g, icd_to_node_dict = make_icd_graph(args.icd_10_map)

    # Add sample counts to ICD-10 graph
    icd_g = annotate_sample_counts(args.samples, icd_g, icd_to_node_dict)

    # Drop nodes smaller than specified minimum number of samples
    icd_g = drop_small_nodes(icd_g, args.minsamp)

    # Write out table of ICD-10 terms with counts of samples
    write_final_table(icd_g, outfile)


if __name__ == '__main__':
    main()
