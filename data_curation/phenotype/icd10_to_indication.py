#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Convert sample + ICD10 pairs into sample + string indication pairs
"""


import csv
import networkx as nx
from re import sub
import argparse
from sys import stdout


def get_children(G, node_id):
    """
    Return all generations of children for a node
    """

    children = list(set([node for node in nx.dfs_preorder_nodes(G, node_id)]))

    return children


def read_icd_list(infile, G, conv, children=True):
    """
    Read a list of ICD-10 codes, and gather all children
    """

    with open(infile) as f:
        terms = f.readlines()
        terms = [x.rstrip() for x in terms]

    children = [get_children(G, conv[t]) for t in terms]
    children = [item for sublist in children for item in sublist]
    children = [G.nodes[t]['icd'] for t in children]

    allterms = list(set(terms + children))

    return allterms


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


def read_icd_indication_dict(icd10map):
    """
    Read UKBB official ICD-10 map into dict of indications
    """

    icd_dict = {}
    junk_words = ['OTHER', 'UNSPECIFIED', 'UNKNOWN']
    split_words = ['AND', 'WITH']
    truncate_words = ['WITHOUT']
    junk_phrases = ['PERSONAL_HISTORY_OF_', 'SITE_NOT_SPECIFIED', 'NOT_SPECIFIED', \
                    'NOT_ELSEWHERE_CLASSIFIED']

    with open(icd10map) as infile:

        reader = csv.reader(infile, delimiter='\t')

        for icd, descrip, node_idx, parent_idx, selectable in reader:

            if icd == 'coding':
                continue

            # Clean description
            if descrip.startswith('Chapter'):
                descrip = descrip.upper().split(' ', 2)[2]
            else:
                descrip = descrip.upper().split(' ', 1)[1]
            descrip = sub('[^A-Za-z0-9\ ]+', '', descrip)
            descrip = ' '.join([w for w in descrip.split(' ') if w not in junk_words])
            for word in split_words:
                descrip = descrip.replace(' ' + word + ' ', ';')
            for word in truncate_words:
                descrip = descrip.split(word)[0].rstrip()
            descrip = sub('[\ ]+', '_', descrip)
            for phrase in junk_phrases:
                descrip = descrip.replace(phrase, '')
            
            if icd not in icd_dict.keys():
                icd_dict[icd] = descrip
            else:
                icd_dict = ';'.join(icd_dict[icd], descrip)

    return icd_dict


def get_parents(G, node_id):
    """
    Return all generations of parents for a node
    """

    parents = list(set(nx.ancestors(G, node_id)))

    return parents


def process_samples(samples, icd_indic_dict, G, icd_node_dict, 
                    whitelist, blacklist, nclist, blacklist_samps,
                    outfile, blacklisted_samps_outfile, 
                    default='NA', report=False):
    """
    Iterate over a file of sample + ICD-10 pairs and convert to indications,
    while applying various whitelists and blacklists
    """

    with open(samples) as infile:

        reader = csv.reader(infile, delimiter='\t')

        for sample, icds in reader:

            # Parse ICD-10 codes
            if icds == '':
                icds = []
            else:
                icds = icds.split(';')

            # Get all parent terms corresponding to original ICD-10 codes
            parents = [get_parents(G, icd_node_dict[i]) for i in icds]
            parents = [item for sublist in parents for item in sublist]
            parents = [G.nodes[t].get('icd', None) for t in parents]
            parents = [t for t in parents if t is not None]

            allterms = list(set(icds + parents))

            # Check for any overlap with sample-level blacklist terms
            # If any hits are found, immediately stop processing sample
            # and record that sample on sample-level blacklist outfile
            sbl_hits = list(set(allterms).intersection(set(blacklist_samps)))
            if len(sbl_hits) > 0:
                reason = 'Matched sample blacklist terms: {0}'
                outline = '{0}\t{1}\n'
                outline = outline.format(sample, 
                                         reason.format(', '.join(sbl_hits)))
                blacklisted_samps_outfile.write(outline)
                continue

            # Restrict to terms on term whitelist
            allterms = [t for t in allterms if t in whitelist]

            # Erase terms that are on term blacklist
            allterms = [t for t in allterms if t not in blacklist]

            # Convert all remaining terms into indications
            if len(allterms) > 0:
                allphenos = [icd_indic_dict.get(t, None) for t in allterms]
                allphenos = [p for p in allphenos if p is not None]
                pheno = ';'.join(allphenos)
            else:
                pheno = default

            outfile.write('{0}\t{1}\n'.format(sample, pheno))


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
    parser.add_argument('--whitelist-terms', help='List of ICD-10 codes to include' +
                        '. Will include the supplied terms and all of their ' + 
                        'descendant terms. [default: include all terms, unless ' + 
                        'blacklisted]', metavar='file')
    parser.add_argument('--blacklist-terms', help='List of ICD-10 codes to exclude' +
                        '. Will exclude the supplied terms and all of their ' + 
                        'descendant terms. [default: exclude no terms]', 
                        metavar='file')
    parser.add_argument('--not-control-terms', help='List of ICD-10 codes to use' +
                        ' when determining which samples represent unaffected ' + 
                        'controls. Any samples without one or more --whitelist-terms ' +
                        'with one or more of these terms will be added to ' +
                        '--blacklisted-samples-outfile. [default: any samples without ' +
                        'any --whitelist-terms are treated as controls]', 
                        metavar='file')
    parser.add_argument('--blacklist-samples', help='List of ICD-10 codes to use' +
                        ' when marking samples for exclusion based on phenotypes. ' + 
                        'Will exclude samples with any of these ICD-10 codes or ' + 
                        'their descendant terms, irrespective of other codes. ' +
                        '[default: exclude no samples]', metavar='file')
    parser.add_argument('--blacklisted-samples-outfile', help='Path to blacklisted ' + 
                        'samples. Only used if --blacklist-samples is specified. ' + 
                        '[default: icd10_sample_blacklist.txt]', metavar='file',
                        default='icd10_sample_blacklist.txt')
    parser.add_argument('-d', '--default', help='String to use if no ICD-10 codes ' +
                        'match. [default: healthy]', default='healthy', 
                        metavar='string')
    parser.add_argument('--report-fails', action='store_true', 
                        help='Print warnings for each failed ICD-10 mapping. ' +
                        '[default: silent]')
    parser.add_argument('-o', '--outfile', help='Path to outfile. ' +
                        '[default: stdout]', metavar='file')
    args = parser.parse_args()

    # Read ICD-10 indication dictionary
    icd_indic_dict = read_icd_indication_dict(args.icd_10_map)

    # Read ICD-10 graph
    icd_g, icd_node_dict = make_icd_graph(args.icd_10_map)

    # Build whitelist terms (+ children)
    if args.whitelist_terms is not None:
        whitelist = read_icd_list(args.whitelist_terms, icd_g, icd_node_dict)
    else:
        whitelist = list(icd_indic_dict.keys())

    # Build blacklist terms (+ children)
    if args.blacklist_terms is not None:
        blacklist = read_icd_list(args.blacklist_terms, icd_g, icd_node_dict)
    else:
        blacklist = []

    # Build non-control terms (+ children)
    if args.not_control_terms is not None:
        nclist = read_icd_list(args.not_control_terms, icd_g, icd_node_dict)
    else:
        nclist = []

    # Build sample blacklist terms (+ children)
    if args.blacklist_samples is not None:
        blacklist_samps = read_icd_list(args.blacklist_samples, 
                                        icd_g, icd_node_dict)
    else:
        blacklist_samps = []

    # Open connection to outfiles
    if args.outfile is None:
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    blacklisted_samps_outfile = open(args.blacklisted_samples_outfile, 'w')
    blacklisted_samps_outfile.write('#sample\treason\n')

    # Process sample + ICD-10 pairings
    process_samples(args.samples, icd_indic_dict, icd_g, icd_node_dict, 
                    whitelist, blacklist, nclist, blacklist_samps, 
                    outfile, blacklisted_samps_outfile, 
                    args.default, args.report_fails)
    

if __name__ == '__main__':
    main()
