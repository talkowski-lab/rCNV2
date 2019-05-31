#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Map a list of phenotype indications to HPO codes

For each row in the input file, all indications will be compared against the 
specified HPO database. All matching HPO terms and their related superterms 
will be returned.
"""


import re
import networkx
import obonet
import csv
import argparse
from sys import stdout


def clean_pheno(string):
    """
    Clean a phenotype string
    """

    keywords_to_remove = ['POSSIBLE_']

    for kw in keywords_to_remove:
        string = string.replace(kw, '')

    string = string.lower()
    string = re.sub('_', ' ', string)
    string = re.sub('[^a-zA-Z0-9 \n\.]', '', string)
    string = string.strip()

    return string


def read_hpo_break_table(break_terms):
    """
    Read --break-terms table and convert to dict
    """

    break_table = {}

    with open(break_terms) as infile:
        reader = csv.reader(infile, delimiter='\t')

        for key, term in reader:
            key = clean_pheno(key)
            if key not in break_table.keys():
                break_table[key] = [term]
            else:
                break_table[key].append(term)

    return break_table


def get_superterms(term, hpo_g):
    """
    Get all HPO terms upstream of a query term
    """

    if hpo_g.nodes.get(term, None) is not None:
        superterms = networkx.descendants(hpo_g, term)
        superterms = [s for s in superterms if s is not None]
    else:
        superterms = []

    return superterms


def make_hpo_lookup_table(hpo_g, supp_terms=None, break_table=None, 
                          eligible_terms=None):
    """
    Convert HPO network graph into simple indication lookup dictionary
    """

    hpo_dict = {}

    # Parse all nodes from HPO graph
    for term, data in hpo_g.nodes(data=True):

        keys = [data.get('name')]

        if data.get('synonym') is not None:
            syns = [s.split('"')[1] for s in data.get('synonym')]
            keys = set(keys + syns)

        keys = [clean_pheno(k) for k in set(keys) if k is not None]

        for key in keys:

            if term in break_table.get(key, []):
                continue
            
            if key not in hpo_dict.keys():
                hpo_dict[key] = [term]
            else:
                if term not in hpo_dict[key]:
                    hpo_dict[key].append(term)

            for sterm in get_superterms(term, hpo_g):
                if sterm not in hpo_dict[key]:
                    hpo_dict[key].append(sterm)

    # Add supplementary terms, as desired
    if supp_terms is not None:
        with open(supp_terms) as supp:
            reader = csv.reader(supp, delimiter='\t')

            for key, terms in reader:
                key = clean_pheno(key)
                terms = list(set(terms.split(';')))

                if key not in hpo_dict.keys():
                    hpo_dict[key] = terms
                else:
                    hpo_dict[key] = hpo_dict[key] + terms

                for term in terms:
                    for sterm in get_superterms(term, hpo_g):
                        if sterm not in hpo_dict[key]:
                            hpo_dict[key].append(sterm)


    # Restrict terms to those appearing in eligibility list, if optioned
    # Remove entries that don't apply to any of the designated eligible terms
    if eligible_terms is not None:
        keys = list(hpo_dict.keys())
        for key in keys:
            terms = [t for t in hpo_dict.get(key) if t in eligible_terms]
            if len(terms) == 0:
                hpo_dict.pop(key)
            else:
                hpo_dict[key] = terms

    return hpo_dict


def parse_indications(indications, split_keys):
    """
    Split a string of indications into a list of terms to be queried vs HPO
    """

    split_ind = re.split('[;/]', indications)

    for split_key in split_keys:
        subsplit = [t.split(split_key) for t in split_ind]
        split_ind = [item for sublist in subsplit for item in sublist]

    split_ind = [s for s in split_ind if s is not '']

    return split_ind


def match_indication(indication, hpo_dict):
    """
    Match a single indication to HPO dictionary, and return a nonredundant list
    of all matching terms and superterms
    """

    query = clean_pheno(indication)

    raw_matches = [hp for key, hp in hpo_dict.items() if query == key]

    raw_matches = set([hp for sublist in raw_matches for hp in sublist])

    matches = list(raw_matches)

    return matches


def parse_phenotypes(phenos, outfile, hpo_dict, default):
    """
    Master function to parse a list of samples and phenotypes, and write 
    the converted HPO terms to output file
    """

    split_keys=['_AND_', '_WITH_', '_OR_', '__', '_PREVIOUS']

    with open(phenos) as infile:

        reader = csv.reader(infile, delimiter='\t')

        for sample, indications in reader:

            matches = []

            for indication in parse_indications(indications, split_keys):
                raw_matches = match_indication(indication, hpo_dict)
                if raw_matches is not None:
                    for match in raw_matches:
                        if match not in matches:
                            matches.append(match)

            # if eligible_terms is not None and len(matches) > 0:
            #     matches = [t for t in matches if t in eligible_terms]

            if len(matches) == 0:
                matches = [default]

            matches.sort()

            outline = '{0}\t{1}\n'.format(sample, ';'.join(matches))

            outfile.write(outline)


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('phenos', help='Test file of indications. One line per ' +
                        'patient. Two tab-delimited columns: patient ID and ' +
                        'string of semicolon-delimited indications.')
    parser.add_argument('--obo', help='Path to HPO .obo file. [default: ' + 
                        'http://purl.obolibrary.org/obo/hp.obo]',
                        default='http://purl.obolibrary.org/obo/hp.obo',
                        metavar='(file|url)')
    parser.add_argument('-s', '--supplementary-terms', help='File specifying ' + 
                        'supplementary terms to consider during HPO matching. ' +
                        'Two tab-delimited columns: keyword to match & HPO term.' +
                        ' [default: no supplementary terms]',
                        metavar='tsv')
    parser.add_argument('-x', '--break-terms', help='File specifying ' + 
                        'HPO-indication pairs to deliberately not consider during ' +
                        'HPO matching. Two tab-delimited columns: indication ' +
                        '& HPO term. [default: use native HPO mapping]',
                        metavar='tsv')
    parser.add_argument('-e', '--eligible-terms', help='List of HPO terms ' +
                        'to retain in final patient classification. [default: ' +
                        'keep all terms]',
                        metavar='tsv')
    parser.add_argument('--no-match-default', help='Default term to use for ' +
                        'samples with no matching terms. [default: "UNKNOWN"]',
                        default='UNKNOWN', dest='default')
    parser.add_argument('-o', '--outfile', help='Path to outfile. ' +
                        '[default: stdout]', metavar='file')
    
    # Process args & options
    args = parser.parse_args()

    if args.outfile is None:
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    
    hpo_g = obonet.read_obo(args.obo)

    if args.break_terms is not None:
        break_table = read_hpo_break_table(args.break_terms)
    else:
        break_table = None

    if args.eligible_terms is not None:
        eligible_terms = [t.strip() for t in open(args.eligible_terms).readlines()]
    else:
        eligible_terms = None

    # Create HPO lookup table
    hpo_dict = make_hpo_lookup_table(hpo_g, 
                                     supp_terms=args.supplementary_terms, 
                                     break_table=break_table,
                                     eligible_terms=eligible_terms)

    # Process phenotypes & write to outfile
    parse_phenotypes(args.phenos, outfile, hpo_dict, args.default)


if __name__ == '__main__':
    main()
