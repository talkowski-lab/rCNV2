#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Determine the minimal set of most specific HPO terms to retain for analysis 
based on a list of patients and corresponding HPO terms
"""

import networkx
import obonet
import csv
import argparse
from sys import stdout


def tally_hpo_counts(pheno_file, ignore_terms):
    """
    Create dictionary of counts of samples per HPO term in input file
    """

    hpo_counts = {}

    with open(pheno_file) as infile:
        reader = csv.reader(infile, delimiter='\t')

        for sample, terms in reader:
            for term in terms.split(';'):

                if ignore_terms is not None:
                    if term in ignore_terms:
                        continue

                if term in hpo_counts.keys():
                    hpo_counts[term] += 1
                else:
                    hpo_counts[term] = 1

    return hpo_counts


def write_raw_counts(hpo_counts, outfile):
    """
    Sort raw counts of samples per HPO term, and write to outfile
    """

    raw_out = open(outfile, 'w')

    sorted_keys = sorted(hpo_counts, key=hpo_counts.get, reverse=True)

    sorted_counts = [(i, hpo_counts[i]) for i in sorted_keys]

    for term, value in sorted_counts:
        raw_out.write('{0}\t{1}\n'.format(term, value))


def prune_small_terms(hpo_counts, min_samples):
    """
    Drop HPO terms with too few samples
    """

    filtered_counts = {}

    for term, count in hpo_counts.items():
        if count >= min_samples:
            filtered_counts[term] = count

    return filtered_counts


def get_related_terms(hpo_counts, hpo_g):
    """
    Define pairs of related terms to be evaluated
    """

    related_terms = []

    for term, count in hpo_counts.items():

        if term in hpo_g.nodes():
            parents = networkx.descendants(hpo_g, term)
            parents = [t for t in parents if t in hpo_counts.keys()]
            children = networkx.ancestors(hpo_g, term)
            children = [t for t in children if t in hpo_counts.keys()]
            relatives = list(set(parents + children))

            for rel in relatives:
                if (term, rel) not in related_terms \
                and (rel, term) not in related_terms:
                    related_terms.append((term, rel))

    return related_terms


def filter_related_term_list(hpo_counts, related_terms, min_diff):
    """
    Retain pairs of related terms with sample size differences < min_diff
    """

    filtered_pairs = []

    for termA, termB in related_terms:
        
        countA = hpo_counts.get(termA, 0)
        countB = hpo_counts.get(termB, 0)
        
        if abs(countA - countB) < min_diff:
            filtered_pairs.append((termA, termB))

    return filtered_pairs


def count_overlapping_samples(related_terms, pheno_file):
    """
    Count number of overlapping samples for pairs of HPO terms
    """

    pairwise_counts = {}
    eligible_terms = []

    for termA, termB in related_terms:
        pairwise_counts[(termA, termB)] = 0

        if termA not in eligible_terms:
            eligible_terms.append(termA)
        if termB not in eligible_terms:
            eligible_terms.append(termB)

    with open(pheno_file) as infile:
        reader = csv.reader(infile, delimiter='\t')

        for sample, terms in reader:
            filtered_terms = [t for t in terms.split(';') if t in eligible_terms]

            if len(filtered_terms) > 1:
                for i in range(0, len(filtered_terms) - 1):
                    for j in range(i + 1, len(filtered_terms)):
                        termA = filtered_terms[i]
                        termB = filtered_terms[j]
                        if (termA, termB) in pairwise_counts.keys():
                            pairwise_counts[(termA, termB)] += 1
                        elif (termB, termA) in pairwise_counts.keys():
                            pairwise_counts[(termB, termA)] += 1

    return pairwise_counts


def write_prune_log(logfile, termA, termB, nameA, nameB, reason):
    """
    Write HPO term pruning message to file
    """

    lmsg = 'Retained "{0}" ({2}) and pruned "{1}" ({3}): ' + reason + '\n'
    logfile.write(lmsg.format(nameA, nameB, termA, termB))


def prune_related_terms(hpo_counts, related_terms, pairwise_counts, hpo_g, 
                        min_diff, filter_log):
    """
    Drop related HPO terms based on sample overlap
    """

    terms_to_prune = []

    if filter_log is not None:
        f_log = open(filter_log, 'w')

    for termA, termB in related_terms:

        countA = hpo_counts.get(termA, 0)
        countB = hpo_counts.get(termB, 0)

        if (termA, termB) in pairwise_counts.keys():
            countAB = pairwise_counts.get((termA, termB), 0)
        else:
            countAB = pairwise_counts.get((termB, termA), 0)

        n_diff = max([countA - countAB, countB - countAB])

        if n_diff >= min_diff:
            continue
        
        levelA = len(networkx.descendants(hpo_g, termA))
        levelB = len(networkx.descendants(hpo_g, termB))

        if levelA < levelB:
            terms_to_prune.append(termB)
            if filter_log is not None:
                write_prune_log(f_log, termA, termB,
                                hpo_g.nodes(data=True)[termA]['name'],
                                hpo_g.nodes(data=True)[termB]['name'],
                                '{3} is higher-order term than {2}')
        elif levelA > levelB:
            terms_to_prune.append(termA)
            if filter_log is not None:
                write_prune_log(f_log, termB, termA,
                                hpo_g.nodes(data=True)[termB]['name'],
                                hpo_g.nodes(data=True)[termA]['name'],
                                '{2} is higher-order term than {3}')
        else:
            if countA >= countB:
                terms_to_prune.append(termB)
                if filter_log is not None:
                    write_prune_log(f_log, termA, termB,
                                    hpo_g.nodes(data=True)[termA]['name'],
                                    hpo_g.nodes(data=True)[termB]['name'],
                                    'both terms are same order, but {2} has ' + \
                                    'more samples than {3}')
            else:
                terms_to_prune.append(termA)
                if filter_log is not None:
                    write_prune_log(f_log, termB, termA,
                                    hpo_g.nodes(data=True)[termB]['name'],
                                    hpo_g.nodes(data=True)[termA]['name'],
                                    'both terms are same order, but {2} has ' + \
                                    'more samples than {3}')

    hpo_counts = { k : v for k, v in hpo_counts.items() if k not in terms_to_prune }

    return hpo_counts


def write_final_hpo_list(hpo_counts, hpo_g, outfile):
    """
    Format final set of HPO terms and write to outfile

    Output format: tsv, with the following columns
        1. HPO code
        2. Plain-text description
        3. Count of samples
        4. Tier
        5. Parent term(s)
        6. Child term(s)
    """

    # Gather metadata to write to file
    terms = list(hpo_counts.keys())
    counts = list(hpo_counts.values())
    tiers = []
    descrips = []
    parents = []
    children = []
    for term, count in hpo_counts.items():

        if term in hpo_g.nodes():
            tiers.append(len(networkx.descendants(hpo_g, term)))
            descrips.append(hpo_g.nodes(data=True)[term]['name'])

            l_parents = networkx.descendants(hpo_g, term)
            l_parents = [t for t in l_parents if t in hpo_counts.keys()]
            if len(l_parents) > 0:
                parents.append(';'.join(l_parents))
            else:
                parents.append('NA')

            l_children = networkx.ancestors(hpo_g, term)
            l_children = [t for t in l_children if t in hpo_counts.keys()]
            if len(l_children) > 0:
                children.append(';'.join(l_children))
            else:
                children.append('NA')
            
        else:
            if term == 'HEALTHY_CONTROL':
                tiers.append(1)
                descrips.append('Unaffected control sample')
                parents.append('NA')
            else:
                tiers.append(2)
                descrips.append('NA')
                parents.append('HP:0000118')
            children.append('NA')

    # Determine sort order
    idx_counts = [(i, v) for i, (k, v) in enumerate(hpo_counts.items())]
    order = [i for i, v in sorted(idx_counts, reverse=True, key=lambda x: x[1])]
    # order = range(0, len(terms))

    # Write header to file
    header = '#HPO_term\tdescription\tsamples\tHPO_tier' + \
             '\tparent_terms\tchild_terms\n'
    outfile.write(header)

    # Iterate over HPO terms and write to file
    newline = '\t'.join(['{' + str(k) + '}' for k in range(0, 6)]) + '\n'
    for i in order:
        outfile.write(newline.format(terms[i], descrips[i], str(counts[i]), 
                                     str(tiers[i]), parents[i], children[i]))


def main():
    """
    Main block
    """

    # Parse arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('phenos', help='Test file of phenotypes. One line per ' +
                        'patient. Two tab-delimited columns: patient ID and ' +
                        'string of semicolon-delimited HPO terms.')
    parser.add_argument('--obo', help='Path to HPO .obo file. [default: ' + 
                        'http://purl.obolibrary.org/obo/hp.obo]',
                        default='http://purl.obolibrary.org/obo/hp.obo',
                        metavar='(file|url)')
    parser.add_argument('--ignore-term', action='append', help='Entirely ' + \
                        'ignore specified term.', metavar='string')
    parser.add_argument('-o', '--outfile', help='Path to outfile. ' +
                        '[default: stdout]', metavar='file')
    parser.add_argument('--raw-counts', help='Path to optional output file ' +
                        'containing all raw sample counts per HPO term. ' + 
                        '[default: do not report raw counts]',
                        metavar='path')
    parser.add_argument('--filter-log', help='Path to optional output file ' +
                        'containing logs for filtered HPO terms. ' + 
                        '[default: do not log filtering steps]',
                        metavar='path')
    parser.add_argument('--min-samples', help='Minimum number of samples per ' + 
                        'HPO term. [default: 1,000]', default=1000, type=int, 
                        metavar='int')
    parser.add_argument('--min-diff', help='Minimum number of samples that ' + 
                        'must differ between HPO terms before pruning one or ' +
                        'the other. [default: 1,000]', default=1000, type=int, 
                        metavar='int')
    args = parser.parse_args()

    # Open connection to obo
    hpo_g = obonet.read_obo(args.obo)

    # Open connection to outfile
    if args.outfile is None:
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Count number of samples per HPO term
    hpo_counts = tally_hpo_counts(args.phenos, args.ignore_term)

    # Write out raw counts, if optioned
    if args.raw_counts is not None:
        write_raw_counts(hpo_counts, args.raw_counts)

    # Prune HPO terms based on minimum number of samples
    hpo_counts = prune_small_terms(hpo_counts, args.min_samples)

    # Prune HPO terms based on relatedness and sample overlap
    related_terms = get_related_terms(hpo_counts, hpo_g)
    related_terms = filter_related_term_list(hpo_counts, related_terms, 
                                             args.min_diff)
    pairwise_counts = count_overlapping_samples(related_terms, args.phenos)
    hpo_counts = prune_related_terms(hpo_counts, related_terms, pairwise_counts, 
                                     hpo_g, args.min_diff, args.filter_log)

    # Format & write final list of pruned HPO terms for analysis
    write_final_hpo_list(hpo_counts, hpo_g, outfile)


if __name__ == '__main__':
    main()
