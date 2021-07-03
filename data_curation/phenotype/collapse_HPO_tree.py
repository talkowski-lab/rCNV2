#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Determine the minimal set of most specific HPO terms to retain for analysis 
based on a list of patients and corresponding HPO terms
"""


import gzip
import networkx
import obonet
import csv
import argparse
from sys import stdout


def load_precomp_pairs_single_cohort(pairs_tsv_in):
    """
    Load precomputed HPO pairs for a single cohorts, if optioned
    Returns a dict of sample counts keyed on sorted HPO term pairs
    """

    counts = {}

    if pairs_tsv_in.endswith('.gz'):
        fin = gzip.open(pairs_tsv_in, 'rt')
    else:
        fin = open(pairs_tsv_in)

    for term1, term2, n1, n2, n12 in csv.reader(fin, delimiter='\t'):

        if term1.startswith('#'):
            continue

        pair_key = tuple(sorted((term1, term2)))
        if pair_key not in counts.keys():
            counts[pair_key] = int(n12)

    return counts


def load_precomp_pairs(hpo_pair_cohorts_tsv):
    """
    Wraps load_precomp_pairs_single_cohort() for all cohorts in --hpo-pair-cohorts
    """

    precomp_pairs = {}

    with open(hpo_pair_cohorts_tsv) as fin:
        for cohort, tsv_path in csv.reader(fin, delimiter='\t'):
            precomp_pairs[cohort] = load_precomp_pairs_single_cohort(tsv_path)

    return precomp_pairs


def tally_hpo_counts(pheno_file, ignore_terms, precomp_pairs=None):
    """
    Create dictionary of counts of samples per HPO term in input file
    """

    hpo_counts = {}

    with open(pheno_file) as infile:
        reader = csv.reader(infile, delimiter='\t')

        for sample, terms in reader:
            # Make sure samples with base term also have HP:0000118
            sterms = terms.split(';')
            if 'HP:0000001' in sterms and 'HP:0000118' not in sterms:
                sterms.append('HP:0000118')

            for term in sterms:

                if ignore_terms is not None:
                    if term in ignore_terms:
                        continue

                if term in hpo_counts.keys():
                    hpo_counts[term] += 1
                else:
                    hpo_counts[term] = 1

    if precomp_pairs is not None:
        for cohort_dict in precomp_pairs.values():

            cohort_pairs = cohort_dict.keys()
            cohort_terms = set([term for pair in cohort_pairs for term in pair])

            for term in cohort_terms:

                if ignore_terms is not None:
                    if term in ignore_terms:
                        continue

                N = cohort_dict.get((term, term), 0)

                if term in hpo_counts.keys():
                    hpo_counts[term] += N
                else:
                    hpo_counts[term] = N

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


def pair_all_terms(terms):
    """
    Make a list of tuples for all combinations of terms
    """

    all_term_pairs = []

    for a in range(len(terms)):
        termA = terms[a]

        for b in range(len(terms)):
            termB = terms[b]
            
            if a >= b:
                continue

            else:
                pair = tuple(sorted([termA, termB]))
                all_term_pairs.append(pair)

    return all_term_pairs


def prune_jaccard_terms(hpo_counts, pairwise_counts, hpo_g, max_jac, filter_log):
    """
    Prune all terms (related or not) based on Jaccard index of samples shared
    between each pair of terms
    """

    dropped_terms = []

    for termA, termB in pairwise_counts.keys():

        if termA in dropped_terms or termB in dropped_terms:
            continue

        if termA == termB:
            continue

        countA = hpo_counts.get(termA, 0)
        countB = hpo_counts.get(termB, 0)

        if (termA, termB) in pairwise_counts.keys():
            countAB = pairwise_counts.get((termA, termB), 0)
        else:
            countAB = pairwise_counts.get((termB, termA), 0)

        jac = countAB / (countA + countB - countAB)

        if jac > max_jac:
            if countA >= countB:
                hpo_counts.pop(termB)
                dropped_terms.append(termB)
                write_prune_log(filter_log, termA, termB,
                                hpo_g.nodes(data=True)[termA]['name'],
                                hpo_g.nodes(data=True)[termB]['name'],
                                '{3} and {2} have Jaccard similarity of ' + \
                                str(round(jac, 2)) + ' and {2} is larger.')
            else:
                hpo_counts.pop(termA)
                dropped_terms.append(termA)
                write_prune_log(filter_log, termB, termA,
                                hpo_g.nodes(data=True)[termB]['name'],
                                hpo_g.nodes(data=True)[termA]['name'],
                                '{3} and {2} have Jaccard similarity of ' + \
                                str(round(jac, 2)) + ' and {2} is larger.')

    return hpo_counts
            


def calc_term_overlap(termsA, termsB):
    """
    Calculate the reciprocal overlap between two lists of HPO terms
    """

    nA = len(termsA)
    nB = len(termsB)

    if nA == 0 or nB == 0:
        ro = 0
    else:
        nOvr = len(set(termsA).intersection(set(termsB)))
        oA = nOvr / nA
        oB = nOvr / nB
        ro = min([oA, oB])

    return ro


def get_related_terms(hpo_counts, hpo_g, min_sib_frac=0.5):
    """
    Define pairs of related terms to be evaluated
    """

    related_terms = []

    parents_d = {}

    siblings_d = {}

    # Get all parents
    for term, count in hpo_counts.items():
        if term in hpo_g.nodes():
            parents = networkx.descendants(hpo_g, term)
            parents = [t for t in parents if t in hpo_counts.keys()]
            parents_d[term] = parents

    # Get all siblings (based on reciprocal sharing of at least X% of parent terms)
    for term, count in hpo_counts.items():
        if term in hpo_g.nodes():
            parents = parents_d[term]
            siblings = []
            for termB, parentsB in parents_d.items():
                if termB != term:
                    if calc_term_overlap(parents, parentsB) >= min_sib_frac:
                        siblings.append(termB)
            siblings_d[term] = siblings

    # Get all children and filter nodes
    for term, count in hpo_counts.items():
        if term in hpo_g.nodes():
            parents = parents_d[term]
            siblings = siblings_d[term]
            children = networkx.ancestors(hpo_g, term)
            children = [t for t in children if t in hpo_counts.keys()]
            relatives = list(set(parents + children + siblings))

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


def count_overlapping_samples(related_terms, pheno_file, precomp_pairs=None):
    """
    Count number of overlapping samples for pairs of HPO terms
    """

    pairwise_counts = {}
    eligible_terms = []

    # Populate list of all eligible terms to be considered for pairwise overlap
    for termA, termB in related_terms:
        term_pair = tuple(sorted([termA, termB]))
        pairwise_counts[term_pair] = 0

        if termA not in eligible_terms:
            eligible_terms.append(termA)
        if termB not in eligible_terms:
            eligible_terms.append(termB)

    # if precomp_pairs is not None:
    #     for cohort_dict in precomp_pairs.values():
    #         for termA, termB in cohort_dict.keys():
    #             if termA == termB:
    #                 continue

    #             if termA in eligible_terms and termB in eligible_terms:
    #                 term_pair = tuple(sorted([termA, termB]))
    #                 if term_pair not in pairwise_counts.keys():
    #                     pairwise_counts[term_pair] = 0

    #             if termA not in eligible_terms:
    #                 eligible_terms.append(termA)
    #             if termB not in eligible_terms:
    #                 eligible_terms.append(termB)

    # Count pairwise overlap from phenos input
    with open(pheno_file) as infile:
        reader = csv.reader(infile, delimiter='\t')

        for sample, terms in reader:
            filtered_terms = [t for t in terms.split(';') if t in eligible_terms]

            if len(filtered_terms) > 1:
                for i in range(0, len(filtered_terms) - 1):
                    for j in range(i + 1, len(filtered_terms)):
                        termA = filtered_terms[i]
                        termB = filtered_terms[j]
                        term_pair = tuple(sorted([termA, termB]))
                        if term_pair in pairwise_counts.keys():
                            pairwise_counts[term_pair] += 1

    # Add pairwise overlaps from precomp_pairs, if optioned
    if precomp_pairs is not None:
        for cohort_dict in precomp_pairs.values():

            for termA, termB in cohort_dict.keys():
                if termA == termB:
                    continue

                term_pair = tuple(sorted([termA, termB]))
                if term_pair in pairwise_counts.keys():
                    pairwise_counts[term_pair] += cohort_dict.get(term_pair, 0)

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

    for termA, termB in related_terms:

        if termA in terms_to_prune:
            if filter_log is not None:
                msg = 'Skipping candidate pair "{0}" ({1}) & "{2}" ({3}): ' + \
                      '{1} has already been pruned.\n'
                filter_log.write(msg.format(hpo_g.nodes(data=True)[termA]['name'], termA, 
                                       hpo_g.nodes(data=True)[termB]['name'], termB))
            continue
        if termB in terms_to_prune:
            if filter_log is not None:
                msg = 'Skipping candidate pair "{0}" ({1}) & "{2}" ({3}): ' + \
                      '{3} has already been pruned.\n'
                filter_log.write(msg.format(hpo_g.nodes(data=True)[termA]['name'], termA, 
                                       hpo_g.nodes(data=True)[termB]['name'], termB))
            continue

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
                write_prune_log(filter_log, termA, termB,
                                hpo_g.nodes(data=True)[termA]['name'],
                                hpo_g.nodes(data=True)[termB]['name'],
                                '{3} is higher-order term than {2}')
        elif levelA > levelB:
            terms_to_prune.append(termA)
            if filter_log is not None:
                write_prune_log(filter_log, termB, termA,
                                hpo_g.nodes(data=True)[termB]['name'],
                                hpo_g.nodes(data=True)[termA]['name'],
                                '{2} is higher-order term than {3}')
        else:
            if countA >= countB:
                terms_to_prune.append(termB)
                if filter_log is not None:
                    write_prune_log(filter_log, termA, termB,
                                    hpo_g.nodes(data=True)[termA]['name'],
                                    hpo_g.nodes(data=True)[termB]['name'],
                                    'both terms are same order, but {2} has ' + \
                                    'at least as many samples as {3}')
            else:
                terms_to_prune.append(termA)
                if filter_log is not None:
                    write_prune_log(filter_log, termB, termA,
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
    parser.add_argument('--hpo-pair-cohorts', help='Two-column .tsv of cohort ' +
                        'name and path to pairwise HPO counts for cohorts where ' +
                        'necessary')
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
    parser.add_argument('--max-jac', help='Maximum Jaccard index for samples ' +
                        'belonging to any two HPO terms before pruning one or ' +
                        'the other. [default: 0.8]', default=0.8, type=float, 
                        metavar='float')
    args = parser.parse_args()

    # Open connection to obo
    hpo_g = obonet.read_obo(args.obo)

    # Open connection to outfile
    if args.outfile is None:
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Open connection to filter log, if optioned
    if args.filter_log is not None:
        filter_log = open(args.filter_log, 'w')
    else:
        filter_log = None

    # Load counts from cohorts with precomputed HPO pairs, if any
    if args.hpo_pair_cohorts is not None:
        precomp_pairs = load_precomp_pairs(args.hpo_pair_cohorts)
    else:
        precomp_pairs = None

    # Count number of samples per HPO term
    hpo_counts = tally_hpo_counts(args.phenos, args.ignore_term, precomp_pairs)

    # Write out raw counts, if optioned
    if args.raw_counts is not None:
        write_raw_counts(hpo_counts, args.raw_counts)

    # Prune HPO terms based on minimum number of samples
    hpo_counts = prune_small_terms(hpo_counts, args.min_samples)

    # Prune HPO terms based on Jaccard index
    all_term_pairs = pair_all_terms(list(hpo_counts.keys()))
    pairwise_counts = count_overlapping_samples(all_term_pairs, args.phenos, 
                                                precomp_pairs)
    hpo_counts = prune_jaccard_terms(hpo_counts, pairwise_counts, hpo_g,
                                     args.max_jac, filter_log)

    # Prune HPO terms based on relatedness and raw sample overlap
    related_terms = get_related_terms(hpo_counts, hpo_g)
    related_terms = filter_related_term_list(hpo_counts, related_terms, 
                                             args.min_diff)
    hpo_counts = prune_related_terms(hpo_counts, related_terms, pairwise_counts, 
                                     hpo_g, args.min_diff, filter_log)

    # Format & write final list of pruned HPO terms for analysis
    write_final_hpo_list(hpo_counts, hpo_g, outfile)


if __name__ == '__main__':
    main()
