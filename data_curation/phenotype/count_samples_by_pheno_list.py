#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Count samples per cohort matching at least one (or none) of an input list of HPO terms
"""


import gzip
import csv
import argparse
from sys import stdout


# Import functions from collapse_HPO_tree.py for convenience
from collapse_HPO_tree import load_precomp_pairs_single_cohort, load_precomp_pairs, tally_hpo_counts


def count_samples_from_tsv(phenofile, keep_hpos, invert=False):
    """
    Count number of samples in phenofile with at least one of keep_hpos
    """

    n = 0

    with open(phenofile) as fin:
        for sid, phenos in csv.reader(fin, delimiter='\t'):
            hits = [h in phenos.split(';') for h in keep_hpos]
            if any(hits):
                if not invert:
                    n += 1
            else:
                if invert:
                    n += 1

    return n


def load_related_terms(hpo_tiers_in):
    """
    Load a dict of HPO term relationships by tier from a precomputed .tsv
    """

    rel_dict = {}

    with open(hpo_tiers_in) as fin:

        reader = csv.reader(fin, delimiter='\t')
        for term, descrip, total_n, tier, parents, children in reader:

            if term.startswith('#'):
                continue

            if parents == 'NA':
                parents = []
            else:
                parents = parents.split(';')

            if children == 'NA':
                children = []
            else:
                children = children.split(';')

            rel_dict[term] = {'parents' : parents, 'children' : children}

    return rel_dict


def estimate_samples_from_pairs(pairs_dict, keep_hpos, invert=False, 
                                consider_controls=False):
    """
    Estimate the number of samples with at least one of keep_hpos from precomputed HPO pairs
    """

    hpo_counts = {h : pairs_dict.get((h, h), 0) for h in keep_hpos}
    hpo_counts = {h : n for h, n in sorted(hpo_counts.items(), key=lambda x: -x[1]) if n > 0}

    used_terms = []
    n = 0

    # Note: this assumes all hpos in keep_hpos are independent, which won't
    # technically be true but should be a reasonably close approximation given
    # that all direct parent-child terms have already been pruned
    for term, nt in hpo_counts.items():
        p_keep = 1
        for oterm in used_terms:
            p_ovr = pairs_dict.get(tuple(sorted((term, oterm))), 0) / nt
            p_keep = p_keep * (1 - p_ovr)
        n += (p_keep * nt)
        used_terms.append(term)

    n = int(round(n, 0))

    if invert:
        n = pairs_dict[('HP:0000118', 'HP:0000118')] - n
        if not consider_controls:
            n += pairs_dict[('HEALTHY_CONTROL', 'HEALTHY_CONTROL')]

    return n


def main():
    """
    Main block
    """

    # Parse arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('phenos', help='Two-column .tsv of cohort name and ' +
                        'path to phenotype file for that cohort.')
    parser.add_argument('hpos', help='List of HPO terms to include')
    parser.add_argument('--hpo-pair-cohorts', help='Two-column .tsv of cohort ' +
                        'name and path to pairwise HPO counts for cohorts where ' +
                        'necessary', metavar='tsv')
    parser.add_argument('--hpo-tier-metadata', help='.tsv of HPO metadata generated ' +
                        'by convert_phenotypes_to_hpo.sh. Required only if ' +
                        '--hpo-pair-cohorts is also provided.')
    parser.add_argument('-v', '--invert', action='store_true', help='Invert sample ' +
                        'inclusion criteria (similar to grep -v) [default: False].')
    parser.add_argument('--consider-controls', action='store_true', help='Consider ' +
                        'control samples when counting [default: exclude controls].')
    parser.add_argument('-o', '--outfile', help='Path to outfile. ' +
                        '[default: stdout]', metavar='tsv')
    args = parser.parse_args()

    # Open connection to outfile
    if args.outfile is None:
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    outfile.write('#cohort\tN_samples\n')

    # Load list of HPOs to retain
    keep_hpos = [x.rstrip() for x in open(args.hpos).readlines()]
    if args.invert and not args.consider_controls:
        keep_hpos.append('HEALTHY_CONTROL')

    # Remove redundant child HPO terms from list of terms to keep for precomputed matrixes
    if args.hpo_tier_metadata is not None:
        rel_dict = load_related_terms(args.hpo_tier_metadata)
        redundant_hpos = []
        for term, rels in rel_dict.items():
            if any([p in keep_hpos for p in rels['parents']]):
                redundant_hpos.append(term)
        keep_hpos = [h for h in keep_hpos if h not in redundant_hpos]

    # Compute counts per cohort in args.phenos
    with open(args.phenos) as tsvin:
        for cohort, phenofile in csv.reader(tsvin, delimiter='\t'):
            n = count_samples_from_tsv(phenofile, keep_hpos, args.invert)
            outfile.write('\t'.join([cohort, str(n)]) + '\n')

    # Estimate counts per cohort from precomputed HPO pairs, if any
    if args.hpo_pair_cohorts is not None:

        if args.hpo_tier_metadata is None:
            exit('Must specify --hpo-tier-metadata with --hpo-pair-cohorts')

        precomp_pairs = load_precomp_pairs(args.hpo_pair_cohorts)

        for cohort, pairs_dict in precomp_pairs.items():
            n = estimate_samples_from_pairs(pairs_dict, keep_hpos, args.invert, 
                                            args.consider_controls)
            outfile.write('\t'.join([cohort, str(n)]) + '\n')

    # Close outfile to clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

