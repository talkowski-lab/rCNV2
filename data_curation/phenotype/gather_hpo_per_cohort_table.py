#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Compute summary table of counts of samples per HPO term per cohort
"""


from collapse_HPO_tree import load_precomp_pairs, load_precomp_pairs_single_cohort
from os import path
import csv
import locale
import argparse
from sys import stdout


locale.setlocale(locale.LC_ALL, 'en_US.utf8')


def gather_counts(cohort_table, hpo_metadata, hpo_dir, precomp_counts={}):
    """
    Gather nested dictionary of sample counts per HPO per study
    """

    with open(hpo_metadata) as f:
        phenos = [l.split('\t')[0] for l in f.readlines()]
        phenos = [p for p in phenos if not p.startswith('#') \
                                    and not p == 'HEALTHY_CONTROL']

    counts = {}

    with open(cohort_table) as infile:
        reader = csv.reader(infile, delimiter='\t')

        for cohort, total_n, case_n, ctrl_n in reader:

            if cohort.startswith('#'):
                continue

            if cohort not in counts.keys():
                counts[cohort] = { 'HEALTHY_CONTROL' : int(ctrl_n) }

            for pheno in phenos:
                if pheno not in counts[cohort].keys():
                    counts[cohort][pheno] = 0

            pheno_fpath = hpo_dir + cohort + '.cleaned_phenos.txt'

            if path.isfile(pheno_fpath):
                with open(pheno_fpath) as pfin:
                    pheno_csv = csv.reader(pfin, delimiter='\t')

                    for sample, terms in pheno_csv:
                        terms = [t for t in terms.split(';') if t in phenos]

                        for pheno in terms:
                            counts[cohort][pheno] += 1

            elif cohort in precomp_counts.keys():
                
                for pheno in phenos:
                    counts[cohort][pheno] = precomp_counts[cohort].get((pheno, pheno), 0)

            else:
                err = 'ERROR: could not locate phenotype information for cohort {}. Exiting.'
                from sys import exit
                exit(err.format(cohort))

    return counts


def apply_pairwise_curation(counts, tsvin):
    """
    Apply manual pairwise similarity pruning
    """

    with open(tsvin) as f:
        for keep, discard in csv.reader(f, delimiter='\t'):
            if keep.startswith('#'):
                continue

            for cohort in counts.keys():
                if discard in counts[cohort].keys():
                    counts[cohort].pop(discard)

    return counts


def gather_metacounts(counts, metalist):
    """
    Combine sample counts across individual cohorts within each metacohort
    """

    metacounts = {}

    with open(metalist) as infile:
        reader = csv.reader(infile, delimiter='\t')

        for name, cohorts in reader:
            if name.startswith('#'):
                continue
            else:
                metacounts[name] = {}

            for cohort in cohorts.split(';'):
                if cohort in counts.keys():
                    for pheno in counts[cohort].keys():
                        if pheno not in metacounts[name].keys():
                            metacounts[name][pheno] = int(counts[cohort][pheno])
                        else:
                            metacounts[name][pheno] += int(counts[cohort][pheno])

    return metacounts


def write_table_header(outfile, counts, html=False):
    """
    Write a formatted header to outfile
    """

    column_names = ['HPO', 'description', 'Total'] + list(counts.keys())

    if html:
        outfile.write('| ' + ' | '.join(column_names) + ' |  \n')
        hr_line = [':---'] * 2 + ['---:'] * (len(counts.keys()) + 1)
        outfile.write('| ' + ' | '.join(hr_line) + ' |  \n')
    else:
        outfile.write("#" + '\t'.join(column_names) + '\n')


def metacohort_abundance_filter(counts, metacounts, min_samples, min_metas):
    """
    Prune counts and metacounts based on minimum sample abundance
    """

    metacohorts = list(metacounts.keys())
    metacohorts = [m for m in metacohorts if m not in 'mega total'.split()]
    hpos = list(metacounts[metacohorts[0]].keys())

    passing_hpos = {}
    for hpo in hpos:
        passing_hpos[hpo] = len([m for m in metacohorts if metacounts[m][hpo] >= min_samples])
    keep_hpos = [m for m, k in passing_hpos.items() if k >= min_metas]

    for cohort in counts.keys():
        for pheno in hpos:
            if pheno not in keep_hpos:
                counts[cohort].pop(pheno)

    for meta in metacounts.keys():
        for pheno in hpos:
            if pheno not in keep_hpos:
                metacounts[meta].pop(pheno)

    return counts, metacounts


def write_table(counts, hpo_metadata, outfile, html=False):
    """
    Master function to process input data and write summary table
    """

    # Write header to outfile
    write_table_header(outfile, counts, html)

    # Add one row per phenotype
    with open(hpo_metadata) as infile:
        reader = csv.reader(infile, delimiter='\t')

        for term, descrip, total_n, tier, parents, children in reader:

            if term.startswith('#'):
                continue

            if term not in counts[list(counts.keys())[0]].keys():
                continue

            if term == 'HEALTHY_CONTROL':
                total_n = 0
                for cohort in counts.keys():
                    if cohort != 'mega':
                        total_n += counts[cohort]['HEALTHY_CONTROL']

            if html:
                total_n = str(locale.format_string("%d", int(total_n), grouping=True))
            else:
                total_n = str(total_n)

            counts_per = [counts[cohort][term] for cohort in counts.keys()]
            if html:
                counts_per = [locale.format_string("%d", int(k), grouping=True) for k in counts_per]
            counts_per = [str(k) for k in counts_per]

            outvals = [term, descrip, total_n] + counts_per

            if html:
                outfile.write('| ' + ' | '.join(outvals) + ' |  \n')
            else:
                outfile.write('\t'.join(outvals) + '\n')


def main():
    """
    Main block
    """

    # Parse arguments & options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('hpo_metadata', help='Table of HPO metadata.')
    parser.add_argument('cohort_table', help='Table of cohorts and total samples.')
    parser.add_argument('hpo_dir', help='Path to directory with cleaned ' + 
                        'sample-level HPO assignments.')
    parser.add_argument('--hpo-pair-cohorts', help='Two-column .tsv of cohort ' +
                        'name and path to pairwise HPO counts for cohorts where ' +
                        'necessary')
    parser.add_argument('--meta-cohorts', help='Path to tsv specifying ' + 
                        'metacohorts to consider. Two columns: metacohort ID ' + 
                        'and semicolon-separated cohorts to include. ' + 
                        '[default: no metacohort output]', dest='metalist')
    parser.add_argument('--meta-out', help='Path to output file for counts ' + 
                        'summarized per metacohort. [default: no metacohort ' + 
                        'output]', dest='meta_out')
    parser.add_argument('--min-per-metacohort', type=int, default=100,
                        help='Minimum number of samples required per metacohort ' +
                        'for at least --min-metacohort total metacohorts. [default 100]')
    parser.add_argument('--min-metacohorts', type=int, default=1,
                        help='Minimum number of metacohorts required to have ' +
                        'at least --min-per-metacohort samples each. [default 1]')
    parser.add_argument('--pairwise-curation', help='Two-column .tsv specifying ' +
                        'manual pairs to prune. First column = term to keep; ' +
                        'second column = term to drop')
    parser.add_argument('-o', '--outfile', help='Path to output file. ' +
                        '[default: stdout]')
    parser.add_argument('--html', action='store_true', help='Output tables in ' + 
                        ' HTML format. [default: flat tsv output]')
    args = parser.parse_args()

    # Clean hpo directory specification
    if args.hpo_dir.endswith('/'):
        hpo_dir = args.hpo_dir
    else:
        hpo_dir = args.hpo_dir + '/'

    # Load counts from cohorts with precomputed HPO pairs, if any
    if args.hpo_pair_cohorts is not None:
        precomp_pairs = load_precomp_pairs(args.hpo_pair_cohorts)
    else:
        precomp_pairs = None

    # Gather table of HPO codes per cohort
    counts = gather_counts(args.cohort_table, args.hpo_metadata, hpo_dir, precomp_pairs)

    # Manually prune pairs that are too similar, if optioned
    if args.pairwise_curation is not None:
        counts = apply_pairwise_curation(counts, args.pairwise_curation)

    # Open connection to outfile
    if args.outfile is None:
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Format output table for metacohorts, if optioned
    if args.metalist is not None \
    and args.meta_out is not None:
        metaoutfile = open(args.meta_out, 'w')
        metacounts = gather_metacounts(counts, args.metalist)

        # Filter HPO terms based on minimum metacohort abundance, if optioned
        counts, metacounts = metacohort_abundance_filter(counts, metacounts,
                                                         args.min_per_metacohort,
                                                         args.min_metacohorts)
    # Format output tables & write to file
    write_table(counts, args.hpo_metadata, outfile, args.html)
    if args.metalist is not None \
    and args.meta_out is not None:
        write_table(metacounts, args.hpo_metadata, metaoutfile, args.html)

if __name__ == '__main__':
    main()

