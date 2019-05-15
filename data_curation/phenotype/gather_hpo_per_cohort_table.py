#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Compute summary table of counts of samples per HPO term per cohort
"""


from os import path
import csv
import locale
import argparse
from sys import stdout


locale.setlocale(locale.LC_ALL, 'en_us')


def gather_counts(cohort_table, hpo_metadata, hpo_dir):
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
                counts[cohort] = { 'HEALTHY_CONTROL' : ctrl_n }

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

    return counts


def write_table_header(outfile, counts):
    """
    Write a HTML-formatted header to outfile
    """

    column_names = ['HPO', 'description', '**Total**'] + list(counts.keys())
    outfile.write('| ' + ' | '.join(column_names) + ' |  \n')

    hr_line = [':---'] * 2 + ['---:'] * (len(counts.keys()) + 1)
    outfile.write('| ' + ' | '.join(hr_line) + ' |  \n')


def write_table(counts, hpo_metadata, outfile):
    """
    Master function to process input data and write summary table
    """

    # Write header to outfile
    write_table_header(outfile, counts)

    # Add one row per phenotype
    with open(hpo_metadata) as infile:
        reader = csv.reader(infile, delimiter='\t')

        for term, descrip, total_n, tier, parents, children in reader:

            if term.startswith('#'):
                continue

            total_n = str(locale.format("%d", total_n, grouping=True))

            counts_per = [counts[cohort][term] for cohort in counts.keys()]
            counts_per = [locale.format("%d", k, grouping=True) for k in counts_per]
            counts_per = [str(k) for k in counts_per]

            outvals = [term, descrip, total_n] + counts_per

            outfile.write('| ' + ' | '.join(outvals) + ' |  \n')


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
    parser.add_argument('hpo_dir', help='Path to directory with cleaned ' + \
                        'sample-level HPO assignments.')
    parser.add_argument('-o', '--outfile', help='Path to output file. ' +
                        '[default: stdout]')
    args = parser.parse_args()

    # Clean hpo directory specification
    if args.hpo_dir.endswith('/'):
        hpo_dir = args.hpo_dir
    else:
        hpo_dir = args.hpo_dir + '/'

    # Gather table of HPO codes per cohort
    counts = gather_counts(args.cohort_table, args.hpo_metadata, hpo_dir)

    # Open connection to outfile
    if args.outfile is None:
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Format output table & write to file
    write_table(counts, args.hpo_metadata, outfile)


if __name__ == '__main__':
    main()

