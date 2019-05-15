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


locale.setlocale(locale.LC_ALL, 'en_US.utf8')


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


def write_table_header(outfile, counts):
    """
    Write a HTML-formatted header to outfile
    """

    column_names = ['HPO', 'description', 'Total'] + list(counts.keys())
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

            if term == 'HEALTHY_CONTROL':
                total_n = 0
                for cohort in counts.keys():
                    total_n += counts[cohort]['HEALTHY_CONTROL']

            total_n = str(locale.format_string("%d", int(total_n), grouping=True))

            counts_per = [counts[cohort][term] for cohort in counts.keys()]
            counts_per = [locale.format_string("%d", int(k), grouping=True) for k in counts_per]
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
    parser.add_argument('--meta-cohorts', help='Path to tsv specifying ' + \
                        'metacohorts to consider. Two columns: metacohort ID ' + \
                        'and semicolon-separated cohorts to include. ' + \
                        '[default: no metacohort output]', dest='metalist')
    parser.add_argument('--meta-out', help='Path to output file for counts ' + \
                        'summarized per metacohort. [default: no metacohort ' + \
                        'output]', dest='meta_out')
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

    # Format output table for metacohorts, if optioned
    if args.metalist is not None \
    and args.meta_out is not None:
        metaoutfile = open(args.meta_out, 'w')
        metacounts = gather_metacounts(counts, args.metalist)
        write_table(metacounts, args.hpo_metadata, metaoutfile)


if __name__ == '__main__':
    main()

