#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Extract significant genes from a single phenotype's meta-analysis BED file
Successor script to get_significant_genes.R
"""


from finemap_genes import format_stat, get_sig_label
from os import path
import gzip
import csv
import argparse
from sys import stdout


def parse_stats(stats_in, primary_p_cutoff, secondary_p_cutoff=0.05, 
                n_nominal_cutoff=2, secondary_or_nominal=True, fdr_q_cutoff=0.05, 
                secondary_for_fdr=False, ew_sig_only=False, p_is_neg_log10=True):
    """
    Input: .bed file of meta-analysis association stats
    Output: set of significant genes
    """

    if path.splitext(stats_in)[1] in '.gz .bgz .bgzip'.split():
        csvin = gzip.open(stats_in, 'rt')
    else:
        csvin = open(stats_in)
    reader = csv.reader(csvin, delimiter='\t')

    sig_genes = set()

    for chrom, start, end, gene, n_nominal, top_cohort, excluded_cohorts, \
        case_freq, control_freq, lnOR, lnOR_lower, lnOR_upper, zscore, \
        primary_p, primary_q, secondary_lnOR, secondary_lnOR_lower, \
        secondary_lnOR_upper, secondary_zscore, secondary_p, secondary_q \
        in reader:

        # Skip header line
        if chrom.startswith('#'):
            continue

        # Determine gene significance
        primary_p = format_stat(primary_p, p_is_neg_log10, 1)
        primary_q = format_stat(primary_q, p_is_neg_log10, 1)
        secondary_p = format_stat(secondary_p, p_is_neg_log10, 1)
        n_nominal = int(n_nominal)
        sig_label = get_sig_label(primary_p, secondary_p, n_nominal, primary_q, 
                                  primary_p_cutoff, secondary_p_cutoff, 
                                  n_nominal_cutoff, secondary_or_nominal, 
                                  fdr_q_cutoff, secondary_for_fdr)

        # Add gene to sig_genes if significant
        if ew_sig_only:
            if sig_label == 'EWS':
                sig_genes.add(gene)
        else:
            if sig_label in 'EWS FDR'.split():
                sig_genes.add(gene)

    csvin.close()

    return sig_genes


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('statslist', help='.tsv of metadata per phenotype. Three ' +
                        'required columns: HPO, path to meta-analysis stats, ' +
                        'and primary P-value cutoff.')
    parser.add_argument('--secondary-p-cutoff', help='Maximum secondary P-value to ' + 
                        'consider as exome-wide significant. [default: 1]', 
                        default=1, type=float)
    parser.add_argument('--min-nominal', help='Minimum number of individual cohorts ' + 
                        'required to be nominally significant to consider ' +
                        'exome-wide significant. [default: 1]', default=1, type=int)
    parser.add_argument('--secondary-or-nominal', dest='secondary_or_nom', 
                        help='Allow genes to meet either --secondary-p-cutoff ' +
                        'or --min-nominal, but do not require both for exome-wide ' +
                        'significance. [default: require both]', default=False, 
                        action='store_true')
    parser.add_argument('--fdr-q-cutoff', help='Maximum FDR Q-value to ' + 
                        'consider as FDR significant. [default: 0.05]', 
                        default=0.05, type=float)
    parser.add_argument('--secondary-for-fdr', action='store_true', 
                        default=False, help='Apply sample secondary and/or min. ' +
                        'nominal criteria when evaluating FDR significance. ' +
                        '[default: apply no secondary criteria to FDR genes]')
    parser.add_argument('-o', '--outfile', help='Output .tsv of final list of ' +
                        'significant genes. [default: stdout]')
    args = parser.parse_args()

    # Open connection to output file
    if args.outfile is None \
    or args.outfile in 'stdout -'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    outfile.write('#gene\thpo\n')

    # Iterate over all HPOs
    with open(args.statslist) as tsvin:
        for hpo, sspath, primary_p_cutoff in csv.reader(tsvin, delimiter='\t'):
            # Get significant genes per HPO
            sig_genes = parse_stats(sspath, float(primary_p_cutoff), 
                                    args.secondary_p_cutoff, args.min_nominal, 
                                    args.secondary_or_nom, args.fdr_q_cutoff, 
                                    args.secondary_for_fdr)
            # Write gene, HPO pairs to outfile
            for gene in sig_genes:
                outfile.write('\t'.join([gene, hpo]) + '\n')

    outfile.close()


if __name__ == '__main__':
    main()

