#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Compute summary table of CNV stats per cohort
"""

from os import path
import csv
import pybedtools as pbt
from numpy import median
import argparse
from sys import stdout


def load_cohort_meta(meta_in):
    """
    Load cohort metadata into dictionary
    """

    meta = {}

    with open(meta_in) as infile:
        reader = csv.reader(infile, delimiter='\t')

        for cohort, total_n, case_n, ctrl_n, bed in reader:

            if cohort not in meta.keys():
                meta[cohort] = { 'samples' : int(total_n),
                                 'cases' : int(case_n),
                                 'controls' : int(ctrl_n),
                                 'bed' : bed }

    return meta


def cnv_stats_single(cnvs):
    """
    Collect basic CNV stats for an input pbt.BedTool
    """

    cnvs_case = cnvs.filter(lambda x: x[5] != 'HEALTHY_CONTROL').saveas()
    cnvs_ctrl = cnvs.filter(lambda x: x[5] == 'HEALTHY_CONTROL').saveas()

    n_cnv_all = len(cnvs)
    med_cnv_all = median([x.length for x in cnvs])

    n_cnv_case = len(cnvs_case)
    if len(cnvs_case) > 0:
        med_cnv_case = median([x.length for x in cnvs_case])
    else:
        med_cnv_case = 0

    n_cnv_ctrl = len(cnvs_ctrl)
    if len(cnvs_ctrl) > 0:    
        med_cnv_ctrl = median([x.length for x in cnvs_ctrl])
    else:
        med_cnv_ctrl = 0

    cnv_stats = { 'n_all' : n_cnv_all,
                  'med_all' : med_cnv_all,
                  'n_case' : n_cnv_case,
                  'med_case' : med_cnv_case,
                  'n_ctrl' : n_cnv_ctrl,
                  'med_ctrl' : med_cnv_ctrl }

    return cnv_stats


def calc_cnv_stats(cohort, meta):
    """
    Calculate CNV metadata for a single cohort
    """

    # Load CNVs
    cnvs = pbt.BedTool(meta[cohort]['bed'])
    dels = cnvs.filter(lambda x: x[4] == 'DEL').saveas()
    dups = cnvs.filter(lambda x: x[4] == 'DUP').saveas()

    # Get stats
    cnv_stats = cnv_stats_single(cnvs)
    del_stats = cnv_stats_single(dels)
    dup_stats = cnv_stats_single(dups)

    res = {'CNV' : cnv_stats,
           'DEL' : del_stats,
           'DUP' : dup_stats }

    return res


def write_html_table(meta, cnv_stats, htmlout):
    """
    Create html-formatted table of CNV stats for README
    """

    # Write header
    colnames = ['Dataset', 'N Cases', 'Case CNVs', 'CNVs /Case', 
                'Case Median Size', 'Case DEL:DUP', 'N Ctrls', 'Ctrl CNVs', 
                'CNVs /Ctrl', 'Ctrl Median Size', 'Ctrl DEL:DUP']
    header = '| ' + ' | '.join(colnames) + ' |  '
    htmlout.write(header + '\n')
    align = ['---', '---:', '---:', '---:', '---:', '---:', '---:', '---:', 
             '---:', '---:', '---:']
    subheader = '| ' + ' | '.join(align) + ' |  '
    htmlout.write(subheader + '\n')

    # Print one line per cohort
    for cohort in meta.keys():
        # Add case data
        n_case = meta[cohort]['cases']
        n_cnv = cnv_stats[cohort]['CNV']['n_case']
        newline = '{} {:,} {:,}'.format(cohort, n_case, n_cnv).split()
        if n_case > 0:
            medsize = cnv_stats[cohort]['CNV']['med_case'] / 1000
            n_del = cnv_stats[cohort]['DEL']['n_case']
            n_dup = cnv_stats[cohort]['DUP']['n_case']
            ratio = n_del / n_dup
            if ratio >= 1:
                ratio_str = '{:.2f}:1'.format(ratio)
            else:
                ratio_str = '1:{:.2f}'.format(1 / ratio)
            newline = newline + [str(round(n_cnv / n_case, 2)),
                                 str(round(medsize, 0)) + ' kb',
                                 ratio_str]
        else:
            newline = newline + ['0', 'NA', 'NA']
        # Add control data
        n_ctrl = meta[cohort]['controls']
        n_cnv = cnv_stats[cohort]['CNV']['n_ctrl']
        newline = newline + '{:,} {:,}'.format(n_ctrl, n_cnv).split()
        if n_ctrl > 0:
            medsize = cnv_stats[cohort]['CNV']['med_ctrl'] / 1000
            n_del = cnv_stats[cohort]['DEL']['n_ctrl']
            n_dup = cnv_stats[cohort]['DUP']['n_ctrl']
            ratio = n_del / n_dup
            if ratio >= 1:
                ratio_str = '{:.2f}:1'.format(ratio)
            else:
                ratio_str = '1:{:.2f}'.format(1 / ratio)
            newline = newline + [str(round(n_cnv / n_ctrl, 2)),
                                 str(round(medsize, 0)) + ' kb', 
                                 ratio_str]
        else:
            newline = newline + ['0', 'NA', 'NA']
        # Write out
        htmlout.write('| ' + ' | '.join(newline) + ' |  \n')


def write_plot_table(meta, cnv_stats, outfile):
    """
    Create table of CNV stats for plotting
    """

    cnvs = 'DEL DUP'.split()

    # Write header
    cols = '{1}_per_{0}\tmed_{0}_{1}_size'
    header = 'cohort'
    for pheno in 'case ctrl'.split():
        header = header + '\tn_' + pheno
        for cnv in cnvs:
            header = '\t'.join([header, cols.format(pheno, cnv)])
        header = header + '\t' + '{0}_DEL_DUP_ratio'.format(pheno)
    outfile.write(header + '\n')

    # Print one line per cohort
    for cohort in meta.keys():
        # Add case data
        n_case = meta[cohort]['cases']
        newline = [cohort, n_case]
        n_del = cnv_stats[cohort]['DEL']['n_case']
        n_dup = cnv_stats[cohort]['DUP']['n_case']
        if n_case > 0:
            newline = newline + [n_del / n_case,
                                 cnv_stats[cohort]['DEL']['med_case'],
                                 n_dup / n_case,
                                 cnv_stats[cohort]['DUP']['med_case'],
                                 n_del / n_dup]
        else:
            newline = newline + [0, 'NA', 0, 'NA', 'NA']
        # Add control data
        n_ctrl = meta[cohort]['controls']
        newline.append(n_ctrl)
        n_del = cnv_stats[cohort]['DEL']['n_ctrl']
        n_dup = cnv_stats[cohort]['DUP']['n_ctrl']
        if n_ctrl > 0:
            newline = newline + [n_del / n_ctrl,
                                 cnv_stats[cohort]['DEL']['med_ctrl'],
                                 n_dup / n_ctrl,
                                 cnv_stats[cohort]['DUP']['med_ctrl'],
                                 n_del / n_dup]
        else:
            newline = newline + [0, 'NA', 0, 'NA', 'NA']
        # Write out
        outfile.write('\t'.join([str(x) for x in newline]) + '\n')


def main():
    """
    Main block
    """

    # Parse arguments & options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cohort_metadata', help='Table of cohort metadata')
    parser.add_argument('--html', help='Path to output table in HTML format. ' +
                        '[default: stdout]', default='stdout')
    parser.add_argument('--tsv', help='Path to output .tsv file.')
    args = parser.parse_args()

    # Open connection to outfile(s)
    if args.html in 'stdout -'.split():
        htmlout = stdout
    else:
        htmlout = open(args.html, 'w')
    if args.tsv is not None:
        tsvout = open(args.tsv, 'w')

    # Read cohort metadata
    meta = load_cohort_meta(args.cohort_metadata)
    
    # Annotate CNV stats for each cohort
    cnv_stats = {}
    for cohort in meta.keys():
        cnv_stats[cohort] = calc_cnv_stats(cohort, meta)

    # Build html table
    if args.html is not None:
        write_html_table(meta, cnv_stats, htmlout)

    # Build tsv table for plotting, if optioned
    write_plot_table(meta, cnv_stats, tsvout)


if __name__ == '__main__':
    main()

