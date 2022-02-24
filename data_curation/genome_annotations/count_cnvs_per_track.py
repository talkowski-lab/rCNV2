#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Count CNVs per annotation track, split by CNV type and metacohort
"""


import csv
import pybedtools as pbt
import pandas as pd
import numpy as np
import argparse
from sys import stdout
from os import path
import subprocess


cnv_types = 'DEL DUP'.split()
phenotypes = 'case control'.split()


def load_cnvs_single(cnvpath, hpos, cnv_type):
    """
    Loads & filters CNVs to a single CNV type for a list of HPOs
    Returns : pbt.BedTool
    """

    cnvs = pbt.BedTool(cnvpath)

    return cnvs.filter(lambda x: len(set(x[-1].split(';')).\
                                     intersection(set(hpos))) > 0).\
                filter(lambda x: x[4] == cnv_type).\
                saveas()


def load_cnvs(cohorts_in, case_hpos, control_hpos):
    """
    Load CNVs and split by phenotype and CNV type
    Returns: 
        1. cnvs : multi-level dict of {cohort : {phenotype : {cnv_type : pbt.BedTool}}}
        2. cnv_counts : multi-level dict of {cohort : {phenotype : {cnv_type : int}}}
        3. cohort_n : multi-level dict of {cohort : {phenotype : int}}
    """

    cnvs = {}
    cnv_counts = {}
    cohort_n = {}
    pdict = {phenotypes[0] : case_hpos, phenotypes[1] : control_hpos}

    with open(cohorts_in) as fin:
        for cohort, n_case, n_control, cnvpath in csv.reader(fin, delimiter='\t'):
            cnvs[cohort] = {ptype : {ctype : load_cnvs_single(cnvpath, pdict[ptype], ctype) \
                                     for ctype in cnv_types} for ptype in phenotypes}
            cnv_counts[cohort] = {ptype : {ctype : len(cnvs[cohort][ptype][ctype]) \
                                  for ctype in cnv_types} for ptype in phenotypes}
            cohort_n[cohort] = {phenotypes[0] : int(n_case), 
                                phenotypes[1] : int(n_control)}

    return cnvs, cnv_counts, cohort_n


def load_counts_by_context(c_by_c_in, cohorts, case_hpo='DEVELOPMENTAL', 
                           control_hpo='HEALTHY_CONTROL'):
    """
    Load counts of CNVs per cohort, CNV type, and phenotype based on genic context
    Returns dict
    """

    # Load & filter table
    if c_by_c_in is not None:
        c_by_c_df = pd.read_csv(c_by_c_in, sep='\t').\
                       rename(columns={'#cohort' : 'cohort'})
        keepers = (c_by_c_df.hpo.isin([case_hpo, control_hpo])) & \
                  (c_by_c_df.gset == 'not_unconstrained')
        c_by_c_df = c_by_c_df.loc[keepers, :]
    else:
        c_by_c_df = None

    # Build nested dict of counts
    counts_by_context = {cnv : {} for cnv in 'DEL DUP'.split()}
    for cnv in counts_by_context.keys():
        for ptype in phenotypes:
            hpo = {'case' : case_hpo, 'control' : control_hpo}[ptype]
            counts_by_context[cnv][ptype] = {}
            for cohort in cohorts:
                if c_by_c_df is None:
                    counts_by_context[cnv][ptype][cohort] = 0
                else:
                    ridx = (c_by_c_df.cohort == cohort) & \
                           (c_by_c_df.hpo == hpo) & \
                           (c_by_c_df.cnv == cnv)
                    counts_by_context[cnv][ptype][cohort] = \
                        c_by_c_df.loc[ridx, 'cnvs'].values[0]

    return counts_by_context


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cohorts', required=True, help='Four-column tsv with ' +
                        'cohort name, n_cases, n_controls, and path to CNV BED.')
    parser.add_argument('--hpo-list', help='.txt file listing HPO terms to ' +
                        'consider for case samples.', required=True)
    parser.add_argument('--track-stats', required=True, help='Tsv of tracks ' + 
                        'with stats as produced by curate_track.py. Assumes first ' + 
                        'column is track name and last column is track path.')
    parser.add_argument('-F', '--frac-overlap', help='Minimum fraction of element ' +
                        'that must be overlapped by a CNV to be counted', 
                        type=float, default=10e-10)
    parser.add_argument('--control-hpo', default='HEALTHY_CONTROL', help='HPO code ' +
                        'to use for control CNV counts. [default: HEALTHY_CONTROL]')
    parser.add_argument('--norm-by-samplesize', action='store_true', help='Return ' +
                        'the number of samples with no CNVs as "_ref" category. ' +
                        '[default: returns number of CNVs not intersecting track]')
    parser.add_argument('--conditional-exclusion', help='BED of cohorts to ' + 
                        'exclude per 200kb window.')
    parser.add_argument('--counts-by-gene-context', help='.tsv of counts of ' +
                        'CNVs per phenotype by genic context. If provided, will ' +
                        'subtract coding CNV carriers from denominators to ' +
                        'better calibrate burden testing.')
    parser.add_argument('-o', '--outfile', help='Path to output BED file. ' +
                        '[default: stdout]', default='stdout')
    parser.add_argument('-z', '--gzip', action='store_true', help='Compress ' +
                        'output tsv with gzip.')
    args = parser.parse_args()

    # Open connections to output file
    if args.outfile in '- stdout /dev/stdout'.split():
        outfile = stdout
        gzip_out = False
    else:
        if path.splitext(args.outfile)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outfile_path = path.splitext(args.outfile)[0]
            gzip_out = True
        else:
            outfile_path = args.outfile
            gzip_out = args.gzip
        outfile = open(outfile_path, 'w')

    # Load list of case HPOs
    with open(args.hpo_list) as fin:
        case_hpos = [x.rstrip() for x in fin.readlines()]
        
    # Load CNVs per cohort and split by CNV type and phenotype
    cnvs, cnv_counts, cohort_n = load_cnvs(args.cohorts, case_hpos, [args.control_hpo])

    # Load conditional exclusion BED per cohort
    if args.conditional_exclusion is not None:
        xbt = pbt.BedTool(args.conditional_exclusion)
        xlist = {c : xbt.filter(lambda x: c in x.name).cut(range(3)).merge() for c in cnvs.keys()}
    else:
        xlist = {c : pbt.BedTool('', from_string=True) for c in cnvs.keys()}

    # Load samples to be subtracted from denominators
    counts_by_context = load_counts_by_context(args.counts_by_gene_context,
                                               cnvs.keys(), 'DEVELOPMENTAL',
                                               args.control_hpo)

    # Iterate over tracks and count CNVs for each metacohort
    with open(args.track_stats) as tsin:
        for track_vals in [x.rstrip().split('\t') for x in tsin.readlines()]:
            # Process header
            if '#' in track_vals[0]:
                out_header = track_vals[:-1] + ['cnv']
                for cohort in cnvs.keys():
                    for ptype in phenotypes:
                        out_header.append('_'.join([cohort, ptype, 'alt']))
                        out_header.append('_'.join([cohort, ptype, 'ref']))
                outfile.write('\t'.join(out_header) + '\n')
                continue

            # Process each track
            tname = str(track_vals[0])
            tbt = pbt.BedTool(str(track_vals[-1])).cut(range(3)).saveas()
            tvals_out_base = track_vals[:-1]
            for ctype in cnv_types:
                tvals_out_sub = tvals_out_base + [ctype]
                for cohort in cnvs.keys():
                    for ptype in phenotypes:
                        cbt = cnvs[cohort][ptype][ctype]
                        xbt = xlist[cohort]
                        tfilt = tbt.intersect(xbt, v=True)
                        hits = len(cbt.intersect(tfilt, u=True, F=args.frac_overlap))
                        tvals_out_sub.append(str(hits))
                        if args.norm_by_samplesize:
                            n_sub = counts_by_context[ctype][ptype][cohort]
                            nohits = cohort_n[cohort][ptype] - hits - n_sub
                            tvals_out_sub.append(str(nohits))
                        else:
                            nohits = cnv_counts[cohort][ptype][ctype] - hits
                            tvals_out_sub.append(str(nohits))
                outfile.write('\t'.join([str(x) for x in tvals_out_sub]) + '\n')

    # Close & compress output file (if optioned)
    outfile.close()
    if gzip_out:
        subprocess.run(['gzip', '-f', outfile_path])


if __name__ == '__main__':
    main()
