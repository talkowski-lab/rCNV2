#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
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


def load_cnvs_single(cnvpath, hpo, cnv_type):
    """
    Loads & filters CNVs for a single HPO & CNV type
    Returns : pbt.BedTool
    """

    cnvs = pbt.BedTool(cnvpath)

    def _hpo_filter(feature, hpo):
        return hpo in feature[-1]

    def _cnv_filter(feature, cnv_type):
        return feature[4] == cnv_type

    return cnvs.filter(_hpo_filter, hpo).\
                filter(_cnv_filter, cnv_type).\
                saveas()


def load_cnvs(cohorts_in, case_hpo, control_hpo):
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
    pdict = {phenotypes[0] : case_hpo, phenotypes[1] : control_hpo}

    with open(cohorts_in) as fin:
        for cohort, n_case, n_control, cnvpath in csv.reader(fin, delimiter='\t'):
            cnvs[cohort] = {ptype : {ctype : load_cnvs_single(cnvpath, pdict[ptype], ctype) \
                                     for ctype in cnv_types} for ptype in phenotypes}
            cnv_counts[cohort] = {ptype : {ctype : len(cnvs[cohort][ptype][ctype]) \
                                  for ctype in cnv_types} for ptype in phenotypes}
            cohort_n[cohort] = {phenotypes[0] : int(n_case), 
                                phenotypes[1] : int(n_control)}

    return cnvs, cnv_counts, cohort_n


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
    parser.add_argument('--track-stats', required=True, help='Tsv of tracks ' + 
                        'with stats as produced by curate_track.py. Assumes first ' + 
                        'column is track name and last column is track path.')
    parser.add_argument('-F', '--frac-overlap', help='Minimum fraction of element ' +
                        'that must be overlapped by a CNV to be counted', 
                        type=float, default=10e-10)
    parser.add_argument('--case-hpo', help='HPO term to consider for case samples. ' +
                        '[default: HP:0000118]', default='HP:0000118')
    parser.add_argument('--control-hpo', default='HEALTHY_CONTROL', help='HPO code ' +
                        'to use for control CNV counts. [default: HEALTHY_CONTROL]')
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
        

    # Load CNVs per cohort and split by CNV type and phenotype
    cnvs, cnv_counts, cohort_n = load_cnvs(args.cohorts, args.case_hpo, args.control_hpo)

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
            tpath = str(track_vals[-1])
            tvals_out_base = track_vals[:-1]
            for ctype in cnv_types:
                tvals_out_sub = tvals_out_base + [ctype]
                for cohort in cnvs.keys():
                    for ptype in phenotypes:
                        cbt = cnvs[cohort][ptype][ctype]
                        hits = len(cbt.intersect(tpath, u=True, F=args.frac_overlap))
                        tvals_out_sub.append(str(hits))
                        nohits = cohort_n[cohort][ptype] - hits
                        tvals_out_sub.append(str(nohits))
                outfile.write('\t'.join([str(x) for x in tvals_out_sub]) + '\n')

    # Close & compress output file (if optioned)
    outfile.close()
    if gzip_out:
        subprocess.run(['gzip', '-f', outfile_path])


if __name__ == '__main__':
    main()
