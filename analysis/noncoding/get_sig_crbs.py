#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Extract significant CRBs from meta-analysis summary statistics
"""


import csv
import pandas as pd
import numpy as np
import argparse
from sys import stdout
from os import path
import pybedtools as pbt
import subprocess


def read_sumstats_as_dict(ss_tsv):
    """
    Load sumstats tsv input as {hpo : str}
    """

    ss = {}

    with open(ss_tsv) as fin:
        for hpo, ssbed in csv.reader(fin, delimiter='\t'):
            if '#' in hpo:
                continue
            else:
                ss[hpo] = ssbed

    return ss


def process_ss_dict(ss_dict, primary_p=10e-8, secondary_p=0.05, n_nom=1, 
                    secondary_or_nom=False):
    """
    Loads & processes sumstats per phenotype to extract significant CRBs
    Returns: {crb : [hpos]} of significant CRBs
    """

    sig_crbs = {}

    for hpo, bedpath in ss_dict.items():

        ss = pd.read_csv(bedpath, sep='\t')

        sig_primary = ss['meta_phred_p'] >= -np.log10(primary_p)
        sig_secondary = ss['meta_phred_p_secondary'] >= -np.log10(secondary_p)
        sig_nom = ss['n_nominal_cohorts'] >= n_nom

        if secondary_or_nom:
            sig = (sig_primary & (sig_secondary | sig_nom))
        else:
            sig = (sig_primary & sig_secondary & sig_nom)

        if sum(sig) > 0:
            for crb in ss['crb_id'][sig].tolist():
                if crb in sig_crbs.keys():
                    sig_crbs[crb].append(hpo)
                else:
                    sig_crbs[crb] = [hpo]

    return sig_crbs


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sumstats', required=True, help='Two-column tsv of ' +
                        'HPO & path to meta-analysis summary statistics BED file.')
    parser.add_argument('--coding-sumstats', help='Two-column tsv of HPO & path ' +
                        'to meta-analysis summary statistics BED file for analysis ' +
                        'from all CNVs (including coding CNVs).')
    parser.add_argument('--primary-p', default=1e-8, type=float, 
                        help='Minimum primary P-value to consider significant.')
    parser.add_argument('--coding-p', default=1e-8, type=float, 
                        help='Minimum primary P-value from --coding-sumstats to ' +
                        'consider significant.')
    parser.add_argument('--secondary-p', default=0.05, type=float, 
                        help='Minimum secondary P-value to consider significant.')
    parser.add_argument('--n-nominal', default=1, type=int, 
                        help='Minimum number of nominally significant cohorts.')
    parser.add_argument('--secondary-or-nominal', action='store_true', 
                        help='Require significant CRBs to satisfy either ' +
                        '--secondary-p or --n-nominal.')
    parser.add_argument('--cnv', default='N/S', help='CNV type.')
    parser.add_argument('-o', '--outfile', default='stdout', help='Path to ' +
                        'output BED.')
    parser.add_argument('-z', '--bgzip', action='store_true', help='Compress ' +
                        'output BED with bgzip.')
    args = parser.parse_args()

    # Open connection to outfile
    if args.outfile in 'stdout /dev/stdout -'.split():
        outfile = stdout
        bgzip = False
    else:
        if path.splitext(args.outfile)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outfile_path = path.splitext(args.outfile)[0]
            bgzip = True
        else:
            outfile_path = args.outfile
            bgzip = args.bgzip
        outfile = open(outfile_path, 'w')

    # Read sumstats
    ss = read_sumstats_as_dict(args.sumstats)
    if args.coding_sumstats is not None:
        css = read_sumstats_as_dict(args.coding_sumstats)
    else:
        css = None

    # Load CRB coordinates as pd.DataFrame
    crb_df = pbt.BedTool(list(ss.values())[0]).cut(range(4)).\
                 to_dataframe(names='chrom start end crb_id'.split())

    # Get list of significant CRBs
    sig_crbs = process_ss_dict(ss, args.primary_p, args.secondary_p, 
                               args.n_nominal, args.secondary_or_nominal)
    if css is not None:
        coding_sig_crbs = process_ss_dict(css, args.coding_p, 1, 0, True)
        sig_ids = list(sig_crbs.keys())
        for crb in sig_ids:
            new_hpos = [h for h in sig_crbs[crb] if h in coding_sig_crbs.get(crb, [])]
            if len(new_hpos) == 0:
                sig_crbs.pop(crb)
            else:
                sig_crbs[crb] = new_hpos

    # Write sig CRBs to file
    outfile.write('\t'.join('#chrom start end crb_id cnv hpos'.split()) + '\n')
    for crb, hpos in sig_crbs.items():
        coords = crb_df[crb_df['crb_id'] == crb].values.tolist()[0]
        out_vals = [str(x) for x in coords] + [args.cnv, ';'.join(hpos)]
        outfile.write('\t'.join(out_vals) + '\n')
    outfile.close()
    if bgzip:
        subprocess.run(['bgzip', '-f', outfile_path])


if __name__ == '__main__':
	main()


