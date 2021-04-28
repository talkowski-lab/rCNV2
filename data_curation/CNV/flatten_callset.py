#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Flatten a single cohort's CNV callset
"""


import gzip
import csv
import pandas as pd
import random
import argparse
from sys import stdout
from os.path import splitext
from athena.utils.misc import bgzip


def load_hpo_pairs(tsv_in):
    """
    Load precomputed HPO pairs .tsv
    """

    elig_hpos = set()
    overlap_probs = {}

    if tsv_in.endswith('.gz'):
        reader = csv.reader(gzip.open(tsv_in, mode='rt'), delimiter='\t')
    else:
        reader = csv.reader(open(tsv_in, mode='rt'), delimiter='\t')

    for hpo1, hpo2, n1, n2, n12 in reader:
        if hpo1.startswith('#'):
            continue
        
        # Add terms to list of all eligible terms
        elig_hpos.add(hpo1)
        elig_hpos.add(hpo2)

        # Compute overlap probability
        ovr = int(n12) / int(n1)
        if hpo1 not in overlap_probs.keys():
            overlap_probs[hpo1] = {}
        if hpo2 not in overlap_probs[hpo1].keys():
            overlap_probs[hpo1][hpo2] = 0
        overlap_probs[hpo1][hpo2] = ovr

    # Return values as nested dict
    return {'hpos' : elig_hpos, 'overlap_probs' : overlap_probs}


def load_input_bed(infile, elig_hpos, cnv_type_col=4):
    """
    Load input BED file as pd.dataframe
    """

    # Load & sort BED file
    cnv_df = pd.read_csv(infile, sep='\t')
    cnv_df.sort_values(cnv_df.columns[:3].tolist(), inplace=True)

    # Subset HPO columns to only those present also in --hpo-pairs
    hpos_in_header = list(cnv_df.columns[cnv_type_col:])
    for hpo in hpos_in_header:
        if hpo == 'TOTAL':
            continue
        if hpo not in elig_hpos:
            cnv_df.drop(hpo, axis=1, inplace=True)

    return cnv_df


def get_hpo_probability(new_hpo, existing_hpos):
    """
    Get the probability of observing 
    """

    if len(existing_hpos) == 0:
        return 1
    else:
        import pdb; pdb.set_trace()


def assign_case_hpos(hpo_counts, overlap_probs, case_n):
    """
    Randomizes synthetic assignment of case HPOs
    """

    cases = [set() for i in range(case_n)]

    hpo_counts = {k : v for k, v in sorted(hpo_counts.items(), reverse=True, 
                                           key=lambda item: item[1])}

    for hpo, n in hpo_counts.items():
        while n > 0:
            for i in range(case_n):
                if hpo in cases[i]:
                    continue
                else:
                    p_add = get_hpo_probability(hpo, cases[i])
                    if random.random() < p_add:
                        cases[i].add(hpo)
                        n -= 1

    import pdb; pdb.set_trace()


def flatten_cnv(vals, overlap_probs, prefix, k, control_hpo='HEALTHY_CONTROL', 
                cnv_type_col=4):
    """
    Flattens a single row from input BED file
    Returns a BED-styled string of flattened CNV calls
    """

    cnvs = []
    cnv_fmt = '{}\t{}\t{}\t{}_{}\t{}\t{}'

    coords = vals[:3].tolist()
    cnvtype = vals[cnv_type_col-1]

    total_n = vals.TOTAL
    control_n = vals[control_hpo]
    case_n = total_n - control_n

    for i in range(control_n):
        k += 1
        newcnv = cnv_fmt.format(*coords, prefix, k, cnvtype, control_hpo)
        cnvs.append(newcnv)

    hpo_counts = {k : v for k, v in vals[cnv_type_col:].to_dict().items() \
                  if k not in ['TOTAL', control_hpo] and v > 0}
    case_hpos = assign_case_hpos(hpo_counts, overlap_probs, case_n)
    for i in range(case_n):
        k += 1
        newcnv = cnv_fmt.format(*coords, prefix, k, cnvtype, case_hpos[i])
        cnvs.append(newcnv)

    return cnvs


def main():
    """
    Main block
    """

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infile', help='Input BED (supports "stdin")')
    parser.add_argument('-o', '--outfile', help='Output BED [default: stdout]',
                        default='stdout')
    parser.add_argument('-p', '--hpo-pairs', help='.tsv of precomputed sample ' +
                        'counts for all HPO pairs', required=True)
    parser.add_argument('-c', '--cohort', required=True, help='Cohort name')
    parser.add_argument('--control-hpo', default='HEALTHY_CONTROL',
                        help='HPO term to consider healthy control [default: "HEALTHY_CONTROL"]')
    parser.add_argument('--cnv-type-column', type=int, default=4,
                        help='Column number specifying CNV type [default: 4]',)
    parser.add_argument('--seed', type=int, default=2021,
                        help='Random seed [default: 2021]',)
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output BED with bgzip.')
    args = parser.parse_args()

    # Load HPO pairs & compute co-occurrence probabilities
    hpo_dict = load_hpo_pairs(args.hpo_pairs)

    # Load input BED
    cnv_df = load_input_bed(args.infile, hpo_dict['hpos'], args.cnv_type_column)

    # Open connection to outfile & write header
    if args.outfile in '- stdout /dev/stdout'.split():
        outpath = 'stdout'
        outfile = stdout
    else:
        if args.outfile.endswith('gz'):
            outpath = splitext(args.outfile)[0]
        else:
            outpath = args.outfile
        outfile = open(outpath, 'w')
    header = '\t'.join('#chr start end name cnv pheno'.split()) + '\n'
    outfile.write(header)

    # Process each row in cnv_df
    k = 0
    random.seed(args.seed)
    for row in cnv_df.iterrows():
        newcnvs = flatten_cnv(row[1], hpo_dict['overlap_probs'], args.cohort, k,
                              args.control_hpo, args.cnv_type_column)
        n_new = len(newcnvs)
        k += n_new
        if n_new > 0:
            outfile.write('\n'.join(newcnvs) + '\n')
    outfile.close()

    # Bgzip --outfile, if optioned
    if args.bgzip and outpath != 'stdout':
        bgzip(outpath)


if __name__ == '__main__':
    main()

