#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Calculate inverse variance-weighted effect sizes and max p-values per phenotype per segment
"""


import pandas as pd
import re
import pybedtools as pbt
import gzip
import numpy as np
import warnings
from pysam import TabixFile
from scipy.stats import norm
import argparse
from sys import stdout


np.seterr(all='raise')


def load_segs(bed_in):
    """
    Loads all rCNV-associated segments
    Returns:
        1. pd.DataFrame of all segment info
        2. dict mapping segment IDs to pbt.BedTool of credible intervals
    """

    segs_df = pd.read_csv(bed_in, sep='\t').rename(columns={'#chr' : 'chr'})
    seg_ids = segs_df.iloc[:, 3].tolist()

    cs_dict = {}
    for seg_id in seg_ids:
        ints = segs_df.cred_interval_coords[segs_df.iloc[:, 3] == seg_id].tolist()[0]
        ints_bt = pbt.BedTool(re.sub(':|-', '\t', ints.replace(';', '\n')), 
                              from_string=True)
        cs_dict[seg_id] = ints_bt

    return segs_df, cs_dict


def tabix_region(bedpath, querybt):
    """
    Uses tabix to extract all windows spanning in querybt
    Returns: pbt.BedTool of all windows in query
    """

    # Format query
    chrom = querybt[0].chrom
    start = str(np.nanmin(querybt.cut(range(3)).to_dataframe().start))
    end = str(np.nanmin(querybt.cut(range(3)).to_dataframe().end))
    region = '{}:{}-{}'.format(chrom, start, end)

    # Extract all windows
    return pbt.BedTool('\n'.join([x for x in TabixFile(bedpath).fetch(region)]),
                       from_string=True)


def ci2se(ci):
    """
    Converts a tuple of (lower, upper) confidence interval bounds to standard error
    """

    ci = sorted(ci)

    return (ci[1] - ci[0]) / (2 * 1.96)


def iv_mean(values, variances, conf=0.95):
    """
    Returns inverse-variance weighted mean of values and conf% confidence interval
    """

    weights = [1 / v for v in variances]

    wsum = np.nansum(weights)

    numerator = np.nansum([x / v for x, v in zip(values, variances)])

    ivm = numerator / wsum

    pooled_se = np.sqrt(1 / wsum)

    ci_dist = norm.ppf(conf) * pooled_se

    return ivm, (ivm - ci_dist, ivm + ci_dist)


def calc_all_effects(cs_dict, phenos, ssdir):
    """
    Compute inverse variance-weighted mean effect size and max p-value for each phenotype for each segment
    Returns: flattened pd.DataFrame (one row per segment-phenotype-CNV pair)
    """

    outcols = 'region_id hpo cnv lnor lnor_lower lnor_upper ' + \
              'pvalue pvalue_secondary peak_lnor'
    fx_df = pd.DataFrame(columns=outcols.split())

    for seg_id, seg_bt in cs_dict.items():
        for pheno, hpo in phenos:
            for cnv in 'DEL DUP'.split():
                # Extracts sumstats for windows completely overlapped by credible interval(s)
                ssbedpath = \
                    '{}/{}.rCNV.{}.sliding_window.meta_analysis.stats.bed.gz'.\
                        format(ssdir, pheno, cnv)
                ss_colnames = \
                    gzip.open(ssbedpath, 'rt').readline().rstrip().replace('#', '').split('\t')
                ssbt = tabix_region(ssbedpath, seg_bt)
                wdf = ssbt.intersect(seg_bt, f=0.99, u=True, wa=True).\
                          to_dataframe(names=ss_colnames)

                # Computes inverse variance-weighted mean of effect sizes
                lnors = wdf.meta_lnOR.to_numpy()
                if np.all(np.isnan(lnors)):
                    iv_lnor, iv_lnor_lower, iv_lnor_upper = 0, np.nan, np.nan
                else:
                    lnor_cis = wdf.loc[:, ['meta_lnOR_lower', 'meta_lnOR_upper']].\
                                   to_records(index=False)
                    lnor_vars = [ci2se(ci) for ci in lnor_cis]
                    iv_lnor, (iv_lnor_lower, iv_lnor_upper) = iv_mean(lnors, lnor_vars)

                # Computes peak effect size
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore', 'All-NaN slice encountered')
                    best_lnor = np.nanmax(lnors)

                # Gathers peak primary and secondary p-values
                if np.all(np.isnan(wdf.meta_neg_log10_p)):
                    best_p_primary = 0
                else:
                    best_p_primary = np.nanmax(wdf.meta_neg_log10_p)
                if np.all(np.isnan(wdf.meta_neg_log10_p_secondary)):
                    best_p_secondary = 0
                else:
                    best_p_secondary = np.nanmax(wdf.meta_neg_log10_p_secondary)

                # Add OR estimate and peak pvalues to fx_df
                newvals = [seg_id, hpo, cnv, iv_lnor, iv_lnor_lower, iv_lnor_upper,
                           best_p_primary, best_p_secondary, best_lnor]
                fx_df = fx_df.append(pd.Series(newvals, index=fx_df.columns),
                                               ignore_index=True)

    return fx_df


def main():
    """
    Main command-line block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('loci', help='BED of final significant loci.')
    parser.add_argument('phenos', help='tsv of prefix, hpo pairs for all phenotypes.')
    parser.add_argument('sumstats_dir', help='path to directory containing ' +
                        'meta-analysis summary statistics for all phenotypes and ' +
                        'CNV types.')
    parser.add_argument('-o', '--outfile', help='Output tsv of summary stats per ' +
                        'segment per phenotype per CNV. [default: stdout]',
                        default='stdout')
    args = parser.parse_args()
    
    # Open connections to output files
    if args.outfile in 'stdout - /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Load segments
    segs_df, cs_dict = load_segs(args.loci)

    # Load phenotypes as (prefix, hpo) tuple pairs
    phenos = [tuple(x.rstrip().split('\t')) for x in open(args.phenos, 'r').readlines()]

    # Calculate effects for each segment for each phenotype
    effects_df = calc_all_effects(cs_dict, phenos, args.sumstats_dir)

    # Write results to outfile
    effects_df.to_csv(outfile, sep='\t', na_rep='NA', index=False)


if __name__ == '__main__':
    main()

