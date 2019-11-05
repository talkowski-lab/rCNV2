#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Refine significant sliding window associations
"""


from os import path
import gzip
import pandas
from numpy import nanmin, isnan, array, median, floor, ceil
import pybedtools as pbt
import csv
import statsmodels.api as sm
from scipy.stats import norm
import argparse


def get_bed_header(bedpath):
    """
    Get header from BED file
    """

    if path.splitext(bedpath)[1] in '.gz .bgz .bgzip'.split():
        header = gzip.GzipFile(bedpath).readline().decode('utf-8').rstrip('\n')
    else:
        header = open(bedpath).readline().rstrip('\n')

    return header


def load_pvalues(pvalues_bed, cnv_type, p_is_phred=False):
    """
    Import p-values matrix as dataframe
    """

    # Read pvalues_bed header to get phenotype keys
    pvals_header = get_bed_header(pvalues_bed)
    pval_cols = pvals_header.rstrip().split('\t')
    pval_cols = [x.replace('.' + cnv_type, '').replace('#', '') for x in pval_cols]

    # Format & transform pvalue matrix
    pvals_bt = pandas.read_table(pvalues_bed, names=pval_cols, comment='#')
    coords = pvals_bt.iloc[:,0:3]
    def _make_window_id(coords_i):
        return 'pval_window_{0}_{1}_{2}'.format(str(coords_i['chr']), 
                                                str(coords_i['start']), 
                                                str(coords_i['end']))
    window_ids = coords.apply(_make_window_id, axis=1)
    pvals = pvals_bt.iloc[:,3:]
    if p_is_phred:
        pvals = pvals.apply(lambda x: 10 ** -x, axis=0)

    pvals_out = pandas.concat([coords, window_ids, pvals], axis=1, ignore_index=True)
    pvals_out.columns = ['chr', 'start', 'end', 'id'] + pval_cols[3:]
    return pvals_out


def get_pvalues(region, pvals, reg_pval_int, cutoff):
    """
    Extract peak p-values and significant phenotypes for a query region
    """
    
    # Extract p-value windows within query region
    hits = reg_pval_int.filter(lambda x: x.chrom == region.chrom and \
                                         x.start == region.start and \
                                         x.end == region.end)
    hit_ids = [x[6] for x in hits]
    pval_hits = pvals[pvals['id'].isin(hit_ids)]

    # Find peak p-value for each phenotype
    p_peak = {}
    signif_hpos = []
    for pheno in list(pvals.columns[4:]):
        hpo = pheno.replace('HP', 'HP:')
        minp = nanmin(pval_hits[pheno])
        peak_ids = list(pval_hits[pval_hits[pheno] == minp]['id'])
        p_peak[hpo] = (minp, peak_ids)
        if not isnan(minp) and minp <= cutoff:
            signif_hpos.append(hpo)
    
    return p_peak, signif_hpos


def get_sentinel_window(p_peaks, pvals):
    """
    Parse dict of (p, [window_id]) tuples to find most significant window
    Selects left-most window in the case of ties
    """
    
    best_p = nanmin([p for (p, ids) in p_peaks.values()])
    sentinel_ids = [ids for (p, ids) in p_peaks.values() if p == best_p]
    sentinel_ids = list(set([i for s in sentinel_ids for i in s]))
    if len(sentinel_ids) > 1:
        sent_windows = pvals[pvals['id'].isin(sentinel_ids)]
        sentinel = list(sent_windows.sort_values(by='start')['id'])[0]
    else:
        sentinel = str(sentinel_ids[0])

    return sentinel
    

def load_regions(regions_in, pvals, samples_df, cnv_type, cutoff, prefix):
    """
    Read BED of regions to refine, and annotate with p-values
    """

    regions = {}

    regions_bt = pbt.BedTool(regions_in)

    pvals_bt = pbt.BedTool.from_dataframe(pvals.iloc[:,0:4])

    reg_pvals_int = regions_bt.intersect(pvals_bt, wa=True, wb=True, F=1.0).saveas()

    for r in regions_bt:
        rname = '{0}_region_{1}:{2}-{3}'.format(prefix, str(r.chrom),
                                                str(r.start), str(r.end))

        p_peaks, signif_hpos = get_pvalues(r, pvals, reg_pvals_int, cutoff)

        sentinel_id = get_sentinel_window(p_peaks, pvals)

        sentinel_pvals = pvals[pvals['id'].isin([sentinel_id])].iloc[:1, 4:].values.flatten().tolist()
        sent_pval_tuples = tuple(zip(pvals.columns[4:].tolist(), sentinel_pvals))
        sent_hpos = [hpo.replace('HP', 'HP:') for hpo, p in sent_pval_tuples if p <= cutoff]

        s_counts = samples_df[samples_df['hpo'].isin(sent_hpos)].loc[:, ['hpo', 'total']]
        s_counts.sort_values(by='total', inplace=True, ascending=False)
        sentinel_pheno = s_counts.iloc[0]['hpo']
        
        if rname not in regions.keys():
            regions[rname] = {'chr' : r.chrom, 
                              'start' : r.start, 
                              'end' : r.end, 
                              'region_id' : rname, 
                              'signif_hpos' : signif_hpos, 
                              'p_peaks' : p_peaks,
                              'sentinel_id' : sentinel_id, 
                              'sentinel_hpo' : sentinel_pheno}

    return regions


def intersect_cnvs(query_bt, cnvs, cnv_type, case_hpos, control_hpo, 
                   rename=None, frac=0):
    """
    Intersect a query region with CNVs from a subset of phenotypes
    Returns case CNVs and control CNVs separately
    """

    cnv_hits = pbt.BedTool(cnvs).intersect(query_bt, wa=True, F=frac)
    
    if cnv_type is not 'CNV':
        cnv_hits = cnv_hits.filter(lambda x: x[4] == cnv_type).saveas()

    if rename is not None:
        def _rename_cnv(cnv, rename):
            cnv.name = '_'.join([rename, cnv.name])
            return cnv
        cnv_hits = cnv_hits.each(_rename_cnv, rename=rename).saveas()
    
    def _pheno_filter(cnv, keep_phenos):
        if len([p for p in cnv[5].split(';') if p in keep_phenos]) > 0:
            return True
        else:
            return False
    
    case_cnvs = cnv_hits.filter(_pheno_filter, keep_phenos=case_hpos).saveas()
    control_cnvs = cnv_hits.filter(_pheno_filter, keep_phenos=[control_hpo]).saveas()

    return case_cnvs, control_cnvs


def sort_cnvs(cnv_bt, sent_mid):
    """
    Order case CNVs based on distance to sentinel midpoint
    Input as bt, output as pandas.dataframe
    """

    dists = cnv_bt.loc[:, 'start':'end'].apply(lambda x: abs(x - sent_mid), axis=1)
    cnv_bt = pandas.concat([cnv_bt, dists.apply(max, axis=1)], axis=1, ignore_index=True)
    cnv_bt.columns = ['chr', 'start', 'end', 'cnv_id', 'cnv', 'phenos', 'abs_min_dist']
    cnv_bt.sort_values(by='abs_min_dist', inplace=True)
    
    return cnv_bt


def sweeting_correction(n_control_ref, n_case_ref, n_control_alt, n_case_alt, cc_sum):
    """
    Apply empirical continuity correction to case & control CNV data
    Per Sweeting et al., Stat. Med., 2004 (section 3.3)
    """

    # Count number of carriers and non-carriers
    n_alt = sum(array(n_control_alt + n_case_alt))
    n_alt_bycohort = list(array(n_control_alt) + array(n_case_alt))
    n_ref = sum(array(n_control_ref + n_case_ref))

    # Require at least one CNV to be observed
    if n_alt > 0:
        nt = n_alt
        R = n_ref / n_alt

        # Pooled odds ratio estimate of all non-zero studies
        nz_idx = [i for i in range(0, len(n_case_ref)) if n_alt_bycohort[i] > 0]
        nz_case_alt = sum([n_case_alt[i] for i in nz_idx])
        nz_case_ref = sum([n_case_ref[i] for i in nz_idx])
        nz_case_odds = nz_case_alt / nz_case_ref
        nz_control_alt = sum([n_control_alt[i] for i in nz_idx])
        nz_control_ref = sum([n_control_ref[i] for i in nz_idx])
        nz_control_odds = nz_control_alt / nz_control_ref
        if nz_control_odds > 0:
            ohat = nz_case_odds / nz_control_odds
        else:
            nz_case_alt += (0.5 * len(n_case_ref))
            nz_case_ref += (0.5 * len(n_case_ref))
            nz_case_odds = nz_case_alt / nz_case_ref
            nz_control_alt += (0.5 * len(n_control_ref))
            nz_control_ref += (0.5 * len(n_control_ref))
            nz_control_odds = nz_control_alt / nz_control_ref
            ohat = nz_case_odds / nz_control_odds

        # Solve for kc & kt
        kc = R / (R + ohat)
        kt = ohat / (R + ohat)

        # Compute continuity corrections
        cor_case_alt = cc_sum * kt
        cor_case_ref = cc_sum * kc
        cor_control_alt = cc_sum * (nt + kt)
        cor_control_ref = cc_sum * ((nt * R) + kc)

        # Apply continuity corrections
        n_control_ref = [x + cor_control_ref for x in n_control_ref]
        n_case_ref = [x + cor_case_ref for x in n_case_ref]
        n_control_alt = [x + cor_control_alt for x in n_control_alt]
        n_case_alt = [x + cor_case_alt for x in n_case_alt]

    else:
        exit('No CNV carriers fed to Sweeting correction')

    return n_control_ref, n_case_ref, n_control_alt, n_case_alt


def mh_test(cohorts, n_control, n_case, n_control_alt, n_case_alt,
            empirical_continuity=True, cc_sum=0.01):
    """
    Calculate P-value using Mantel-Haenszel meta-analysis method
    """

    # Apply empirical continuity correction
    n_control_ref = list(array(n_control) - array(n_control_alt))
    n_case_ref = list(array(n_case) - array(n_case_alt))
    if empirical_continuity:
        n_control_ref, n_case_ref, n_control_alt, n_case_alt = \
            sweeting_correction(n_control_ref, n_case_ref, n_control_alt, 
                                n_case_alt, cc_sum)

    # Create list of 2x2 tables
    cohort_idxs = list(range(0, len(cohorts)))
    def _make_cont_table(idx):
        tb = array([[n_control_ref[idx], n_case_ref[idx]], 
                    [n_control_alt[idx], n_case_alt[idx]]])
        return tb
    meta_tabs = [_make_cont_table(i) for i in cohort_idxs]
    strat_tab = sm.stats.StratifiedTable(meta_tabs)

    # Calculate odds ratio, chisq, and p-value
    oddsratio = strat_tab.oddsratio_pooled
    chisq = strat_tab.test_null_odds().statistic
    pval = strat_tab.test_null_odds().pvalue
    if oddsratio < 1:
        pval = 1 - pval

    return oddsratio, chisq, pval


def get_min_case_cnvs(sentinel, case_cnvs, control_cnvs, cohorts, samples_df, 
                      case_hpo, control_hpo, cutoff):
    """
    Get minimum number of case CNVs required to achieve a significant result
    """

    def _count_total_samples_by_cohort(samples_df, hpo, cohorts):
        counts = [int(samples_df[samples_df['hpo'] == hpo][cohort].values) \
                 for cohort in cohorts]
        return counts

    def _count_cnvs_by_cohort(cnvs, cohorts):
        cnv_ids = [x.split('_')[0] for x in cnvs['cnv_id'].values]
        counts = [len([x for x in cnv_ids if x == cohort]) for cohort in cohorts]
        return counts

    # Get all control count data
    n_control = _count_total_samples_by_cohort(samples_df, control_hpo, cohorts)
    n_control_alt = _count_cnvs_by_cohort(control_cnvs, cohorts)

    # Get baseline case count data
    n_case = _count_total_samples_by_cohort(samples_df, case_hpo, cohorts)
    n_case_alt_max = _count_cnvs_by_cohort(case_cnvs, cohorts)
    sum_n_case_alt_max = sum(array(n_case_alt_max))

    # Increment case CNVs one at a time until significant p-value is reached
    k = 0
    best_p = 1
    while best_p > cutoff:
        if k < sum_n_case_alt_max:
            k += 1
            n_case_alt = _count_cnvs_by_cohort(case_cnvs.iloc[0:k, :], cohorts)
            oddsratio, chisq, new_p = mh_test(cohorts, n_control, n_case, 
                                              n_control_alt, n_case_alt)
            best_p = min([best_p, new_p])
            # Helpful debugging print out
            # print('\t'.join([str(k), str(oddsratio), str(chisq), str(new_p), str(best_p)]))
        else:
            err = 'All case CNVs exhausted to refine {0}, and genome-wide ' + \
                  'significance was never achieved.'
            exit(err.format(sentinel))
    
    return k


def refine_sentinel(regions, rid, pvals, samples_df, cnvs_table, 
                    cnv_type, control_hpo, cutoff, resolution=10000):
    """
    Define minimum critical region for a single sentiel window
    """

    # Get sentinel information
    sentinel = regions[rid]['sentinel_id']
    sent_w = pvals[pvals['id'].isin([sentinel])]
    sent_bt = pbt.BedTool.from_dataframe(sent_w.iloc[:,0:3])
    sent_mid = (sent_bt[0].start + sent_bt[0].end) / 2
    sent_hpo = regions[rid]['sentinel_hpo']

    # Intersect sentinel window with CNVs from sentinel phenotype
    cohorts = [s.split('\t')[0] for s in open(cnvs_table).readlines()]
    phenos = regions[rid]['signif_hpos']
    cnvs = {}
    with open(cnvs_table) as cnv_fin:
        reader = csv.reader(cnv_fin, delimiter='\t')
        for cname, cnvpath in reader:
            case_bt, control_bt = intersect_cnvs(sent_bt, cnvpath, cnv_type, 
                                                 [sent_hpo], control_hpo,
                                                 rename=cname, frac=0.5)
            cnvs[cname] = {'case_cnvs' : case_bt,
                           'control_cnvs' : control_bt}

    # Pool and sort all CNVs across metacohorts for cases and controls
    case_cnvs_l = [x['case_cnvs'] for x in cnvs.values()]
    case_cnvs = pandas.concat([bt.to_dataframe() for bt in case_cnvs_l if len(bt) > 0], 
                              axis=0, ignore_index=True)
    case_cnvs = sort_cnvs(case_cnvs, sent_mid)
    control_cnvs_l = [x['control_cnvs'] for x in cnvs.values()]
    control_cnvs = pandas.concat([bt.to_dataframe() for bt in control_cnvs_l if len(bt) > 0], 
                                 axis=0, ignore_index=True)
    control_cnvs = sort_cnvs(control_cnvs, sent_mid)

    # Get minimum set of case CNVs required for genome-wide significant result
    k_case_cnvs = get_min_case_cnvs(sentinel, case_cnvs, control_cnvs, cohorts, 
                                    samples_df, sent_hpo, control_hpo, cutoff)
    min_case_cnvs = case_cnvs.iloc[0:k_case_cnvs, :]
    def _round_to(x, to):
        to * round(x / to)
    refined_chrom = regions[rid]['chr']
    refined_start = int(resolution * floor(median(min_case_cnvs['start']) / resolution))
    refined_end = int(resolution * ceil(median(min_case_cnvs['end']) / resolution))
    refined_size = int(refined_end - refined_start)

    # Print result
    msg = 'Region {0}, originally {1}kb, was refined to {2}kb critical region ' + \
          '({3}:{4}-{5}) using {6} of {7} case CNVs for {8}.\n'
    print(msg.format(rid, str(int(regions[rid]['end']-regions[rid]['start'])/1000),
                     str(refined_size/1000), refined_chrom, str(refined_start), 
                     str(refined_end), str(k_case_cnvs), str(len(case_cnvs)), 
                     sent_hpo))


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('regions', help='BED file of regions to refine.')
    parser.add_argument('cnvs_table', help='Two-column TSV file specifying ' +
                        'name and path to CNV BED file for each cohort.')
    parser.add_argument('pvalues', help='BED file of original sliding windows ' +
                         'and p-values for each phenotype.')
    parser.add_argument('pheno_table', help='table with counts of samples per ' +
                        'HPO term per cohort.')
    parser.add_argument('--cnv-type', help='type of CNVs to evaluate. [default: ' + 
                        'use all CNVs]', choices=['CNV', 'DEL', 'DUP'], 
                        default='CNV')
    parser.add_argument('--cutoff', help='P-value of significance threshold. ' + 
                        '[default: 10e-8]', default=10e-8, type=float)
    parser.add_argument('--p-is-phred', help='supplied P-values are Phred-scaled ' +
                        '(-log10[P]). [default: False]', default=False, 
                        action='store_true')
    parser.add_argument('--control-hpo', default='HEALTHY_CONTROL', help='HPO code ' +
                        'to use for control CNV counts. [default: HEALTHY_CONTROL]')
    parser.add_argument('--resolution', help='Round refined region coordinates ' + 
                        'up to the nearest X bp. [default: 10000]', default=10000, 
                        type=int)
    parser.add_argument('--prefix', help='string to append to all refined regions. ' +
                        '[default: CNV type (specified above)]', default=None)

    args = parser.parse_args()

    # Reset args
    if args.prefix is None:
        args.prefix = args.cnv_type

    # Load p-value matrix
    pvals = load_pvalues(args.pvalues, args.cnv_type, args.p_is_phred)

    # Load case/control sample counts
    samples_df = pandas.read_table(args.pheno_table, header=0, comment=None)
    samples_df.columns = [x.replace('#', '').lower() for x in samples_df.columns]

    # Read original list of regions to refine
    regions = load_regions(args.regions, pvals, samples_df, args.cnv_type, 
                           args.cutoff, args.prefix)

    # Refine regions one at a time
    for rid in regions.keys():
        refine_sentinel(regions, rid, pvals, samples_df, args.cnvs_table, 
                        args.cnv_type, args.control_hpo, args.cutoff, 
                        args.resolution)


if __name__ == '__main__':
    main()
