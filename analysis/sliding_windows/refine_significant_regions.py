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
import pandas as pd
import pybedtools as pbt
import csv
import numpy as np
import statsmodels.api as sm
import scipy.stats as sstats
import tempfile
import subprocess
import argparse
from sys import stdout
from operator import itemgetter


def get_bed_header(bedpath):
    """
    Get header from BED file
    """

    if path.splitext(bedpath)[1] in '.gz .bgz .bgzip'.split():
        header = gzip.GzipFile(bedpath).readline().decode('utf-8').rstrip('\n')
    else:
        header = open(bedpath).readline().rstrip('\n')

    return header


def _make_window_id(coords):
    """
    Create window ID based on coordinates
    """

    return 'window_{0}_{1}_{2}'.format(str(coords['chr']), 
                                       str(coords['start']), 
                                       str(coords['end']))


def load_sig_df(sig_df, cnv_type):
    """
    Load a matrix of True/False significance indications per window per phenotype
    """

    # Read sig_df header to get phenotype keys
    sig_header = get_bed_header(sig_df)
    sig_cols = sig_header.rstrip().split('\t')
    sig_cols = [x.replace('.' + cnv_type, '').replace('#', '').replace('HP', 'HP:') \
                for x in sig_cols]

    # Read sig matrix and add window IDs
    sig_df = pd.read_table(sig_df, names=sig_cols, comment='#')
    coords = sig_df.iloc[:, 0:3]
    window_ids = coords.apply(_make_window_id, axis=1)
    sig_labs = sig_df.iloc[:, 3:]
    sig_out = pd.concat([coords, window_ids, sig_labs], axis=1, ignore_index=True)
    sig_out.columns = ['chr', 'start', 'end', 'id'] + sig_cols[3:]

    # Subset to BedTool of windows significant in at least one phenotype
    sig_windows = sig_out[sig_out.iloc[:, 4:].apply(any, axis=1)].iloc[:, 0:4]
    sig_bt = pbt.BedTool.from_dataframe(sig_windows)

    return sig_out, sig_bt


def load_pvalues(pvalues_bed, sig_bt, cnv_type, p_is_phred=False):
    """
    Import p-values matrix as dataframe
    """

    # Read pvalues_bed header to get phenotype keys
    pvals_header = get_bed_header(pvalues_bed)
    pval_cols = pvals_header.rstrip().split('\t')
    pval_cols = [x.replace('.' + cnv_type, '').replace('#', '').replace('HP', 'HP:') \
                 for x in pval_cols]

    # Format & transform pvalue matrix
    pvals_bt = pbt.BedTool(pvalues_bed).intersect(sig_bt, wa=True, u=True)
    pvals_df = pvals_bt.to_dataframe(names=pval_cols)
    coords = pvals_df.iloc[:,0:3]
    window_ids = coords.apply(_make_window_id, axis=1)
    pvals = pvals_df.iloc[:,3:]
    if p_is_phred:
        pvals = pvals.apply(lambda x: 10 ** -x, axis=0)

    pvals_out = pd.concat([coords, window_ids, pvals], axis=1, ignore_index=True)
    pvals_out.columns = ['chr', 'start', 'end', 'id'] + pval_cols[3:]

    # Restrict to p-values of significant windows
    sig_ids = list(sig_bt.to_dataframe().iloc[:, 3].values.flatten())

    return pvals_out[pvals_out['id'].isin(sig_ids)]


def load_cnvs(cnv_bed, cnv_type, whitelist=None, rename=None):
    """
    Load & filter CNVs from BED file into pbt.BedTool
    """

    cnvs = pbt.BedTool(cnv_bed)
    
    if cnv_type is not 'CNV':
        cnvs = cnvs.filter(lambda x: x[4] == cnv_type)

    if whitelist is not None:
        cnvs = cnvs.intersect(whitelist, wa=True, u=True)

    if rename is not None:
        def _rename_cnv(cnv, rename):
            cnv.name = '_'.join([rename, cnv.name])
            return cnv
        cnvs = cnvs.each(_rename_cnv, rename=rename)

    cnvs = cnvs.saveas()

    return cnvs


def load_hpos(hpo_file):
    """
    Returns list of lists of HPO terms per sample
    """

    hpos = []

    with open(hpo_file) as fin:
        reader = csv.reader(fin, delimiter='\t')
        for sample, terms in reader:
            hpos.append(terms.split(';'))

    return hpos


def load_cohort_info(info_file, regions_bt, cnv_type):
    """
    Read metadata per cohort and save as dict
    """

    cohort_info = {}

    with open(info_file) as fin:
        reader = csv.reader(fin, delimiter='\t')

        for cohort, cnv_bed, hpo_file in reader:
            cnvs = load_cnvs(cnv_bed, cnv_type, regions_bt, cohort)
            hpos = load_hpos(hpo_file)
            cohort_info[cohort] = {'cnvs' : cnvs,
                                   'hpos' : hpos}

    return cohort_info


def load_p_cutoffs(p_cutoffs_in):
    """
    Read a dictionary of p-value cutoffs per phenotype
    """

    p_cutoffs = {}

    with open(p_cutoffs_in) as tsvin:
        reader = csv.reader(tsvin, delimiter='\t')
        for hpo, cutoff in reader:
            if hpo.startswith('#'):
                continue
            else:
                p_cutoffs[hpo.replace('HP', 'HP:')] = float(cutoff)

    return p_cutoffs


def get_sentinel_window(sig_windows):
    """
    Parse list of dicts of window : [(hpo1, p1), (hpo2, p2)] tuples to find
    the most significant window. Selects left-most window in the case of ties.
    """
    
    best_p = np.nanmin([p for (p, ids) in p_peaks.values()])
    sentinel_ids = [ids for (p, ids) in p_peaks.values() if p == best_p]
    sentinel_ids = list(set([i for s in sentinel_ids for i in s]))
    if len(sentinel_ids) > 1:
        sent_windows = pvals[pvals['id'].isin(sentinel_ids)]
        sentinel = list(sent_windows.sort_values(by='start')['id'])[0]
    else:
        sentinel = str(sentinel_ids[0])

    return sentinel


def get_sig_windows(rbt, sig_df, sig_bt, pvals, pad_sentinel):
    """
    Returns a dict of significant windows overlapping a query BedTool (rbt)

    Annotates each window with all significant phenotypes (if any) and their
    associated p-values
    """

    all_hpos = list(sig_df.columns[4:])

    sig_windows = {}

    # Get window IDs of significant windows overlapping rbt
    sig_hit_bt = sig_bt.intersect(rbt, wa=True, u=True).saveas()
    if len(sig_hit_bt) > 0:
        sig_hit_df = sig_bt.intersect(rbt, wa=True, u=True).to_dataframe()
        sig_wids = list(sig_hit_df['name'].values.flatten())

        # Annotate each window with tuples of sigificant hpos and pvalues
        for wid in sig_wids:
            w_sig_bool = list(sig_df[sig_df['id'] == wid].values.flatten()[4:])
            w_hpo_sig = [hpo for hpo, sig in tuple(zip(all_hpos, w_sig_bool)) if sig]
            w_hpo_sig_p = list(pvals[pvals['id'] == wid][w_hpo_sig].values.flatten())
            sig_windows[wid] = tuple(zip(w_hpo_sig, w_hpo_sig_p))
        
        # Get list of all significant hpos
        all_sig_hpos = list(set([hpo for item in sig_windows.values() \
                                 for hpo, p in item]))

        # Get window ID and p-value of most significant bin
        # In the case of multiple, equally significant bins, take the first encountered
        # (Note: this will always be the left-most in coordinate space)
        best_p = np.nanmin([p for item in sig_windows.values() for hpo, p in item])
        for wdat in sig_windows.items():
            p_matches = [hpo for hpo, p in wdat[1] if p == best_p]
            if len(p_matches) > 0:
                best_wid = wdat[0]
                best_hpo = p_matches[0]
                break

        # Get list of all hpos with significant association within ± pad_sentinel
        best_w_bt = pbt.BedTool.from_dataframe(sig_df[sig_df['id'] == best_wid])
        def _pad_window(feature, pad):
            feature['start'] -= pad
            feature['end'] += pad
            return feature
        best_w_bt = best_w_bt.each(_pad_window, pad=pad_sentinel)
        other_best_ws = sig_bt.intersect(best_w_bt, wa=True, u=True).to_dataframe()
        other_best_wids = (other_best_ws.iloc[:, 3].values.flatten())
        all_best_hpos = []
        for wdat in sig_windows.items():
            if wdat[0] in other_best_wids:
                all_best_hpos = all_best_hpos + [hpo for hpo, p in wdat[1]]

        sentinel = {'window_id' : best_wid, 
                    'best_hpo' : best_hpo, 
                    'all_hpos' : list(set(all_best_hpos)),
                    'pvalue' : best_p}

    else:
        all_sig_hpos = []
        sentinel = None

    return sig_windows, all_sig_hpos, sentinel


def load_regions(regions_in, sig_df, sig_bt, pvals, cnv_type, pad_sentinel, prefix):
    """
    Read BED of regions to be refined, and annotate with metadata
    """

    regions = {}
    rname_fmt = '{0}_query_region_{1}:{2}-{3}'

    regions_bt = pbt.BedTool(regions_in)

    for region in regions_bt:
        rname = rname_fmt.format(prefix, str(region.chrom), str(region.start), 
                                 str(region.end))
        rstr = '{0}\t{1}\t{2}\n'.format(region.chrom, region.start, region.end)
        rbt = pbt.BedTool(rstr, from_string=True)

        # Get dict keyed on significant window ids overlapping region, 
        # along with significant associated hpos and their p-values per window
        # Also return a list of all significant hpos 
        # and information on the current sentinel window
        sig_windows, sig_hpos, sentinel = get_sig_windows(rbt, sig_df, sig_bt, 
                                                          pvals, pad_sentinel)
        
        if rname not in regions.keys():
            regions[rname] = {'chr' : region.chrom, 
                              'start' : region.start, 
                              'end' : region.end, 
                              'region_id' : rname, 
                              'sig_windows' : sig_windows,
                              'sig_hpos' : sig_hpos, 
                              'sentinel' : sentinel}

    return regions, regions_bt


def print_region_info(regions, rid, logfile):
    """
    Print basic region info prior to refinement
    """

    print('\nBeginning refinement of {0}...'.format(rid), file=logfile)
    rsize = regions[rid]['end'] - regions[rid]['start']
    print('  * Query region size: {:,} kb'.format(int(round(rsize/1000))), 
          file=logfile)
    n_phenos = len(regions[rid]['sig_hpos'])
    pheno_str = ', '.join(regions[rid]['sig_hpos'])
    print('  * Associated with {0} phenotypes ({1})'.format(str(n_phenos), 
                                                              pheno_str), 
          file=logfile)


def get_cc_counts(cohorts, case_hpos, control_hpos):
    """
    Parses per-cohort sample phenotypes to count effective cases and controls
    """

    def _count_samples(samples, hpos):
        k = 0
        for sample in samples:
            for hpo in sample:
                if hpo in hpos:
                    k += 1
                    break
        return k

    cc_counts = {}

    for cohort in cohorts.items():
        n_case = _count_samples(cohort[1]['hpos'], case_hpos)
        n_control = _count_samples(cohort[1]['hpos'], control_hpos)
        cc_counts[cohort[0]] = {'n_case' : n_case,
                                'n_control' : n_control}

    return cc_counts


def filter_cnvs(query_bt, cnv_bt, hpos, max_size=10000000, 
                frac=0.5, exclude_cnvs=None, include_cnvs=None):
    """
    Intersect a query region with CNVs from a subset of phenotypes
    """

    cnv_hits = cnv_bt.intersect(query_bt, wa=True, F=frac)
    
    def _pheno_filter(cnv, keep_phenos):
        if len([p for p in cnv[5].split(';') if p in keep_phenos]) > 0:
            return True
        else:
            return False

    if exclude_cnvs is not None:
        cnv_hits = cnv_hits.filter(lambda x: x.name not in exclude_cnvs)

    if include_cnvs is not None:
        cnv_hits = cnv_hits.filter(lambda x: x.name in include_cnvs)

    cnv_hits = cnv_hits.filter(lambda x: x.length <= max_size)
    
    return cnv_hits.filter(_pheno_filter, keep_phenos=hpos).saveas()


def sort_cnvs(cnv_bt, sent_mid):
    """
    Order case CNVs based on distance to sentinel midpoint
    Input as bt, output as pd.dataframe
    """

    dists = cnv_bt.loc[:, 'start':'end'].apply(lambda x: abs(x - sent_mid), axis=1)
    cnv_bt = pd.concat([cnv_bt, dists.apply(max, axis=1)], axis=1, ignore_index=True)
    cnv_bt.columns = ['chr', 'start', 'end', 'cnv_id', 'cnv', 'phenos', 'abs_min_dist']
    cnv_bt.sort_values(by='abs_min_dist', inplace=True)
    
    return cnv_bt


def sweeting_correction(n_control_ref, n_case_ref, n_control_alt, n_case_alt, cc_sum):
    """
    Apply empirical continuity correction to case & control CNV data
    Per Sweeting et al., Stat. Med., 2004 (section 3.3)
    """

    # Count number of carriers and non-carriers
    n_alt = sum(np.array(n_control_alt + n_case_alt))
    n_alt_bycohort = list(np.array(n_control_alt) + np.array(n_case_alt))
    n_ref = sum(np.array(n_control_ref + n_case_ref))

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
            min_p=1, min_or_lower=0, min_nominal=0, empirical_continuity=True, 
            cc_sum=0.01):
    """
    Conduct Mantel-Haenszel meta-analysis

    Note: this function is deprecated in favor of calling exact same Rscript
          with subprocess() to ensure mathematical consistency
    """

    # Only run if any CNV carriers are observed
    if sum(np.array(n_control_alt + n_case_alt)) > 0:

        # Apply empirical continuity correction
        n_control_ref = list(np.array(n_control) - np.array(n_control_alt))
        n_case_ref = list(np.array(n_case) - np.array(n_case_alt))
        if empirical_continuity :
            n_control_ref, n_case_ref, n_control_alt, n_case_alt = \
                sweeting_correction(n_control_ref, n_case_ref, n_control_alt, 
                                    n_case_alt, cc_sum)

        # Create list of 2x2 tables
        cohort_idxs = list(range(0, len(cohorts)))
        def _make_cont_table(idx):
            tb = np.array([[n_control_ref[idx], n_case_ref[idx]], 
                        [n_control_alt[idx], n_case_alt[idx]]])
            return tb
        meta_tabs = [_make_cont_table(i) for i in cohort_idxs]

        # Compute number of cohorts at nominal significance
        fisher_p = [sstats.fisher_exact(x, 'greater')[1] for x in meta_tabs]
        n_nom = len([p for p in fisher_p if p <= 0.05])

        # Calculate odds ratio, chisq, and p-value
        strat_tab = sm.stats.StratifiedTable(meta_tabs)
        np.seterr(divide = 'ignore')
        oddsratio = strat_tab.oddsratio_pooled
        oddsratio_ci = strat_tab.oddsratio_pooled_confint()
        chisq = strat_tab.test_null_odds().statistic
        pval = strat_tab.test_null_odds().pvalue
        if oddsratio < 1:
            pval = 1 - pval

    else:
        n_nom = 0
        oddsratio = 1
        oddsratio_ci = (np.NaN, np.NaN)
        chisq = 0
        pval = 1

    # Make determination on overall test significance
    if pval <= min_p \
    and oddsratio_ci[0] >= min_or_lower \
    and n_nom >= min_nominal:
        sig = True
    else:
        sig = False

    # # HELPFUL DEBUG:
    # msg = '\t'.join([str(x) for x in [pval, oddsratio, oddsratio_ci[0], chisq, n_nom]])

    return oddsratio, oddsratio_ci, chisq, pval, n_nom, sig


def meta_analysis(model, cohorts, n_control, n_case, n_control_alt, n_case_alt,
                  min_p=1, min_secondary_p=1, min_or_lower=0, min_nominal=0, 
                  empirical_continuity=True, cc_sum=0.01):
    """
    Perform random-effects or Mantel-Haenszel meta-analysis

    Wraps subprocess.call() to the Rscript used for discovery meta-analysis
    This is necessary to ensure mathematical consistency between discovery & refinement
    """

    n_case_ref = list(np.array(n_case) - np.array(n_case_alt))
    n_control_ref = list(np.array(n_control) - np.array(n_control_alt))
    cohort_idxs = list(range(0, len(cohorts)))

    # Compute number of cohorts at nominal significance
    def _make_cont_table(idx):
        tb = np.array([[n_control_ref[idx], n_case_ref[idx]], 
                    [n_control_alt[idx], n_case_alt[idx]]])
        return tb
    meta_tabs = [_make_cont_table(i) for i in cohort_idxs]
    fisher_p = [sstats.fisher_exact(x, 'greater')[1] for x in meta_tabs]
    n_nom = len([p for p in fisher_p if p <= 0.05])

    # Wraps Rscript call in context of temporary directory for easy cleanup
    colnames = '#chr start end case_alt case_ref control_alt control_ref fisher_phred_p'
    in_header = '\t'.join(colnames.split())
    with tempfile.TemporaryDirectory() as tmpdir:
        # Writes input files for meta-analysis script
        m_fn = tmpdir + '/meta.in.tsv'
        m_fin = open(m_fn, 'w')    
        for i in cohort_idxs:
            tmp_fn = '{}/{}.tmp'.format(tmpdir, cohorts[i])
            tmp_fin = open(tmp_fn, 'w')
            tmp_fin.write(in_header + '\n')
            counts = [n_case_alt[i], n_case_ref[i], 
                      n_control_alt[i], n_control_ref[i]]
            dummy = '1 1 2'.split() + [str(x) for x in counts] + ['1']
            tmp_fin.write('\t'.join(dummy) + '\n')
            tmp_fin.close()
            m_fin.write('{}\t{}\n'.format(cohorts[i], tmp_fn))
        m_fin.close()

        # Calls Rscript
        script = '/'.join([path.dirname(path.realpath(__file__)),
                           'window_meta_analysis.R'])
        m_fo = tmpdir + '/meta.out.tsv'
        cstr = '{} --model {} --p-is-phred {} {}'
        subprocess.call(cstr.format(script, model, m_fn, m_fo), shell=True)

        # Read results back into python
        res = pd.read_csv(m_fo, delim_whitespace=True)
        pval = 10 ** -float(res['meta_phred_p'])
        secondary_pval = 10 ** -float(res['meta_phred_p_secondary'])
        oddsratio = float(res['meta_lnOR'])
        oddsratio_ci = [float(res['meta_lnOR_lower']), 
                        float(res['meta_lnOR_upper'])]
        chisq = float(res.iloc[:, -2])
        if pval <= min_p \
        and secondary_pval <= min_secondary_p \
        and oddsratio_ci[0] >= min_or_lower \
        and n_nom >= min_nominal:
            sig = True
        else:
            sig = False

    # # HELPFUL DEBUG:
    # msg = '\t'.join([str(x) for x in [pval, oddsratio, oddsratio_ci[0], chisq, n_nom]])
    # print(msg + '\n')

    return oddsratio, oddsratio_ci, chisq, pval, secondary_pval, n_nom, sig


def _count_cnvs_by_cohort(cnvs, cohorts):
    """
    Count CNVs per cohort from an input pd.dataframe 
    """

    cnv_ids = [x.split('_')[0] for x in cnvs['cnv_id'].values]
    counts = [len([x for x in cnv_ids if x == cohort]) for cohort in cohorts]
    
    return counts


def get_min_case_cnvs(case_cnvs, control_cnvs, cc_counts, min_p, min_secondary_p,
                      min_or_lower, min_nominal, model):
    """
    Get minimum number of case CNVs required to achieve a significant result
    """

    cohorts = list(cc_counts.keys())

    # Get all control count data
    n_control = [v['n_control'] for (k, v) in cc_counts.items()]
    n_control_alt = _count_cnvs_by_cohort(control_cnvs, cohorts)

    # Get baseline case count data
    n_case = [v['n_case'] for (k, v) in cc_counts.items()]
    n_case_alt_max = _count_cnvs_by_cohort(case_cnvs, cohorts)
    sum_n_case_alt_max = sum(np.array(n_case_alt_max))

    # Increment case CNVs one at a time until significance
    k = 0
    sig = False
    if sum_n_case_alt_max > 0:
        while sig is False:
            k += 1
            n_case_alt = _count_cnvs_by_cohort(case_cnvs.iloc[0:k, :], cohorts)
            oddsratio, or_ci, chisq, new_p, new_secondary_p, n_nom, sig \
                = meta_analysis(model, cohorts, n_control, n_case, n_control_alt, 
                                n_case_alt, min_p, min_secondary_p, min_or_lower, 
                                min_nominal)
            if k == sum_n_case_alt_max:
                break
    
    return k, sig


def calc_credible_interval(cnvs_df, method='density', resolution=10000, cred=0.8):
    """
    Given a set of CNVs (as pd.dataframe), calculate the credible interval
    that defines their minimum region of overlap, either based on breakpoint
    positions or CNV density
    """

    if method not in 'density breakpoint'.split():
        exit('Credible interval method \'{}\' is not recognized.'.format(method))

    tail_width = (1 - cred) / 2
    q_tails = [tail_width, 1 - tail_width]

    if method == 'breakpoint':
        cred_start_exact = np.quantile(cnvs_df['start'], q_tails[0])
        cred_end_exact = np.quantile(cnvs_df['end'], q_tails[1])
        cred_start = int(resolution * np.floor(cred_start_exact / resolution))
        cred_end = int(resolution * np.ceil(cred_end_exact / resolution))

    if method == 'density':
        cnvs_bt = pbt.BedTool.from_dataframe(cnvs_df.iloc[:, :3])
        chrom = cnvs_bt[0].chrom
        bin_min = int(resolution * np.floor(np.nanmin(cnvs_df['start']) / resolution))
        bin_max = int(resolution * np.ceil(np.nanmax(cnvs_df['end']) / resolution))
        binrange_bt = pbt.BedTool('\t'.join([chrom, str(bin_min), str(bin_max)]), 
                                  from_string=True)
        bins_bt = pbt.BedTool().window_maker(b=binrange_bt, w=resolution)
        bin_counts = [int(x[3]) for x in bins_bt.coverage(cnvs_bt, c=True)]
        bin_cdf = np.cumsum(bin_counts)
        cdf_max = np.nanmax(bin_cdf)
        cred_cdf_min = q_tails[0] * np.nanmax(bin_cdf)
        cred_cdf_max = q_tails[1] * np.nanmax(bin_cdf)
        bins_df = bins_bt.to_dataframe()
        bins_df['cdf'] = bin_cdf
        cred_start = np.nanmax(np.array(bins_df[bins_df['cdf'] <= cred_cdf_min]['start']))
        cred_end = np.nanmin(bins_df[bins_df['cdf'] >= cred_cdf_max]['end'])

    return cred_start, cred_end


def refine_sentinel(regions, rid, sig_df, sig_bt, pvals, cohorts, cnv_type, 
                    control_hpo, p_cutoffs, p_cutoff_ladder, min_secondary_p,
                    min_or_lower, min_nominal, model, exclude_cnvs, logfile, 
                    max_size=10000000, min_case_cnvs=3, resolution=10000, cred=0.8, ci_method='density'):
    """
    Define minimum critical region for a single sentinel window
    """

    rbt = pbt.BedTool('\t'.join([regions[rid]['chr'], 
                                 str(regions[rid]['start']), 
                                 str(regions[rid]['end'])]),
                      from_string=True)

    # Get sentinel information
    sentinel = regions[rid]['sentinel']
    sent_wid = sentinel['window_id']
    sent_hpos = sentinel['all_hpos']
    region_hpos = regions[rid]['sig_hpos']
    sent_bt = pbt.BedTool.from_dataframe(sig_df[sig_df['id'].isin([sent_wid])].iloc[:, 0:3])
    sent_mid = (sent_bt[0].start + sent_bt[0].end) / 2

    # Print sentinel information for logging
    print('  * Sentinel window: {0}'.format(sent_wid), file=logfile)
    print('  * Strongest association: P = {:.2E} in {}'.format(sentinel['pvalue'], 
                                                               sentinel['best_hpo']), 
          file=logfile)
    print('  * Sentinel window associated with {0} phenotypes ({1})'.\
          format(str(len(sent_hpos)), ', '.join(sent_hpos)), file=logfile)

    # Get case/control sample counts for sentinel phenotypes
    # case_hpos = regions[rid]['sig_hpos']
    cc_counts = get_cc_counts(cohorts, sent_hpos, [control_hpo])
    total_n_case = sum([x['n_case'] for x in cc_counts.values()])
    total_n_control = sum([x['n_control'] for x in cc_counts.values()])

    # Set p-value threshold according to number of cases
    if p_cutoff_ladder is not None:
        min_p = np.nanmax(p_cutoff_ladder['min_p'][p_cutoff_ladder['n_cases'] <= total_n_case])
    else:
        min_p = np.nanmax([p_cutoffs[hpo] for hpo in sent_hpos])

    # Intersect region and sentinel window with CNVs from sentinel phenotype
    # Note: do *not* filter control CNVs on excluded CNVs or size
    cnvs = {}
    for cohort in cohorts.items():
        region_case_cnvs = filter_cnvs(rbt, cohort[1]['cnvs'], region_hpos, 
                                       frac=0.0001, exclude_cnvs=exclude_cnvs)
        region_control_cnvs = filter_cnvs(rbt, cohort[1]['cnvs'], [control_hpo],
                                          frac=0.0001)
        all_case_cnvs = filter_cnvs(sent_bt, cohort[1]['cnvs'], sent_hpos, 
                                    frac=0.0001, exclude_cnvs=exclude_cnvs)
        all_control_cnvs = filter_cnvs(sent_bt, cohort[1]['cnvs'], [control_hpo],
                                       frac=0.0001)
        case_cnvs = filter_cnvs(sent_bt, cohort[1]['cnvs'], sent_hpos, 
                                max_size=max_size, frac=0.5, 
                                exclude_cnvs=exclude_cnvs)
        control_cnvs = filter_cnvs(sent_bt, cohort[1]['cnvs'], [control_hpo],
                                   frac=0.5)
        cnvs[cohort[0]] = {'region_case_cnvs' : region_case_cnvs,
                           'region_control_cnvs' : region_control_cnvs,
                           'all_case_cnvs' : all_case_cnvs,
                           'all_control_cnvs' : all_control_cnvs,
                           'case_cnvs' : case_cnvs,
                           'control_cnvs' : control_cnvs}
    
    # Pool and sort all CNVs across metacohorts for cases and controls
    case_cnvs_l = [x['case_cnvs'] for x in cnvs.values()]
    n_case_cnvs = np.sum([len(x) for x in case_cnvs_l])
    if n_case_cnvs > 0:
        case_cnvs = pd.concat([bt.to_dataframe() for bt in case_cnvs_l if len(bt) > 0], 
                              axis=0, ignore_index=True)
        case_cnvs = sort_cnvs(case_cnvs, sent_mid)
    else:
        case_cnvs = pd.DataFrame(columns='chrom start end cnv_id cnv phenos abs_min_dist'.split())
    control_cnvs_l = [x['control_cnvs'] for x in cnvs.values()]
    control_cnvs = pd.concat([bt.to_dataframe() for bt in control_cnvs_l if len(bt) > 0], 
                                 axis=0, ignore_index=True)
    control_cnvs = sort_cnvs(control_cnvs, sent_mid)

    # Print sentinel sample & CNV data
    print('  * Effective case count: {:,}'.format(total_n_case), 
          file=logfile)
    print('  * Found {:,} qualifying case {}s'.format(len(case_cnvs), cnv_type), 
          file=logfile)
    print('  * Effective control count: {:,}'.format(total_n_control), 
          file=logfile)
    print('  * Found {:,} qualifying control {}s'.format(len(control_cnvs), cnv_type), 
          file=logfile)
    print('  * Effective genome-wide threshold: {:.2E}'.format(min_p), 
          file=logfile)

    # Get minimum set of case CNVs required for genome-wide significance
    k_case_cnvs, sig = get_min_case_cnvs(case_cnvs, control_cnvs, cc_counts, 
                                         min_p, min_secondary_p, min_or_lower, 
                                         min_nominal, model)

    if sig:
        min_case_cnvs = case_cnvs.iloc[0:k_case_cnvs, :]
        min_case_cnv_ids = list(min_case_cnvs['cnv_id'].values.flatten())
        n_eff_case_cnvs = len(case_cnvs)

    # If qualifying, window-overlapping CNVs don't reach significance, try again
    # with all window-overlapping CNVs (not just those that passed overlap filter)
    else:
        # Pool and sort all CNVs across metacohorts for cases and controls
        all_case_cnvs_l = [x['all_case_cnvs'] for x in cnvs.values()]
        all_case_cnvs = pd.concat([bt.to_dataframe() for bt in all_case_cnvs_l \
                                       if len(bt) > 0], axis=0, ignore_index=True)
        all_case_cnvs = sort_cnvs(all_case_cnvs, sent_mid)
        all_control_cnvs_l = [x['all_control_cnvs'] for x in cnvs.values()]
        all_control_cnvs = pd.concat([bt.to_dataframe() for bt in \
                                          all_control_cnvs_l if len(bt) > 0], axis=0, 
                                         ignore_index=True)
        all_control_cnvs = sort_cnvs(all_control_cnvs, sent_mid)

        # Rerun incremental analysis 
        msg = '  * The {:,} qualifying case {}s were insufficient to establish ' + \
              'significance. Re-trying with all {:,} {}s overlapping window.'
        print(msg.format(len(case_cnvs), cnv_type, len(all_case_cnvs), cnv_type), 
              file=logfile)

        k_case_cnvs, sig = get_min_case_cnvs(all_case_cnvs, all_control_cnvs, 
                                             cc_counts, min_p, min_secondary_p, 
                                             min_or_lower, min_nominal, model)
        if sig:
            min_case_cnvs = all_case_cnvs.iloc[0:k_case_cnvs, :]
            min_case_cnv_ids = list(min_case_cnvs['cnv_id'].values.flatten())
            n_eff_case_cnvs = len(all_case_cnvs)

        # If all window-overlapping CNVs can't reach significance, stop here
        # and use the full set of all CNVs overlapping the window
        else:
            msg = '  * The larger set of {:,} case {}s were insufficient to establish ' + \
                  'significance. Using all CNVs to define minimal credible region.'
            print(msg.format(len(all_case_cnvs), cnv_type), file=logfile)
            min_case_cnvs = all_case_cnvs
            min_case_cnv_ids = list(min_case_cnvs['cnv_id'].values.flatten())
            n_eff_case_cnvs = len(all_case_cnvs)

    # Calculate credible set of coordinates to define the refined region
    refined_chrom = regions[rid]['chr']
    refined_start, refined_end = calc_credible_interval(min_case_cnvs, ci_method, 
                                                        resolution, cred)
    refined_size = int(refined_end - refined_start)

    # Print result to log
    msg = '  * Successfully refined to {:,} kb minimal credible region ' + \
          '({}:{:,}-{:,}) using {:,} of {:,} case CNVs'
    print(msg.format(int(np.floor(refined_size / 1000)), refined_chrom, 
                     refined_start, refined_end, k_case_cnvs, n_eff_case_cnvs), 
          file=logfile)

    # Get any other case CNVs matching description of critical region but not 
    # included in genome-wide significance minimal set or significance test
    refined_bt = pbt.BedTool('{}\t{}\t{}'.format(refined_chrom, 
                                                 str(refined_start),
                                                 str(refined_end)),
                             from_string=True)
    supp_case_ids = []
    for cnvbt in [x['region_case_cnvs'] for x in cnvs.values()]:
        ids = [x.name for x in cnvbt.intersect(refined_bt, f=0.50, u=True)]
        supp_case_ids = supp_case_ids + ids

    all_case_cnv_ids = list(set(min_case_cnv_ids + supp_case_ids))

    # Print number of CNVs excluded to log
    msg = '  * Excluding {:,} case {}s for subsequent conditional analyses due ' + \
          'to at least 50% coverage by minimal credible region.'
    print(msg.format(len(all_case_cnv_ids), cnv_type), file=logfile)

    # Return refined region
    new_region = {'chr' : refined_chrom, 
                  'start' : refined_start, 
                  'end' : refined_end, 
                  'size' : refined_size, 
                  'sentinel' : sentinel,
                  'case_CNV_ids' : all_case_cnv_ids}

    return new_region


def get_cnvs_by_id(cohorts, ids):
    """
    Extract a BedTool of CNVs across cohorts based on a list of CNV IDs
    """

    filt_cnvs = []

    for cohort in cohorts.values():
        cnvs = cohort['cnvs']
        filt_cnvs.append(cnvs.filter(lambda x: x.name in ids).saveas())

    cnvs_df = pd.concat([bt.to_dataframe() for bt in filt_cnvs if len(bt) > 0], 
                            axis=0, ignore_index=True)
    cnvs_df.columns = ['chr', 'start', 'end', 'cnv_id', 'cnv', 'phenos']

    cnvs_bt = pbt.BedTool.from_dataframe(cnvs_df)

    return cnvs_bt


def update_pvalues(window_ids, pvals, sig_df, sig_bt, cohorts, exclude_cnvs, 
                   exclude_cnv_bt, p_cutoffs, min_secondary_p, min_nominal, 
                   min_or_lower, model, max_size, min_case_cnvs, control_hpo):
    """
    Recalculate p-values and significance matrix for significant hpos from provided windows
    """

    # Iterate over all provided windows
    for wid in window_ids:

        wbt = sig_bt.filter(lambda x: x.name == wid).saveas()
        wpval_orig = pvals[pvals['id'] == wid].iloc[:, 4:].values.flatten()

        # Only run if any excluded CNVs overlap the window
        # If no excluded CNVs overlap the window, there's no point in retesting
        # as P-value cannot have changed
        if len(exclude_cnv_bt.intersect(wbt, u=True)) == 0:
            continue

        # Iterate over all HPOs and recalculate p-values
        new_pvals = []
        for hpo in pvals.columns[4:]:

            # Only recompute p-value if hpo was originally significant
            if sig_df.loc[sig_df['id'] == wid, hpo].values[0]:

                # Get case/control sample counts for phenotype
                cc_counts = get_cc_counts(cohorts, [hpo], [control_hpo])
                n_control = [v['n_control'] for (k, v) in cc_counts.items()]
                n_case = [v['n_case'] for (k, v) in cc_counts.items()]

                # Intersect sentinel window with CNVs from phenotype
                # Note: do *not* exclude control CNVs based on used cnvs
                cnvs = {}
                for cohort in cohorts.items():
                    case_cnvs = filter_cnvs(wbt, cohort[1]['cnvs'], [hpo], 
                                            max_size=max_size, 
                                            exclude_cnvs=exclude_cnvs)
                    control_cnvs = filter_cnvs(wbt, cohort[1]['cnvs'], 
                                               [control_hpo])
                    cnvs[cohort[0]] = {'case_cnvs' : case_cnvs,
                                       'control_cnvs' : control_cnvs}
                
                # Pool all CNVs across metacohorts for cases and controls
                case_cnvs_l = [x['case_cnvs'] for x in cnvs.values()]
                if len([x for x in case_cnvs_l if len(x) > 0]) > 0:
                    case_cnvs = pd.concat([bt.to_dataframe() for bt in \
                                               case_cnvs_l if len(bt) > 0], 
                                              axis=0, ignore_index=True)
                    case_cnvs.columns = ['chr', 'start', 'end', 'cnv_id', \
                                         'cnv', 'phenos']
                    n_case_alt = _count_cnvs_by_cohort(case_cnvs, cohorts)
                else:
                    n_case_alt = [0] * len(cohorts)
                control_cnvs_l = [x['control_cnvs'] for x in cnvs.values()]
                if len([x for x in control_cnvs_l if len(x) > 0]) > 0:
                    control_cnvs = pd.concat([bt.to_dataframe() for bt in \
                                                  control_cnvs_l if len(bt) > 0], 
                                                 axis=0, ignore_index=True)
                    control_cnvs.columns = ['chr', 'start', 'end', 'cnv_id', \
                                            'cnv', 'phenos']
                    n_control_alt = _count_cnvs_by_cohort(control_cnvs, cohorts)
                else:
                    n_control_alt = [0] * len(cohorts)

                # Get HPO-specific p-value cutoff
                min_p = p_cutoffs[hpo]

                # Automatically set p-value to 1 if not enough case CNVs are observed
                # (This will cause this window/hpo pair to be skipped in secondary
                #  refinement analyses)
                if len(case_cnvs) < min_case_cnvs:
                    new_pvals.append(1)
                    sig_df.loc[sig_df['id'] == wid, hpo] = False
                
                else:
                    # Run meta-analysis & update p-value
                    oddsratio, or_ci, chisq, new_p, new_secondary_p, n_nom, sig \
                        = meta_analysis(model, list(cohorts.keys()), n_control, 
                                        n_case, n_control_alt, n_case_alt, min_p, 
                                        min_secondary_p, min_or_lower, min_nominal)
                    new_pvals.append(new_p)

                    # Update significance label if HPO is no longer significant
                    if not sig:
                        sig_df.loc[sig_df['id'] == wid, hpo] = False

            # If hpo was not originally significant, keep old pvalue
            else:
                new_pvals.append(pvals.loc[pvals['id'] == wid, hpo].values[0])
        
        # Assign new p-values
        pvals.loc[pvals['id'] == wid, pvals.columns[4:]] = new_pvals


def update_regions(rids, regions, sig_df, orig_sig_df, sig_bt, pvals, 
                   pad_sentinel):
    """
    Update sentinel & significance info of region(s)
    """

    for rid in rids:
        rstr = '{0}\t{1}\t{2}\n'.format(regions[rid]['chr'], 
                                        str(regions[rid]['start']), 
                                        str(regions[rid]['end']))
        rbt = pbt.BedTool(rstr, from_string=True)


        sig_windows, sig_hpos, sentinel = get_sig_windows(rbt, sig_df, sig_bt, 
                                                          pvals, pad_sentinel)

        regions[rid]['sig_windows'] = sig_windows
        regions[rid]['sig_hpos'] = sig_hpos
        regions[rid]['sentinel'] = sentinel


def assign_cnvs_to_mcrs(mcrs, cohorts):
    """
    Assign CNVs to MCRs
    """

    # Make BedTool of MCRs
    mcrs_str = ''
    for mcr in mcrs:
        mcr_s = '\t'.join([str(x) for x in [mcr['chr'], mcr['start'], 
                                            mcr['end'], mcr['id']]])
        mcrs_str = mcrs_str + mcr_s + '\n'
    mcrs_bt = pbt.BedTool(mcrs_str, from_string=True)

    # Intersect CNVs with MCRs
    cnvs_l = [x['cnvs'].intersect(mcrs_bt, wo=True) for x in cohorts.values()]
    def _get_pair_info(x):
        aob = int(x[10]) / x.length
        boa = int(x[10]) / (int(x[8]) - int(x[7]))
        res = {'cnv_id' : x.name, 'mcr_id' : x[9], 
               'frac_cnv_overlapped' : aob, 
               'frac_mcr_overlapped' : boa, 
               'min_overlap' : np.nanmin([aob, boa])}
        return res
    pair_info = [_get_pair_info(x) for cnvs in cnvs_l for x in cnvs]
    cnv_ids = list(set([x['cnv_id'] for x in pair_info]))

    # Assign each CNV to the MCR it overlaps best
    mcr_cnvs = {}
    for cnv_id in cnv_ids:
        hits = [(x['mcr_id'], x['min_overlap']) for x in pair_info 
                if x['cnv_id'] == cnv_id]
        best_ovr = np.nanmax([ovr for mcr, ovr in hits])
        best_hit = [mcr for mcr, ovr in hits if ovr == best_ovr][0]
        if best_hit not in mcr_cnvs.keys():
            mcr_cnvs[best_hit] = [cnv_id]
        else:
            mcr_cnvs[best_hit].append(cnv_id)

    return mcr_cnvs


def get_final_stats(mcrs, cohorts, orig_sig_df, orig_pvals, control_hpo, model):
    """
    Calculate final association statistics for all minimal credible regions
    """

    all_cnvs = pbt.BedTool.cat(*[x['cnvs'] for x in cohorts.values()], 
                               postmerge=False)
    orig_sig_bt = pbt.BedTool.from_dataframe(orig_sig_df.iloc[:, :4])
    mcr_stats = {}

    # Assign CNVs to MCRs
    mcr_cnvs = assign_cnvs_to_mcrs(mcrs, cohorts)

    # Iterate over MCRs
    for mcr in mcrs:

        mcr_id = mcr['id']
        mcr_stats[mcr_id] = {'chr' : mcr['chr'], 
                             'start' : mcr['start'], 
                             'end' : mcr['end'], 
                             'size' : mcr['size'],
                             'id' : mcr_id}

        mcr_bt_s = '\t'.join([str(x) for x in [mcr['chr'], mcr['start'], mcr['end']]])
        mcr_bt = pbt.BedTool(mcr_bt_s, from_string=True)

        # Get original phenotypes with associations in MCR
        orig_sig_windows, orig_sig_hpos, orig_sentinel \
            = get_sig_windows(mcr_bt, orig_sig_df, orig_sig_bt, orig_pvals, 
                              pad_sentinel=100)
        mcr_stats[mcr_id]['hpos'] = orig_sig_hpos
        orig_sig_wids = [x.name for x in orig_sig_bt.intersect(mcr_bt, u=True)]

        # Iterate over each originally significant phenotype
        for hpo in orig_sig_hpos:

            # Get case/control sample counts for phenotype
            cc_counts = get_cc_counts(cohorts, [hpo], [control_hpo])
            n_control = [v['n_control'] for (k, v) in cc_counts.items()]
            n_case = [v['n_case'] for (k, v) in cc_counts.items()]

            # Extract original sentinel ID & p-value
            hpo_sig_wids = orig_sig_df[orig_sig_df[hpo]]['id']
            hpo_mcr_sig_wids = [x for x in hpo_sig_wids if x in orig_sig_wids]
            best_p = np.nanmin(orig_pvals[orig_pvals['id'].isin(hpo_mcr_sig_wids)][hpo])
            best_wid = orig_pvals[orig_pvals[hpo] == best_p]['id'].values[0]
            best_w_bt = orig_sig_bt.filter(lambda x: x.name == best_wid).saveas()

            # Compute odds ratio for original sentinel stats
            cnv_df_names = ['chr', 'start', 'end', 'cnv_id', 'cnv', 'pheno']
            control_cnvs = filter_cnvs(best_w_bt, all_cnvs, hpos=[control_hpo])
            control_cnvs_df = control_cnvs.to_dataframe(names=cnv_df_names)
            n_control_alt = _count_cnvs_by_cohort(control_cnvs_df, cohorts.keys())
            case_cnvs = filter_cnvs(best_w_bt, all_cnvs, hpos=[hpo])
            case_cnvs_df = case_cnvs.to_dataframe(names=cnv_df_names)
            n_case_alt = _count_cnvs_by_cohort(case_cnvs_df, cohorts.keys())
            oddsratio, or_ci, chisq, new_p, new_secondary_p, n_nom, sig \
                = meta_analysis(model, list(cohorts.keys()), n_control, n_case, n_control_alt, n_case_alt)
            hpo_sent_stats = {'chr' : mcr['chr'],
                              'start' : best_w_bt[0].start,
                              'end' : best_w_bt[0].end,
                              'wid' : best_wid,
                              'oddsratio' : oddsratio,
                              'or_ci' : or_ci, 
                              'pvalue' : best_p}

            # Compute conditional stats for MCR
            mcr_control_cnvs = filter_cnvs(mcr_bt, all_cnvs, hpos=[control_hpo], 
                                           frac=10e-10, include_cnvs=mcr_cnvs[mcr_id])
            mcr_control_cnvs_df = mcr_control_cnvs.to_dataframe(names=cnv_df_names)
            mcr_n_control_alt = _count_cnvs_by_cohort(mcr_control_cnvs_df, cohorts.keys())
            mcr_case_cnvs = filter_cnvs(mcr_bt, all_cnvs, hpos=[hpo], frac=10e-10, 
                                        include_cnvs=mcr_cnvs[mcr_id])
            mcr_case_cnvs_df = mcr_case_cnvs.to_dataframe(names=cnv_df_names)
            mcr_n_case_alt = _count_cnvs_by_cohort(mcr_case_cnvs_df, cohorts.keys())
            mcr_oddsratio, mcr_or_ci, mcr_chisq, mcr_pval, mcr_secondary_pval, mcr_n_nom, mcr_sig \
                = meta_analysis(model, list(cohorts.keys()), n_control, n_case, mcr_n_control_alt, mcr_n_case_alt)
            hpo_mcr_stats = {'oddsratio' : mcr_oddsratio,
                             'or_ci' : mcr_or_ci, 
                             'pvalue' : mcr_pval}

            # Add stats for hpo to MCR stats dict
            mcr_stats[mcr_id][hpo] = {'sentinel' : hpo_sent_stats,
                                      'mcr' : hpo_mcr_stats}

        # Compute pooled stats across all associated phenotypes for region-level output
        cc_counts = get_cc_counts(cohorts, orig_sig_hpos, [control_hpo])
        n_control = [v['n_control'] for (k, v) in cc_counts.items()]
        n_case = [v['n_case'] for (k, v) in cc_counts.items()]
        mcr_control_cnvs = filter_cnvs(mcr_bt, all_cnvs, hpos=[control_hpo], 
                                       frac=10e-10, include_cnvs=mcr_cnvs[mcr_id])
        mcr_control_cnvs_df = mcr_control_cnvs.to_dataframe(names=cnv_df_names)
        mcr_n_control_alt = _count_cnvs_by_cohort(mcr_control_cnvs_df, cohorts.keys())
        mcr_case_cnvs = filter_cnvs(mcr_bt, all_cnvs, hpos=orig_sig_hpos, frac=10e-10, 
                                    include_cnvs=mcr_cnvs[mcr_id])
        mcr_case_cnvs_df = mcr_case_cnvs.to_dataframe(names=cnv_df_names)
        mcr_n_case_alt = _count_cnvs_by_cohort(mcr_case_cnvs_df, cohorts.keys())
        mcr_oddsratio, mcr_or_ci, mcr_chisq, mcr_pval, mcr_secondary_pval, mcr_n_nom, mcr_sig \
            = meta_analysis(model, list(cohorts.keys()), n_control, n_case, mcr_n_control_alt, mcr_n_case_alt)

        # Add pooled stats to MCR stats dict
        mcr_stats[mcr_id]['pooled_oddsratio'] = mcr_oddsratio
        mcr_stats[mcr_id]['pooled_or_ci'] = mcr_or_ci
        mcr_stats[mcr_id]['pooled_pvalue'] = mcr_pval

    return mcr_stats


def write_associations(mcr_stats, fout, cnv_type):
    """
    Write association-level stats to outfile
    """

    # Write header
    header_cols = '#chr start end region_id size cnv hpo or or_ci_lower or_ci_upper ' + \
                  'pvalue sentinel_start sentinel_end sentinel_id sentinel_or ' + \
                  'sentinel_or_lower sentinel_or_upper sentinel_pvalue'
    fout.write('\t'.join(header_cols.split()) + '\n')

    # Iterate over HPOs for each MCR
    for mcr in mcr_stats.values():
        base_mcr_info = '\t'.join([mcr['chr'], str(mcr['start']), str(mcr['end']), 
                                   mcr['id'], str(mcr['size']), cnv_type])

        for hpo in mcr['hpos']:
            hpo_mcr_stats = mcr[hpo]['mcr']
            hpo_mcr_info = '\t'.join([hpo, 
                                      str(round(hpo_mcr_stats['oddsratio'], 2)),
                                      str(round(hpo_mcr_stats['or_ci'][0], 2)),
                                      str(round(hpo_mcr_stats['or_ci'][1], 2)),
                                      '{:.2e}'.format(hpo_mcr_stats['pvalue'])])

            hpo_sent_stats = mcr[hpo]['sentinel']
            hpo_sent_info = '\t'.join([hpo_sent_stats['wid'],
                                       str(hpo_sent_stats['start']),
                                       str(hpo_sent_stats['end']),
                                       str(round(hpo_sent_stats['oddsratio'], 2)),
                                       str(round(hpo_sent_stats['or_ci'][0], 2)),
                                       str(round(hpo_sent_stats['or_ci'][1], 2)),
                                       '{:.2e}'.format(hpo_sent_stats['pvalue'])])

            fout.write('\t'.join([base_mcr_info, hpo_mcr_info, hpo_sent_info]) + '\n')


def write_regions(mcr_stats, fout, cnv_type):
    """
    Write region-level stats to outfile
    """

    # Write header
    header_cols = '#chr start end region_id size cnv pooled_or pooled_or_ci_lower ' + \
                  'pooled_or_ci_upper pooled_pvalue hpos'
    fout.write('\t'.join(header_cols.split()) + '\n')

    # Iterate over MCRs
    for mcr in mcr_stats.values():
        base_mcr_info = '\t'.join([mcr['chr'], str(mcr['start']), str(mcr['end']), 
                                  mcr['id'], str(mcr['size']), cnv_type])
        assoc_stats = '\t'.join([str(round(mcr['pooled_oddsratio'], 2)),
                                 str(round(mcr['pooled_or_ci'][0], 2)),
                                 str(round(mcr['pooled_or_ci'][1], 2)),
                                 '{:.2e}'.format(mcr['pooled_pvalue'])])
        hpos = ';'.join(mcr['hpos'])

        fout.write('\t'.join([base_mcr_info, assoc_stats, hpos]) + '\n')


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('regions', help='BED file of regions to refine.')
    parser.add_argument('cohort_info', help='Three-column TSV file specifying ' +
                        'cohort name, path to CNV BED file, and path to sample-' +
                        'level HPO terms.')
    parser.add_argument('pvalues', help='BED file of original sliding windows ' +
                         'and p-values for each phenotype.')
    parser.add_argument('sig_df', help='BED file of original sliding windows ' +
                         'with TRUE/FALSE label of significance for each phenotype.')
    parser.add_argument('associations_out', help='Path to output BED with HPO-level ' +
                        'association stats (one line per phenotype-region pair.')
    parser.add_argument('regions_out', help='Path to output BED with region-level ' +
                        'association stats (one line per region.')
    parser.add_argument('--cnv-type', help='Type of CNVs to evaluate. [default: ' + 
                        'use all CNVs]', choices=['CNV', 'DEL', 'DUP'], 
                        default='CNV')
    parser.add_argument('--model', help='Meta-analysis model to use. [default: density]', 
                        default='mh', choices=['mh', 're'])
    parser.add_argument('--hpo-p-cutoffs', help='.tsv of p-value cutoffs per phenotype. ' + 
                        '[default: 10e-8 for all phenotypes]')
    parser.add_argument('--p-cutoff-ladder', help='.tsv of p-value cutoffs for a ' + 
                        'range of case counts. Effective case count for sentinel ' +
                        'window will round down to nearest case sample size, if ' +
                        'provided. [default: use P-values from --hpo-p-cutoffs]')
    parser.add_argument('--p-is-phred', help='Supplied P-values are Phred-scaled ' +
                        '(-log10[P]). [default: False]', default=False, 
                        action='store_true')
    parser.add_argument('--secondary-p-cutoff', help='Maximum secondary P-value to ' + 
                        'consider as significant. [default: 1]', default=1, type=float)
    parser.add_argument('--min-or-lower', help='Minimum lower bound of 95% confidence ' + 
                        'interval for odds ratio. Supply as untransformed OR. [default: 1]', 
                        default=1, type=float)
    parser.add_argument('--retest-min-or-lower', help='Minimum lower bound of 95% confidence ' + 
                        'interval for odds ratio used when evaluating secondary sentinels. ' + 
                        'Supply as untransformed OR. [default: 1]', default=1, type=float)
    parser.add_argument('--max-cnv-size', help='Maximum size of CNVs to include when ' +
                        'attempting to refine minimal credible regions [default: 3Mb]', 
                        default=10000000, type=int)
    parser.add_argument('--min-nominal', help='Minimum number of individual cohorts ' + 
                        'required to be nominally significant to report an independent ' +
                        'minimal credible region. [default: 1]', 
                        default=1, type=int)
    parser.add_argument('--min-case-cnvs', help='Minimum count of CNVs required ' +
                        'per independent minimal credible region [default: 1]', 
                        default=1, type=int)
    parser.add_argument('--control-hpo', default='HEALTHY_CONTROL', help='HPO code ' +
                        'to use for control CNV counts. [default: HEALTHY_CONTROL]')
    parser.add_argument('--resolution', help='Round refined region coordinates ' + 
                        'up to the nearest X bp. [default: 10000]', default=10000, 
                        type=int)
    parser.add_argument('--credible-interval', help='Specify credible interval ' + 
                        'to use when reporting refined regions [default: 0.9]', 
                        default=0.9, type=float)
    parser.add_argument('--ci-method', help='Method to use when defining minimal ' +
                        'credible intervals. [default: density]', 
                        default='density', choices=['density', 'breakpoint'])
    parser.add_argument('--pad-sentinel', help='Distance from sentinel window to ' + 
                        'look for other significant associations. [default: 100000]', 
                        default=100000, type=int)
    parser.add_argument('--prefix', help='String to append to all refined regions. ' +
                        '[default: CNV type (specified above)]', default=None)
    parser.add_argument('--log', help='Output file for all progress logging. ' +
                        '[default: stdout]', default='stdout')

    args = parser.parse_args()

    # Postprocess args as necessary
    if args.prefix is None:
        args.prefix = args.cnv_type
    min_or_lower = np.log(args.min_or_lower)
    retest_min_or_lower = np.log(args.retest_min_or_lower)

    # Open connections to outfiles
    assocs_fout = open(args.associations_out, 'w')
    regions_fout = open(args.regions_out, 'w')
    if args.log in 'stdout - /dev/stdout'.split():
        logfile = stdout
    else:
        logfile = open(args.log, 'w')

    # Load p-values, significance matrix & BedTool of significant windows
    sig_df, sig_bt = load_sig_df(args.sig_df, args.cnv_type)
    pvals = load_pvalues(args.pvalues, sig_bt, args.cnv_type, args.p_is_phred)

    # Save a copy of original p-values and sig matrix for final reporting
    orig_pvals = pvals.copy(deep=True)
    orig_sig_df = sig_df.copy(deep=True)

    # Read original list of regions to refine
    regions, regions_bt = load_regions(args.regions, sig_df, sig_bt, pvals, 
                                       args.cnv_type, args.pad_sentinel, 
                                       args.prefix)
    orig_regions = regions.copy()
    msg = '* Identified {:,} query regions to be refined. Beginning refinement:\n'
    print(msg.format(len(regions)), file=logfile)

    # Load cohort info, including CNVs & HPOs
    cohorts = load_cohort_info(args.cohort_info, regions_bt, args.cnv_type)
    if args.hpo_p_cutoffs is None:
        hpo_p_cutoffs = {hpo : 1e-8 for hpo in pvals.columns[4:]}
    else:
        hpo_p_cutoffs = load_p_cutoffs(args.hpo_p_cutoffs)

    # Load p-value cutoff ladder, if optioned
    if args.p_cutoff_ladder is not None:
        p_cutoff_ladder = pd.read_csv(args.p_cutoff_ladder, sep='\t', comment='#', 
                                      names='n_cases min_p'.split())
    else:
        p_cutoff_ladder = None

    # Refine regions one at a time
    k = 0
    all_mcrs = []
    used_cnvs = []
    used_cnv_bt = None
    for rid in regions.keys():
        print_region_info(regions, rid, logfile)

        i = 0
        # Sequentially process all sentinels
        while len(regions[rid]['sig_hpos']) > 0:
            i += 1
            k += 1
            print('--Refining sentinel window #{0}'.format(str(i)), file=logfile)
            new_region = refine_sentinel(regions, rid, sig_df, sig_bt, pvals, 
                                         cohorts, args.cnv_type, args.control_hpo, 
                                         hpo_p_cutoffs, p_cutoff_ladder, 
                                         args.secondary_p_cutoff, min_or_lower, 
                                         args.min_nominal, args.model, used_cnvs, 
                                         logfile, args.max_cnv_size, args.min_case_cnvs, 
                                         args.resolution, args.credible_interval, 
                                         args.ci_method)
            all_mcrs.append(new_region)

            # Add used CNVs to blacklist
            used_cnvs = list(set(used_cnvs + new_region['case_CNV_ids']))
            used_cnv_bt = get_cnvs_by_id(cohorts, used_cnvs)

            # Update pvalues and significance for windows in current region
            wids_to_update = list(regions[rid]['sig_windows'].keys())
            update_pvalues(wids_to_update, pvals, sig_df, sig_bt, cohorts, 
                           used_cnvs, used_cnv_bt, hpo_p_cutoffs, args.secondary_p_cutoff, 
                           args.min_nominal, retest_min_or_lower, args.model, 
                           args.max_cnv_size, args.min_case_cnvs, args.control_hpo)
            sig_windows = sig_df[sig_df.iloc[:, 4:].apply(any, axis=1)].iloc[:, 0:4]
            sig_bt = pbt.BedTool.from_dataframe(sig_windows)
            update_regions([rid], regions, sig_df, orig_sig_df, sig_bt, pvals, 
                           args.pad_sentinel)

        print('--No additional sentinel windows identified; region fully refined.\n',
              file=logfile)

    # Sort & name final MCRs
    all_mcrs = sorted(all_mcrs, key=lambda x: (int(x['chr']), int(x['start']), int(x['end'])))
    for i in range(0, k): 
        all_mcrs[i]['id'] = '_'.join([args.prefix, 'min_credible_region', str(all_mcrs[i]['chr']), str(i + 1)])

    # Calculate final association statistics for minimal critical regions
    msg = '* Finished refining {:,} minimal credible regions from {:,} query regions. ' + \
          'Now calculating conditional association statistics.\n'
    print(msg.format(len(all_mcrs), len(regions)), file=logfile)
    final_stats = get_final_stats(all_mcrs, cohorts, orig_sig_df, orig_pvals,
                                  args.control_hpo, args.model)

    # Write final results to output files
    write_associations(final_stats, assocs_fout, args.cnv_type)
    write_regions(final_stats, regions_fout, args.cnv_type)


if __name__ == '__main__':
    main()
