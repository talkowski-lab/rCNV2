#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Identify, cluster, and refine all significant segments per HPO from rCNV sliding window analysis
"""


from os import path
import gzip
import csv
import pybedtools as pbt
import numpy as np
import pandas as pd
from scipy.stats import norm
import string
import networkx as nx
from sklearn.cluster import KMeans
from operator import itemgetter
from itertools import combinations
from collections import Counter
import argparse
from sys import stdout


def format_stat(x, is_neg_log10=False, na_val=0):
    """
    Helper function to format & convert p-values as needed
    """

    if x == 'NA':
        return float(na_val)
    else:
        if is_neg_log10:
            return 10 ** -float(x)
        else:
            return float(x)


def get_sig_label(primary_p, secondary_p, n_nominal, primary_q, primary_p_cutoff, 
                  secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                  secondary_or_nominal=True, fdr_q_cutoff=0.05,
                  secondary_for_fdr=False):
    """
    Checks if a window should be considered exome-wide or FDR significant
    """

    # Run all comparisons
    primary_p_is_sig = (primary_p < primary_p_cutoff)
    secondary_p_is_sig = (secondary_p < secondary_p_cutoff)
    n_nominal_is_sig = (n_nominal >= n_nominal_cutoff)
    fdr_q_is_sig = (primary_q < fdr_q_cutoff)

    # Determine secondary criteria
    if secondary_p_is_sig and n_nominal_is_sig:
        secondary_is_sig = True
    elif secondary_or_nominal and (secondary_p_is_sig or n_nominal_is_sig):
        secondary_is_sig = True
    else:
        secondary_is_sig = False

    # First consider genome-wide significance
    if primary_p_is_sig and secondary_is_sig:
        return 'GWS'
    
    # Second consider FDR significance
    elif fdr_q_is_sig and secondary_for_fdr and secondary_is_sig:
        return 'FDR'
    elif fdr_q_is_sig and not secondary_for_fdr:
        return 'FDR'

    # Otherwise, non-significant
    else:
        return 'NS'


def ci2se(ci):
    """
    Converts a tuple of (lower, upper) confidence interval bounds to standard error
    """

    ci = sorted(ci)

    return (ci[1] - ci[0]) / (2 * 1.96)


def make_cs(refine_res, cs_val=0.95, cs_merge_buffer=200000, 
            cov_df=None, jac_cutoff=0.8, sig_wids=[]):
    """
    Merge windows into credible set based on sum of ranked PIPs
    Pads CS with ±cs_merge_buffer; recommended to set to size of one full bin
    """

    cs_windows = []

    vals = pd.DataFrame.from_dict(refine_res, orient='index')
    vals = vals.sort_values(by='PIP', ascending=False)

    cs_sum = 0
    for i in range(len(vals)):
        cs_windows.append(vals.index[i])
        cs_sum += vals.PIP[i]

        if cs_sum >= cs_val:
            break

    # Convert window IDs back to pbt.BedTool and merge to nonredundant intervals
    # Include all windows with sufficient covariance with any significant window from credset
    if cov_df is not None:
        cs_sig_windows = set(sig_wids).intersection(set(cs_windows))
        best_cov = cov_df.loc[cov_df.index.isin(cs_sig_windows), :].max()
        cs_windows = list(set(cs_windows).union(set(best_cov.index[best_cov >= jac_cutoff])))
        credset_bt = pbt.BedTool('\n'.join([x.replace('_', '\t') for x in cs_windows]), 
                                 from_string=True).sort().merge()

    # Otherwise, use fixed merge distance
    else:
        credset_bt = pbt.BedTool('\n'.join([x.replace('_', '\t') for x in cs_windows]), 
                                 from_string=True).sort().merge(d=cs_merge_buffer)
    
    credset_coords = [[x.chrom, x.start, x.end] for x in credset_bt]
    
    return credset_coords, credset_bt, cs_windows


def refine(window_priors, window_info, null_variance=0.42 ** 2, cs_val=0.95, 
           cs_merge_buffer=200000, cov_df=None, jac_cutoff=0.8):
    """
    Refines a set of windows given their prior probability and effect size estimates
    Inputs:
        window_priors : dict of windows with prior probabilities
        window_info : dict of window association stats as processed by process_hpo()
        null_se : float, variance of OR estimates under the null. By default,
                  this value is set to a 5% chance that ORs are > 2, as suggested
                  by Wakefield 2009
    """

    windows = list(window_priors.keys())

    if len(windows) == 1:
        refine_res = {windows[0] : {'ABF' : None, 'PIP' : 1}}
    else:
        refine_res = {}

        # Compute ABF per window
        for window in windows:
            # From Wakefield, 2009
            theta = window_info[window]['lnOR']
            se = ci2se((window_info[window]['lnOR_lower'], window_info[window]['lnOR_upper']))
            V = se ** 2
            if V > 0:
                zsq = (theta ** 2) / V
                W = null_variance
                ABF = np.sqrt((V+W) / V) * np.exp((-zsq / 2) * (W / (V+W)))

                # Wakefield 2009 formulates BF relative to H0. We need to invert to 
                # obtain evidence & posterior for H1 (i.e., true non-zero effect)
                # However, we also are only testing for a _positive_ effect in cases,
                # so we will only invert ABF if theta >= 0. This is necessary to
                # prevent sites with enrichments in controls being scored with high ABF
                if theta >= 0:
                    ABF = 1 / ABF
    
            else:
                ABF = 0

            refine_res[window] = {'ABF' : ABF}

        # Compute PIP per window as fraction of total BF, adjusted for prior probs
        # As per Mahajan 2018 (T2D refining paper, Nat. Genet.)
        posteriors = {}
        for window in windows:
            posteriors[window] = refine_res[window]['ABF'] * window_priors[window]
        posterior_sum = np.sum(list(posteriors.values()))
        for window in windows:
            refine_res[window]['PIP'] = posteriors[window] / posterior_sum

    # Define credible set as sum of ranked PIPs >= credset_val
    sig_wids = [k for k, v in window_info.items() if v['gw_sig'] or v['fdr_sig']]
    credset_coords, credset_bt, credset_windows \
        = make_cs(refine_res, cs_val, cs_merge_buffer, cov_df, jac_cutoff, sig_wids)

    return refine_res, credset_coords, credset_bt, credset_windows


def parse_stats(stats_in, primary_p_cutoff, p_is_neg_log10=True, 
                secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                secondary_or_nominal=True, fdr_q_cutoff=0.05,
                secondary_for_fdr=False, sig_only=False, keep_windows=None,
                refine_secondary=False):
    """
    Parse all association stats for a single phenotype
    """

    stats_dict = {}    

    if path.splitext(stats_in)[1] in '.gz .bgz .bgzip'.split():
        csvin = gzip.open(stats_in, 'rt')
    else:
        csvin = open(stats_in)
    reader = csv.reader(csvin, delimiter='\t')


    for chrom, start, end, n_nominal, top_cohort, excluded_cohorts, case_freq, \
        control_freq, lnOR, lnOR_lower, lnOR_upper, zscore, primary_p, primary_q, \
        secondary_lnOR, secondary_lnOR_lower, secondary_lnOR_upper, \
        secondary_zscore, secondary_p, secondary_q \
        in reader:

        # Skip header line
        if chrom.startswith('#'):
            continue

        # Set window id
        window = '_'.join([chrom, start, end])

        # If optioned, restrict on window name
        if keep_windows is not None:
            if window not in keep_windows:
                continue

        # Clean up window data
        primary_p = format_stat(primary_p, p_is_neg_log10, 1)
        primary_q = format_stat(primary_q, p_is_neg_log10, 1)
        secondary_p = format_stat(secondary_p, p_is_neg_log10, 1)
        n_nominal = int(n_nominal)
        sig_label = get_sig_label(primary_p, secondary_p, n_nominal, primary_q, 
                                  primary_p_cutoff, secondary_p_cutoff, 
                                  n_nominal_cutoff, secondary_or_nominal, 
                                  fdr_q_cutoff, secondary_for_fdr)
        gw_sig = False
        fdr_sig = False
        if sig_label == 'GWS':
            gw_sig = True
            fdr_sig = True
        elif sig_label == 'FDR':
            fdr_sig = True
        case_freq = format_stat(case_freq)
        control_freq = format_stat(control_freq)
        if refine_secondary:
            use_lnOR = format_stat(secondary_lnOR)
            use_lnOR_lower = format_stat(secondary_lnOR_lower)
            use_lnOR_upper = format_stat(secondary_lnOR_upper)
            use_zscore = format_stat(secondary_zscore)
        else:
            use_lnOR = format_stat(lnOR)
            use_lnOR_lower = format_stat(lnOR_lower)
            use_lnOR_upper = format_stat(lnOR_upper)
            use_zscore = format_stat(zscore)
        window_stats = {'case_freq' : case_freq, 'control_freq' : control_freq,
                        'lnOR' : use_lnOR, 'lnOR_lower' : use_lnOR_lower,
                        'lnOR_upper' : use_lnOR_upper, 'zscore' : use_zscore,
                        'primary_p' : primary_p, 'primary_q' : primary_q, 
                        'secondary_p' : secondary_p, 'n_nominal' : n_nominal, 
                        'gw_sig' : gw_sig, 'fdr_sig' : fdr_sig}

        # Store window association stats
        if sig_only:
            if sig_label in 'GWS FDR'.split():
                window_bt = pbt.BedTool('\t'.join([chrom, start, end, window]), 
                                      from_string=True)
                window_stats['window_bt'] = window_bt
                stats_dict[window] = window_stats
        else:
            window_bt = pbt.BedTool('\t'.join([chrom, start, end, window]), 
                                  from_string=True)
            window_stats['window_bt'] = window_bt
            stats_dict[window] = window_stats

    csvin.close()

    return stats_dict


def process_hpo(hpo, stats_in, primary_p_cutoff, p_is_neg_log10=True, 
                secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                secondary_or_nominal=True, fdr_q_cutoff=0.05, 
                secondary_for_fdr=False, block_merge_dist=1000000, 
                block_prefix='window_block', null_variance=0.42 ** 2, 
                refine_secondary=False, cs_val=0.95):
    """
    Loads & processes all necessary data for a single phenotype
    Returns a dict with the following entries:
        sig_windows : dict of sig windows, with each entry corresponding to stats
                      for a single significant window
        all_windows : dict of sig windows + all windows within block_merge_dist with
                      stats per window
        blocks : pbt.BedTool of all clustered windows to be refined
    """

    hpo_info = {'blocks' : {}}

    # First pass: parse data for significant windows only
    hpo_info['sig_windows'] = parse_stats(stats_in, primary_p_cutoff, p_is_neg_log10, 
                                          secondary_p_cutoff, n_nominal_cutoff, 
                                          secondary_or_nominal, fdr_q_cutoff,
                                          secondary_for_fdr, sig_only=True, 
                                          refine_secondary=refine_secondary)

    # Second pass: parse data for all windows within block_merge_dist of sig_windows
    if len(hpo_info['sig_windows']) > 0:
        # Make bt of significant windows
        sig_window_bts = [g['window_bt'] for g in hpo_info['sig_windows'].values()]
        if len(sig_window_bts) > 1:
            sig_windows_bt = sig_window_bts[0].cat(*sig_window_bts[1:], postmerge=False).sort()
        else:
            sig_windows_bt = sig_window_bts[0]

        # Intersect sig windows with all windows
        all_windows_bt = pbt.BedTool(stats_in).cut(range(3)).sort()
        nearby_windows_bt = all_windows_bt.closest(sig_windows_bt.sort(), d=True).\
                               filter(lambda x: int(x[-1]) > -1 and \
                                                int(x[-1]) <= block_merge_dist).\
                               saveas()
        nearby_windows = ['_'.join([x.chrom, str(x.start), str(x.end)]) for x in nearby_windows_bt]

        # Gather window stats
        hpo_info['all_windows'] = parse_stats(stats_in, primary_p_cutoff, p_is_neg_log10, 
                                              secondary_p_cutoff, n_nominal_cutoff, 
                                              secondary_or_nominal, fdr_q_cutoff,
                                              secondary_for_fdr, sig_only=False,
                                              keep_windows=nearby_windows, 
                                              refine_secondary=refine_secondary)

        # Cluster significant windows into blocks to be refined
        window_bts = [g['window_bt'] for g in hpo_info['all_windows'].values()]
        if len(window_bts) > 1:
            windows_bt = window_bts[0].cat(*window_bts[1:], postmerge=False).sort()
        else:
            windows_bt = window_bts[0]
        blocks = windows_bt.merge(d=block_merge_dist, c=4, o='distinct')

        # Assign block IDs
        blocks_wids = {'_'.join([hpo, block_prefix, str(k)]) : block \
                       for k, block in enumerate(blocks)}

        # Perform initial refinment of each block with flat prior
        for block_id, block in blocks_wids.items():
            windows = block[3].split(',')
            window_priors = {window : 1 / len(windows) for window in windows}
            refine_res, credset_coords, credset_bt, credset_windows \
                = refine(window_priors, hpo_info['all_windows'], null_variance,
                         cs_val, block_merge_dist)
            hpo_info['blocks'][block_id] = {'coords' : block,
                                            'refine_res' : refine_res,
                                            'credset_coords' : credset_coords,
                                            'credset_bt' : credset_bt,
                                            'credset_windows' : credset_windows}

    # If no windows are significant, add empty placeholder dict for all windows
    else:
        hpo_info['all_windows'] = {}

    return hpo_info


def load_all_hpos(statslist, secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                  secondary_or_nominal=True, fdr_q_cutoff=0.05, 
                  secondary_for_fdr=False, block_merge_dist=200000, 
                  block_prefix='window_block', refine_secondary=False, 
                  cs_val=0.95):
    """
    Wrapper function to process each HPO with process_hpo()
    Returns a dict with one entry per HPO
    """

    hpo_data = {}

    with open(statslist) as infile:
        reader = csv.reader(infile, delimiter='\t')
        for hpo, stats_in, pval, in reader:
            print('Loading data from {}...'.format(hpo))
            primary_p_cutoff = float(pval)
            hpo_data[hpo] = process_hpo(hpo, stats_in, primary_p_cutoff, 
                                        p_is_neg_log10=True, 
                                        secondary_p_cutoff=secondary_p_cutoff, 
                                        n_nominal_cutoff=n_nominal_cutoff, 
                                        secondary_or_nominal=secondary_or_nominal,
                                        fdr_q_cutoff=fdr_q_cutoff,
                                        secondary_for_fdr=secondary_for_fdr,
                                        block_merge_dist=block_merge_dist,
                                        block_prefix=block_prefix,
                                        refine_secondary=refine_secondary,
                                        cs_val=cs_val)

    return hpo_data


def calc_cnv_cov(cnvbed, hpo_data, cnv, frac=0.5, max_search_dist=20000000):
    """
    Compute a CNV covariance matrix for cases from each HPO per chromosome
    """

    cnv_cov = {}
    cnvbt_orig = pbt.BedTool(cnvbed)
    contigs = set([x.chrom for x in cnvbt_orig])

    # Iterate over each HPO
    for hpo, hdat in hpo_data.items():
        print('Computing covariance matrixes for {}...'.format(hpo))

        # Make single bedtool of all windows per contig
        wbt_dict = {contig : {'all_wids' : set()} for contig in contigs}
        for wid in hdat['all_windows'].keys():
            contig = wid.split('_')[0]
            wbt_dict[contig]['all_wids'].add(wid)
        for contig in contigs:
            wbt_str = ''
            for wid in wbt_dict[contig]['all_wids']:
                wbt_str += '\t'.join(wid.split('_') + [wid]) + '\n'
            wbt_dict[contig]['wbt'] = pbt.BedTool(wbt_str, from_string=True)

        # Filter CNVs by HPO and CNV type
        cnvbt = cnvbt_orig.filter(lambda x: hpo in x[5])
        if cnv != 'NS':
            cnvbt = cnvbt.filter(lambda x: x[4] == cnv).saveas()

        # Make covariance matrix of all by all windows per chromosome
        cov_dfs = {}
        for contig in contigs:
            
            # Filter CNVs and windows to contig of interest
            cnvbt_contig = cnvbt.filter(lambda x: x.chrom == contig)
            wbt_contig = wbt_dict[contig]['wbt']
            all_contig_wids = wbt_dict[contig]['all_wids']

            # Make dict mapping window ID to dict of set(CNV ids)
            cnvs_per_window = {wid : set() for wid in all_contig_wids}
            for hit in cnvbt_contig.intersect(wbt_contig, wa=True, wb=True, F=frac):
                cnvid = hit[3]
                wid = hit[-1]
                cnvs_per_window[wid].add(cnvid)            
            
            # Compute covarance for all pairs of windows
            cov_dfs[contig] = pd.DataFrame(columns=all_contig_wids)
            for wid_a in all_contig_wids:
                jac_l = []
                
                # If first window has no CNVs, Jaccard index = 0 for all mates
                cnvs_a = cnvs_per_window[wid_a]
                if len(cnvs_a) == 0:
                    cov_dfs[contig].loc[wid_a] = [0.0] * len(all_contig_wids)
                    continue

                for wid_b in all_contig_wids:
                    # If the Jaccard index has already been computed, 
                    # can copy value across matrix diagonal
                    if wid_b in cov_dfs[contig].index:
                        jac_l.append(cov_dfs[contig].loc[wid_b, wid_a])
                        continue

                    # If second window has no CNVs, Jaccard index = 0
                    cnvs_b = cnvs_per_window[wid_b]
                    if len(cnvs_b) == 0:
                        jac_l.append(0.0)
                        continue

                    # Otherwise, compute Jaccard index as long as windows are
                    # closer than max_search_dist apart
                    mid_a = np.mean([int(x) for x in wid_a.split('_')[1:]])
                    mid_b = np.mean([int(x) for x in wid_b.split('_')[1:]])
                    if np.abs(mid_b - mid_a) > max_search_dist:
                        jac_l.append(0.0)
                    else:
                        jac_l.append(len(cnvs_a.intersection(cnvs_b)) / len(cnvs_a.union(cnvs_b)))

                cov_dfs[contig].loc[wid_a] = jac_l

        cnv_cov[hpo] = cov_dfs

    return cnv_cov


def merge_blocks_by_cov(hpo_data, cnv_cov, jac_cutoff=0.8, 
                        block_merge_dist=200000, cs_val=0.95, 
                        block_prefix='window_block'):
    """
    Merge blocks from the same HPO by CNV covariance
    """

    for hpo, hdat in hpo_data.items():

        all_bids = list(hdat['blocks'].keys())

        # Construct graph of all blocks
        G = nx.Graph()
        for bid in all_bids:
            G.add_node(bid)

        # Add edges between pairs of nodes where at least one pair of windows from 
        # their credible intervals shares CNV cov >= jac_cutoff
        for bid_a in all_bids:
            chrom_a = hdat['blocks'][bid_a]['credset_coords'][0][0]
            
            for bid_b in all_bids:
                chrom_b = hdat['blocks'][bid_b]['credset_coords'][0][0]
                
                # Only process nonredundant block pairs on the same chromosome
                if bid_a == bid_b or chrom_a != chrom_b:
                    continue

                cov_df = cnv_cov[hpo][chrom_a]
                wids_a = hdat['blocks'][bid_a]['credset_windows']
                wids_b = hdat['blocks'][bid_b]['credset_windows']
                cov_df = cov_df.loc[cov_df.index.isin(wids_a), 
                                    cov_df.columns.isin(wids_b)]
                best_jac = cov_df.max().max()
                if best_jac >= jac_cutoff:
                    G.add_edge(bid_a, bid_b)

        # Collapse all subgraphs of two or more nodes
        k = 0
        for cluster in nx.connected_components(G):
            if len(cluster) > 1:
                k += 1
                new_bid = '_'.join([hpo, block_prefix, 'merged', str(k)])

                # Take union of all windows
                windows = set()
                for bid in cluster:
                    windows.update(hdat['blocks'][bid]['refine_res'].keys())

                # Update refinement
                window_priors = {window : 1 / len(windows) for window in windows}
                refine_res, credset_coords, credset_bt, credset_windows = \
                    refine(window_priors, hdat['all_windows'], cs_val=cs_val, 
                           cs_merge_buffer=block_merge_dist)

                # Determine maximum significance level of any window in credible set
                if any([hdat['all_windows'][wid]['gw_sig'] for wid in credset_windows]):
                    credset_max_sig = 'genome_wide'
                elif any([hdat['all_windows'][wid]['fdr_sig'] for wid in credset_windows]):
                    credset_max_sig = 'FDR'
                else:
                    credset_max_sig = 'not_significant'

                # Add new merged block to hpo_data
                block = credset_bt.merge(d=int(10e10))
                hpo_data[hpo]['blocks'][new_bid] = \
                    {'coords' : block,
                    'refine_res' : refine_res,
                    'credset_coords' : credset_coords,
                    'credset_bt' : credset_bt,
                    'credset_windows' : credset_windows,
                    'credset_max_sig' : credset_max_sig}

                # Remove all original member blocks
                for bid in cluster:
                    hpo_data[hpo]['blocks'].pop(bid)

    return hpo_data


def clump_windows(cov_df, sig_wids, jac_cutoff=0.2):
    """
    inputs:
        cov_df : an input matrix of CNV covariance for pairs of windows
        jac_cutoff : minimum Jaccard index to treat windows as non-independent
    outputs:
        a list of clumps of window IDs
    """

    # Make graph of all windows
    all_wids = list(cov_df.columns)
    wg = nx.Graph()
    for wid in all_wids:
        wg.add_node(wid)

    # Annotate edges with Jaccard index if >= jac_cutoff
    for sig_wid in sig_wids:
        for other_wid in all_wids:
            jac = cov_df.loc[sig_wid, other_wid]
            if jac >= jac_cutoff:
                wg.add_edge(sig_wid, other_wid)
                wg.edges[sig_wid, other_wid]['jac'] = jac

    clumps = []
    for subg in nx.connected_components(wg):
        if len(set(sig_wids).intersection(set(subg))) > 1:
            clumps.append(list(subg))
    
    return clumps


def split_blocks_by_cov(hpo_data, cnv_cov, jac_cutoff=0.2, block_prefix='window_block', 
                        null_variance=0.42 ** 2, cs_val=0.95):
    """
    Split blocks based on CNV covariance
    """

    for hpo, hdat in hpo_data.items():

        kn = 0
        sig_wids = list(hdat['sig_windows'].keys())
        orig_bids = list(hdat['blocks'].keys())

        for block_id in orig_bids:
            bdat = hdat['blocks'][block_id]
            block_wids = list(bdat['refine_res'].keys())
            sig_block_wids = list(set(sig_wids).intersection(set(block_wids)))
            chrom = bdat['credset_coords'][0][0]
            cov_df = cnv_cov[hpo][chrom]
            # Note: only consider covariance between sig windows when assessing independence
            # to avoid chaining of single-linkage between intermediate blocks
            cov_df = cov_df.loc[cov_df.index.isin(sig_block_wids), 
                                cov_df.columns.isin(sig_block_wids)]
            window_clumps = clump_windows(cov_df, sig_block_wids, jac_cutoff)

            # If multiple independent CNV clumps are found, split into independent blocks
            if len(window_clumps) > 1:

                for i in range(len(window_clumps)):

                    # Collect data corresponding to split block
                    windows = window_clumps[i]
                    windows_bt = pbt.BedTool('\n'.join(['\t'.join(x.split('_') + [x]) for x in windows]), 
                                             from_string=True).sort()
                    coords = windows_bt.merge(c=4, o='distinct')
                    window_priors = {window : 1 / len(windows) for window in windows}
                    refine_res, credset_coords, credset_bt, credset_windows \
                        = refine(window_priors, hdat['all_windows'], null_variance,
                                 cs_val, cs_merge_buffer=0)
                    kn += 1
                    new_block_id = '_'.join([hpo, block_prefix, 'split', str(kn)])

                    # Determine maximum significance level of any window in credible set
                    if any([hdat['all_windows'][wid]['gw_sig'] for wid in credset_windows]):
                        credset_max_sig = 'genome_wide'
                    elif any([hdat['all_windows'][wid]['fdr_sig'] for wid in credset_windows]):
                        credset_max_sig = 'FDR'
                    else:
                        credset_max_sig = 'not_significant'

                    # Update hpo_data
                    hpo_data[hpo]['blocks'][new_block_id] = \
                        {'coords' : coords,
                         'refine_res' : refine_res,
                         'credset_coords' : credset_coords,
                         'credset_bt' : credset_bt,
                         'credset_windows' : credset_windows,
                         'credset_max_sig' : credset_max_sig}

                # Remove original block after splitting
                hpo_data[hpo]['blocks'].pop(block_id)

    return hpo_data


def assign_or_quantiles(hpo_data, n_or_bins=1):
    """
    Assign all associations into quantiles based on effect size
    """

    # Gather effect size estimate of most significant window per block
    lnors = {}
    for hpo, hdat in hpo_data.items():
        for bid, bdat in hdat['blocks'].items():
            best_p = 1
            best_wid = None
            best_lnor = 0
            for wid in list(bdat['refine_res'].keys()):
                w_p = hdat['all_windows'][wid]['primary_p']
                if w_p < best_p:
                    best_p = w_p
                    best_wid = wid
                    best_lnor = hdat['all_windows'][wid]['lnOR']
            lnors[bid] = {'lnor' : best_lnor, 'hpo' : hpo}

    # Assign blocks to quantiles
    quants = np.floor(n_or_bins * np.argsort([x['lnor'] for x in lnors.values()]) / len(lnors))
    qdict = {a : int(b) for a, b in zip(lnors.keys(), quants)}
    for bid in qdict.keys():
        hpo_data[lnors[bid]['hpo']]['blocks'][bid]['lnor_quantile'] = qdict[bid]

    return hpo_data


def estimate_null_variance_basic(hpo_data, Wsq, dev_hpos=[], n_or_bins=1, 
                                 split_gw_fdr=False):
    """
    Estimates null variance from average of all significant windows and best window per block
    """

    # Compute 2 null variance estimates for refining:
    # 1. Mean of all significant windows
    # 2. Mean of all top windows (one per block)
    vardict_sig = {hpo : {'gw' : {i : [] for i in range(n_or_bins + 1)},
                          'fdr' : {i : [] for i in range(n_or_bins + 1)}} \
                   for hpo in hpo_data.keys()}
    vardict_best = {hpo : {'gw' : {i : [] for i in range(n_or_bins + 1)},
                           'fdr' : {i : [] for i in range(n_or_bins + 1)}} \
                    for hpo in hpo_data.keys()}

    # Collect variance estimates from each significant block
    for hpo, dat in hpo_data.items():
        for bdat in dat['blocks'].values():

            # Get block-wide odds ratio quantile and significance level
            lnor_q = bdat['lnor_quantile']
            gw_sig = any([dat['all_windows'][wid]['gw_sig'] for wid \
                          in bdat['refine_res'].keys()])
            if gw_sig:
                w_sig = 'gw'
            else:
                w_sig = 'fdr'
            if not split_gw_fdr:
                gw_sig = True

            # Get odds ratio for each significant window within the block 
            sig_wids = [wid for wid in bdat['refine_res'].keys() \
                        if any([dat['all_windows'][wid]['gw_sig'],
                                dat['all_windows'][wid]['fdr_sig']])]
            sig_lnors = {wid: dat['all_windows'][wid]['lnOR'] for wid in sig_wids}
            sig_vars = {wid : (float(lnor) / 1.96) ** 2 for wid, lnor in sig_lnors.items()}

            # Get odds ratio of single window with strongest P-value
            bpvals = [(window, dat['all_windows'][window]['primary_p']) \
                      for window in sig_wids]
            best_wid = sorted(bpvals, key=lambda x: x[1])[0][0]
            best_var = sig_vars[best_wid]

            # Update overall variance estimate collectors
            vardict_sig[hpo][w_sig][lnor_q] += list(sig_vars.values())
            vardict_best[hpo][w_sig][lnor_q].append(best_var)

    # Summarize variance estimates across strata
    for sig in 'gw fdr'.split():
        for i in range(n_or_bins):

            # Developmental HPOs
            var_sig_dev_all = [vardict_sig[h][sig][i] for h in dev_hpos]
            var_sig_dev = np.nanmean([v for sub in var_sig_dev_all for v in sub])
            var_best_dev_all = [vardict_best[h][sig][i] for h in dev_hpos]
            var_best_dev = np.nanmean([v for sub in var_best_dev_all for v in sub])
            for hpo in dev_hpos:
                Wsq[hpo][sig][i] += [var_sig_dev, var_best_dev]

            # All other HPOs
            other_hpos = list(set(hpo_data.keys()).difference(set(dev_hpos)))
            var_sig_other_all = [vardict_sig[h][sig][i] for h in other_hpos]
            var_sig_other = np.nanmean([v for sub in var_sig_other_all for v in sub])
            var_best_other_all = [vardict_best[h][sig][i] for h in other_hpos]
            var_best_other = np.nanmean([v for sub in var_best_other_all for v in sub])
            for hpo in other_hpos:
                Wsq[hpo][sig][i] += [var_sig_other, var_best_other]

    return Wsq


def estimate_null_variance_gs(gs_lists, statslist, Wsq, single_gs_hpo=False, 
                              n_or_bins=1):
    """
    Estimates null variance from the average of a list of known causal windows
    """


    statspaths = {h : p for h, p in [x.rstrip().split('\t')[:2] \
                  for x in open(statslist).readlines()]}
    with gzip.open(list(statspaths.values())[0], 'rt') as ex_statfile:
        statscols = ex_statfile.readline().rstrip().split('\t')

    # Estimate null variance for each entry in gs_lists
    for gspath in gs_lists:
        for hpo, statspath in statspaths.items():
            # Intersect sumstats for phenotype with GS regions
            gsdf = pbt.BedTool(statspath).\
                       intersect(pbt.BedTool(gspath), u=True, f=1.0).\
                       to_dataframe(names=statscols)
            gsdf['window'] = gsdf[['#chr', 'start', 'end']].astype(str).\
                                 aggregate('_'.join, axis=1)

            # Read effect sizes per window and convert to mean variance
            stats = gsdf.loc[:, 'window meta_lnOR'.split()].\
                         rename(columns={'meta_lnOR' : 'lnOR'})
            gs_var = np.nanmean((stats.lnOR.astype(float) / 1.96) ** 2)

            # Update Wsq estimates for all sig. and effect size quantiles
            if single_gs_hpo:
                for hpo in Wsq.keys():
                    for sig in 'gw fdr'.split():
                        for i in range(n_or_bins):
                            Wsq[hpo][sig][i].append(gs_var)
                break
            else:
                for sig in 'gw fdr'.split():
                    for i in range(n_or_bins):
                        Wsq[hpo][sig][i].append(gs_var)

    return Wsq


def estimate_null_variance(hpo_data, gs_list, statslist, dev_hpos=[], 
                           n_or_bins=1, split_gw_fdr=False, single_gs_hpo=False):
    """
    Estimate null variance for fine-mapping, optionally splitting on several covariates
    """

    # Build null variance hierarchy
    Wsq = {hpo : {'gw' : {i : [] for i in range(n_or_bins)}, 
                  'fdr' : {i : [] for i in range(n_or_bins)}} \
           for hpo in hpo_data.keys()}

    # Estimate null variance using 2+ approaches
    #   1. all significant windows
    #   2. most significant window per block
    #   3+. known causal regions (optional; can be multiple lists)
    Wsq = estimate_null_variance_basic(hpo_data, Wsq, dev_hpos, n_or_bins, split_gw_fdr)
    if gs_list is not None:
        with open(gs_list) as gsf:
            Wsq = estimate_null_variance_gs([x.rstrip() for x in gsf.readlines()], 
                                            statslist, Wsq, single_gs_hpo, n_or_bins)

    # Prune nans from all null variance estimates
    for hpo in Wsq.keys():
        for sig in 'gw fdr'.split():
            for i in range(n_or_bins):
                var_list = Wsq[hpo][sig][i]
                Wsq[hpo][sig][i] = list(np.array(var_list)[~np.isnan(var_list)])

    return Wsq


def output_null_var_tsv(Wsq, outfile):
    """
    Write table of null variance estimates to .tsv file
    """

    with open(outfile, 'w') as fout:
        header = '\t'.join('#HPO sig odds_ratio_quantile null_variance'.split())
        fout.write(header + '\n')

        for hpo in Wsq.keys():
            for sig in 'gw fdr'.split():
                for i, var_ests in Wsq[hpo][sig].items():
                    out_fmt = '\t'.join([hpo, sig, str(i), '{}']) + '\n'
                    for v in var_ests:
                        fout.write(out_fmt.format(v))


def update_refine(hpo_data, Wsq, split_gw_fdr=False, cs_val=0.95, 
                  block_merge_dist=200000, cnv_cov=None, jac_cutoff=0.8):
    """
    Update initial refinement results (flat prior) with null variances
    If multiple Wsqs are provided, uses Bayesian model averaging across them
    """

    for hpo, hdat in hpo_data.items():
        for block_id, bdat in hdat['blocks'].items():
            windows = list(bdat['refine_res'].keys())
            sig_wids = [w for w in windows if w in hdat['sig_windows'].keys()]
            window_priors = {window : 1 / len(windows) for window in windows}

            # Get block-wide odds ratio quantile and significance level
            lnor_q = bdat['lnor_quantile']
            gw_sig = any([hdat['all_windows'][wid]['gw_sig'] for wid in windows])
            if gw_sig or not split_gw_fdr:
                b_sig = 'gw'
            else:
                b_sig = 'fdr'

            # Subset to covariance matrix from this chromosome, if optioned
            if cnv_cov is not None:
                chrom = bdat['credset_coords'][0][0]
                cov_df = cnv_cov[hpo][chrom]
            else:
                cov_df = None
            
            # Compute ABFs and PIPs at each null variance estimate
            bma_input = []
            for W in Wsq[hpo][b_sig][lnor_q]:
                refine_res, credset_coords, credset_bt, credset_windows \
                    = refine(window_priors, hpo_data[hpo]['all_windows'], W, 
                             cs_val, block_merge_dist, cov_df, jac_cutoff)
                bma_input.append(refine_res)

            # Average ABFs for each window
            refine_res = {w : {} for w in windows}
            for window in windows:
                refine_res[window]['ABF'] \
                    = np.nanmean([s[window]['ABF'] for s in bma_input])

            # Recompute PIPs according to averaged ABFs
            ABF_sum = np.nansum([x['ABF'] for x in refine_res.values()])
            for window in windows:
                refine_res[window]['PIP'] = refine_res[window]['ABF'] / ABF_sum
            
            # Recompute credible set based on BMA PIPs
            credset_coords, credset_bt, credset_windows = \
                make_cs(refine_res, cs_val, cov_df=cov_df, jac_cutoff=jac_cutoff, 
                        sig_wids=sig_wids)

            # Determine maximum significance level of any window in credible set
            if any([hpo_data[hpo]['all_windows'][wid]['gw_sig'] for wid in credset_windows]):
                credset_max_sig = 'genome_wide'
            elif any([hpo_data[hpo]['all_windows'][wid]['fdr_sig'] for wid in credset_windows]):
                credset_max_sig = 'FDR'
            else:
                credset_max_sig = 'not_significant'

            hpo_data[hpo]['blocks'][block_id] = {'refine_res' : refine_res,
                                                 'credset_coords' : credset_coords,
                                                 'credset_bt' : credset_bt,
                                                 'credset_windows' : credset_windows,
                                                 'credset_max_sig' : credset_max_sig}

    return hpo_data


def get_cytobands(bt, cyto_bed):
    """
    Return cytoband nomenclature for all cytobands overlapped by bt
    Note: assumes all entries in bt are on the same chromosome
    """

    chrom = bt[0].chrom
    bands = sorted(list(set([x[-2] for x in bt.intersect(cyto_bed, wb=True)])))
    if len(bands) > 1:
        bandrange = '{}{}-{}'.format(chrom, bands[0], bands[-1])
    else:
        bandrange = chrom + bands[0]

    return bandrange


def rename_blocks(block_dict, cyto_bed, block_prefix):
    """
    Rename segments using cytoband nomenclature
    Input: dict of { old_id : pbt.BedTool of credible intervals }
    Returns a dict mapping { old_id : new_id }
    """

    alpha_map = {i : l for i, l in enumerate(string.ascii_uppercase)}

    rename_dict = {}

    # Map cytobands to all blocks
    cytomap = {}
    for old_id, bt in block_dict.items():
        cytomap[old_id] = get_cytobands(bt, cyto_bed)
    band_counter = {b : 0 for b in set(list(cytomap.values()))}

    # Rename blocks according to new scheme
    for old_id in block_dict.keys():
        bandrange = cytomap[old_id]
        k = Counter(cytomap.values()).get(bandrange, 0)
        if k < 2:
            new_id = '_'.join([block_prefix, bandrange])
        else:
            new_id = '_'.join([block_prefix, bandrange, 
                               alpha_map[band_counter[bandrange]]])
            band_counter[bandrange] += 1
        rename_dict[old_id] = new_id
    
    return rename_dict


def rename_all_blocks(hpo_data, cyto_bed, block_prefix):
    """
    Rename segments using cytoband nomenclature
    """

    for hpo, hdat in hpo_data.items():

        block_dict = {bid : bd['credset_bt'] for bid, bd in hdat['blocks'].items()}

        rename_dict = rename_blocks(block_dict, cyto_bed, '_'.join([hpo, block_prefix]))

        for old_id, new_id in rename_dict.items():
            hpo_data[hpo]['blocks'][new_id] = hpo_data[hpo]['blocks'][old_id]
            hpo_data[hpo]['blocks'].pop(old_id)
    
    return hpo_data


def iv_mean(values, variances, conf=0.95):
    """
    Returns inverse-variance weighted mean of values and conf% confidence interval
    """

    # Note: must restrict to values with positive, non-zero variance
    keep_idxs = np.where(np.array(variances) > 0)[0].tolist()
    values = itemgetter(*keep_idxs)(values)
    if type(values) is not tuple:
        values = (values, )
    variances = itemgetter(*keep_idxs)(variances)
    if type(variances) is not tuple:
        variances = (variances, )

    weights = [1 / v for v in variances]

    wsum = np.nansum(weights)

    numerator = np.nansum([x / v for x, v in zip(values, variances)])

    ivm = numerator / wsum

    pooled_se = np.sqrt(1 / wsum)

    ci_dist = norm.ppf(conf) * pooled_se

    return ivm, (ivm - ci_dist, ivm + ci_dist)


def output_assoc_bed(hpo_data, cyto_bed, outfile, cnv='NS'):
    """
    Format final list of credible sets with summary statistics
    """
    
    cols = 'chr start_min end_max credible_set_id cnv hpo sig_level ' + \
           'cytoband mean_control_freq mean_case_freq ' + \
           'pooled_ln_or pooled_ln_or_ci_lower pooled_ln_or_ci_upper ' + \
           'best_pvalue n_cred_intervals cred_interval_coords cred_intervals_size'
    outfile.write('#' + '\t'.join(cols.split()) + '\n')

    for hpo, hdat in hpo_data.items():
        for block_id, binfo in hpo_data[hpo]['blocks'].items():
            # Get basic credible set info
            n_cred = len(binfo['credset_coords'])
            cred_coords = ['{}:{}-{}'.format(x[0], x[1], x[2]) for x in binfo['credset_coords']]
            cred_size = np.nansum([x.length for x in binfo['credset_bt']])
            chrom = str(binfo['credset_coords'][0][0])
            start = str(np.nanmin(binfo['credset_bt'].to_dataframe().start))
            end = str(np.nanmax(binfo['credset_bt'].to_dataframe().end))
            windows = sorted(list(set([w for w in binfo['credset_windows'] if w in hdat['all_windows']])))
            wdat = hpo_data[hpo]['all_windows']
            control_freq = np.nanmean([wdat[w]['control_freq'] for w in windows])
            case_freq = np.nanmean([wdat[w]['case_freq'] for w in windows])
            cytoband = get_cytobands(binfo['credset_bt'], cyto_bed)
            sig_level = binfo['credset_max_sig']

            # Get pooled effect size as inverse-variance weighted mean of all windows
            n_windows = len(windows)
            if n_windows > 1:
                wors = [wdat[w]['lnOR'] for w in windows]
                wvars = [ci2se((wdat[w]['lnOR_lower'], wdat[w]['lnOR_upper'])) ** 2 \
                         for w in windows]
                lnor, lnor_ci = iv_mean(wors, wvars)
                lnor_lower, lnor_upper = sorted(lnor_ci)

            else:
                lnor = wdat[windows[0]]['lnOR']
                lnor_lower = wdat[windows[0]]['lnOR_lower']
                lnor_upper = wdat[windows[0]]['lnOR_upper']

            best_p = np.nanmin([wdat[w]['primary_p'] for w in windows])
            if best_p == 0:
                best_z = np.nanmax([wdat[w]['zscore'] for w in windows])
                best_p = norm.sf(best_z)

            # Write credset stats to file
            outline = '\t'.join([chrom, start, end, block_id, cnv, hpo, sig_level, cytoband])
            outnums_fmt = '\t{:.3E}\t{:.3E}\t{:.3}\t{:.3}\t{:.3}\t{:.3E}'
            outline += outnums_fmt.format(control_freq, case_freq, lnor, lnor_lower, 
                                          lnor_upper, best_p)
            outline += '\t' + '\t'.join([str(n_cred), ';'.join(cred_coords), 
                                         str(cred_size)]) + '\n'
            if sig_level != 'not_significant':
                outfile.write(outline)

    outfile.close()


def cluster_credsets(hpo_data, block_merge_dist=200000, kmeans_subclustering=False):
    """
    Cluster credible sets across HPOs to collapse overlapping regions
    """

    # Pool credible sets, tagged with block ID
    pooled_creds_str = ''
    for hdat in hpo_data.values():
        for bid, binfo in hdat['blocks'].items():
            for ci in binfo['credset_bt']:
                pooled_creds_str += '\t'.join([ci.chrom, str(ci.start), 
                                               str(ci.end), bid]) + '\n'
    pooled_creds_bt = pbt.BedTool(pooled_creds_str, from_string=True).sort().saveas()
    merged_creds_bt = pooled_creds_bt.merge(c=4, o='distinct', d=block_merge_dist)

    # Build nx.Graph() of credible sets to be clustered
    G = nx.Graph()
    csids = \
        list(set([e for s in [x.split(',') for x in merged_creds_bt.to_dataframe().\
            iloc[:, 3].tolist()] for e in s]))
    G.add_nodes_from(csids)
    for x in merged_creds_bt:
        creds = x[3].split(',')
        if len(creds) > 1:
            G.add_edges_from(list(combinations(creds, 2)))

    # If optioned, attempt to split subgraphs containing multiple blocks from the same HPO
    # This corrects for daisy chaining between distinct clumps identified
    # earlier using split_blocks_by_cov()
    if kmeans_subclustering:
        v1_clusters = list(nx.connected_components(G))
        for g in v1_clusters:

            # Only need to correct clusters containing >1 independent block from the same HPO
            hpo_counts = Counter([x.split('_')[0] for x in g])
            if any(np.array(list(hpo_counts.values())) > 1):

                # Get info on all blocks in cluster
                all_bids = list(g)
                block_bt = pooled_creds_bt.filter(lambda x: x[3] in all_bids).saveas()

                # Compute matrix of coverage between all blocks
                block_cov = pd.DataFrame(columns=all_bids)
                for bid in all_bids:
                    bid_bt = block_bt.filter(lambda x: x[3] == bid).saveas()
                    cvals = {}
                    for x in block_bt.coverage(bid_bt):
                        ovr_bp = int(x[5])
                        union = bid_bt.cat(pbt.BedTool('\t'.join(x[:4]), 
                                           from_string=True))
                        union_bp = np.sum([len(x) for x in union])
                        cvals[x[3]] = ovr_bp / union_bp
                    block_cov.loc[bid] = block_cov.columns.map(cvals)

                # Initial hierarchical clustering such that k clusters are obtained,
                # where k is the max number of independent blocks from any one HPO
                kmeans = KMeans(n_clusters=max(hpo_counts.values()), random_state=2021).fit(block_cov)
                newclusters = {k : set() for k in range(max(hpo_counts.values()))}
                for k, bid in zip(kmeans.labels_, block_cov.index):
                    newclusters[k].add(bid)

                # Prune all edges between new subclusters
                for members in newclusters.values():
                    for bid in members:
                        for other_id in [x for x in all_bids if x not in members]:
                            if G.has_edge(bid, other_id):
                                G.remove_edge(bid, other_id)

    # Format each cluster in the graph of credsets
    clustered_credsets = \
        { 'clustered_region_' + str(i) : x for i, x \
            in enumerate(nx.connected_components(G))}
    
    return clustered_credsets


def define_joint_credints(hpo_data, hpo_dict, cs_val, cnv_cov=None, 
                          jac_cutoff=0.8, ncase_dict={}, return_all=True):
    """
    Function to average PIPs across all windows from two or more HPOs to define a 
    consensus set of credible intervals
    """

    # Take union of all windows evaluated in all phenotypes
    all_wids = [hpo_data[hpo]['blocks'][bid]['refine_res'].keys() \
                for bid, hpo in hpo_dict.items()]
    all_wids = set([w for l in all_wids for w in l])
    pip_df = pd.DataFrame(index=all_wids)

    # Compile list of significant window IDs
    sig_wids = set()
    for hpo in hpo_dict.values():
        for wid in all_wids:
            if wid in hpo_data[hpo]['sig_windows'].keys():
                sig_wids.add(wid)

    # Map PIPs from each HPO onto pip_df
    for bid, hpo in hpo_dict.items():
        pip_map = {k : v['PIP'] for k, v in hpo_data[hpo]['blocks'][bid]['refine_res'].items()}
        pip_df[hpo] = pip_df.index.map(pip_map)

    # Average PIPs across all HPOs
    pip_df.fillna(value=0, inplace=True)
    refine_res = pip_df.mean(axis=1).to_frame(name='PIP').to_dict(orient='index')

    # Take CNV covariance from HPO with largest sample size for assessing CNV overlap
    # If significant window is not present for largest HPO, iterate through 
    # other HPOs ordered by size until all significant windows are added
    if cnv_cov is not None:
        hpos_by_n = [x[0] for x in sorted(ncase_dict.items(), key=lambda x: x[1], reverse=True)]
        chrom = list(all_wids)[0].split('_')[0]
        for i, hpo in enumerate([h for h in hpos_by_n if h in hpo_dict.values()]):
            if i == 0:
                cov_df = cnv_cov[hpo][chrom]
                cov_df = cov_df.loc[cov_df.index.isin(sig_wids), :]
            else:
                sup_df = cnv_cov[hpo][chrom]
                sup_df = sup_df.loc[sup_df.index.isin(sig_wids), :]
                new_rows = ~sup_df.index.isin(cov_df.index)
                new_cols = ~sup_df.columns.isin(cov_df.columns)
                # Three-step update of cov_df
                # First step: update missing columns for existing rows, if any
                if any(new_cols):
                    cov_df = pd.merge(cov_df, sup_df.loc[~new_rows, new_cols], 
                                      how='left', left_index=True, right_index=True)
                # Second step: rbind new rows, if any
                if any(new_rows):
                    cov_df = cov_df.append(sup_df.loc[new_rows, :])
                # Third step: fill any NaNs with sup_df
                cov_df.update(sup_df, overwrite=False)
    else:
        cov_df = None
    cov_df.fillna(0, inplace=True)

    # Redefine credible intervals
    credset_coords, credset_bt, credset_windows = \
        make_cs(refine_res, cs_val, cov_df=cov_df, jac_cutoff=jac_cutoff, 
                sig_wids=sig_wids)

    if return_all:
        return refine_res, credset_coords, credset_bt, credset_windows
    else:
        return credset_bt


def joint_refine_all(hpo_data, final_loci, cs_val=0.95, cnv_cov=None, 
                     jac_cutoff=0.8, ncase_dict={}, cs_merge_buffer=200000):
    """
    Wrapper to call joint_refinement() on all multi-HPO clusters
    """

    for members in final_loci.values():
        n_members = len(members)
        if n_members > 1:
            hpo_dict = {cs : cs.split('_')[0] for cs in members}
            hpos = sorted(list(set(hpo_dict.values())))
            n_hpos = len(hpos)

            # Collapse blocks from the same HPO within the same cluster prior to joint refinement
            hpo_counts = Counter(hpo_dict.values())
            for hpo, n in hpo_counts.items():
                if n > 1:
                    # Collect data corresponding to merged blocks, taking average of PIPs within a single HPO
                    hdat = hpo_data[hpo]
                    hbids = [k for k, v in hpo_dict.items() if v == hpo]
                    windows = set()
                    sig_windows = set()
                    for bid in hbids:
                        bwids = list(hdat['blocks'][bid]['refine_res'].keys())
                        windows.update(bwids)
                        sig_windows.update([w for w in hdat['sig_windows'] if w in bwids])
                    windows = list(windows)
                    sig_windows = list(sig_windows)
                    chrom = windows[0].split('_')[0]
                    pip_df = pd.DataFrame(index=windows)
                    for bid in hbids:
                        pip_map = {k : v['PIP'] for k, v in hdat['blocks'][bid]['refine_res'].items()}
                        pip_df[bid] = pip_df.index.map(pip_map)
                    pip_df.fillna(value=0, inplace=True)
                    refine_res = pip_df.mean(axis=1).to_frame(name='PIP').to_dict(orient='index')
                    credset_coords, credset_bt, credset_windows = \
                        make_cs(refine_res, cs_val, cs_merge_buffer, 
                                cnv_cov[hpo][chrom], jac_cutoff, sig_windows)
                    coords = [chrom, 
                              np.nanmin([x.start for x in credset_bt]), 
                              np.nanmax([x.stop for x in credset_bt])]

                    # Determine maximum significance level of any window in credible set
                    if any([hdat['all_windows'][wid]['gw_sig'] for wid in credset_windows]):
                        credset_max_sig = 'genome_wide'
                    elif any([hdat['all_windows'][wid]['fdr_sig'] for wid in credset_windows]):
                        credset_max_sig = 'FDR'
                    else:
                        credset_max_sig = 'not_significant'

                    # Update hpo_data for first block and remove other blocks
                    for i, bid in enumerate(hbids):
                        if i == 0:
                            hpo_data[hpo]['blocks'][bid] = \
                                {'coords' : coords,
                                 'refine_res' : refine_res,
                                 'credset_coords' : credset_coords,
                                 'credset_bt' : credset_bt,
                                 'credset_windows' : credset_windows,
                                 'credset_max_sig' : credset_max_sig}
                        else:
                            hpo_data[hpo]['blocks'].pop(bid)
                            hpo_dict.pop(bid)                    

            # Joint refinement of blocks after collapsing redundant blocks (see above)
            credints_dict = {cs : hpo_data[hpo]['blocks'][cs]['credset_coords'] \
                             for cs, hpo in hpo_dict.items()}
            credints_bts = [hpo_data[hpo]['blocks'][cs]['credset_bt'] \
                            for cs, hpo in hpo_dict.items()]
            refine_res, cs_coords, cs_bt, cs_windows = \
                define_joint_credints(hpo_data, hpo_dict, cs_val, cnv_cov, \
                                      jac_cutoff, ncase_dict)

            # Update hpo_data with results from joint refinement
            new_coords = cs_bt.merge()
            for bid, hpo in hpo_dict.items():

                # Determine maximum significance level of any window in credible set
                hcs_windows = [w for w in cs_windows if w in hpo_data[hpo]['all_windows'].keys()]
                if any([hpo_data[hpo]['all_windows'][wid]['gw_sig'] for wid in hcs_windows]):
                    credset_max_sig = 'genome_wide'
                elif any([hpo_data[hpo]['all_windows'][wid]['fdr_sig'] for wid in hcs_windows]):
                    credset_max_sig = 'FDR'
                else:
                    credset_max_sig = 'not_significant'

                # Update data
                hpo_data[hpo]['blocks'][bid] = \
                    {'coords' : new_coords,
                     'refine_res' : refine_res,
                     'credset_coords' : cs_coords,
                     'credset_bt' : cs_bt,
                     'credset_windows' : cs_windows,
                     'credset_max_sig' : credset_max_sig}

    return hpo_data


def prune_redundant_blocks(hpo_data):
    """
    Searches for blocks from the same HPO with overlapping credible intervals
    and retains just one
    """

    for hpo, hdat in hpo_data.items():
        if len(hdat['blocks']) == 0:
            continue

        # First step: make a graph of all blocks where edges indicate overlapping credible intervals
        G = nx.Graph()
        G.add_nodes_from(hdat['blocks'].keys())
        cs_bt_strs = []
        for bid, bdat in hdat['blocks'].items():
            cs_bt_strs += ['{}\t{}\t{}\t{}\n'.format(*x, bid) for x in bdat['credset_coords']]
        cs_bt = pbt.BedTool(''.join(cs_bt_strs), from_string=True)
        for hit in cs_bt.sort().merge(c=4, o='distinct'):
            bids = hit[3].split(',')
            if len(bids) > 1:
                for bid_a in bids:
                    for bid_b in bids:
                        if bid_a != bid_b:
                            G.add_edge(bid_a, bid_b)

        # Second step: resolve subgraphs with multiple nodes
        for g in nx.connected_components(G):
            if len(g) > 1:
                # Gather evidence for each block (significance level and size)
                criteria = {bid : (hdat['blocks'][bid]['credset_max_sig'], 
                            np.sum([x.length for x in hdat['blocks'][bid]['credset_bt']])) \
                            for bid in g}
                # Keep blocks with higher significance level (GW over FDR)
                # Break ties by taking larger block
                criteria = {k : v for k, v in sorted(criteria.items(), 
                                                     key=lambda x: x[1][1], 
                                                     reverse=True)}
                criteria = {k : v for k, v in sorted(criteria.items(), 
                                                     key=lambda x: x[1][0].lower(), 
                                                     reverse=True)}
                for i, bid in enumerate(criteria.keys()):
                    if i > 0:
                        hpo_data[hpo]['blocks'].pop(bid)

    return hpo_data


def output_loci_bed(hpo_data, final_loci, cyto_bed, outfile, ncase_dict, cnv='NS', 
                    block_prefix=None, joint_credsets=False, cs_val=0.95,
                    cnv_cov=None, jac_cutoff=0.8):
    """
    Format final list of collapsed credible sets and compute pooled summary statistics
    """

    if block_prefix is None:
        if cnv == 'NS':
            region_id_prefix = 'merged_segment'
        else:
            region_id_prefix = 'merged_{}_segment'.format(cnv)
    else:
        region_id_prefix = 'merged_{}_segment'.format(block_prefix)

    cols = 'chr start_min end_max region_id cnv best_sig_level cytoband ' + \
           'pooled_control_freq pooled_case_freq ' + \
           'pooled_ln_or pooled_ln_or_ci_lower pooled_ln_or_ci_upper ' + \
           'min_ln_or max_ln_or n_hpos hpos n_constituent_assocs constituent_assocs ' + \
           'n_cred_intervals cred_interval_coords cred_intervals_size'
    outfile.write('#' + '\t'.join(cols.split()) + '\n')

    out_presort = {}

    for k, members in enumerate(final_loci.values()):
        # Get basic information for credsets in final cluster
        n_members = len(members)
        hpo_dict = {cs : cs.split('_')[0] for cs in members}
        hpos = sorted(list(set(hpo_dict.values())))
        n_hpos = len(hpos)
        credints_dict = {cs : hpo_data[hpo]['blocks'][cs]['credset_coords'] for cs, hpo in hpo_dict.items()}
        credints_bts = [hpo_data[hpo]['blocks'][cs]['credset_bt'] for cs, hpo in hpo_dict.items()]
        if n_members > 1:
            if joint_credsets:
                credints_bt = \
                    define_joint_credints(hpo_data, hpo_dict, cs_val, cnv_cov, 
                                          jac_cutoff, ncase_dict, return_all=False)
            else:
                credints_bt = credints_bts[0].cat(*credints_bts[1:], postmerge=False).sort().merge()
        else:
            credints_bt = credints_bts[0].sort().merge()
        n_credints = len(credints_bt)
        credints_size = np.nansum([x.length for x in credints_bt])
        credints_coords = ['{}:{}-{}'.format(x.chrom, x.start, x.end) for x in credints_bt]
        
        # Get region-level basic information
        chrom = credints_bt[0].chrom
        start = str(np.nanmin(credints_bt.to_dataframe().start))
        end = str(np.nanmax(credints_bt.to_dataframe().end))
        cytoband = get_cytobands(credints_bt, cyto_bed)
        region_id = '_'.join([region_id_prefix, cytoband, str(k)])

        # Summarize HPO-specific information pooled across all windows from each
        # contributing credset (note: *not* all windows for all merged cred intervals)
        windows_dict = {}
        for bid, hpo in hpo_dict.items():
            windows_dict[hpo] = [w for w in hpo_data[hpo]['blocks'][bid]['credset_windows'] \
                                 if w in hpo_data[hpo]['all_windows'].keys()]
        # Compute pooled control & case frequencies as mean weighted by np.sqrt(N_cases)
        control_freq_dict, case_freq_dict = {}, {}
        for hpo, windows in windows_dict.items():
            control_freq_dict[hpo] = \
                np.nanmean([hpo_data[hpo]['all_windows'][w]['control_freq'] for w in windows])
            case_freq_dict[hpo] = \
                np.nanmean([hpo_data[hpo]['all_windows'][w]['case_freq'] for w in windows])
        control_freq = np.nanmean(list(control_freq_dict.values()))
        case_weights = [np.sqrt(ncase_dict[hpo]) for hpo in case_freq_dict.keys()]
        case_freq = np.average(list(case_freq_dict.values()), weights=case_weights)
        # Compute pooled effect size as inverse-variance weighted average
        lnor_means, lnor_cis = {}, {}
        for hpo in hpos:
            wdat = hpo_data[hpo]['all_windows']
            hlnors = [wdat[w]['lnOR'] for w in windows_dict[hpo]]
            hvars = [ci2se((wdat[w]['lnOR_lower'], wdat[w]['lnOR_upper'])) ** 2 \
                         for w in windows_dict[hpo]]
            hlnor, hlnor_ci = iv_mean(hlnors, hvars)
            lnor_means[hpo] = hlnor
            lnor_cis[hpo] = sorted(hlnor_ci)
        min_lnor = np.nanmin(list(lnor_means.values()))
        max_lnor = np.nanmax(list(lnor_means.values()))
        lnor, lnor_ci = \
            iv_mean(list(lnor_means.values()), 
                    [ci2se(tuple(ci)) ** 2 for ci in lnor_cis.values()])
        # Get best significance level from any window
        sig_levels = [hpo_data[hpo]['blocks'][bid].get('credset_max_sig') \
                      for bid, hpo in hpo_dict.items()]
        if 'genome_wide' in sig_levels:
            best_sig_level = 'genome_wide'
        elif 'FDR' in sig_levels:
            best_sig_level = 'FDR'
        else:
            best_sig_level = 'not_significant'

        # Prepare to write region stats to file
        out_front = '\t'.join([chrom, start, end])
        out_back = '\t'.join([cnv, best_sig_level, cytoband])
        outnums_fmt = '\t{:.3E}\t{:.3E}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}'
        out_back += outnums_fmt.format(control_freq, case_freq, lnor, lnor_ci[0], 
                                       lnor_ci[1], min_lnor, max_lnor)
        out_back += '\t' + '\t'.join([str(n_hpos), ';'.join(hpos), 
                                      str(n_members), ';'.join(sorted(members)),
                                      str(n_credints), ';'.join(credints_coords),
                                      str(credints_size)]) + '\n'
        if best_sig_level != 'not_significant':
            out_presort[(int(chrom), int(start), int(end))] = [out_front, region_id, out_back]

    # Final renaming of blocks with cytoband nomenclature using rename_blocks()
    block_dict = {v[1] : pbt.BedTool(v[0], from_string=True) for v in out_presort.values()}
    rename_dict = rename_blocks(block_dict, cyto_bed, region_id_prefix)

    # Iterate over sorted blocks and write to file
    for i, key in enumerate(sorted(out_presort.keys())):
        outline = '\t'.join([out_presort[key][0], 
                             rename_dict[out_presort[key][1]], 
                             out_presort[key][2]])
        outfile.write(outline)

    outfile.close()


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('statslist', help='tsv of metadata per phenotype. Three ' +
                        'required columns: HPO, path to meta-analysis stats, ' +
                        'and primary P-value cutoff.')
    parser.add_argument('hpos_by_cohort', help='tsv of sample sizes per HPO.')
    parser.add_argument('--secondary-p-cutoff', help='Maximum secondary P-value to ' + 
                        'consider as significant. [default: 1]', default=1, type=float)
    parser.add_argument('--cnv', help='Indicate CNV type. [default: NS]', default='NS')
    parser.add_argument('--min-nominal', help='Minimum number of individual cohorts ' + 
                        'required to be nominally significant to consider ' +
                        'significant. [default: 0]', default=0, type=int)
    parser.add_argument('--secondary-or-nominal', dest='secondary_or_nom', 
                        help='Allow windows to meet either --secondary-p-cutoff ' +
                        'or --min-nominal, but do not require both. ' +
                        '[default: require both]', default=False, action='store_true')
    parser.add_argument('--fdr-q-cutoff', help='Maximum FDR Q-value to ' + 
                        'consider as FDR significant. [default: 0.05]', 
                        default=0.05, type=float)
    parser.add_argument('--secondary-for-fdr', action='store_true', 
                        default=False, help='Apply sample secondary and/or min. ' +
                        'nominal criteria when evaluating FDR significance. ' +
                        '[default: apply no secondary criteria to FDR segments]')
    parser.add_argument('--cnv-bed', help='BED file of all CNVs from all cohorts. ' +
                        'Used for computing CNV covariance between pairs of windows.')
    parser.add_argument('--credible-sets', dest='cs_val', type=float, default=0.95,
                        help='Credible set value. [default: 0.95]')
    parser.add_argument('--joint-credset-definition', action='store_true', 
                        default=False, dest='joint_credsets', help='After ' +
                        'clustering and refinement, add a final step of joint ' +
                        'credible set definition across all HPOs associated at ' +
                        'each locus [default: define credsets per HPO per locus]')
    parser.add_argument('--search-distance', help='Load non-significant windows within ' +
                        'this distance of at least one significant window. [default: 1Mb]', 
                        default=1000000, type=int)
    parser.add_argument('--refine-pad', help='Distance to pad each significant window ' +
                        'prior to refinement. [default: 200kb]', default=200000, 
                        type=int)
    parser.add_argument('--block-jaccard', help='If --cnv-bed is provided, subdivide ' +
                        'blocks of significant windows if their maximum Jaccard index of ' +
                        'overlapping CNVs is not greater than this value. [default: 0.2]', 
                        default=0.2, type=float, dest='block_jac')
    parser.add_argument('--window-jaccard', help='If --cnv-bed is provided, clump ' +
                        'windows into credible set coordinates if their Jaccard index of ' +
                        'overlapping CNVs with a credible set is at least as large ' +
                        'as this value. [default: 0.8]', 
                        default=0.8, type=float, dest='window_jac')
    parser.add_argument('--known-causal-loci-list', help='.tsv list of paths to ' +
                        '.bed lists of known causal loci. Used for estimating null ' +
                        'variance. Can be specified multiple times. [default: ' +
                        'no known causal regions]', dest='gs_list')
    parser.add_argument('--single-gs-hpo', action='store_true', default=False,
                        help='Use first row in --variance-estimation-statslist ' +
                        'to estimate null variance for all HPOs. [default: ' +
                        'estimate HPO-specific null variance]')
    parser.add_argument('--variance-estimation-statslist', help='.tsv of paths to ' +
                        'meta-analysis per phenotype to use with --known-causal-loci-list, ' +
                        'if different from main statslist. [default: use main statslist]',
                        dest='var_est_statslist')
    parser.add_argument('--developmental-hpos', help='.txt file listing which ' +
                        'HPOs to treat as developmental when estimating null ' +
                        'variance. [default: pool all HPOs for null variance]')
    parser.add_argument('--n-effect-size-bins', dest='n_or_bins', type=int, default=1,
                        help='Number of effect size quantiles to evaluate when ' +
                        'estimating null variance [default: do not split ' +
                        'associations by effect size]')
    parser.add_argument('--estimate-variance-by-sig', dest='split_gw_fdr', 
                        action='store_true', help='Estimate null variance separately ' +
                        'for genome-wide significant and FDR-significant segments. ' +
                        '[default: pool all segments when estimating null variance]')
    parser.add_argument('--cytobands', help='BED of chromosome cytobands. Used ' +
                        'for naming regions. Optional. [default: name based on coordinates]')
    parser.add_argument('--refine-secondary', action='store_true', default=False,
                        help='Use secondary association statistics for segment refinement ' + 
                        '[default: use primary association stats]')
    parser.add_argument('--sig-loci-bed', help='Output BED of significant segments ' +
                        'and their overall association statistics.')
    parser.add_argument('--sig-assoc-bed', help='Output BED of significant ' +
                        'segment-phenotype pairs and their corresponding association ' +
                        'statistics.')
    parser.add_argument('--null-variance-estimates-tsv', help='Output .tsv with ' +
                        'full table of estimated null variances.', dest='null_var_out')
    parser.add_argument('--prefix', help='Prefix for naming loci & associations.')
    args = parser.parse_args()

    # Set block prefix
    block_prefix_components = ['segment']
    if args.prefix is not None:
        block_prefix_components.insert(0, args.prefix)
    if args.cnv is not None and args.cnv != 'NS':
        block_prefix_components.insert(0, args.cnv)
    else:
        block_prefix_components.insert(0, 'sig')
    block_prefix = '_'.join(block_prefix_components)

    # Process data per hpo
    hpo_data = load_all_hpos(args.statslist, args.secondary_p_cutoff, 
                             args.min_nominal, args.secondary_or_nom, 
                             args.fdr_q_cutoff, args.secondary_for_fdr,
                             args.search_distance, block_prefix, 
                             args.refine_secondary, args.cs_val)

    # Read list of developmental HPOs:
    if args.developmental_hpos is not None:
        dev_hpos = [x.rstrip() for x in open(args.developmental_hpos).readlines()]
    else:
        dev_hpos = []

    # Compute CNV covariance for all window pairs involving at least one significant window
    if args.cnv_bed is not None:
        cnv_cov = calc_cnv_cov(args.cnv_bed, hpo_data, args.cnv)
    else:
        cnv_cov = None

    # Collapse blocks with high CNV covariance prior to updated refinement
    hpo_data = merge_blocks_by_cov(hpo_data, cnv_cov, args.window_jac, 
                                   args.refine_pad, args.cs_val, block_prefix)

    # Group blocks by effect size quantile
    hpo_data = assign_or_quantiles(hpo_data, args.n_or_bins)

    # Estimate null variance
    if args.var_est_statslist is not None:
        var_est_statslist = args.var_est_statslist
    else:
        var_est_statslist = args.statslist
    Wsq = estimate_null_variance(hpo_data, args.gs_list, var_est_statslist, 
                                 dev_hpos, args.n_or_bins, args.split_gw_fdr,
                                 args.single_gs_hpo)

    # Update original refinement results with BMA of re-estimated null variances
    hpo_data = update_refine(hpo_data, Wsq, args.split_gw_fdr, args.cs_val, 
                             args.refine_pad, cnv_cov, args.window_jac)

    # Split blocks based on CNV covariance, if optioned
    if cnv_cov is not None:
        hpo_data = split_blocks_by_cov(hpo_data, cnv_cov, args.block_jac, block_prefix)

    # Read dict of N_case per HPO
    ncase_df = pd.read_csv(args.hpos_by_cohort, sep='\t').loc[:, '#HPO Total'.split()]
    ncase_df.index = ncase_df.iloc[:, 0]
    ncase_dict = ncase_df.drop(columns='#HPO').transpose().to_dict(orient='records')[0]

    # Joint refinement of overlapping credible intervals across HPOs, if optioned
    if args.joint_credsets:
        # First pass of clustering overlapping credible sets across HPOs
        # While attempting to split clusters such that each cluster has only one block per HPO
        clustered_loci_v1 = cluster_credsets(hpo_data, block_merge_dist=args.refine_pad, 
                                             kmeans_subclustering=True)

        # Joint refinement of clustered blocks
        hpo_data = joint_refine_all(hpo_data, clustered_loci_v1, args.cs_val, 
                                    cnv_cov, args.window_jac, ncase_dict,
                                    args.refine_pad)

        # Remove any lingering redundant blocks after joint refinement
        hpo_data = prune_redundant_blocks(hpo_data)

    # Rename all blocks according to cytobands of credible sets, if optioned
    if args.cytobands is not None:
        hpo_data = rename_all_blocks(hpo_data, args.cytobands, block_prefix)

    # Format & write final table of significant associations
    if args.sig_assoc_bed is not None:
        sig_assoc_bed = open(args.sig_assoc_bed, 'w')
        output_assoc_bed(hpo_data, args.cytobands, sig_assoc_bed, args.cnv)

    # Cluster final associations based on simple overlap between credints
    if args.joint_credsets:
        final_loci = cluster_credsets(hpo_data, block_merge_dist=0)
    else:
        final_loci = cluster_credsets(hpo_data, block_merge_dist=args.refine_pad)

    # Format & write final table of significant regions
    # No joint refinement necessary at this step because it has already been accomplished above
    if args.sig_loci_bed is not None:
        sig_loci_bed = open(args.sig_loci_bed, 'w')
        output_loci_bed(hpo_data, final_loci, args.cytobands, sig_loci_bed, 
                        ncase_dict, args.cnv, block_prefix.replace('_segment', ''),
                        False, args.cs_val, cnv_cov, args.window_jac)


if __name__ == '__main__':
    main()
