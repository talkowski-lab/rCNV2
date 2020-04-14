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
import argparse
from sys import stdout


def format_stat(x, is_phred=False, na_val=0):
    """
    Helper function to format & convert p-values as needed
    """

    if x == 'NA':
        return float(na_val)
    else:
        if is_phred:
            return 10 ** -float(x)
        else:
            return float(x)


def is_window_sig(primary_p, secondary_p, n_nominal, primary_p_cutoff, 
                  secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                  secondary_or_nominal=True):
    """
    Checks if a window should be considered genome-wide significant
    """

    if primary_p >= primary_p_cutoff:
        return False
    else:
        secondary = (secondary_p < secondary_p_cutoff)
        n_nom = (n_nominal >= n_nominal_cutoff)
        if secondary_or_nominal:
            if not any([secondary, n_nom]):
                return False
        elif not all([secondary, n_nom]):
            return False

    return True


def ci2se(ci):
    """
    Converts a tuple of (lower, upper) confidence interval bounds to standard error
    """

    ci = sorted(ci)

    return (ci[1] - ci[0]) / (2 * 1.96)


def make_cs(refine_res, cs_val=0.95):
    """
    Merge windows into credible set based on sum of ranked PIPs
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
    credset_bt = pbt.BedTool('\n'.join([x.replace('_', '\t') for x in cs_windows]), 
                             from_string=True).sort().merge()
    credset_coords = [[x.chrom, x.start, x.end] for x in credset_bt]
    
    return credset_coords, credset_bt, cs_windows


def refine(window_priors, window_info, null_variance=0.42 ** 2, cs_val=0.95):
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
    credset_coords, credset_bt, credset_windows = make_cs(refine_res, cs_val)

    return refine_res, credset_coords, credset_bt, credset_windows


def parse_stats(stats_in, primary_p_cutoff, p_is_phred=True, 
                secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                secondary_or_nominal=True, sig_only=False, keep_windows=None,
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


    for chrom, start, end, n_nominal, top_cohort, case_freq, control_freq, \
        lnOR, lnOR_lower, lnOR_upper, zscore, primary_p, secondary_lnOR, \
        secondary_lnOR_lower, secondary_lnOR_upper, secondary_zscore, secondary_p \
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
        primary_p = format_stat(primary_p, p_is_phred, 1)
        secondary_p = format_stat(secondary_p, p_is_phred, 1)
        n_nominal = int(n_nominal)
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

        # Store window association stats
        if sig_only:
            if is_window_sig(primary_p, secondary_p, n_nominal, primary_p_cutoff,
                           secondary_p_cutoff, n_nominal_cutoff, secondary_or_nominal):
                window_bt = pbt.BedTool('\t'.join([chrom, start, end, window]), 
                                      from_string=True)
                window_stats = {'case_freq' : case_freq, 'control_freq' : control_freq,
                                'lnOR' : use_lnOR, 'lnOR_lower' : use_lnOR_lower,
                                'lnOR_upper' : use_lnOR_upper, 'zscore' : use_zscore,
                                'primary_p' : primary_p, 'secondary_p' : secondary_p,
                                'n_nominal' : n_nominal, 'window_bt' : window_bt}
                stats_dict[window] = window_stats
        else:
            window_bt = pbt.BedTool('\t'.join([chrom, start, end, window]), 
                                  from_string=True)
            window_stats = {'case_freq' : case_freq, 'control_freq' : control_freq,
                            'lnOR' : use_lnOR, 'lnOR_lower' : use_lnOR_lower,
                            'lnOR_upper' : use_lnOR_upper, 'zscore' : use_zscore,
                            'primary_p' : primary_p, 'secondary_p' : secondary_p,
                            'n_nominal' : n_nominal, 'window_bt' : window_bt}
            stats_dict[window] = window_stats

    csvin.close()

    return stats_dict


def process_hpo(hpo, stats_in, primary_p_cutoff, p_is_phred=True, 
                secondary_p_cutoff=0.05, n_nominal_cutoff=2, 
                secondary_or_nominal=True, block_merge_dist=500000, 
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
    se_by_chrom = {}

    # First pass: parse data for significant windows only
    hpo_info['sig_windows'] = parse_stats(stats_in, primary_p_cutoff, p_is_phred, 
                                          secondary_p_cutoff, n_nominal_cutoff, 
                                          secondary_or_nominal, sig_only=True,
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
        hpo_info['all_windows'] = parse_stats(stats_in, primary_p_cutoff, p_is_phred, 
                                              secondary_p_cutoff, n_nominal_cutoff, 
                                              secondary_or_nominal, sig_only=False,
                                              keep_windows=nearby_windows, 
                                              refine_secondary=refine_secondary)

        # Cluster significant windows into blocks to be refined
        window_bts = [g['window_bt'] for g in hpo_info['all_windows'].values()]
        if len(window_bts) > 1:
            windows_bt = window_bts[0].cat(*window_bts[1:], postmerge=False).sort()
        else:
            windows_bt = window_bts[0]
        blocks = windows_bt.merge(d=block_merge_dist, c=4, o='distinct')

        # Perform initial refinment of each block with flat prior
        k = 0
        for block in blocks:
            k += 1
            block_id = '_'.join([hpo, block_prefix, str(k)])
            windows = block[3].split(',')
            window_priors = {window : 1 / len(windows) for window in windows}
            refine_res, credset_coords, credset_bt, credset_windows \
                = refine(window_priors, hpo_info['all_windows'], null_variance,
                         cs_val)
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
                  secondary_or_nominal=True, block_merge_dist=500000,
                  refine_secondary=False, cs_val=0.95):
    """
    Wrapper function to process each HPO with process_hpo()
    Returns a dict with one entry per HPO
    """

    hpo_data = {}

    with open(statslist) as infile:
        reader = csv.reader(infile, delimiter='\t')
        for hpo, stats_in, pval, in reader:
            primary_p_cutoff = float(pval)
            hpo_data[hpo] = process_hpo(hpo, stats_in, primary_p_cutoff, 
                                        p_is_phred=True, 
                                        secondary_p_cutoff=secondary_p_cutoff, 
                                        n_nominal_cutoff=n_nominal_cutoff, 
                                        secondary_or_nominal=secondary_or_nominal,
                                        block_merge_dist=block_merge_dist,
                                        refine_secondary=refine_secondary,
                                        cs_val=cs_val)

    return hpo_data


def estimate_null_variance_basic(hpo_data):
    """
    Estimates null variance per phenotype from average of all significant windows
    """

    vardict_all = {hpo : {} for hpo in hpo_data.keys()}
    vardict_sig = {hpo : [] for hpo in hpo_data.keys()}
    vardict_best = {hpo : [] for hpo in hpo_data.keys()}

    for hpo, dat in hpo_data.items():
        
        for window, gdat in dat['all_windows'].items():
            # Estimate null variance from effect size per Wakefield, AJHG, 2007
            var = (float(gdat['lnOR']) / 1.96) ** 2
            vardict_all[hpo][window] = var
            if window in dat['sig_windows'].keys():
                vardict_sig[hpo].append(var)
        
        for bdat in dat['blocks'].values():
            bpvals = [(window, dat['all_windows'][window]['primary_p']) for window \
                          in bdat['refine_res'].keys()]
            best_window = sorted(bpvals, key=lambda x: x[1])[0][0]
            vardict_best[hpo].append(vardict_all[hpo][best_window])

    # Compute 2 null variance estimates for refining:
    # 1. Mean of all significant windows
    # 2. Mean of all top windows (one per block)

    v1 = np.nanmean([x for l in vardict_sig.values() for x in l if x > 0])
    v2 = np.nanmean([x for l in vardict_best.values() for x in l if x > 0])

    return [v1, v2]


def estimate_null_variance_gs(gs_lists, statslist, refine_secondary=False):
    """
    Estimates null variance for all cases from the average of a list of known causal windows
    """

    var = []

    statspath = open(statslist).readline().rstrip().split('\t')[1]
    statsbed = pbt.BedTool(statspath)
    if path.splitext(statspath)[1] in '.gz .bgz .bgzip'.split():
        statsfin = gzip.open(statspath, 'rt')
    else:
        statsfin = open(statspath)
    statscols = statsfin.readline().rstrip().split('\t')

    for gspath in gs_lists:
        # Intersect sumstats for highest-level phenotype with GS regions
        gsdf = statsbed.intersect(pbt.BedTool(gspath), u=True, f=1.0).\
                   to_dataframe(names=statscols)
        gsdf['window'] = gsdf[['#chr', 'start', 'end']].astype(str).\
                             aggregate('_'.join, axis=1)

        # Read effect sizes per window from highest level phenotype
        if refine_secondary:
            keep_cols = 'window meta_lnOR_secondary'.split()
        else:
            keep_cols = 'window meta_lnOR'.split()
        stats = gsdf.loc[:, keep_cols].rename(columns={'meta_lnOR' : 'lnOR',
                                                       'meta_lnOR_secondary' : 'lnOR'})
        gs_vars = (stats.lnOR.astype(float) / 1.96) ** 2
        var.append(float(np.nanmean(gs_vars)))

    return var


def update_refine(hpo_data, Wsq, cs_val=0.95):
    """
    Update initial refinement results (flat prior) with null variances
    If multiple Wsqs are provided, uses Bayesian model averaging across them
    """

    # TODO: implement BMA over multiple Wsq values
    import pdb; pdb.set_trace()
    for hpo, hdat in hpo_data.items():
        for block_id, bdat in hdat['blocks'].items():
            windows = list(bdat['refine_res'].keys())
            window_priors = {window : 1 / len(windows) for window in windows}
            refine_res, credset_coords, credset_bt, credset_windows \
                = refine(window_priors, hpo_data[hpo]['all_windows'], W, cs_val)
            hpo_data[hpo]['blocks'][block_id] = {'refine_res' : refine_res,
                                                 'credset_coords' : credset_coords,
                                                 'credset_bt' : credset_bt,
                                                 'credset_windows' : credset_windows}

    return hpo_data


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


def output_assoc_bed(hpo_data, outfile, cnv='NS'):
    """
    Format final list of credible sets with summary statistics
    """
    
    cols = 'chr start_min end_max credible_set_id cnv hpo mean_control_freq mean_case_freq ' + \
           'pooled_ln_or pooled_ln_or_ci_lower pooled_ln_or_ci_upper ' + \
           'best_pvalue n_cred_intervals cred_interval_coords cred_intervals_size'
    outfile.write('#' + '\t'.join(cols.split()) + '\n')

    for hpo in hpo_data.keys():
        for block_id, binfo in hpo_data[hpo]['blocks'].items():
            # Get basic credible set info
            n_cred = len(binfo['credset_coords'])
            cred_coords = ['{}:{}-{}'.format(x[0], x[1], x[2]) for x in binfo['credset_coords']]
            cred_size = np.nansum([x.length for x in binfo['credset_bt']])
            chrom = str(binfo['credset_coords'][0][0])
            start = str(np.nanmin(binfo['credset_bt'].to_dataframe().start))
            end = str(np.nanmax(binfo['credset_bt'].to_dataframe().end))
            windows = sorted(list(set(binfo['credset_windows'])))
            wdat = hpo_data[hpo]['all_windows']
            control_freq = np.nanmean([wdat[w]['control_freq'] for w in windows])
            case_freq = np.nanmean([wdat[w]['case_freq'] for w in windows])

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

            # Write window stats to file
            outline = '\t'.join([chrom, start, end, block_id, cnv, hpo])
            outnums_fmt = '\t{:.3E}\t{:.3E}\t{:.3}\t{:.3}\t{:.3}\t{:.3E}'
            outline += outnums_fmt.format(control_freq, case_freq, lnor, lnor_lower, 
                                          lnor_upper, best_p)
            outline += '\t' + '\t'.join([str(n_cred), ';'.join(cred_coords), 
                                         str(cred_size)]) + '\n'
            outfile.write(outline)

    outfile.close()


def cluster_credsets(hpo_data):
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
    merged_creds_bt = pbt.BedTool(pooled_creds_str, from_string=True).sort().\
                          merge(c=4, o='distinct')

    import pdb; pdb.set_trace()


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
    parser.add_argument('--secondary-p-cutoff', help='Maximum secondary P-value to ' + 
                        'consider as significant. [default: 1]', default=1, type=float)
    parser.add_argument('--cnv', help='Indicate CNV type. [default: NS]', default='NS')
    parser.add_argument('--min-nominal', help='Minimum number of individual cohorts ' + 
                        'required to be nominally significant to consider ' +
                        'significant. [default: 1]', default=1, type=int)
    parser.add_argument('--secondary-or-nominal', dest='secondary_or_nom', 
                        help='Allow windows to meet either --secondary-p-cutoff ' +
                        'or --min-nominal, but do not require both. ' +
                        '[default: require both]', default=False, action='store_true')
    parser.add_argument('--credible-sets', dest='cs_val', type=float, default=0.95,
                        help='Credible set value. [default: 0.95]')
    parser.add_argument('--distance', help='Distance to pad each significant window ' +
                        'prior to refinement. [default: 1Mb]', default=1000000, 
                        type=int)
    parser.add_argument('--known-causal-loci-list', help='.tsv list of paths to ' +
                        '.bed lists of known causal loci. Used for estimating null ' +
                        'variance. Can be specified multiple times. [default: ' +
                        'no known causal regions]', dest='gs_list')
    parser.add_argument('--refine-secondary', action='store_true', default=False,
                        help='Use secondary association statistics for segment refinement ' + 
                        '[default: use primary association stats]')
    parser.add_argument('--sig-loci-bed', help='Output BED of significant segments ' +
                        'and their overall association statistics.')
    parser.add_argument('--sig-assoc-bed', help='Output BED of significant ' +
                        'segment-phenotype pairs and their corresponding association ' +
                        'statistics.')
    parser.add_argument('hpos_by_cohort', help='tsv of sample sizes per HPO.')
    args = parser.parse_args()

    # Process data per hpo
    hpo_data = load_all_hpos(args.statslist, args.secondary_p_cutoff, 
                             args.min_nominal, args.secondary_or_nom, 
                             args.distance, args.refine_secondary, args.cs_val)

    # Estimate null variance based on:
    #   1. all significant windows
    #   2. most significant window per block
    #   3. known causal regions (optional; can be multiple lists)
    Wsq = estimate_null_variance_basic(hpo_data)
    if args.gs_list is not None:
        with open(args.gs_list) as gsf:
            Wsq += estimate_null_variance_gs(gsf.read().splitlines(), args.statslist)
    Wsq = sorted(Wsq)
    print('Null variance estimates: ' + ', '.join([str(round(x, 3)) for x in Wsq]))

    # Update original refinement results with BMA of re-estimated null variances
    hpo_data = update_refine(hpo_data, Wsq, args.cs_val)

    # Format & write final table of significant associations
    if args.sig_assoc_bed is not None:
        sig_assoc_bed = open(args.sig_assoc_bed, 'w')
        output_assoc_bed(hpo_data, sig_assoc_bed, args.cnv)

    # # Cluster credible sets across HPOs
    # final_loci = cluster_credsets(hpo_data)

    # # Read dict of N_case per HPO
    # ncase_df = pd.read_csv(args.hpos_by_cohort, sep='\t').loc[:, '#HPO Total'.split()]
    # ncase_df.index = ncase_df.iloc[:, 0]
    # ncase_dict = ncase_df.drop(columns='#HPO').transpose().to_dict(orient='records')[0]

    # # Format & write final table of significant windows
    # if args.sig_loci_bed is not None:
    #     sig_loci_bed = open(args.sig_loci_bed, 'w')
    #     output_loci_bed(hpo_data, sig_loci_bed, ncase_dict, args.cnv)


if __name__ == '__main__':
    main()
