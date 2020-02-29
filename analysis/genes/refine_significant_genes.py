#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Refine significant gene associations
"""


from os import path
import gzip
import pandas as pd
import pybedtools as pbt
import csv
from operator import itemgetter
import numpy as np
import tempfile
import subprocess
from math import sqrt
from scipy.stats import norm
import argparse
from sys import stdout


def get_bed_header(bedpath):
    """
    Get header from BED file
    """

    if path.splitext(bedpath)[1] in '.gz .bgz .bgzip'.split():
        header = gzip.GzipFile(bedpath).readline().decode('utf-8').rstrip('\n')
    else:
        header = open(bedpath).readline().rstrip('\n')

    return header


def load_sig_genes(sig_df, cnv_type):
    """
    Load a dictionary of significant genes and their associated HPOs
    """

    sig_genes = {}

    # Read sig_df header to get phenotype keys
    sig_header = get_bed_header(sig_df)
    sig_cols = sig_header.rstrip().split('\t')
    sig_cols = [x.replace('.' + cnv_type, '').replace('#', '').replace('HP', 'HP:') \
                for x in sig_cols]

    # Read sig matrix and extract genes significant in at least one phenotype
    sig_df = pd.read_table(sig_df, names=sig_cols, comment='#')
    hpos = sig_df.iloc[:, 4:].columns.tolist()
    sig_genes_l = sig_df[sig_df.iloc[:, 4:].apply(any, axis=1)].iloc[:, 3].tolist()
    for gene in sig_genes_l:
        sig_hpos = sig_df.columns[4:][sig_df[sig_df['gene'] == gene].iloc[:, 4:]\
                   .to_numpy().flatten()].tolist()
        if gene not in sig_genes.keys():
            sig_genes[gene] = {'sig_hpos' : sig_hpos}
        else:
            sig_genes[gene]['sig_hpos'] = list(set(sig_genes[gene]['sig_hpos'] + sig_hpos))

    # Subset sig matrix to genes significant at least once
    sig_df = sig_df.loc[sig_df['gene'].isin(sig_genes.keys()), :]

    return sig_genes, sig_df


def load_pvalues(pvalues_bed, sig_genes, cnv_type, p_is_phred=False):
    """
    Import p-values matrix as dataframe and annotate significant genes with P-values
    """

    # Read pvalues_bed header to get phenotype keys
    pvals_header = get_bed_header(pvalues_bed)
    pval_cols = pvals_header.rstrip().split('\t')
    pval_cols = [x.replace('.' + cnv_type, '').replace('#', '').replace('HP', 'HP:') \
                 for x in pval_cols]

    # Format & transform pvalue matrix
    pvals_df = pd.read_table(pvalues_bed, names=pval_cols, comment='#')
    pvals_df = pvals_df.loc[pvals_df['gene'].isin(sig_genes.keys()), :]
    pvals = pvals_df.iloc[:, 4:]
    if p_is_phred:
        pvals = pvals.apply(lambda x: 10 ** -x, axis=0)
    pvals_df = pd.concat([pvals_df.iloc[:, 0:4], pvals], axis=1)
    
    # Annotate each significant gene with original P-value and designate top hpo
    for gene in sig_genes.keys():
        sig_genes[gene]['pvals'] = {}
        best_p = 1
        best_hpo = None
        for hpo in sig_genes[gene]['sig_hpos']:
            pval = pvals_df.loc[pvals_df['gene'] == gene, 
                                pvals_df.columns == hpo].to_numpy()[0][0]
            sig_genes[gene]['pvals'][hpo] = pval
            if pval < best_p:
                best_p = pval
                best_hpo = hpo
        sig_genes[gene]['best_p'] = best_p
        sig_genes[gene]['best_hpo'] = best_hpo

    return pvals_df, sig_genes


def load_gene_blocks(blocks_in, sig_genes, prefix):
    """
    Load BED file of significant gene blocks to refine, and annotate with HPOs & pvalues
    """

    blocks = {}

    if path.splitext(blocks_in)[1] in '.gz .bgz .bgzip'.split():
        fin = gzip.open(blocks_in, 'rt')
    else:
        fin = open(blocks_in)

    i = 0
    reader = csv.reader(fin, delimiter='\t')
    for chrom, start, end, genes in reader:
        i += 1
        block_id = '_'.join([prefix, 'gene_block', str(i)])
        genes = genes.split(',')

        # Determine sentinel gene & stats
        best_p = 1
        best_hpo = None
        sent_hpos = None
        sent_gene = None
        equal_sent_genes = []
        for gene in genes:
            hpo_gene = sig_genes[gene]['best_hpo']
            p_gene = sig_genes[gene]['pvals'][hpo_gene]
            if p_gene < best_p:
                best_p = p_gene
                best_hpo = hpo_gene
                sent_hpos = sig_genes[gene]['sig_hpos']
                sent_gene = gene
            elif p_gene == best_p:
                equal_sent_genes.append(gene)

        all_hpos = [g for l in [sig_genes[gene]['sig_hpos'] for gene in genes] for g in l]

        blocks[block_id] = {'chrom' : chrom, 'start' : start, 'end' : end, 
                            'genes' : genes,
                            'sent_gene' : sent_gene,
                            'equal_sent_genes' : equal_sent_genes,
                            'sent_best_p' : best_p,
                            'sent_best_hpo' : best_hpo,
                            'sent_hpos' : sent_hpos,
                            'block_hpos' : sorted(list(set(all_hpos)))}

    return blocks


def process_gtf(gtf_in, sig_genes, bl_list, logfile, xcov=0.3):
    """
    Read gtf, format entries, subset to significant genes, and compute various metadata
    """

    print('Loading gene models from input gtf...', file=logfile)

    gtfbt = pbt.BedTool(gtf_in)

    # Build lists of eligible gene names and transcript IDs
    genes, transcripts = [], []

    for f in gtfbt:
        if f.fields[2] == 'transcript':
            gname = f.attrs['gene_name']
            if gname in sig_genes.keys():
                tname = f.attrs['transcript_id']
                if gname not in genes:
                    genes.append(gname)
                if tname not in transcripts:
                    transcripts.append(tname)

    # Filter & clean records in gtf
    def _filter_gtf(feature):
        """
        Restrict GTF features to desired elements
        """
        if feature.fields[2] in 'exon transcript'.split() \
        and feature.attrs['gene_name'] in genes \
        and feature.attrs['transcript_id'] in transcripts:
            return True
        else:
            return False

    attrs_to_drop = 'gene_id gene_type gene_status transcript_type ' + \
                    'transcript_status transcript_name protein_id ' + \
                    'tag ccdsid havana_gene havana_transcript'
    attrs_to_drop = attrs_to_drop.split()

    def _clean_feature(feature):
        """
        Clean unnecessary fields & info from GTF features
        """
        for key in attrs_to_drop:
            if key in feature.attrs.keys():
                feature.attrs.pop(key)
        return feature

    gtfbt = gtfbt.filter(_filter_gtf).filter(_clean_feature).saveas()

    def _blacklist_genes(bt, bl, xcov):
        """
        Remove genes based on overlap with blacklist
        """

        txcov = bt.coverage(bl).filter(lambda x: x[2] == 'transcript')
        keepers = [x.attrs['gene_name'] for x in txcov if float(x.fields[-1]) < xcov]
        return bt.filter(lambda x: x.attrs['gene_name'] in keepers)

    # Filter genes based on blacklist overlap
    if bl_list is not None:
        for bl in bl_list:
            gtfbt = _blacklist_genes(gtfbt, bl, xcov).saveas()
    keepers = list(set([x.attrs['gene_name'] for x in gtfbt]))
    genes = [g for g in genes if g in keepers]
    transcripts = [g for g in transcripts if g in keepers]

    # Build dictionary of cds lengths per gene
    cds_dict = {}
    for e in gtfbt.filter(lambda x: x.fields[2] == 'exon'):
        gname = e.attrs['gene_name']
        if gname not in cds_dict.keys():
            cds_dict[gname] = e.length
        else:
            cds_dict[gname] += e.length

    # Make separate BedTools for exons and transcripts
    txbt = gtfbt.filter(lambda x: x.fields[2] == 'transcript').saveas()
    exonbt = gtfbt.filter(lambda x: x.fields[2] == 'exon').saveas()

    # Add coordinate info to sig_genes
    for gene in sig_genes.keys():
        sig_genes[gene]['txbt'] = txbt.filter(lambda x: x.attrs['gene_name'] == gene).saveas()
        sig_genes[gene]['exonbt'] = exonbt.filter(lambda x: x.attrs['gene_name'] == gene).saveas()

    # Report results of loaded genes
    msg = '  * Loaded {:,} {}'
    print(msg.format(len(txbt), 'transcripts'), file=logfile)
    print(msg.format(len(exonbt), 'exons'), file=logfile)
    print(msg.format(len(genes), 'gene symbols') + '\n', file=logfile)

    return gtfbt, txbt, exonbt, genes, transcripts, cds_dict, sig_genes


def overlap_cnvs_exons(cnvbt, exonbt, cds_dict, min_cds_ovr, max_genes):
    """
    Compute restricted overlaps of CNVs vs genes and genes vs CNVs
    """

    cnvs_per_gene = {}
    cnvs_per_gene_raw = {}
    genes_per_cnv = {}
    cnv_cds_sums = {}

    for i in cnvbt.intersect(exonbt, wo=True):
        # Get basic feature-intersection info
        cnvid = i[3]
        exinfo = i[14].strip(';').split(';')
        exinfo_name = [x for x in exinfo if x.startswith('gene_name ')]
        gene = exinfo_name[0].split(' ')[1].replace('"', '')
        ovrlen = int(i[-1])

        # Add CNV ID to list for gene, if necessary
        if gene not in cnvs_per_gene_raw.keys():
            cnvs_per_gene_raw[gene] = []
        if cnvid not in cnvs_per_gene_raw[gene]:
            cnvs_per_gene_raw[gene].append(cnvid)

        # Add bp overlap to list for CNV
        if cnvid not in cnv_cds_sums:
            cnv_cds_sums[cnvid] = {}
        if gene not in cnv_cds_sums[cnvid].keys():
            cnv_cds_sums[cnvid][gene] = ovrlen
        else:
            cnv_cds_sums[cnvid][gene] += ovrlen

    # Scale cds overlap as fraction of total gene per CNV, 
    # and add CNV/gene to list if observed overlap >= min overlap
    for cnvid in cnv_cds_sums.keys():
        for gene, ovrbp in cnv_cds_sums[cnvid].items():
            cds_ovr = ovrbp / cds_dict[gene]
            if cds_ovr >= min_cds_ovr:
                if gene not in cnvs_per_gene.keys():
                    cnvs_per_gene[gene] = [cnvid]
                else:
                    cnvs_per_gene[gene].append(cnvid)
                if cnvid not in genes_per_cnv.keys():
                    genes_per_cnv[cnvid] = [gene]
                else:
                    genes_per_cnv[cnvid].append(gene)

    # Exclude CNVs based on overlapping more than max_genes
    cnvs_to_exclude = [cid for cid, genes in genes_per_cnv.items() \
                       if len(genes) > max_genes]
    for cnvid in cnvs_to_exclude:
        if cnvid in genes_per_cnv:
            genes_per_cnv.pop(cnvid)
    for gene, cnvs in cnvs_per_gene.items():
        cnvs_per_gene[gene] = [c for c in cnvs if c not in cnvs_to_exclude]
    
    return genes_per_cnv, cnvs_per_gene


def load_cnvs(cnv_bed, cnv_type, exonbt, cds_dict, min_cds_ovr, max_genes, 
              control_hpo, rename=None, pad_controls=0):
    """
    Load & filter CNVs from BED file into pbt.BedTool
    """

    cnvs = pbt.BedTool(cnv_bed).intersect(exonbt, wa=True, u=True)
    
    if cnv_type is not 'CNV':
        cnvs = cnvs.filter(lambda x: x[4] == cnv_type)

    def _pad_control_cnv(feature, pad_controls, control_hpo):
        """
        Add a fixed distance to control CNV breakpoints
        """

        if feature[5] == control_hpo:
            feature.start = max([0, feature.start - pad_controls])
            feature.stop = feature.stop + pad_controls
        return feature

    cnvs = cnvs.each(_pad_control_cnv, pad_controls, control_hpo)

    if rename is not None:
        def _rename_cnv(cnv, rename):
            cnv.name = '_'.join([rename, cnv.name])
            return cnv
        cnvs = cnvs.each(_rename_cnv, rename=rename)

    cnvs = cnvs.saveas()

    genes_per_cnv, cnvs_per_gene = overlap_cnvs_exons(cnvs, exonbt, cds_dict, 
                                                      min_cds_ovr, max_genes)

    cnv_hpos = { c.name : c[5].split(';') for c in cnvs }

    return cnvs, genes_per_cnv, cnvs_per_gene, cnv_hpos


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


def load_cohort_info(info_file, cnv_type, exonbt, cds_dict, min_cds_ovr, 
                     max_genes, control_hpo, pad_controls=0):
    """
    Read metadata per cohort and save as dict
    """

    cohort_info = {}

    with open(info_file) as fin:
        reader = csv.reader(fin, delimiter='\t')

        for cohort, cnv_bed, hpo_file in reader:
            cnvs, genes_per_cnv, cnvs_per_gene, cnv_hpos \
                = load_cnvs(cnv_bed, cnv_type, exonbt, cds_dict, min_cds_ovr, 
                            max_genes, control_hpo, cohort, pad_controls)
            hpos = load_hpos(hpo_file)
            cohort_info[cohort] = {'cnvs' : cnvs,
                                   'genes_per_cnv' : genes_per_cnv,
                                   'cnvs_per_gene' : cnvs_per_gene,
                                   'cnv_hpos' : cnv_hpos,
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


def print_block_info(blocks, blockid, sig_genes, logfile):
    """
    Print basic gene block info prior to refinement
    """

    print('Beginning refinement of {0}...'.format(blockid), file=logfile)

    binfo = blocks[blockid]
    ngenes = len(binfo['genes'])
    if ngenes == 1:
        ven = 'gene'
    else:
        ven = 'genes'
    block_size = int(binfo['end']) - int(binfo['start'])
    n_hpos = len(binfo['block_hpos'])
    genes_w_pvals = [(gene, sig_genes[gene]['best_p']) for gene in binfo['genes']]
    sorted_genes = sorted(genes_w_pvals, key=itemgetter(1))
    genes_w_pvals_fmt = ', '.join(['{} ({:.2E})'.format(g, p) for g, p in sorted_genes])
    
    msg = '  * Block of {:,} significant {} spanning {:,} kb ({}:{:,}-{:,})'
    print(msg.format(ngenes, ven, int(round(block_size / 1000)), binfo['chrom'],
                     int(binfo['start']), int(binfo['end'])), file=logfile)

    msg = '  * Associated with {:,} phenotypes ({})'
    print(msg.format(n_hpos, ', '.join(binfo['block_hpos'])), file=logfile)

    msg = '  * Significant genes (ordered by best P-value): {}'
    print(msg.format(genes_w_pvals_fmt), file=logfile)


def get_cohort_cnvs(cohorts, in_genes, ex_genes, in_hpos, cnvs_to_ignore=[]):
    """
    Get list of CNVs per metacohort meeting the following criteria
    1. CNV is not in cnvs_to_ignore
    2. CNV has phenotype in in_hpos
    3. CNV hits in_genes
    4. CNV _does not hit_ ex_genes
    """

    cohort_cnvs = {m : [] for m in cohorts.keys()}
    
    for cohort in cohorts.keys():
        for cnv, cnvgenes in cohorts[cohort]['genes_per_cnv'].items():
            if cnv not in cnvs_to_ignore \
            and len([h for h in cohorts[cohort]['cnv_hpos'][cnv] if h in in_hpos]) > 0 \
            and len([g for g in cnvgenes if g in in_genes]) > 0 \
            and len([g for g in cnvgenes if g in ex_genes]) == 0:
                cohort_cnvs[cohort].append(cnv)

    return cohort_cnvs


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


def meta_analysis(model, cohorts, cnv_dict, n_sample_dict):
    """
    Perform random-effects or Mantel-Haenszel meta-analysis

    Wraps subprocess.call() to the Rscript used for discovery meta-analysis
    This is necessary to ensure mathematical consistency
    """

    # Compute case/control alts
    n_control = [d['n_control'] for d in n_sample_dict.values()]
    n_case = [d['n_case'] for d in n_sample_dict.values()]
    n_control_alt = [d['n_control'] for d in cnv_dict.values()]
    n_case_alt = [d['n_case'] for d in cnv_dict.values()]
    n_case_ref = list(np.array(n_case) - np.array(n_case_alt))
    n_control_ref = list(np.array(n_control) - np.array(n_control_alt))
    
    cohort_idxs = list(range(0, len(cohorts)))

    # Wraps Rscript call in context of temporary directory for easy cleanup
    colnames = '#chr start end gene case_alt case_ref control_alt control_ref fisher_phred_p'
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
            dummy = '1 1 2 dummy'.split() + [str(x) for x in counts] + ['1']
            tmp_fin.write('\t'.join(dummy) + '\n')
            tmp_fin.close()
            m_fin.write('{}\t{}\n'.format(cohorts[i], tmp_fn))
        m_fin.close()

        # Calls Rscript
        script = '/'.join([path.dirname(path.realpath(__file__)),
                           'gene_meta_analysis.R'])
        m_fo = tmpdir + '/meta.out.tsv'
        cstr = '{} --model {} --p-is-phred {} {}'
        subprocess.call(cstr.format(script, model, m_fn, m_fo), shell=True)

        # Read results back into python
        res = pd.read_csv(m_fo, delim_whitespace=True)
        pval = 10 ** -float(res['meta_phred_p'])
        oddsratio = float(res['meta_lnOR'])
        oddsratio_ci = [float(res['meta_lnOR_lower']), 
                        float(res['meta_lnOR_upper'])]
        oddsratio_se = ((oddsratio_ci[1] - oddsratio_ci[0]) / 2) / 1.96
        chisq = float(res.iloc[:, -2])

    return oddsratio, oddsratio_ci, oddsratio_se, chisq, pval


def compare_geneset_oddsratios(cred_genes, other_genes, cred_hpos, cohorts, 
                               meta_model, control_hpo):
    """
    Compare odds ratios of cred_genes vs other_genes for cred_hpos
    """

    # Get case CNVs specific to cred genes and to other genes
    cred_case_cnvs = get_cohort_cnvs(cohorts, cred_genes, other_genes, cred_hpos)
    other_case_cnvs = get_cohort_cnvs(cohorts, other_genes, cred_genes, cred_hpos)

    # # Get case CNVs that hit (but not specific to) cred genes and to other genes
    # cred_case_cnvs = get_cohort_cnvs(cohorts, cred_genes, [], cred_hpos)
    # other_case_cnvs = get_cohort_cnvs(cohorts, other_genes, [], cred_hpos)

    # Get control CNVs that are specific to cred genes and to other genes
    cred_control_cnvs = get_cohort_cnvs(cohorts, cred_genes, other_genes, [control_hpo])
    other_control_cnvs = get_cohort_cnvs(cohorts, other_genes, cred_genes, [control_hpo])

    # # Get control CNVs that hit (but not specific to) cred genes and other genes
    # cred_control_cnvs = get_cohort_cnvs(cohorts, cred_genes, [], [control_hpo])
    # other_control_cnvs = get_cohort_cnvs(cohorts, other_genes, [], [control_hpo])

    # Summarize CNV counts
    cred_cnvs = {}
    other_cnvs = {}
    for cohort in cohorts.keys():
        cred_cnvs[cohort] = {'n_case' : len(cred_case_cnvs[cohort]),
                             'n_control' : len(cred_control_cnvs[cohort])}
        other_cnvs[cohort] = {'n_case' : len(other_case_cnvs[cohort]),
                              'n_control' : len(other_control_cnvs[cohort])}

    # Get effective number of cases and controls per cohort
    cc_counts = get_cc_counts(cohorts, cred_hpos, [control_hpo])

    # Compute meta-analysis odds ratios
    cred_or, cred_or_ci, cred_or_se, cred_chisq, cred_p \
        = meta_analysis(meta_model, list(cohorts.keys()), cred_cnvs, cc_counts)
    other_or, other_or_ci, other_or_se, other_chisq, other_p \
        = meta_analysis(meta_model, list(cohorts.keys()), other_cnvs, cc_counts)

    # Compare ORs and compute p-value
    diff = cred_or - other_or
    pooled_se = sqrt((cred_or_se ** 2) + (other_or_se ** 2))
    diff_z = diff / pooled_se
    diff_p = norm.sf(diff_z)

    if not np.isnan(cred_or):
        import pdb; pdb.set_trace()

    return cred_or, cred_or_se, other_or, other_or_se, diff_z, diff_p


def compare_oddsratio_vs_null(cred_genes, other_genes, cred_hpos, cohorts, 
                               meta_model, control_hpo):
    """
    Compare odds ratios of CNVs specific to cred_genes vs other_genes for cred_hpos
    """

    # Get case & control CNVs specific to cred genes
    cred_case_cnvs = get_cohort_cnvs(cohorts, cred_genes, other_genes, cred_hpos)
    cred_control_cnvs = get_cohort_cnvs(cohorts, cred_genes, other_genes, [control_hpo])

    # Summarize CNV counts
    cred_cnvs = {}
    for cohort in cohorts.keys():
        cred_cnvs[cohort] = {'n_case' : len(cred_case_cnvs[cohort]),
                             'n_control' : len(cred_control_cnvs[cohort])}

    # Get effective number of cases and controls per cohort
    cc_counts = get_cc_counts(cohorts, cred_hpos, [control_hpo])

    # Compute meta-analysis odds ratios
    cred_or, cred_or_ci, cred_or_se, cred_chisq, cred_p \
        = meta_analysis(meta_model, list(cohorts.keys()), cred_cnvs, cc_counts)

    return cred_or, cred_or_se, cred_chisq, cred_p


def refine_block(binfo, blockid, sig_genes, cohorts, meta_model, compare, 
                 control_hpo, logfile):
    """
    Wrapper function to conduct all refinement steps for a single block
    """

    # Order genes by best p-value
    genes_w_pvals = [(gene, sig_genes[gene]['best_p']) for gene in binfo['genes']]
    sorted_genes = sorted(genes_w_pvals, key=itemgetter(1))
    genes = [g for g, p in sorted_genes]
    ngenes = len(genes)

    # Add genes to credible set one at a time until OR is significant after Bonf. 
    # correction for number of tests performed
    p = 1
    k = 1
    threshold = 0.05
    while p >= threshold / k and k < ngenes:
        # Define genes in & out of credible set and other metadata
        cred_genes = genes[0:k]
        other_genes = genes[k:ngenes]
        # cred_hpos = [sig_genes[g]['sig_hpos'] for g in cred_genes]
        # cred_hpos = list(set([h for l in cred_hpos for h in l]))
        cred_hpos = binfo['block_hpos']

        # Compare OR of cred_genes vs other_genes for cred_hpos
        if compare == 'other':
            cred_or, cred_or_se, other_or, other_or_se, diff_z, p \
                = compare_geneset_oddsratios(cred_genes, other_genes, cred_hpos, 
                                             cohorts, meta_model, control_hpo)

            # Print logging results:
            msg = '  * Evaluating {:,}-gene credible set (OR = {:.2} ± {:.2}) vs. ' + \
                  'other {} genes (OR = {:.2} ± {:.2}): Z = {:.2}, P = {:.2E}'
            print(msg.format(len(cred_genes), cred_or, cred_or_se, len(other_genes),
                             other_or, other_or_se, diff_z, p), file=logfile)
        # Otherwise, compare OR of cred_genes vs null of OR = 1
        elif compare == 'null':
            cred_or, cred_or_se, chisq, p \
                = compare_oddsratio_vs_null(cred_genes, other_genes, cred_hpos, 
                                            cohorts, meta_model, control_hpo)

            # Print logging results:
            msg = '  * Evaluating {:,}-gene credible set (OR = {:.2} ± {:.2}) vs. ' + \
                  'null (OR = 1): R-E meta chisq = {:.2}, P = {:.2E}'
            print(msg.format(len(cred_genes), cred_or, cred_or_se, chisq, p), file=logfile)

        if np.isnan(p):
            p = 1

        k += 1

    if k < ngenes:
        print('Success!\n', file=logfile)
    else:
        print('Unable to prune any genes...\n', file=logfile)


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('blocks', help='BED file of blocks of genes to refine.')
    parser.add_argument('gtf', help='GTF of genes to consider.')
    parser.add_argument('cohort_info', help='Three-column TSV file specifying ' +
                        'cohort name, path to CNV BED file, and path to sample-' +
                        'level HPO terms.')
    parser.add_argument('pvalues', help='BED file of original genes and p-values ' +
                        'for each phenotype.')
    parser.add_argument('sig_df', help='BED file of original genes with TRUE/FALSE ' +
                        'label of significance for each phenotype.')
    parser.add_argument('associations_out', help='Path to output BED with HPO-level ' +
                        'association stats (one line per phenotype-gene pair.')
    parser.add_argument('genes_out', help='Path to output BED with gene-level ' +
                        'association stats (one line per gene.')
    parser.add_argument('--cnv-type', help='Type of CNVs to evaluate. [default: ' + 
                        'use all CNVs]', choices=['CNV', 'DEL', 'DUP'], 
                        default='CNV')
    parser.add_argument('--pad-controls', help='Distance to be added to control ' +
                        'breakpoints. [default: 0]',
                        type=float, default=0)
    parser.add_argument('--min-cds-ovr', help='Minimum coding sequence overlap ' +
                        'to consider a CNV gene-overlapping. [default: 0.2]',
                        type=float, default=0.2)
    parser.add_argument('--max-genes', help='Maximum number of genes overlapped by ' + 
                        'a CNV before being the CNV is excluded. [default: 20000]',
                        type=int, default=20000)
    parser.add_argument('-x', '--blacklist', action='append', help='BED file ' +
                        'containing regions to blacklist based on overlap with ' +
                        'genes. May be specified multiple times.')
    parser.add_argument('--xcov', type=float, help='Maximum coverage ' + 
                        'by any blacklist before excluding a gene. [default: 0.3]',
                        default=0.3)
    parser.add_argument('--model', help='Meta-analysis model to use. [default: "re"]', 
                        default='re', choices=['mh', 're'])
    parser.add_argument('--comparison', dest='compare', help='Comparison to perform ' +
                        'for refinement. [default: compare to other genes]', 
                        default='other', choices=['other', 'null'])
    parser.add_argument('--hpo-p-cutoffs', help='.tsv of p-value cutoffs per phenotype. ' + 
                        '[default: 10e-8 for all phenotypes]')
    parser.add_argument('--p-cutoff-ladder', help='.tsv of p-value cutoffs for a ' + 
                        'range of case counts. Effective case count for sentinel ' +
                        'gene will round down to nearest case sample size, if ' +
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
    parser.add_argument('--min-nominal', help='Minimum number of individual cohorts ' + 
                        'required to be nominally significant to report an independent ' +
                        'minimal credible gene. [default: 1]', 
                        default=1, type=int)
    parser.add_argument('--max-cnv-size', help='Maximum size of CNVs to include in ' +
                        'first pass when attempting to refine blocks of genes ' +
                        '[default: 3Mb]', default=3000000, type=int)
    parser.add_argument('--secondary-or-nominal', dest='secondary_or_nom', 
                        help='Allow genes to meet either --secondary-p-cutoff ' +
                        'or --min-nominal, but do not require both. ' +
                        '[default: require both]', default=False, action='store_true')
    parser.add_argument('--min-case-cnvs', help='Minimum count of CNVs required ' +
                        'per independent minimal credible gene [default: 1]', 
                        default=1, type=int)
    parser.add_argument('--control-hpo', default='HEALTHY_CONTROL', help='HPO code ' +
                        'to use for control CNV counts. [default: HEALTHY_CONTROL]')
    parser.add_argument('--pad-sentinel', help='Distance from sentinel gene to ' + 
                        'look for other significant associations. [default: 100000]', 
                        default=100000, type=int)
    parser.add_argument('--prefix', help='String to append to all refined genes. ' +
                        '[default: CNV type (specified above)]', default=None)
    parser.add_argument('--log', help='Output file for all progress logging. ' +
                        '[default: stdout]', default='stdout')

    args = parser.parse_args()

    # Postprocess args as necessary
    if args.prefix is None:
        prefix = args.cnv_type
    else:
        prefix = args.prefix
    min_or_lower = np.log(args.min_or_lower)
    retest_min_or_lower = np.log(args.retest_min_or_lower)

    # Open connections to outfiles
    assocs_fout = open(args.associations_out, 'w')
    genes_fout = open(args.genes_out, 'w')
    if args.log in 'stdout - /dev/stdout'.split():
        logfile = stdout
    else:
        logfile = open(args.log, 'w')

    # Load dictionary of significant genes, and significance/pval matrices
    sig_genes, sig_df = load_sig_genes(args.sig_df, args.cnv_type)
    pvals_df, sig_genes = load_pvalues(args.pvalues, sig_genes, args.cnv_type, 
                                       args.p_is_phred)

    # Save a copy of original p-values and sig matrix for final reporting
    orig_pvals_df = pvals_df.copy(deep=True)
    orig_sig_df = sig_df.copy(deep=True)

    # Load gene blocks to be refined
    blocks = load_gene_blocks(args.blocks, sig_genes, prefix)

    msg = 'Identified {:,} gene blocks to be refined.\n'
    print(msg.format(len(blocks)), file=logfile)

    # Extract relevant data for significant genes from input GTF
    gtfbt, txbt, exonbt, genes, transcripts, cds_dict, sig_genes \
        = process_gtf(args.gtf, sig_genes, args.blacklist, logfile, args.xcov)

    # Load cohort info, including CNVs & HPOs
    cohorts = load_cohort_info(args.cohort_info, args.cnv_type, exonbt, cds_dict,
                               args.min_cds_ovr, args.max_genes,
                               args.control_hpo, args.pad_controls)
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

    # Refine blocks one at a time
    final_blocks = {}
    used_cnvs = []
    used_cnv_bt = None
    for blockid in blocks.keys():
        print_block_info(blocks, blockid, sig_genes, logfile)
        new_block = refine_block(blocks[blockid], blockid, sig_genes, cohorts,
                                 args.model, args.compare, args.control_hpo, 
                                 logfile)

    # Clean up & flush buffer
    logfile.close()


if __name__ == '__main__':
    main()

