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
import numpy as np
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
        sig_genes[gene]['best_hpo'] = best_hpo

    return pvals_df, sig_genes


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

    # Report results of loaded genes
    msg = '  * Loaded {:,} {}'
    print(msg.format(len(txbt), 'transcripts'), file=logfile)
    print(msg.format(len(exonbt), 'exons'), file=logfile)
    print(msg.format(len(genes), 'gene symbols') + '\n', file=logfile)

    return gtfbt, txbt, exonbt, genes, transcripts, cds_dict


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

    import pdb; pdb.set_trace()

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

    return cnvs, genes_per_cnv, cnvs_per_gene


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
            cnvs, genes_per_cnv, cnvs_per_gene \
                = load_cnvs(cnv_bed, cnv_type, exonbt, cds_dict, min_cds_ovr, 
                            max_genes, control_hpo, cohort, pad_controls)
            hpos = load_hpos(hpo_file)
            cohort_info[cohort] = {'cnvs' : cnvs,
                                   'genes_per_cnv' : genes_per_cnv,
                                   'cnvs_per_gene' : cnvs_per_gene,
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


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('genes', help='BED file of blocks of genes to refine.')
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
        args.prefix = args.cnv_type
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

    # Extract relevant data for significant genes from input GTF
    gtfbt, txbt, exonbt, genes, transcripts, cds_dict \
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

    import pdb; pdb.set_trace()



if __name__ == '__main__':
    main()

