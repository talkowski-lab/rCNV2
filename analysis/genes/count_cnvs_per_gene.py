#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Count & weight CNVs overlapping genes
"""


import pybedtools as pbt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
import argparse
from sys import stdout
from os import path
import subprocess
import gzip


def process_cnvs(bedpath, pad_controls, case_hpos, control_hpo, max_size=None):
    """
    Read CNVs & extend control CNV breakpoints by a specified distance
    """

    cnvbt = pbt.BedTool(bedpath)

    def _pad_control_cnv(feature, pad_controls, control_hpo):
        """
        Add a fixed distance to control CNV breakpoints
        """

        if feature[5] == control_hpo:
            feature.start = max([0, feature.start-pad_controls])
            feature.stop = feature.stop+pad_controls
        return feature

    cnvbt = cnvbt.each(_pad_control_cnv, pad_controls, control_hpo)

    # Filter on max size, if optioned
    if max_size is not None:
        cnvbt = cnvbt.filter(lambda x: x.length <= max_size)

    def _filter_phenos(feature, case_hpos, control_hpo):
        """
        Filter CNVs based on phenotype
        """
        if 'KEEP_ALL_SAMPLES' in case_hpos:
            # Keep everything if no case hpo is specified
            return True
        else:
            keep_hpos = case_hpos + [control_hpo]
            cnv_hpos = feature[5].split(';')
            if len(set(keep_hpos).intersection(set(cnv_hpos))) > 0:
                return True
            else:
                return False

    cnvbt = cnvbt.filter(_filter_phenos, case_hpos, control_hpo).saveas()
    
    return cnvbt


def process_gtf(gtf_in, bl_list, xcov=0.3):
    """
    Read gtf, format entries, and compute various metadata
    """

    gtfbt = pbt.BedTool(gtf_in)

    # Build lists of eligible gene names and transcript IDs
    genes, transcripts = [], []

    for f in gtfbt:
        if f.fields[2] == 'transcript':
            gname = f.attrs['gene_name']
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

    return gtfbt, txbt, exonbt, genes, transcripts, cds_dict


def overlap_cnvs_exons(cnvbt, exonbt, cds_dict, weight_mode, min_cds_ovr, max_genes):
    """
    Compute fraction of CDS overlapped per gene for each CNV
    """

    cnvs_per_gene = {}
    cnv_cds_sums = {}

    for i in cnvbt.intersect(exonbt, wo=True):
        # Get basic feature-intersection info
        cnvid = i[3]
        exinfo = i[14].strip(';').split(';')
        exinfo_name = [x for x in exinfo if x.startswith('gene_name ')]
        gene = exinfo_name[0].split(' ')[1].replace('"', '')
        ovrlen = int(i[-1])

        # Add CNV ID to list for gene, if necessary
        if gene not in cnvs_per_gene.keys():
            cnvs_per_gene[gene] = []
        if cnvid not in cnvs_per_gene[gene]:
            cnvs_per_gene[gene].append(cnvid)

        # Add bp overlap to list for CNV
        if cnvid not in cnv_cds_sums:
            cnv_cds_sums[cnvid] = {}
        if gene not in cnv_cds_sums[cnvid].keys():
            cnv_cds_sums[cnvid][gene] = ovrlen
        else:
            cnv_cds_sums[cnvid][gene] += ovrlen

    # Scale cds overlap as fraction of total gene per CNV
    cnv_weights = {}
    for cnvid in cnv_cds_sums.keys():
        if cnvid not in cnv_weights.keys():
            cnv_weights[cnvid] = {}
        for gene, ovrbp in cnv_cds_sums[cnvid].items():
            cds_ovr = ovrbp / cds_dict[gene]
            if cds_ovr >= min_cds_ovr:
                cnv_weights[cnvid][gene] = cds_ovr
            # If overlap is less than required, remove cnv from cnvs_per_gene
            else:
                cnvs_per_gene[gene] = [c for c in cnvs_per_gene[gene] if c != cnvid]

    # Weight CNVs
    for cnvid in cnv_cds_sums.keys():
        if len(cnv_weights[cnvid]) > 0:
            wsum = sum(cnv_weights[cnvid].values())

            # Skip to next CNV if no weighting is specified
            if weight_mode is None:
                continue

            # Otherwise, weight genes by CDS overlap
            if weight_mode == 'weak' \
            or wsum < 1:
                cnv_weights[cnvid] = {gene: w / wsum for gene, w 
                                      in cnv_weights[cnvid].items()}
            elif weight_mode == 'strong':
                cnv_weights[cnvid] = {gene: w / (2 ** (wsum - 1)) for gene, w 
                                      in cnv_weights[cnvid].items()}
    
    # Exclude CNVs based on overlapping more than max_genes
    cnvs_to_exclude = [cid for cid in cnv_cds_sums.keys() \
                       if len(cnv_weights.get(cid, [])) > max_genes]

    # Collapse counts per gene
    raw_counts = {gene: len([x for x in cnvs if x not in cnvs_to_exclude]) \
                  for gene, cnvs in cnvs_per_gene.items()}
    weighted_counts = {}
    for cnv, weights in cnv_weights.items():
        if len(weights) > 0 and cnv not in cnvs_to_exclude:
            for gene, w in weights.items():
                if gene not in weighted_counts.keys():
                    weighted_counts[gene] = w
                else:
                    weighted_counts[gene] += w

    return raw_counts, weighted_counts, cnv_weights


def get_bayes_weights(case_cnvbt, case_cnv_weights, control_cnvbt, control_cnv_weights):
    """
    Compute weighted pseudocounts as a function of size & number of genes
    """

    # Prepare data
    case_vals = [(1, x.length, np.sum(list(case_cnv_weights.get(x.name, {}).values())), x.name)
                 for x in case_cnvbt]
    control_vals = [(0, x.length, np.sum(list(control_cnv_weights.get(x.name, {}).values())), x.name)
                 for x in control_cnvbt]
    df = pd.DataFrame(case_vals + control_vals, columns='pheno size genes cnvid'.split())

    # Fit logistic regression
    x = df.loc[:, 'size genes'.split()]
    y = df['pheno']
    model = LogisticRegression(solver='saga', random_state=0, penalty='none',
                               max_iter=1000).fit(x, y)
    probs = pd.DataFrame(model.predict_proba(x), columns='p_control p_case'.split())
    df = pd.concat([df, probs], axis=1)

    # Apply pseudocount weights to genes for each CNV
    # Pseudocounts normalized to median = 1 independently for cases & controls
    case_scalar = np.median(2 * (1 - df.loc[df['pheno'] == 1, 'p_case']))
    for cnvid in case_cnv_weights.keys():
        pcase = np.float(df.loc[df['cnvid'] == cnvid, 'p_case'])
        pcount = (2 * (1 - pcase)) / case_scalar
        for gene, frac in case_cnv_weights[cnvid].items():
            case_cnv_weights[cnvid][gene] = frac * pcount
    control_scalar = np.median(2 * df.loc[df['pheno'] == 0, 'p_case'])
    for cnvid in control_cnv_weights.keys():
        pcase = np.float(df.loc[df['cnvid'] == cnvid, 'p_case'])
        pcount = (2 * pcase) / control_scalar
        for gene, frac in control_cnv_weights[cnvid].items():
            control_cnv_weights[cnvid][gene] = frac * pcount

    # Sum weights per gene
    case_weights = {}
    for cnv, weights in case_cnv_weights.items():
        if len(weights) > 0:
            for gene, w in weights.items():
                if gene not in case_weights.keys():
                    case_weights[gene] = w
                else:
                    case_weights[gene] += w
    control_weights = {}
    for cnv, weights in control_cnv_weights.items():
        if len(weights) > 0:
            for gene, w in weights.items():
                if gene not in control_weights.keys():
                    control_weights[gene] = w
                else:
                    control_weights[gene] += w

    return case_weights, case_cnv_weights, control_weights, control_cnv_weights


def make_output_table(outbed, txbt, genes, cds_dict, control_counts,
                      control_weights, case_counts, case_weights):
    """
    Format master table of counts per gene and write to outbed
    """

    h_cols = '#chr start end gene cds control_cnvs control_cnvs_weighted ' + \
             'case_cnvs case_cnvs_weighted'
    header = '\t'.join(h_cols.split())
    outbed.write(header + '\n')

    for i in txbt:
        gene_id = i.attrs['gene_name']
        gline = [i.chrom, i.start, i.end, gene_id, cds_dict.get(gene_id, 'NA'),
                 control_counts.get(gene_id, 0), control_weights.get(gene_id, 0),
                 case_counts.get(gene_id, 0), case_weights.get(gene_id, 0)]
        gline_str = '\t'.join([str(x) for x in gline]) + '\n'
        outbed.write(gline_str)

    # Must close output file to flush buffer (end of chr22 sometimes gets clipped)
    outbed.close()


def write_annotated_cnvs(cnvbt, cnvs_out, case_cnv_weights, control_cnv_weights, 
                         control_hpo='HEALTHY_CONTROL', annotate_cds=False, 
                         max_float=6):
    """
    Format & write BED file of CNVs annotated with gene overlaps
    """

    # Write header to output file
    hcols = '#chr start end cnvid cnv phenos ngenes total_weight genes'
    if annotate_cds:
        hcols += ' cds_per_gene'
    cnvs_out.write('\t'.join(hcols.split()) + '\n')

    for cnv in cnvbt:
        cnvid = cnv.name
        if cnv.fields[-1] == control_hpo:
            cnv_dict = control_cnv_weights
        else:
            cnv_dict = case_cnv_weights
        hits = cnv_dict.get(cnvid, {})
        ngenes = len(hits.keys())
        sumweight = np.nansum(list(hits.values()))
        sorted_genes = sorted(hits.keys())
        genes = ';'.join(sorted_genes)
        outfields = cnv.fields + [str(ngenes), str(sumweight), genes]
        if annotate_cds:
            sorted_cds = [str(round(hits[g], max_float)) for g in sorted_genes]
            outfields.append(';'.join(sorted_cds))
        cnvs_out.write('\t'.join(outfields) + '\n')


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cnvs', help='CNV BED file to compare vs GTF.')
    parser.add_argument('gtf', help='GTF of genes to consider.')
    parser.add_argument('--pad-controls', help='Distance to be added to control ' +
                        'breakpoints. [default: 0]',
                        type=float, default=0)
    parser.add_argument('--weight-mode', help='Specify behavior for distributing ' +
                        'weight for multi-gene CNVs. [default: no weighting]',
                        choices=['weak', 'strong', 'bayesian'], default=None)
    parser.add_argument('--min-cds-ovr', help='Minimum coding sequence overlap ' +
                        'to consider a CNV gene-overlapping. [default: 0.2]',
                        type=float, default=0.2)
    parser.add_argument('--max-cnv-size', help='Maximum CNV size to be included. ' +
                        '[default: No limit]', type=int, default=None)
    parser.add_argument('--max-genes', help='Maximum number of genes overlapped by ' + 
                        'a CNV before being the CNV is excluded. [default: 20000]',
                        type=int, default=20000)
    parser.add_argument('--max-genes-in-cases-only', action='store_true',
                        help='Only apply --max-genes to case CNVs. [default: ' +
                        'apply to both case & control CNVs]')
    parser.add_argument('-t', '--type', help='Type of CNV to include (DEL/DUP). ' +
                        '[default: all]')
    parser.add_argument('--hpo', help='HPO term to consider for case samples. ' +
                        'If no --hpo is supplied, will count all CNVs ' +
                        'irrespective of phenotype. Can also provide a semicolon-' +
                        'delimited list of HPOs to consider at once.')
    parser.add_argument('--control-hpo', default='HEALTHY_CONTROL', help='HPO code ' +
                        'to use for control CNV counts. [default: HEALTHY_CONTROL]',
                        dest='control_hpo')
    parser.add_argument('-x', '--blacklist', action='append', help='BED file ' +
                        'containing regions to blacklist based on overlap with ' +
                        'genes. May be specified multiple times.')
    parser.add_argument('--xcov', type=float, help='Maximum coverage ' + 
                        'by any blacklist before excluding a gene. [default: 0.3]',
                        default=0.3)
    parser.add_argument('-o', '--outbed', help='Path to output file. ' +
                        '[default: stdout]')
    parser.add_argument('--cnvs-out', help='Path to output BED file for CNVs ' + 
                        'annotated with genes disrupted.')
    parser.add_argument('--annotate-cds-per-gene', default=False, action='store_true',
                        help='Append CDS overlapped per gene to CNV BED output file.')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output BED with bgzip.')
    parser.add_argument('--bgzip-cnvs-out', dest='bgzip_cnvs_out', action='store_true', 
                        help='Compress --cnvs-out BED with bgzip.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print ' + 
                        'diagnostics.')
    
    args = parser.parse_args()

    # Open connection to output file
    if args.outbed is None \
    or args.outbed in 'stdout -'.split():
        outbed = stdout
    else:
        if path.splitext(args.outbed)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outbed_path = path.splitext(args.outbed)[0]
        else:
            outbed_path = args.outbed
        outbed = open(outbed_path, 'w')

    if args.hpo is None:
        case_hpos = ['KEEP_ALL_SAMPLES']
    else:
        case_hpos = args.hpo.split(';')

    # Process input CNV BED file
    cnvbt = process_cnvs(args.cnvs, args.pad_controls, case_hpos, args.control_hpo,
                         args.max_cnv_size)

    if args.type is not None:
        if args.type != 'CNV':
            cnvbt = cnvbt.filter(lambda x: args.type in x.fields).saveas()

    # Extract relevant data from input GTF
    gtfbt, txbt, exonbt, genes, transcripts, cds_dict \
        = process_gtf(args.gtf, args.blacklist, args.xcov)
    if args.verbose:
        msg = 'Loaded {:,} {} from input gtf'
        print(msg.format(len(txbt), 'transcripts'))
        print(msg.format(len(exonbt), 'exons'))
        print(msg.format(len(genes), 'gene symbols'))

    # Intersect CNVs with exons
    case_cnvbt = cnvbt.filter(lambda x: args.control_hpo not in x[5].split(';')).saveas()
    case_counts, case_weights, case_cnv_weights \
        = overlap_cnvs_exons(case_cnvbt, exonbt, cds_dict, args.weight_mode, 
                             args.min_cds_ovr, args.max_genes)
    max_genes_controls = args.max_genes
    if args.max_genes_in_cases_only:
        max_genes_controls = 20000
    control_cnvbt = cnvbt.filter(lambda x: args.control_hpo in x[5].split(';')).saveas()
    control_counts, control_weights, control_cnv_weights \
        = overlap_cnvs_exons(control_cnvbt, exonbt, cds_dict, args.weight_mode, 
                             args.min_cds_ovr, max_genes_controls)

    # Compute Bayesian weights, if optioned
    if args.weight_mode == 'bayesian' \
    and len(case_cnv_weights) > 0 \
    and len(control_cnv_weights) > 0:
        case_weights, case_cnv_weights, control_weights, control_cnv_weights \
            = get_bayes_weights(case_cnvbt, case_cnv_weights, control_cnvbt, control_cnv_weights)
    
    # Format output main counts table and write to outfile
    make_output_table(outbed, txbt, genes, cds_dict, control_counts, 
                      control_weights, case_counts, case_weights)
    if args.outbed is not None \
    and args.outbed not in 'stdout -'.split() \
    and args.bgzip:
        subprocess.run(['bgzip', '-f', outbed_path])

    # If optioned, format annotated CNV BED file and write to args.cnvs_out
    if args.cnvs_out is not None:
        if path.splitext(args.cnvs_out)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            cnvs_outpath = path.splitext(args.cnvs_out)[0]
            bgzip_cnvs_out = True
        else:
            cnvs_outpath = args.cnvs_out
            bgzip_cnvs_out = args.bgzip_cnvs_out
        cnvs_out = open(cnvs_outpath, 'w')
        write_annotated_cnvs(cnvbt, cnvs_out, case_cnv_weights, 
                             control_cnv_weights, args.control_hpo,
                             args.annotate_cds_per_gene)
        if bgzip_cnvs_out:
            subprocess.run(['bgzip', '-f', cnvs_outpath])


if __name__ == '__main__':
    main()
