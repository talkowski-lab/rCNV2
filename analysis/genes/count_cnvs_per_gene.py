#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Count & weight CNVs overlapping genes
"""


import pybedtools
import argparse
from sys import stdout
from os import path
import subprocess
import gzip


def process_cnvs(bedpath, pad_controls, case_hpo, control_hpo):
    """
    Read CNVs & extend control CNV breakpoints by a specified distance
    """

    cnvbt = pybedtools.BedTool(bedpath)

    def _pad_control_cnv(feature, pad_controls, control_hpo):
        """
        Add a fixed distance to control CNV breakpoints
        """

        if feature[5] == control_hpo:
            feature.start = max([0, feature.start-pad_controls])
            feature.stop = feature.stop+pad_controls
        return feature

    cnvbt = cnvbt.each(_pad_control_cnv, pad_controls, control_hpo)

    def _filter_phenos(feature, case_hpo, control_hpo):
        """
        Filter CNVs based on phenotype
        """
        if case_hpo is 'KEEP_ALL_SAMPLES':
            # Keep everything if no case hpo is specified
            return True
        else:
            keep_hpos = [case_hpo, control_hpo]
            cnv_hpos = feature[5].split(';')
            if len(set(keep_hpos).intersection(set(cnv_hpos))) > 0:
                return True
            else:
                return False

    cnvbt = cnvbt.filter(_filter_phenos, case_hpo, control_hpo).saveas()

    return cnvbt


def process_gtf(gtf_in):
    """
    Read gtf, format entries, and compute various metadata
    """

    gtfbt = pybedtools.BedTool(gtf_in)

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


def overlap_cnvs_exons(cnvbt, exonbt, cds_dict):
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

    # Calculate weighted counts per gene for each CNV
    cnv_weights = {}
    for cnvid in cnv_cds_sums.keys():
        if cnvid not in cnv_weights.keys():
            cnv_weights[cnvid] = {}
        for gene, ovrbp in cnv_cds_sums[cnvid].items():
            cnv_weights[cnvid][gene] = ovrbp / cds_dict[gene]
        wsum = sum(cnv_weights[cnvid].values())
        cnv_weights[cnvid] = {gene: w / wsum for gene, w in cnv_weights[cnvid].items()}

    # Collapse counts per gene
    raw_counts = {gene: len(cnvs) for gene, cnvs in cnvs_per_gene.items()}
    weighted_counts = {}
    for cnv, weights in cnv_weights.items():
        for gene, w in weights.items():
            if gene not in weighted_counts.keys():
                weighted_counts[gene] = w
            else:
                weighted_counts[gene] += w

    return raw_counts, weighted_counts


def make_output_table(outbed, txbt, genes, cds_dict, control_counts,
                      control_weights, case_counts, case_weights):
    """
    Format master table of counts per gene and write to outbed
    """

    h_cols = '#chr start end gene cds control_count control_count_weighted ' + \
             'case_count case_count_weighted'
    header = '\t'.join(h_cols.split())
    outbed.write(header + '\n')

    for i in txbt:
        gene = i.attrs['gene_name']
        gline = [i.chrom, i.start, i.end, gene, cds_dict.get(gene, 'NA'),
                 control_counts.get(gene, 0), control_weights.get(gene, 0),
                 case_counts.get(gene, 0), case_weights.get(gene, 0)]
        outbed.write('\t'.join([str(x) for x in gline]) + '\n')


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cnvs', help='CNV BED file to compare vs windows.')
    parser.add_argument('gtf', help='GTF of genes to consider.')
    parser.add_argument('--pad-controls', help='Distance to be added to control ' +
                        'breakpoints. [default: 0]',
                        type=float, default=0.75)
    parser.add_argument('-t', '--type', help='Type of CNV to include (DEL/DUP). ' +
                        '[default: all]')
    parser.add_argument('--hpo', help='HPO term to consider for case samples. ' +
                        'If no --hpo is supplied, will count all CNVs ' +
                        'irrespective of phenotype.')
    parser.add_argument('--control-hpo', default='HEALTHY_CONTROL', help='HPO code ' +
                        'to use for control CNV counts. [default: HEALTHY_CONTROL]',
                        dest='control_hpo')
    parser.add_argument('-o', '--outbed', help='Path to output file. ' +
                        '[default: stdout]')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output BED with bgzip.')
    
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
        case_hpo = 'KEEP_ALL_SAMPLES'
    else:
        case_hpo = args.hpo

    # Process input CNV BED file
    cnvbt = process_cnvs(args.cnvs, args.pad_controls, case_hpo, args.control_hpo)

    if args.type is not None:
        if args.type != 'CNV':
            cnvbt = cnvbt.filter(lambda x: args.type in x.fields).saveas()

    # Extract relevant data from input GTF
    gtfbt, txbt, exonbt, genes, transcripts, cds_dict = process_gtf(args.gtf)

    # Intersect CNVs with exons
    case_cnvbt = cnvbt.filter(lambda x: args.control_hpo not in x[5].split(';'))
    case_counts, case_weights = overlap_cnvs_exons(case_cnvbt, exonbt, cds_dict)
    control_cnvbt = cnvbt.filter(lambda x: args.control_hpo in x[5].split(';'))
    control_counts, control_weights = overlap_cnvs_exons(control_cnvbt, exonbt, cds_dict)
    
    # Format output table and write to outfile
    make_output_table(outbed, txbt, genes, cds_dict, control_counts, 
                      control_weights, case_counts, case_weights)
    if args.outbed is not None \
    and args.outbed not in 'stdout -'.split() \
    and args.bgzip:
        subprocess.run(['bgzip', '-f', outbed_path])


if __name__ == '__main__':
    main()
