#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Count CNVs overlapping cis-regulatory blocks (CRBs)
"""


import pybedtools as pbt
import argparse
from sys import stdout
from os import path
import subprocess


def process_cnvs(bedpath, pad_controls, case_hpo, control_hpo, max_size=None):
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


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cnvs', required=True, help='CNV BED file to compare vs CRBs.')
    parser.add_argument('--crbs', required=True, help='CRB BED file.')
    parser.add_argument('--elements', required=True, help='BED file of elements per CRB.')
    parser.add_argument('--pad-controls', help='Distance to be added to control ' +
                        'breakpoints. [default: 0]',
                        type=float, default=0)
    parser.add_argument('--min-element-ovr', help='Minimum fraction of each element ' +
                        'that must be overlapped by a CNV before counting.',
                        type=float, default=1e-10)
    parser.add_argument('--min-frac-all-elements', help='Minimum fraction of all elements ' +
                        'per CRB that must be overlapped by a single CNV before counting.',
                        type=float, default=1e-10)
    parser.add_argument('--max-cnv-size', help='Maximum CNV size to be included. ' +
                        '[default: No limit]', type=int, default=None)
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
        bgzip = False
    else:
        if path.splitext(args.outbed)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outbed_path = path.splitext(args.outbed)[0]
            bgzip = True
        else:
            outbed_path = args.outbed
            bgzip = args.bgzip
        outbed = open(outbed_path, 'w')

    if args.hpo is None:
        case_hpo = 'KEEP_ALL_SAMPLES'
    else:
        case_hpo = args.hpo

    # Process input CNV BED file
    cnvbt = process_cnvs(args.cnvs, args.pad_controls, case_hpo, args.control_hpo,
                         args.max_cnv_size)

    if args.type is not None:
        if args.type != 'CNV':
            cnvbt = cnvbt.filter(lambda x: args.type in x.fields).saveas()

    # # Extract relevant data from input GTF
    # gtfbt, txbt, exonbt, genes, transcripts, cds_dict \
    #     = process_gtf(args.gtf, args.blacklist, args.xcov)
    # if args.verbose:
    #     msg = 'Loaded {:,} {} from input gtf'
    #     print(msg.format(len(txbt), 'transcripts'))
    #     print(msg.format(len(exonbt), 'exons'))
    #     print(msg.format(len(genes), 'gene symbols'))

    # # Intersect CNVs with exons
    # case_cnvbt = cnvbt.filter(lambda x: args.control_hpo not in x[5].split(';')).saveas()
    # case_counts, case_weights, case_cnv_weights \
    #     = overlap_cnvs_exons(case_cnvbt, exonbt, cds_dict, args.weight_mode, 
    #                          args.min_cds_ovr, args.max_genes)
    # max_genes_controls = args.max_genes
    # if args.max_genes_in_cases_only:
    #     max_genes_controls = 20000
    # control_cnvbt = cnvbt.filter(lambda x: args.control_hpo in x[5].split(';')).saveas()
    # control_counts, control_weights, control_cnv_weights \
    #     = overlap_cnvs_exons(control_cnvbt, exonbt, cds_dict, args.weight_mode, 
    #                          args.min_cds_ovr, max_genes_controls)

    # # Compute Bayesian weights, if optioned
    # if args.weight_mode == 'bayesian' \
    # and len(case_cnv_weights) > 0 \
    # and len(control_cnv_weights) > 0:
    #     case_weights, case_cnv_weights, control_weights, control_cnv_weights \
    #         = get_bayes_weights(case_cnvbt, case_cnv_weights, control_cnvbt, control_cnv_weights)
    
    # # Format output main counts table and write to outfile
    # make_output_table(outbed, txbt, genes, cds_dict, control_counts, 
    #                   control_weights, case_counts, case_weights)
    # if bgzip:
    #     subprocess.run(['bgzip', '-f', outbed_path])


if __name__ == '__main__':
    main()
