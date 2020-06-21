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


def count_cnvs_per_crb(cnvbt, crb_path, element_path, min_element_ovr=1e-10,
                       min_frac_all_elements=1e-10):
    """
    Count CNVs per CRB based on intersection with CRB elements
    """

    cnv_ids = [x.name for x in cnvbt]

    crbs = {x.name : {'cnvs' : [], 'n_ele' : int(x[-1])} for x in pbt.BedTool(crb_path)}

    # Count number of elements hit per CRB per CNV
    hit_dict = {cnv_id : {} for cnv_id in cnv_ids}
    element_hits = cnvbt.intersect(element_path, wb=True, F=min_element_ovr)
    for cnv_id, crb_id in [(x[3], x[-1]) for x in element_hits]:
        if crb_id not in hit_dict[cnv_id].keys():
            hit_dict[cnv_id][crb_id] = 1
        else:
            hit_dict[cnv_id][crb_id] += 1

    # Count CNVs per CRB only if they hit >= (min_element_ovr * n_elements) total elements
    for cnv_id, ecounts in hit_dict.items():
        for crb_id, n_hit in ecounts.items():
            n_total = crbs[crb_id]['n_ele']
            if n_hit / n_total >= min_frac_all_elements:
                crbs[crb_id]['cnvs'].append(cnv_id)

    for crb_id in crbs.keys():
        crbs[crb_id].pop('n_ele')

    return crbs


def make_output_table(outbed, cnvbt, crb_counts, crb_path, case_hpo, control_hpo):
    """
    Format master table of counts per CRB and write to outbed
    """

    cnv_hpos = {x.name : x[-1].split(';') for x in cnvbt}

    h_cols = '#chr start end crb_id control_cnvs case_cnvs'
    header = '\t'.join(h_cols.split())
    outbed.write(header + '\n')

    for i in pbt.BedTool(crb_path):
        crb_id = i.name
        crb_cnvs = crb_counts[crb_id]['cnvs']
        k_control = len([c for c in crb_cnvs if control_hpo in cnv_hpos.get(c)])
        k_case = len([c for c in crb_cnvs if case_hpo in cnv_hpos.get(c)])
        crbline = [i.chrom, i.start, i.end, i.name, k_control, k_case]
        crbline_str = '\t'.join([str(x) for x in crbline]) + '\n'
        outbed.write(crbline_str)

    # Must close output file to flush buffer (end of chr22 sometimes gets clipped)
    outbed.close()


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

    # Intersect CNVs with CRBs
    crb_counts = \
        count_cnvs_per_crb(cnvbt, args.crbs, args.elements, args.min_element_ovr, 
                               args.min_frac_all_elements)

    # Format output main counts table and write to outfile
    make_output_table(outbed, cnvbt, crb_counts, args.crbs, args.hpo, args.control_hpo)
    if bgzip:
        subprocess.run(['bgzip', '-f', outbed_path])


if __name__ == '__main__':
    main()
