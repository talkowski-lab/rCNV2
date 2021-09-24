#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Curates counts of functional SVs from gnomAD-SV VCF
"""

import pandas as pd
import numpy as np
import argparse
from os.path import splitext
import pysam
import subprocess


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gnomad-sv-vcf', help='gnomad-SV VCF. Required.', required=True)
    parser.add_argument('--genes', help='list of gene symbols to consider. Required.',
                        required=True)
    parser.add_argument('-o', '--outfile', help='Path to output tsv file. [default: ' +
                        'stdout]', default='stdout')
    parser.add_argument('-z', '--gzip', action='store_true', help='Compress ' + 
                        'output tsv file with gzip. [Default: do not compress]')
    args = parser.parse_args()

    gz_suf = 'gz gzip'.split()

    # Open connection to outfile
    if args.outfile in '- stdout'.split():
        outfile = stdout
    else:
        if splitext(args.outfile)[-1].replace('.' ,'') in gz_suf:
            outfile_path = splitext(args.outfile)[0]
        else:
            outfile_path = args.outfile
        outfile = open(outfile_path, 'w')

    # Build dict of genes & counts
    zeros = {'lof_any' : 0, 'lof_del' : 0, 'lof_other' : 0, 'cg' : 0, 'ied' : 0}
    gdict = {g.rstrip() : zeros.copy() for g in open(args.genes).readlines()}

    # Parse VCF
    for sv in pysam.VariantFile(args.gnomad_sv_vcf).fetch():
        if 'PASS' not in sv.filter:
            continue 

        if 'PROTEIN_CODING__LOF' in sv.info.keys():
            for gene in sv.info.get('PROTEIN_CODING__LOF', []):
                if gene in gdict.keys():
                    gdict[gene]['lof_any'] += 1
                    if sv.info['SVTYPE'] == 'DEL':
                        gdict[gene]['lof_del'] += 1
                    else:
                        gdict[gene]['lof_other'] += 1

        if 'PROTEIN_CODING__COPY_GAIN' in sv.info.keys():
            for gene in sv.info.get('PROTEIN_CODING__COPY_GAIN', []):
                if gene in gdict.keys():
                    gdict[gene]['cg'] += 1

        if 'PROTEIN_CODING__DUP_LOF' in sv.info.keys():
            for gene in sv.info.get('PROTEIN_CODING__DUP_LOF', []):
                if gene in gdict.keys():
                    gdict[gene]['ied'] += 1

    # Write to outfile
    header = '\t'.join('#gene lof_any lof_del lof_other cg ied'.split())
    outfile.write(header + '\n')
    for gene, gdat in gdict.items():
         gline = gene + '\t' + '\t'.join([str(x) for x in gdat.values()])
         outfile.write(gline + '\n')
    outfile.close()

    # Gzip output, if optioned
    if args.gzip:
        subprocess.run(['gzip', '-f', outfile_path])


if __name__ == '__main__':
    main()

