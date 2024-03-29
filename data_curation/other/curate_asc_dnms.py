#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Curates de novo mutations from ASC supplemental tables
Expects the equivalent of Supplementary Table 2.3 from Satterstrom et al., Cell, 2020
Can handle the supplement from Fo et al., medRxiv, 2021, with --version fu_2021
"""

import pandas as pd
import numpy as np
import argparse
from os.path import splitext
import subprocess


lof_csqs = 'frameshift_variant stop_gained splice_donor_variant splice_acceptor_variant'.split()


def load_dnms(dnms_in, controls=False, version='satterstrom_2020'):
    """
    Load DNMs and extract columns & rows of interest
    Return: list of (gene, consequence) tuples
    """

    # Load DNMs as pd.DataFrame
    dnm_df = pd.read_csv(dnms_in, sep='\t')

    # Set code for affectedness status
    if controls:
        affected_key = 1
    else:
        affected_key = 2

    # Parse table based on version of ASC analysis specified
    if version == 'satterstrom_2020':

        # Restrict to high-confidence coding variants
        keep_rows = (dnm_df.Affected_Status == affected_key) \
                    & dnm_df.Coding \
                    & (dnm_df.Confidence == 'HIGH')
        dnm_df = dnm_df.loc[keep_rows, :]

        # Restrict to columns of interest
        keep_cols = 'GENE_NAME VEP_functional_class_canonical_simplified'.split()

    elif version == 'fu_2021':

        # Restrict to high-confidence coding variants
        dnm_df = dnm_df[dnm_df.Affected_Status == affected_key]

        # Restrict to columns of interest
        keep_cols = 'Gene Simplified_csq'.split()

    # Return data as list of (gene, consequence) tuples
    return list(dnm_df.loc[:, keep_cols].itertuples(index=False, name=None))


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dnm-tsv', help='tsv of de novo mutations from Satterstrom ' +
                        'et al., Cell, 2020. Required.', required=True)
    parser.add_argument('--genes', help='list of gene symbols to consider. Required.',
                        required=True)
    parser.add_argument('--controls', help='output DNMs in controls. [default: cases]',
                        action='store_true')
    parser.add_argument('-v', '--version', default='satterstrom_2020', help='Specify ' +
                        'version of ASC analysis to expect. [default: satterstrom_2020]',
                        choices=['satterstrom_2020', 'fu_2021'])
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
    zeros = {'lof' : 0, 'mis' : 0, 'syn' : 0}
    gdict = {g.rstrip() : zeros.copy() for g in open(args.genes).readlines()}

    # Load & filter DNMs
    dnms = load_dnms(args.dnm_tsv, args.controls, args.version)

    # Count DNMs per gene
    for gene, csq in dnms:

        if gene not in gdict.keys():
            continue

        if csq in lof_csqs:
            gdict[gene]['lof'] += 1
        elif 'missense' in csq:
            gdict[gene]['mis'] += 1
        elif 'synonymous' in csq:
            gdict[gene]['syn'] += 1

    # Write to outfile
    header = '\t'.join('#gene lof mis syn'.split())
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

