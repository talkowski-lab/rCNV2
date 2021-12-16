#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Preprocess protein features downloaded from UniProt/SwissProt
"""


import argparse
import pandas as pd


def subset_by_genes(pdat, elig_genes):
    """
    Subset pd.DataFrame based on elig_genes
    """

    # First, find direct matches between primary gene symbol and our gene list
    primary_hits = pdat.gene.isin(elig_genes)

    # Next, try to match into alternative gene symbols
    noprim = pdat[~primary_hits & ~pdat['Gene names'].isna()]
    def _match_secondaries(gstr, elig_genes):
        match = pd.NA
        for gene in gstr.split():
            if gene in elig_genes:
                match = gene
                break
        return match
    noprim.gene = noprim['Gene names'].apply(_match_secondaries, elig_genes=elig_genes)
    rescued = noprim[~noprim.gene.isna()]

    return pd.concat([pdat[primary_hits], rescued], axis=0)


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tsv', help='.tsv downloaded from UniProt')
    parser.add_argument('eligible_genes', help='list of eligible gene symbols')
    parser.add_argument('outfile', help='path to output .tsv')
    args = parser.parse_args()

    # Read UniProt data as pd.DataFrame
    pdat = pd.read_csv(args.tsv, delimiter='\t')
    pdat.rename(columns={'Gene names  (primary )' : 'gene'}, inplace=True)

    # Subset to eligible genes
    elig_genes = [g.rstrip() for g in open(args.eligible_genes).readlines()]
    pdat = subset_by_genes(pdat, elig_genes)

    # Cast numeric column types
    pdat.Length = pdat.Length.astype(int)
    pdat.Mass = pdat.Mass.apply(lambda x: int(x.replace(',', '')))

    # Gather number of interactors
    def _count_interactors(istr):
        if pd.isna(istr):
            return 0
        else:
            return len(istr.split('; '))
    pdat['ppi_degree'] = pdat['Interacts with'].apply(_count_interactors)
    
    # Generate counts of protein features
    pdat['intramembrane_domains'] = \
        pdat['Intramembrane'].apply(lambda x: str(x).count('INTRAMEM'))
    pdat['transmembrane_domains'] = \
        pdat['Transmembrane'].apply(lambda x: str(x).count('TRANSMEM'))
    pdat['alpha_helixes'] = \
        pdat['Helix'].apply(lambda x: str(x).count('HELIX'))
    pdat['beta_sheets'] = \
        pdat['Beta strand'].apply(lambda x: str(x).count('STRAND'))
    pdat['beta_sheets'] = \
        pdat['Beta strand'].apply(lambda x: str(x).count('STRAND'))
    pdat['metal_binding_sites'] = \
        pdat['Metal binding'].apply(lambda x: str(x).count('METAL '))
    pdat['nucleotide_binding_sites'] = \
        pdat['Nucleotide binding'].apply(lambda x: str(x).count('NP_BIND'))
    pdat['active_sites'] = \
        pdat['Active site'].apply(lambda x: str(x).count('ACT_SITE'))

    # Rename selected columns
    pdat.rename(columns={'Length' : 'aa_length', 'Mass' : 'protein_mass'}, inplace=True)

    # Deduplicate redundant gene symbols
    pdat = pdat[~pdat.gene.duplicated()]

    # Subset to columns of interest and write to outfile
    cols_to_keep = 'gene aa_length protein_mass ppi_degree intramembrane_domains ' + \
                   'transmembrane_domains alpha_helixes beta_sheets ' + \
                   'metal_binding_sites nucleotide_binding_sites active_sites'
    out_df = pdat[cols_to_keep.split()]
    out_df.to_csv(args.outfile, sep='\t', na_rep='NA', index=False)


if __name__ == '__main__':
    main()

