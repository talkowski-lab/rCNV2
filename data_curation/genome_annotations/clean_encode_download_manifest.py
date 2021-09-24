#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Filters and cleans an ENCODE experiment metadata matrix
More info: https://www.encodeproject.org/help/batch-download/
"""


import pandas as pd
import argparse


def load_manifest(manifest_url):
    """
    Download, import, and clean manifest
    """

    # Read manifest direct from URL
    df = pd.read_csv(manifest_url, sep='\t')

    # Restrict to:
    #   1. hg19-aligned
    #   2. BED-like
    #   3. released files only
    #   4. ENCODE data audit compliant
    #   5. No ENCODE data audit errors
    keep_rows = \
        (df['File assembly'] == 'hg19') & \
        (df['File type'].str.contains('bed', case=False)) & \
        (df['File Status'] == 'released') & \
        (df['Audit NOT_COMPLIANT'].astype(str) == 'nan') & \
        (df['Audit ERROR'].astype(str) == 'nan') & \
        (df['Biosample treatments'].astype(str) == 'nan') & \
        (~df['Output type'].str.contains('background', case=False))

    # Filter to informative columns
    keep_cols = \
        ['File accession', 'Output type', 'Experiment accession', 'Assay',
         'Biosample term name', 'Biosample type', 'Experiment target', 
         'Experiment date released', 'Project', 'Lab', 'File download URL']
    
    return df.loc[keep_rows, keep_cols]


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('files', help='files.txt from ENCODE batch download. ' +
                        'Assumes manifest url is first line.')
    parser.add_argument('--outfile', required=True, help='Output tsv for cleaned manifest. ' +
                        'Outfiles ending in .gz will be automatically compressed.')
    parser.add_argument('--tracklist', required=True, help='Output tracklist.')
    args = parser.parse_args()

    # Download & clean manifest
    with open(args.files) as fin:
        manifest = load_manifest(fin.readline().rstrip())

    # Write cleaned manifest to output file
    manifest.to_csv(args.outfile, sep='\t', na_rep='NA', header=True, index=False)

    # Write tracklist to output file
    pd.DataFrame(manifest['File download URL']).\
       to_csv(args.tracklist, header=False, index=False)


if __name__ == '__main__':
    main()

