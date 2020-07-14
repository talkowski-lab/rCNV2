#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Preprocess GTEx gene X sample matrix and distill various summary matrices
"""


import pandas as pd
import numpy as np
from os import path
import gzip
from scipy.stats import median_absolute_deviation as mad
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import subprocess
import argparse


def get_genelist(infile):
    """
    Extract list of genes and transcripts to consider
    """

    df = pd.read_csv(infile, delimiter='\t', header=None)

    ensg_ids = [x.split('.')[0] for x in df.iloc[:, 0].tolist()]
    # enst_ids = df.iloc[:, 1].tolist()
    genes = df.iloc[:, 4].tolist()

    return genes, ensg_ids


def load_gtex_samples(manifest_in):
    """
    Load GTEx sample manifest and return dict of samples per tissue
    """

    tissues = {}

    df = pd.read_csv(manifest_in, delimiter='\t')

    for tissue in sorted(list(set(df.SMTSD.tolist()))):
        samples = sorted(list(set(df.SAMPID[df.SMTSD == tissue].tolist())))
        tissues[tissue] = samples
    
    return tissues


def scale_expression(values):
    """
    Apply log10(x+1) scaling to a list of expression values
    """

    values = [float(x) for x in values]
    nonzero = [x for x in values if x > 0]
    if len(nonzero) > 0:
        # noise = np.min(nonzero) / 10
        return np.log10([x + 1 for x in values])
    else:
        return np.array([0 for x in values])


def get_gene_stats(xvals, col_idxs, tissues):
    """
    Compute summary stats across all samples for a given gene & tissue
    """

    xmin, xq1, xmed, xmean, xq3, xmax, xsd, xmad = [], [], [], [], [], [], [], []

    for tissue in tissues.keys():
        tidx = [col_idxs[s] for s in tissues[tissue] if s in col_idxs.keys()]
        if len(tidx) > 0:
            tvals = xvals[tidx]
            xmin.append(np.nanmin(tvals))
            xq1.append(np.nanquantile(tvals, q=0.25))
            xmed.append(np.nanmedian(tvals))
            xmean.append(np.nanmean(tvals))
            xq3.append(np.nanquantile(tvals, q=0.75))
            xmax.append(np.nanmax(tvals))
            xsd.append(np.nanstd(tvals))
            xmad.append(mad(tvals))
        else:
            xmin.append(np.nan)
            xq1.append(np.nan)
            xmed.append(np.nan)
            xmean.append(np.nan)
            xq3.append(np.nan)
            xmax.append(np.nan)
            xsd.append(np.nan)
            xmad.append(np.nan)
    
    return xmin, xq1, xmed, xmean, xq3, xmax, xsd, xmad


def load_gtex_matrix(gtex_in, genes, ensg_ids, tissues, quiet=False):
    """
    Reads GTEx expression matrix and filters to passing genes
    """

    # Note: for memory purposes, each line is processed one at a time

    df_cols = ['gene'] + list(tissues.keys())
    gtex_min = pd.DataFrame(columns=df_cols)
    gtex_q1 = pd.DataFrame(columns=df_cols)
    gtex_median = pd.DataFrame(columns=df_cols)
    gtex_mean = pd.DataFrame(columns=df_cols)
    gtex_q3 = pd.DataFrame(columns=df_cols)
    gtex_max = pd.DataFrame(columns=df_cols)
    gtex_sd = pd.DataFrame(columns=df_cols)
    gtex_mad = pd.DataFrame(columns=df_cols)

    if path.splitext(gtex_in)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
        infile = gzip.open(gtex_in, 'rt')
    else:
        infile = open(gtex_in)

    col_idxs = {}

    r = 0
    for k, line in enumerate(infile):
        # Skip descriptive header lines
        if len(line.split('\t')) < 3:
            continue

        vals = line.rstrip().split('\t')

        # Store indexes corresponding to sample IDs
        if vals[0] == 'Name':
            for i, sid in enumerate(vals):
                if i > 1:
                    col_idxs[sid] = i - 2

        # Process each gene if canonical
        else:
            if not any([vals[0].split('.')[0] in ensg_ids, vals[1] in genes]):
                continue

            # DEBUG: print line number
            if not quiet:
                r += 1
                print('\t'.join([str(r), str(k), vals[1]]))

            # Compute expression stats
            xvals = scale_expression(vals[2:])
            xmin, xq1, xmedian, xmean, xq3, xmax, xsd, xmad = \
                get_gene_stats(xvals, col_idxs, tissues)

            # Add gene stats to each summary df
            gtex_min = gtex_min.append(pd.Series([vals[1]] + xmin, index=gtex_min.columns),
                                       ignore_index=True)
            gtex_q1 = gtex_q1.append(pd.Series([vals[1]] + xq1, index=gtex_q1.columns),
                                     ignore_index=True)
            gtex_median = gtex_median.append(pd.Series([vals[1]] + xmedian, 
                                                       index=gtex_median.columns),
                                             ignore_index=True)
            gtex_mean = gtex_mean.append(pd.Series([vals[1]] + xmean, 
                                                   index=gtex_mean.columns),
                                         ignore_index=True)
            gtex_q3 = gtex_q3.append(pd.Series([vals[1]] + xq3, 
                                                index=gtex_q3.columns),
                                     ignore_index=True)
            gtex_max = gtex_max.append(pd.Series([vals[1]] + xmax, 
                                                 index=gtex_max.columns),
                                       ignore_index=True)
            gtex_sd = gtex_sd.append(pd.Series([vals[1]] + xsd, 
                                                 index=gtex_sd.columns),
                                       ignore_index=True)
            gtex_mad = gtex_mad.append(pd.Series([vals[1]] + xmad, 
                                                 index=gtex_mad.columns),
                                       ignore_index=True)

    return gtex_min, gtex_q1, gtex_median, gtex_mean, gtex_q3, gtex_max, gtex_sd, gtex_mad


def pca_expression(matrix, n_pcs=20):
    """
    Reduce an expression matrix to its n_pcs top principal components
    """

    # Clean input matrix
    X = StandardScaler().fit_transform(matrix.drop(columns='gene').fillna(value=0))

    # PCA
    pca = PCA(n_components=n_pcs).fit_transform(X)
    pc_names = ['expression_component_' + str(i) for i in range(1, n_pcs+1)]
    pcadf = pd.DataFrame(pca, columns=pc_names)

    # Convert back to dataframe with genes as first column
    return pd.concat([matrix.loc[:, 'gene'], pcadf], axis=1)


def write_matrix(matrix, filename, gzip):
    """
    Write matrix to file and gzip if optioned
    """

    matrix.rename(columns={'gene' : '#gene'}, inplace=True)
    
    matrix.to_csv(filename, sep='\t', index=False, na_rep='NA')

    if gzip:
        subprocess.run(['gzip', '-f', filename])


def main():
    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('canonical', help='tsv of canonical transcripts')
    parser.add_argument('gtex_matrix', help='GTEx sample X gene expression matrix')
    parser.add_argument('gtex_sample_manifest', help='GTEx sample manifest tsv')
    parser.add_argument('-p', '--prefix', default='GTEx_summary', help='Prefix ' +
                        'for output matrices')
    parser.add_argument('--n-pcs', default=20, help='Number of principal components ' +
                        'to retain', type=int)
    parser.add_argument('-z', '--gzip', action='store_true', help='Gzip outputs')
    args = parser.parse_args()

    # Load list of genes & transcripts to allow
    genes, ensg_ids = get_genelist(args.canonical)

    # Load GTEx sample manifest
    tissues = load_gtex_samples(args.gtex_sample_manifest)

    # Load sample X gene expression matrix
    gtex_min, gtex_q1, gtex_median, gtex_mean, gtex_q3, gtex_max, gtex_sd, gtex_mad = \
        load_gtex_matrix(args.gtex_matrix, genes, ensg_ids, tissues)

    # PCA of gene-level average expression levels across tissues
    gtex_pca = pca_expression(gtex_mean, args.n_pcs)

    # Write matrices to output files
    write_matrix(gtex_min, args.prefix + '.min.tsv', args.gzip)
    write_matrix(gtex_q1, args.prefix + '.q1.tsv', args.gzip)
    write_matrix(gtex_median, args.prefix + '.median.tsv', args.gzip)
    write_matrix(gtex_mean, args.prefix + '.mean.tsv', args.gzip)
    write_matrix(gtex_q3, args.prefix + '.q3.tsv', args.gzip)
    write_matrix(gtex_max, args.prefix + '.max.tsv', args.gzip)
    write_matrix(gtex_sd, args.prefix + '.sd.tsv', args.gzip)
    write_matrix(gtex_mad, args.prefix + '.mad.tsv', args.gzip)
    write_matrix(gtex_pca, args.prefix + '.pca.tsv', args.gzip)


if __name__ == '__main__':
    main()
