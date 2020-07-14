#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Gather lists of genes with evidence of expression variability in GTEx
"""


import pandas as pd
import numpy as np
from os import path
import gzip
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


def process_gtex(gtex_in, genes, ensg_ids, tissues, min_mean_tpm=1, min_sample_tpm=1, 
                 min_tissues=1, min_variable_cv=0.5, max_invariant_cv=0.3, 
                 min_variable_prop_low=0.005, max_invariant_prop_low=0.001, 
                 min_variable_prop_high=0.05, max_invariant_prop_high=0.01, 
                 quiet=False):
    """
    Identifies variable expressors from GTEx master sample X gene quantification matrix
    """

    # Instantiate collectors
    low_variable, low_invariant, high_variable, high_invariant = set(), set(), set(), set()

    # # DEBUG
    # dbgout = open('mean_props.test.tsv', 'w')

    # Note: for memory purposes, each line is processed one at a time
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
        gene = vals[1]

        # Store indexes corresponding to sample IDs
        if vals[0] == 'Name':
            for i, sid in enumerate(vals):
                if i > 1:
                    col_idxs[sid] = i

        # Process each gene if canonical
        else:
            if not any([vals[0].split('.')[0] in ensg_ids, vals[1] in genes]):
                continue

            # DEBUG: print line number
            if not quiet:
                r += 1
                print('\t'.join([str(r), str(k), vals[1]]))

            # Compute number of expression outliers per tissue
            prop_lows, prop_highs, cvs, tmins, tmaxs = [], [], [], [], []
            n_expressed_tissues = 0
            for tissue, samples in tissues.items():

                # Get expression values for tissue
                tidxs = [col_idxs[s] for s in samples if s in col_idxs.keys()]
                tvals = [float(vals[i]) for i in tidxs]
                tmean = np.nanmean(tvals)
                if len(tvals) > 0:
                    tq1 = np.quantile(tvals, 0.25)
                else:
                    tq1 = 0

                # Only process gene & tissue pairs with average expression > min_mean_tpm
                # and first quartile of expression > 0
                if tmean < min_mean_tpm and tq1 > 0:
                    continue
                else:
                    n_expressed_tissues += 1

                # Compute expression statistics per tissue:
                # 1. Proportion of low outlier samples per tissue
                # 2. Proportion of high outlier samples per tissue
                # 3. Min expression value per tissue
                # 4. Max expression value per tissue
                # Outlier defined as 1.5IQR above/below Q1 or Q3
                nsamp = len(tvals)
                if nsamp > 0:

                    tmins.append(np.nanmin(tvals))
                    tmaxs.append(np.nanmax(tvals))
                    
                    cv = np.nanstd(tvals) / tmean
                    cvs.append(cv)
                    
                    tq13 = np.quantile(tvals, [0.25, 0.75])
                    tiqr = tq13[1] - tq13[0]
                    tlow_cut = tq13[0] - (1.5 * tiqr)
                    thigh_cut = tq13[1] + (1.5 * tiqr)
                    
                    nlow = len([v for v in tvals if v <= tlow_cut])
                    prop_low = nlow / nsamp
                    prop_lows.append(prop_low)
                    
                    nhigh = len([v for v in tvals if v >= thigh_cut])
                    prop_high = nhigh / nsamp
                    prop_highs.append(prop_high)

            # Add gene to lists depending on expression properties
            if n_expressed_tissues >= min_tissues:
                avg_prop_low = np.nanmean(prop_lows)
                avg_prop_high = np.nanmean(prop_highs)
                avg_cv = np.nanmean(cvs)
                global_tmin = np.nanmin(tmins)

                # Low variable vs invariant
                if avg_cv >= min_variable_cv \
                and avg_prop_low >= min_variable_prop_low:
                    low_variable.add(gene)
                elif global_tmin >= min_sample_tpm \
                and avg_cv <= max_invariant_cv \
                and avg_prop_low <= max_invariant_prop_low:
                    low_invariant.add(gene)

                # High variable vs invariant
                if avg_cv >= min_variable_cv \
                and avg_prop_low >= min_variable_prop_high:
                    high_variable.add(gene)
                elif avg_cv <= max_invariant_cv \
                and avg_prop_high <= max_invariant_prop_high:
                    high_invariant.add(gene)

                # dbgdat = [gene, n_expressed_tissues, avg_prop_low, avg_prop_high, avg_cv, global_tmin]
                # dbgout.write('\t'.join([str(x) for x in dbgdat]) + '\n')

    # dbgout.close()

    return low_variable, low_invariant, high_variable, high_invariant


def write_genelist(genes, prefix, suffix):
    """
    Write gene list to file
    """

    outname = '{}.{}.genes.list'.format(prefix, suffix)

    fout = open(outname, 'w')

    for gene in genes:
        fout.write(gene + '\n')

    fout.close()



def main():
    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('canonical', help='tsv of canonical transcripts')
    parser.add_argument('gtex_matrix', help='GTEx sample X gene expression matrix')
    parser.add_argument('gtex_sample_manifest', help='GTEx sample manifest tsv')
    parser.add_argument('--min-mean-tpm', default=1, type=float, help='Minimum ' +
                        'mean TPM to consider a gene as expressed.')
    parser.add_argument('--min-sample-tpm', default=1, type=float, help='Minimum ' +
                        'TPM across all samples to consider a gene as expressed. ' +
                        'Note: applied only for the definition of low invariant genes.')
    parser.add_argument('--min-tissues', default=1, type=int, help='Minimum ' +
                        'number of tissues with expression > --min-tpm before ' +
                        'being considered.')
    parser.add_argument('--variable-tpm-cv', default=0.5, type=float, help='Minimum ' +
                        'coefficient of variation in expression levels among ' +
                        'samples to consider a gene variable.')
    parser.add_argument('--invariant-tpm-cv', default=0.3, type=float, help='Maximum ' +
                        'coefficient of variation in expression levels among ' +
                        'samples to consider a gene invariant.')
    parser.add_argument('--min-variable-prop-low', default=0.005, type=float, 
                        help='Minimum average proportion of low expression ' + 
                        'outliers per tissue to consider a gene variable.')
    parser.add_argument('--max-invariant-prop-low', default=0.001, type=float, 
                        help='Maximum average proportion of low expression ' + 
                        'outliers per tissue to consider a gene invariant.')
    parser.add_argument('--min-variable-prop-high', default=0.05, type=float, 
                        help='Minimum average proportion of high expression ' + 
                        'outliers per tissue to consider a gene variable.')
    parser.add_argument('--max-invariant-prop-high', default=0.01, type=float, 
                        help='Maximum average proportion of high expression ' + 
                        'outliers per tissue to consider a gene invariant.')
    parser.add_argument('-p', '--prefix', default='GTEx_summary', help='Prefix ' +
                        'for output gene lists')
    args = parser.parse_args()

    # Load list of genes & transcripts to allow
    genes, ensg_ids = get_genelist(args.canonical)

    # Load GTEx sample manifest
    tissues = load_gtex_samples(args.gtex_sample_manifest)

    # Process expression data for all genes
    low_variable, low_invariant, high_variable, high_invariant = \
        process_gtex(args.gtex_matrix, genes, ensg_ids, tissues, args.min_mean_tpm, 
                     args.min_sample_tpm, args.min_tissues, 
                     args.variable_tpm_cv, args.invariant_tpm_cv, 
                     args.min_variable_prop_low, args.max_invariant_prop_low,
                     args.min_variable_prop_high, args.max_invariant_prop_high)


    # Write gene lists to files
    write_genelist(list(low_variable), args.prefix, 'low_expression_variable')
    write_genelist(list(low_invariant), args.prefix, 'low_expression_invariant')
    write_genelist(list(high_variable), args.prefix, 'high_expression_variable')
    write_genelist(list(high_invariant), args.prefix, 'high_expression_invariant')


if __name__ == '__main__':
    main()

