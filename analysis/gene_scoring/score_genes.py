#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Fit model & score dosage sensitivity for all genes based on BFDP
"""


from os import path
import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod import families
from scipy.stats import norm
# from statsmodels.stats.multitest import multipletests as padjust
import argparse
from sys import stdout


def load_stats(stats_in):
    """
    Load & format stats per gene
    """

    # Read data
    keep_cols = 'gene chrom bfdp'.split()
    ss = pd.read_csv(stats_in, delimiter='\t').loc[:, keep_cols]
    ss.set_axis(ss.gene, axis=0, inplace=True)
    ss.drop(labels='gene', axis=1, inplace=True)

    # Convert bfdp to numeric
    ss['bfdp'] = pd.to_numeric(ss['bfdp'], errors='coerce')

    return ss


def pair_chroms(sumstats, chroms, quiet=False):
    """
    Divide all chromosomes into pairs based on total number of genes
    """

    # Count genes per chromosome and sort by number of genes
    gpc = {chrom : sum(sumstats.chrom == chrom) for chrom in chroms}
    gpc = {c : g for c, g in sorted(gpc.items(), key=lambda item: item[1], reverse=True)}
    sortchroms = list(gpc.keys())

    chrompairs = []
    for i in range(0, int(len(chroms) / 2) ):
        topchrom = sortchroms[i]
        bottomchrom = sortchroms[-(i+1)]
        chrompairs.append((topchrom, bottomchrom))

    if not quiet:
        print('\nChromosome pairings:\n\nchrA\tchrB\tgenes')
        for pair in chrompairs:
            ngenes = gpc[pair[0]] + gpc[pair[1]]
            print('{}\t{}\t{:,}'.format(pair[0], pair[1], ngenes))

    return chrompairs


def load_features(features_in, sumstats):
    """
    Load gene features, standardize, and subset to genes present in any gene block
    """

    # Read data & drop coordinates
    features = pd.read_csv(features_in, delimiter='\t')
    for colname in '#chr chr start end'.split():
        if colname in features.columns:
            features.drop(labels=colname, axis=1, inplace=True)

    # Fill missing values with column-wise means and normalize
    for col in features.columns.tolist()[1:]:
        missing = (features[col] == '.')
        if missing.sum() > 0:
            fmean = np.nanmean(features.loc[~missing, col].astype(float))
            features.loc[missing, col] = fmean
    features.iloc[:, 1:] = scale(features.iloc[:, 1:])

    # Move gene name to row index
    features.set_axis(features.gene, axis=0, inplace=True)
    features.drop(labels='gene', axis=1, inplace=True)

    return features


def rmse(pairs):
    """
    Compute root mean-squared error for a list of numeric tuples
    """

    mse = [(a - b) ** 2 for a, b in pairs]
    
    return np.sqrt(np.sum(mse) / len(pairs))


def fit_model(features, sumstats, train_pairs, test_pairs, logit_alpha, l1_l2_mix):
    """
    Fit glm to genes on train_chroms and calculate RMSE on test_chroms
    Note: train_pairs and test_pairs are tuples of chromosome names
    """

    # Format chromosomes
    if type(train_pairs) is list:
        trainchroms = sorted(list(sum(train_pairs, ())))
    else:
        trainchroms = sorted(list(train_pairs))
    if type(test_pairs) is list:
        testchroms = sorted(list(sum(test_pairs, ())))
    else:
        testchroms = sorted(list(test_pairs))
    allchroms = sorted(trainchroms + testchroms)

    # Join sumstats with features for logistic regression, subset to 
    # chromosomes of interest, and drop genes with NaN bfdps
    logit_df = sumstats.merge(features, how='left', left_index=True, 
                              right_index=True)
    logit_df = logit_df.loc[logit_df.chrom.isin(allchroms), :].dropna()
    train_df = logit_df.loc[logit_df.chrom.isin(trainchroms), :].\
                   drop(labels='chrom', axis=1)
    test_df = logit_df.loc[logit_df.chrom.isin(testchroms), :].\
                  drop(labels='chrom', axis=1)

    # Fit logit GLM & predict on test set
    glm = GLM(train_df.bfdp, train_df.drop(labels='bfdp', axis=1),
              family=families.Binomial())
    if logit_alpha is None:
        logit = glm.fit()
    else:
        logit = glm.fit_regularized(L1_wt = 1 - l1_l2_mix, alpha=logit_alpha)
    test_bfdps = logit.predict(test_df.drop(labels='bfdp', axis=1))
    test_bfdps.name = 'pred'

    # Compute RMSE of bfdps for test set
    test_vals = test_df.merge(test_bfdps, left_index=True, right_index=True).\
                    loc[:, 'bfdp pred'.split()]
    test_rmse = rmse(test_vals.to_records(index=False))

    return logit, test_rmse


def fit_model_cv(features, sumstats, all_pairs, logit_alpha, l1_l2_mix):
    """
    Fit glm with N-fold CV 
    """

    models = []
    for i in range(len(all_pairs)):
        train_pairs = [x for k, x in enumerate(all_pairs) if k != i]
        test_pairs = all_pairs[i]
        model_fit, test_rmse = \
            fit_model(features, sumstats, train_pairs, test_pairs, logit_alpha, 
                      l1_l2_mix)
        models.append((test_rmse, model_fit))

    best_rmse, best_model = \
        [(r, m) for r, m in sorted(models, key=lambda x: x[0])][0]

    return best_model, best_rmse


def predict_bfdps(features, sumstats, chrompairs, logit_alpha, l1_l2_mix, quiet=False):
    """
    Predicts bfdps for all genes with N-fold CV
    """

    pred_bfdps = pd.DataFrame(columns=['pred_bfdp'], index=features.index)

    # Generate predictions for each pair of chromosomes in serial
    for i in range(len(chrompairs)):
        train_pairs = [x for k, x in enumerate(chrompairs) if k != i]
        pred_chroms = sorted(list(chrompairs[i]))
        pred_genes = list(sumstats.index[sumstats.chrom.isin(pred_chroms)])

        if not quiet:
            msg = 'Now scoring {:,} genes from chromsomes {} & {}'
            print(msg.format(len(pred_genes), pred_chroms[0], pred_chroms[1]))
        
        # Fit model
        model, rmse = fit_model_cv(features, sumstats, train_pairs, logit_alpha, 
                                   l1_l2_mix)

        # Predict bfdps
        pred_df = features.loc[features.index.isin(pred_genes), :]
        pred_vals = model.predict(pred_df)
        pred_vals.name = 'pred_bfdp'
        pred_bfdps.update(pred_vals)

    # Remove genes with no predictions (due to being excluded from feature collection)
    pred_bfdps = pred_bfdps.dropna()

    # Compute final gene score as pnorm from Z-score of pred_bfdps
    pred_bfdps['score'] = norm.sf(scale(pred_bfdps.pred_bfdp))
    pred_bfdps['quantile'] = 100 * pred_bfdps.pred_bfdp.transform('rank') / len(pred_bfdps)

    # # Compute Q-value as B-H FDR adjusted score
    # pred_bfdps['fdr'] = padjust(1 - pred_bfdps['score'], method='fdr_bh')[1]
    
    return pred_bfdps


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('stats', help='.tsv of stats with bfdp per gene.')
    parser.add_argument('features', help='.bed.gz of functional features per gene.')
    parser.add_argument('--regularization-alpha', dest='logit_alpha', type=float,
                        help='Regularization penalty weight for logit glm. Must ' +
                        'be in ~ [0, 1]. [default: no regularization]', default=0.1)
    parser.add_argument('--regularization-l1-l2-mix', dest='l1_l2_mix', type=float,
                        help='Regularization parameter (elastic net alpha) for ' +
                        'logit glm. 0 = L1, 1 = L2, (0, 1) = elastic net. ' +
                        '[default: L2 regularization]', default=1)
    parser.add_argument('-o', '--outfile', default='stdout', help='Output .tsv of ' +
                        'scores per gene. [default: stdout]')
    args = parser.parse_args()

    # Open connections to output files
    if args.outfile in 'stdout - /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Import gene stats
    sumstats = load_stats(args.stats)
    chroms = sorted(np.unique(sumstats.chrom.values))

    # Pair all chromosomes
    chrompairs = pair_chroms(sumstats, chroms)

    # Read gene features
    features = load_features(args.features, sumstats)

    # Predict bfdps for all genes
    pred_bfdps = predict_bfdps(features, sumstats, chrompairs, args.logit_alpha, 
                               args.l1_l2_mix)

    # DEV: write pred_bfdps out to file to decide how to handle quantile norm
    pred_bfdps.to_csv(outfile, sep='\t', index=True, na_rep='NA')


if __name__ == '__main__':
    main()


