#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Fit model & score dosage sensitivity for all genes based on BFDP
"""


model_options = 'logit svm randomforest lda naivebayes sgd neuralnet'.split()


from os import path
import pybedtools as pbt
import pandas as pd
import numpy as np
from sklearn.preprocessing import scale, StandardScaler
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod import families
# from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDAC
from sklearn.naive_bayes import GaussianNB as GNBC
from sklearn.linear_model import SGDClassifier as SGDC
from sklearn.neural_network import MLPClassifier as MLPC
from scipy.stats import norm
from random import seed, shuffle
import argparse
from sys import stdout


def get_chrom_arms(features_in, centromeres_in):
    """
    Assign genes to chromosome arms based on centromere coordinates
    Returns: dict of gene : chromosome arm
    """

    genes_bt = pbt.BedTool(features_in).cut(range(4)).saveas()
    chroms = list(set([x.chrom for x in genes_bt]))
    cen_bt = pbt.BedTool(centromeres_in).\
                 filter(lambda x: x[3] == 'centromere' and x.chrom in chroms).\
                 saveas()

    arm_dict = {}

    for x in genes_bt.closest(cen_bt, D='b'):
        gene = x[3]
        if int(x[-1]) >= 0:
            arm_dict[gene] = x.chrom + 'p'
        else:
            arm_dict[gene] = x.chrom + 'q'

    return arm_dict


def load_stats(stats_in, arm_dict=None):
    """
    Load & format stats per gene
    """

    # Read data
    keep_cols = 'gene chrom bfdp'.split()
    ss = pd.read_csv(stats_in, delimiter='\t').loc[:, keep_cols]
    ss.set_axis(ss.gene, axis=0, inplace=True)
    ss.drop(labels='gene', axis=1, inplace=True)

    # Rewrite gene chromosomes, if optioned
    if arm_dict is not None:
        ss['chrom'] = ss.index.map(arm_dict)

    # Convert bfdp to numeric
    ss['bfdp'] = pd.to_numeric(ss['bfdp'], errors='coerce')

    return ss


def pair_chroms(sumstats, chroms, xbed=None, n_pairs=11, quiet=False):
    """
    Divide all chromosomes into pairs based on total number of genes
    Or, if --centromeres is provided, groups chromosome arms into 11 groups
    Removes genes in xbed from calculations, if provided
    """

    # Load blacklist, if any provided
    if xbed is not None:
        xgenes = pd.read_csv(xbed, sep='\t').iloc[:, 3].values.tolist()
    else:
        xgenes = []

    # Count genes per chromosome and sort by number of genes
    gpc = {chrom : sum((sumstats.chrom == chrom) & ~(sumstats.index.isin(xgenes))) for chrom in chroms}
    xgpc = {chrom : sum((sumstats.chrom == chrom) & (sumstats.index.isin(xgenes))) for chrom in chroms}
    gpc = {c : g for c, g in sorted(gpc.items(), key=lambda item: item[1], reverse=True)}
    sortchroms = list(gpc.keys())

    chrompairs_dict = {x : {'chroms' : ( ), 'ngenes' : 0} for x in range(n_pairs)}
    for chrom, ngenes in gpc.items():
        target_pair = list(chrompairs_dict.keys())[0]
        chrompairs_dict[target_pair]['chroms'] += (chrom, )
        chrompairs_dict[target_pair]['ngenes'] += ngenes
        chrompairs_dict = {c : n for c, n in sorted(chrompairs_dict.items(), key=lambda x: x[1]['ngenes'])}

    chrompairs = [x['chroms'] for x in chrompairs_dict.values()]

    if not quiet:
        print('\nChromosome pairings:\n\ngenes\tbl_genes\tchroms')
        for pair in chrompairs:
            ngenes = np.nansum([gpc[c] for c in pair])
            nxgenes = np.nansum([xgpc[c] for c in pair])
            print('{:,}\t{:,}\t{}'.format(ngenes, nxgenes, ', '.join(list(pair))))

    return chrompairs


def load_features(features_in):
    """
    Load & standardize gene features
    """

    # Read data & drop coordinates
    features = pd.read_csv(features_in, delimiter='\t')
    for colname in '#chr chr start end'.split():
        if colname in features.columns:
            features.drop(labels=colname, axis=1, inplace=True)

    # Fill missing values with column-wise means
    for col in features.columns.tolist()[1:]:
        missing = (features[col] == '.')
        if missing.sum() > 0:
            fmean = np.nanmean(features.loc[~missing, col].astype(float))
            features.loc[missing, col] = fmean

    # Move gene name to row index
    features.set_axis(features.gene, axis=0, inplace=True)
    features.drop(labels='gene', axis=1, inplace=True)

    # Apply sklearn standard normalization to features
    scaler = StandardScaler().fit(features)
    scaled_features = pd.DataFrame(scaler.transform(features),
                                   columns=features.columns,
                                   index=features.index)

    return features


def rmse(pairs):
    """
    Compute root mean-squared error for a list of numeric tuples
    """

    mse = [(a - b) ** 2 for a, b in pairs]
    
    return np.sqrt(np.sum(mse) / len(pairs))


def fit_model(features, sumstats, train_genes, test_genes,  
              model='logit', logit_alpha=0.1, l1_l2_mix=1):
    """
    Fit classifier to train_genes and calculate RMSE on test_genes
    """

    all_genes = train_genes + test_genes

    # Join sumstats with features for logistic regression, subset to 
    # genes of interest, and drop genes with NaN BFDPs 
    full_df = sumstats.merge(features, how='left', left_index=True, 
                              right_index=True)
    full_df = full_df.loc[full_df.index.isin(all_genes), :].dropna()
    train_df = full_df.loc[full_df.index.isin(train_genes), :].\
                   drop(labels='chrom', axis=1)
    test_df = full_df.loc[full_df.index.isin(test_genes), :].\
                  drop(labels='chrom', axis=1)

    # Fit according to specified model
    if model == 'logit':
        # Fit logit GLM & predict on test set
        glm = GLM(train_df.bfdp, train_df.drop(labels='bfdp', axis=1),
                  family=families.Binomial())
        if logit_alpha is None:
            fitted_model = glm.fit()
        else:
            fitted_model = glm.fit_regularized(L1_wt = 1 - l1_l2_mix, alpha=logit_alpha)
        test_bfdps = fitted_model.predict(test_df.drop(labels='bfdp', axis=1))
        test_bfdps.name = 'pred'

    else:
        # Instantiate classifier, depending on model
        if model == 'svm':
            classifier = SVC(C=logit_alpha, random_state=0, max_iter=1000, 
                             probability=True, kernel='linear', break_ties=True)
        elif model == 'randomforest':
            classifier = RFC(random_state=0)
        elif model == 'lda':
            classifier = LDAC()
        elif model == 'naivebayes':
            classifier = GNBC()
        elif model == 'sgd':
            classifier = SGDC(loss='modified_huber', penalty='elasticnet', 
                              alpha=logit_alpha, l1_ratio=1 - l1_l2_mix, 
                              random_state=0)
        elif model == 'neuralnet':
            classifier = MLPC(hidden_layer_sizes=(10, 5, 2), activation='logistic',
                              solver='adam', alpha=logit_alpha, random_state=0)

        # Fit sklearn model & predict on test set
        fitted_model = classifier.fit(train_df.drop(labels='bfdp', axis=1), 
                                      np.round(train_df.bfdp))
        test_bfdps = pd.Series(fitted_model.predict_proba(test_df.drop(labels='bfdp', axis=1))[:, 1],
                               name='pred', index=test_df.index)

    # Compute RMSE of bfdps for test set
    test_vals = test_df.merge(test_bfdps, left_index=True, right_index=True).\
                    loc[:, 'bfdp pred'.split()]
    test_rmse = rmse(test_vals.to_records(index=False))

    return fitted_model, test_rmse


def split_genes(genes, n_splits=10, random_seed=2020):
    """
    Randomly partition a list of genes into n equal groups
    """

    target = int(np.ceil(len(genes) / float(n_splits)))

    seed(random_seed)
    shuffle(genes)

    return [genes[int(i * target):int((i + 1) * target)] for i in range(n_splits)]


def fit_model_cv(features, sumstats, all_pairs, xgenes=[], model='logit', 
                 logit_alpha=0.1, l1_l2_mix=1, random=True):
    """
    Fit model with N-fold CV 
    
    """

    models = []

    # Randomly subsample genes across all_pairs for train/test sets, if specified
    if random:
        all_chroms = [c for s in all_pairs for c in s]
        all_genes = list(sumstats.index[(sumstats.chrom.isin(all_chroms)) & \
                                        ~(sumstats.index.isin(xgenes))])
        gene_splits = split_genes(all_genes)
        for i in range(10):
            train_splits = [s for k, s in enumerate(gene_splits) if k != i]
            train_genes = [g for s in train_splits for g in s]
            test_genes = gene_splits[i]
            model_fit, test_rmse = \
                fit_model(features, sumstats, train_genes, test_genes, model, 
                          logit_alpha, l1_l2_mix)
            models.append((test_rmse, model_fit))

    # Otherwise, partition train/test sets based on chromosome pairs
    else:
        for i in range(len(all_pairs)):
            train_pairs = [x for k, x in enumerate(all_pairs) if k != i]
            train_chroms = [c for s in train_pairs for c in s]
            train_genes = list(sumstats.index[(sumstats.chrom.isin(train_chroms)) & \
                                              ~(sumstats.index.isin(xgenes))])
            test_genes = list(sumstats.index[(sumstats.chrom.isin(list(all_pairs[i]))) & \
                                             ~(sumstats.index.isin(xgenes))])
            model_fit, test_rmse = \
                fit_model(features, sumstats, train_genes, test_genes, model, 
                          logit_alpha, l1_l2_mix)
            models.append((test_rmse, model_fit))

    best_rmse, best_model = \
        [(r, m) for r, m in sorted(models, key=lambda x: x[0])][0]

    return best_model, best_rmse


def predict_bfdps(features, sumstats, chrompairs, xbed=None, model='logit', 
                  logit_alpha=0.1, l1_l2_mix=1, random=True, quiet=False):
    """
    Predicts bfdps for all genes with N-fold CV
    """

    pred_bfdps = pd.DataFrame(columns=['pred_bfdp'], index=features.index)

    # Load gene blacklist, if any provided
    if xbed is not None:
        xgenes = pd.read_csv(xbed, sep='\t').iloc[:, 3].values.tolist()
    else:
        xgenes = []

    # Generate predictions for each pair of chromosomes in serial
    for i in range(len(chrompairs)):
        train_pairs = [x for k, x in enumerate(chrompairs) if k != i]
        pred_chroms = sorted(list(chrompairs[i]))
        pred_genes = list(sumstats.index[sumstats.chrom.isin(pred_chroms)])

        if not quiet:
            msg = 'Now scoring {:,} genes from chromsomes {}'
            print(msg.format(len(pred_genes), ', '.join(pred_chroms)))
        
        # Fit model
        fitted_model, rmse = \
            fit_model_cv(features, sumstats, train_pairs, xgenes, model, 
                         logit_alpha, l1_l2_mix, random)

        # Predict bfdps
        pred_df = features.loc[features.index.isin(pred_genes), :]
        if model == 'logit':
            pred_vals = fitted_model.predict(pred_df)
        else:
            pred_class_probs = fitted_model.predict_proba(pred_df)
            pred_vals = pd.Series(pred_class_probs[:, 1], name='pred', index=pred_df.index)
        pred_vals.name = 'pred_bfdp'
        pred_bfdps.update(pred_vals)

    # Remove genes with no predictions (due to being excluded from feature collection)
    pred_bfdps = pred_bfdps.dropna()

    # Compute final gene score as pnorm from Z-score of pred_bfdps
    pred_bfdps['score'] = norm.sf(scale(pred_bfdps.pred_bfdp))
    pred_bfdps['quantile'] = 100 * pred_bfdps.pred_bfdp.transform('rank') / len(pred_bfdps)
    
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
    parser.add_argument('-c', '--centromeres', help='Centromeres BED.')
    parser.add_argument('-x', '--blacklist', help='Training blacklist BED.')
    parser.add_argument('-m', '--model', choices=model_options, default='logit',
                        help='Choice of classifier. [default: logit]')
    parser.add_argument('--chromsplit', action='store_true', default=False,
                        help='Use chromsome partitions for train/test sets. ' +
                        '[default: randomly sample from chromosome partitions]')
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

    # If centromeres BED file is specified, rewrite all gene chromosomes to also
    # include p/q arm designation for train/test balancing
    if args.centromeres is not None:
        arm_dict = get_chrom_arms(args.features, args.centromeres)
    else:
        arm_dict = None

    # Import gene stats
    sumstats = load_stats(args.stats, arm_dict)
    chroms = sorted(np.unique(sumstats.chrom.values))

    # Pair all chromosomes (or group chromosome arms, if optioned)
    chrompairs = pair_chroms(sumstats, chroms, args.blacklist)

    # Read gene features
    features = load_features(args.features)

    # Predict bfdps for all genes
    pred_bfdps = predict_bfdps(features, sumstats, chrompairs, args.blacklist,
                               args.model, args.logit_alpha, args.l1_l2_mix,
                               random=(not args.chromsplit))

    # Write predicted BFDPs to outfile, along with scaled score & quantile
    pred_bfdps.to_csv(outfile, sep='\t', index=True, na_rep='NA')


if __name__ == '__main__':
    main()


