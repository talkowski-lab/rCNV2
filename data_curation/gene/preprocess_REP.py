#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Preprocess Roadmap Epigenomics ChromHMM data and distill various summary matrices
"""


import pybedtools as pbt
import csv
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import median_absolute_deviation as mad
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import subprocess
import argparse


def load_gtf(gtf_in):
    """
    Read gtf & filter to minimal info required
    """

    gtfbt = pbt.BedTool(gtf_in)

    # Build lists of eligible gene names and transcript IDs
    genes, transcripts = [], []

    for f in gtfbt:
        if f.fields[2] == 'transcript':
            gname = f.attrs['gene_name']
            tname = f.attrs['transcript_id']
            if gname not in genes:
                genes.append(gname)
            if tname not in transcripts:
                transcripts.append(tname)

    # Filter & clean records in gtf
    def _filter_gtf(feature):
        """
        Restrict GTF features to desired elements
        """
        if feature.attrs['gene_name'] in genes \
        and feature.attrs['transcript_id'] in transcripts \
        and feature.fields[2] == 'transcript':
            return True
        else:
            return False

    gtfbt_str = ''
    for x in gtfbt.filter(_filter_gtf):
        gstr = '\t'.join([str(v) for v in [x.chrom, x.start, x.end, x.attrs['gene_name']]]) + '\n'
        gtfbt_str += gstr

    gtfbt = pbt.BedTool(gtfbt_str, from_string=True).saveas()

    return gtfbt, genes, transcripts


def load_state_manifest(state_manifest_in):
    """
    Load dict of Roadmap ChromHMM states
    """

    states = []

    with open(state_manifest_in) as infile:
        reader = csv.reader(infile, delimiter='\t')
        for state, name in reader:
            code = '_'.join([state.replace('E', ''), name])
            states.append(code)

    return states


def load_sample_manifest(sample_manifest_in):
    """
    Load dict of Roadmap epigenome samples
    """

    samples = {}

    with open(sample_manifest_in) as infile:
        reader = csv.reader(infile, delimiter='\t')
        for eid, mnemonic, name, anatomy, biotype in reader:
            if '#' in eid:
                continue
            samples[eid] = {'mnemonic' : mnemonic, 'name' : name,
                            'anatomy' : anatomy, 'biotype' : biotype}

    return samples


def parse_chromhmm_bed(path, states):
    """
    Splits a ChromHMM BED by state
    Returns a dict of {state : pbt.BedTool}
    """

    bts = {}

    # Load BED as pd.DataFrame
    cdf = pd.read_csv(path, sep='\t', names='chrom start end state'.split())
    cdf['chrom'] = cdf.chrom.apply(lambda x: x.replace('chr', ''))

    # Split by state & convert to pbt.BedTool
    for state in states:
        bts[state] = pbt.BedTool.\
                         from_dataframe(cdf.loc[cdf.state == state, :].iloc[:, 0:3]).\
                         sort().merge()

    return bts


def load_chromhmm_beds(chromdir, samples, states, bed_suffix):
    """
    Loads all ChromHMM state BEDs able to be found in chromdir
    Splits each BED by state
    Returns a dict of {eid : {state : pbt.BedTool}}
    """

    chrombeds = {}

    for eid in samples.keys():

        # Check for ChromHMM BED, and load if found
        bedpath = chromdir + '/' + eid + bed_suffix
        if Path(bedpath).exists():
            print('Loading ChromHMM BED for {} ({})'.\
                  format(samples[eid]['name'], bedpath))
            chrombeds[eid] = parse_chromhmm_bed(bedpath, states)
        else:
            print('Warning: unable to locate ChromHMM BED for {} ({})'.\
                  format(samples[eid]['name'], bedpath))
            continue

    print('Finished loading ChromHMM for {} samples'.format(len(chrombeds)))

    return chrombeds


def calc_chrom_coverage(genes_bt, chrombeds, state, samples):
    """
    Compute coverage matrix of gene X tissue for a single ChromHMM state
    """

    gdf = genes_bt.to_dataframe(names='chrom start end gene'.split())

    for eid, cbts in  chrombeds.items():
        print('Annotating {} in {}...'.format(state, samples[eid]['name']))
        cov_bt = genes_bt.coverage(cbts[state])
        vals = pd.Series([x[-1] for x in cov_bt]).astype(float)
        gdf['_'.join([eid, state])] = vals

    # Reformat dataframe as matrix of values with gene name as index
    gdf.index = gdf.gene
    gdf.drop(columns='chrom start end gene'.split(), inplace=True)

    return gdf


def compute_summary_stats(chromcovs, states):
    """
    Compute summary stats of genes X states
    Multi-return of pd.DataFrame
    """

    # Compute summary dataframes
    mins = pd.concat([cdf.apply(np.nanmin, axis=1) for cdf in chromcovs.values()], axis=1)
    q1s = pd.concat([cdf.apply(np.nanquantile, q=0.25, axis=1) for cdf in chromcovs.values()], axis=1)
    means = pd.concat([cdf.apply(np.nanmean, axis=1) for cdf in chromcovs.values()], axis=1)
    medians = pd.concat([cdf.apply(np.nanmedian, axis=1) for cdf in chromcovs.values()], axis=1)
    q3s = pd.concat([cdf.apply(np.nanquantile, q=0.75, axis=1) for cdf in chromcovs.values()], axis=1)
    maxs = pd.concat([cdf.apply(np.nanmax, axis=1) for cdf in chromcovs.values()], axis=1)
    sds = pd.concat([cdf.apply(np.nanstd, axis=1) for cdf in chromcovs.values()], axis=1)
    mads = pd.concat([cdf.apply(mad, nan_policy='omit', axis=1) for cdf in chromcovs.values()], axis=1)
    
    # Rename columns
    for df in [mins, q1s, means, medians, q3s, maxs, sds, mads]:
        df.columns = states
    
    return mins, q1s, means, medians, q3s, maxs, sds, mads


def pca_states(chromcovs, n_pcs=20):
    """
    Reduce an expression matrix to its n_pcs top principal components
    """

    # Column-wise join of all states across all tissues
    matrix = pd.concat(list(chromcovs.values()), axis=1)

    # Clean input matrix
    X = StandardScaler().fit_transform(matrix.fillna(value=0))

    # PCA
    pca = PCA(n_components=n_pcs).fit_transform(X)
    pc_names = ['chromatin_component_' + str(i) for i in range(1, n_pcs+1)]
    pcadf = pd.DataFrame(pca, columns=pc_names)
    pcadf.index = matrix.index

    # Convert back to dataframe with genes as index
    return pcadf


def write_matrix(matrix, filename, gzip):
    """
    Write matrix to file and gzip if optioned
    """

    matrix.index.rename(name='#gene', inplace=True)
    
    matrix.to_csv(filename, sep='\t', index=True, na_rep='NA')

    if gzip:
        subprocess.run(['gzip', '-f', filename])


def main():
    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gtf', help='gtf of canonical transcripts')
    parser.add_argument('chromdir', help='directory containing ChromHMM BEDs')
    parser.add_argument('--state-manifest', required=True, help='ChromHMM state manifest tsv')
    parser.add_argument('--sample-manifest', required=True, help='REP sample manifest tsv')
    parser.add_argument('--bed-suffix', help='Suffix for all ChromHMM BEDs.',
                        default='_18_core_K27ac_mnemonics.bed.gz')
    parser.add_argument('-p', '--prefix', default='REP_summary', help='Prefix ' +
                        'for output matrices')
    parser.add_argument('--n-pcs', default=20, help='Number of principal components ' +
                        'to retain', type=int)
    parser.add_argument('-z', '--gzip', action='store_true', help='Gzip outputs')
    args = parser.parse_args()

    # Load GTF & extract transcript coordinates
    genes_bt, genes, transcripts = load_gtf(args.gtf)

    # Load Roadmap manifests
    states = load_state_manifest(args.state_manifest)
    samples = load_sample_manifest(args.sample_manifest)

    # Load Roadmap ChromHMM BEDs & split by state
    chrombeds = load_chromhmm_beds(args.chromdir, samples, states, args.bed_suffix)

    # Compute coverage stats per sample per state
    chromcovs = {state : calc_chrom_coverage(genes_bt, chrombeds, state, samples) for state in states}

    # Compute summary stats across all tissues per gene per state
    rep_mins, rep_q1s, rep_means, rep_medians, rep_q3s, rep_maxs, rep_sds, rep_mads = \
        compute_summary_stats(chromcovs, states)

    # PCA of all states across tissues
    rep_pca = pca_states(chromcovs, args.n_pcs)

    # Write matrices to output files
    write_matrix(rep_mins, args.prefix + '.min.tsv', args.gzip)
    write_matrix(rep_q1s, args.prefix + '.q1.tsv', args.gzip)
    write_matrix(rep_medians, args.prefix + '.median.tsv', args.gzip)
    write_matrix(rep_means, args.prefix + '.mean.tsv', args.gzip)
    write_matrix(rep_q3s, args.prefix + '.q3.tsv', args.gzip)
    write_matrix(rep_maxs, args.prefix + '.max.tsv', args.gzip)
    write_matrix(rep_sds, args.prefix + '.sd.tsv', args.gzip)
    write_matrix(rep_mads, args.prefix + '.mad.tsv', args.gzip)
    write_matrix(rep_pca, args.prefix + '.pca.tsv', args.gzip)


if __name__ == '__main__':
    main()
