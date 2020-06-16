#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Collect features per gene for an input gtf
"""


import pybedtools as pbt
import numpy as np
import pandas as pd
from pysam import faidx
import csv
from athena.mutrate import add_local_track
import argparse
from sys import stdout
from os import path
import subprocess
import gzip


def process_gtf(gtf_in):
    """
    Read gtf & filter to minimal info required
    """

    gtfbt = pbt.BedTool(gtf_in)

    # Build lists of eligible gene names and ensembl IDs
    genes, ensg_ids, transcripts = [], [], []
    ensg_to_gene, gene_to_ensg = {}, {}

    for f in gtfbt:
        if f.fields[2] == 'transcript':
            gname = f.attrs['gene_name']
            ensg_id = f.attrs['gene_id']
            tname = f.attrs['transcript_id']
            if gname not in genes:
                genes.append(gname)
            if ensg_id not in ensg_ids:
                ensg_ids.append(ensg_id)
            if tname not in transcripts:
                transcripts.append(tname)
            if ensg_id not in ensg_to_gene.keys():
                ensg_to_gene[ensg_id] = gname
            if gname not in gene_to_ensg.keys():
                gene_to_ensg[gname] = ensg_id

    # Filter & clean records in gtf
    def _filter_gtf(feature):
        """
        Restrict GTF features to desired elements
        """
        if feature.fields[2] in 'exon transcript'.split() \
        and feature.attrs['gene_name'] in genes \
        and feature.attrs['transcript_id'] in transcripts:
            return True
        else:
            return False

    attrs_to_drop = 'gene_id gene_type gene_status transcript_type ' + \
                    'transcript_status transcript_name protein_id ' + \
                    'tag ccdsid havana_gene havana_transcript'
    attrs_to_drop = attrs_to_drop.split()

    def _clean_feature(feature):
        """
        Clean unnecessary fields & info from GTF features
        """
        for key in attrs_to_drop:
            if key in feature.attrs.keys():
                feature.attrs.pop(key)
        return feature

    gtfbt = gtfbt.filter(_filter_gtf).filter(_clean_feature).saveas()

    # Make separate BedTools for exons and transcripts
    txbt = gtfbt.filter(lambda x: x.fields[2] == 'transcript').saveas()
    exonbt = gtfbt.filter(lambda x: x.fields[2] == 'exon').saveas()

    return gtfbt, txbt, exonbt, genes, ensg_ids, transcripts, ensg_to_gene, gene_to_ensg


def load_cens_tels(chrom_stats_bed):
    """
    Read BED of centromere & telomere coordinates into dictionary
    """

    chrom_stats_bt = pbt.BedTool(chrom_stats_bed)

    chrom_stats = {}

    for x in chrom_stats_bt:
        chrom = x.chrom
        if chrom not in chrom_stats.keys():
            chrom_stats[chrom] = {}
        ftype = x[3]
        fbt = pbt.BedTool('\t'.join([x.chrom, str(x.start), str(x.end)]) + '\n',
                          from_string=True)
        chrom_stats[chrom][ftype] = fbt
        if ftype == 'centromere':
            chrom_stats[chrom]['p_length'] = x.start
        if ftype == 'q_telomere':
            chrom_stats[chrom]['chrom_length'] = x.end

    for chrom in chrom_stats.keys():
        chrom_len = chrom_stats[chrom]['chrom_length']
        cen_end = chrom_stats[chrom]['centromere'][0].end
        chrom_stats[chrom]['q_length'] = chrom_len - cen_end

    return chrom_stats


def get_neighbor_dists(tx, txdf):
    """
    Get absolute distances for one query transcript and all other transcripts
    Dev note: necessary to avoid unexplained bottlenecks on some chromosomes 
              using pbt.absolute_distance()
    """

    txdf_sub = txdf[(txdf['chr'] == int(tx.chrom)) & \
                    (txdf['end'] != tx.end) & \
                    (txdf['start'] != tx.start + 1)]
    txints = txdf_sub.iloc[:, 1:].to_numpy()

    def _get_dist(r1, r2):
        """
        Calculate absolute distance between two intervals
        Adapted from: https://stackoverflow.com/a/16843530
        """
        x, y = sorted((r1, r2))
        if x[0] <= x[1] < y[0] and all(y[0] <= y[1] for y in (r1, r2)):
            return y[0] - x[1]
        return 0

    return [_get_dist(list(x), [tx.start + 1, tx.end]) for x in txints]


def get_tx_stats(genes, txbt, max_dist=1000000):
    """
    Collect dict of lengths & relations of transcripts
    """

    txdf = txbt.to_dataframe().iloc[:, [0, 3, 4]]
    txdf.columns = ['chr', 'start', 'end']

    tx_stats = {}

    for tx in txbt:
        gene_name = tx.attrs['gene_name']
        txlen = tx.length
        txcoords = '\t'.join([tx.chrom, str(tx.start), str(tx.end)]) + '\n'
        txstrand = tx.strand
        dists = get_neighbor_dists(tx, txdf)
        mindist = max([np.nanmin(dists), 1])
        n_nearby = len([i for i in dists if i <= max_dist])
        if gene_name not in tx_stats.keys():
            tx_stats[gene_name] = {'tx_coords' : txcoords,
                                   'tx_len' : [txlen],
                                   'tx_strand' : txstrand,
                                   'nearest_gene' : mindist,
                                   'genes_within_1mb' : n_nearby}
        else:
            tx_stats[gene_name]['tx_coords'] \
                = tx_stats[gene_name]['tx_coords'] + txcoords
            tx_stats[gene_name]['tx_len'].append(txlen)
            tx_stats[gene_name]['tx_strand'].append(txstrand)
            tx_stats[gene_name]['nearest_gene'].append(mindist)
            tx_stats[gene_name]['genes_within_1mb'].append(n_nearby)

    for gene in genes:
        if gene in tx_stats.keys():
            tx_stats[gene]['tx_coords'] = pbt.BedTool(tx_stats[gene]['tx_coords'],
                                                      from_string=True)
            tx_stats[gene]['tx_len'] = np.nanmedian(tx_stats[gene]['tx_len'])
            tx_stats[gene]['nearest_gene'] = np.nanmedian(tx_stats[gene]['nearest_gene'])
            tx_stats[gene]['genes_within_1mb'] = np.nanmedian(tx_stats[gene]['genes_within_1mb'])
        else:
            tx_stats[gene]['tx_coords'] = pbt.BedTool('\n', from_string=True)
            tx_stats[gene]['tx_len'] = np.nan
            tx_stats[gene]['nearest_gene'] = np.nan
            tx_stats[gene]['genes_within_1mb'] = np.nan

    return tx_stats


def calc_interval_stats(bedt, default=1):
    """
    Calculate basic statistics for a BedTool of intervals
    """

    ni = len(bedt)
    if ni > 0:
        minsize = int(np.nanmin([x.length for x in bedt]))
        maxsize = int(np.nanmax([x.length for x in bedt]))
        medsize = np.nanmedian([x.length for x in bedt])
        sumsize = int(np.nansum([x.length for x in bedt.merge()]))
    else:
        minsize = default
        maxsize = default
        medsize = default
        sumsize = default

    istats = {'n' : ni,
              'min_size' : minsize,
              'median_size' : medsize,
              'max_size' : maxsize,
              'summed_size' : sumsize}

    return istats


def get_exon_intron_stats(genes, tx_stats, exonbt, min_intron_size):
    """
    Collect exon & intron stats
    """

    exon_stats = {}
    intron_stats = {}

    # Get BedTool of exons per gene
    for ex in exonbt:
        gene_name = ex.attrs['gene_name']
        excoords = '\t'.join([ex.chrom, str(ex.start), str(ex.end)]) + '\n'
        if gene_name not in exon_stats.keys():
            exon_stats[gene_name] = {'exon_coords' : excoords}
        else:
            exon_stats[gene_name]['exon_coords'] \
                = exon_stats[gene_name]['exon_coords'] + excoords
    for gene in genes:
        if gene in exon_stats.keys():
            exon_stats[gene]['exon_coords'] = pbt.BedTool(exon_stats[gene]['exon_coords'],
                                                          from_string=True)
        else:
            exon_stats[gene]['exon_coords'] = pbt.BedTool('\n', from_string=True)

    # Get BedTool of introns per gene
    # Define introns as non-exon sections of canonical transcript that start
    # after the first base of the first exon and end before the last base of the
    # last exon, and whose length is >= min_intron_size
    for gene in genes:
        tx = tx_stats[gene]['tx_coords']
        exons = exon_stats[gene]['exon_coords']
        ex_min = np.nanmin(exons.to_dataframe()['start'])
        ex_max = np.nanmax(exons.to_dataframe()['end'])
        introns = tx.subtract(exons)\
                    .filter(lambda x: x.length >= min_intron_size and \
                                      x.start >= ex_min and \
                                      x.end <= ex_max).saveas()
        intron_stats[gene] = {'intron_coords' : introns}

    # Compute exon & intron stats per gene
    for gene in genes:
        exon_stats[gene]['stats'] = calc_interval_stats(exon_stats[gene]['exon_coords'])
        intron_stats[gene]['stats'] = calc_interval_stats(intron_stats[gene]['intron_coords'])

    return exon_stats, intron_stats


def get_chrom_pos_stats(genes, tx_stats, chrom_stats):
    """
    Calculate position of each gene relative to centromere/telomere position and
    chromosome sizes
    """

    chrom_pos_stats = {}

    for gene in genes:
        tx_bed = tx_stats[gene]['tx_coords']
        chrom = tx_bed[0].chrom

        if tx_bed[0].end <= chrom_stats[chrom]['centromere'][0].start:
            arm = 'p'
        else:
            arm = 'q'

        arm_len = chrom_stats[chrom]['{}_length'.format(arm)]
        cen_dist = tx_bed.closest(chrom_stats[chrom]['centromere'], d=True)[0][-1]
        cen_dist = np.abs(int(cen_dist))
        tel_dist = tx_bed.closest(chrom_stats[chrom]['_'.join([arm, 'telomere'])], d=True)[0][-1]
        tel_dist = np.abs(int(tel_dist))

        chrom_pos_stats[gene] = {'cen_dist' : cen_dist,
                                 'cen_dist_norm' : cen_dist / arm_len,
                                 'tel_dist' : tel_dist,
                                 'tel_dist_norm' : tel_dist / arm_len}

    return chrom_pos_stats


def calc_gc(chrom, start, end, ref_fasta, cpg=False):
    """
    Calculate GC content of an interval
    Dev note: ran into tmp file memory issues with pybedtools nucleotide_content()
    """

    seq = faidx(ref_fasta, '{}:{}-{}'.format(chrom, str(start), str(end)))\
          .replace('\n', '').upper().replace('N', '')
    n_gc = seq.count('G') + seq.count('G')

    if cpg:
        n_cpg = seq.count('CG')
        return n_gc / len(seq), n_cpg

    else:
        return n_gc / len(seq)


def get_genomic_features(genes, txbt, exonbt, tx_stats, min_intron_size=4, 
                         chrom_stats=None, ref_fasta=None, athena_tracks=None, 
                         no_scaling=False):
    """
    Collect various genomic features per gene
    """

    # Exon & intron stats
    exon_stats, intron_stats = get_exon_intron_stats(genes, tx_stats, exonbt,
                                                     min_intron_size)

    # Compile feature headers for output file
    header_cols = 'gene_length nearest_gene genes_within_1mb cds_length n_exons ' + \
                  'min_exon_size med_exon_size max_exon_size min_intron_size ' + \
                  'med_intron_size max_intron_size'
    header_cols = header_cols.split()

    # Extract basic genomic features per gene
    gfeats_tmp = {}
    for gene in genes:
        # Scale features unless specified otherwise
        if not no_scaling:
            gfeats = [np.log10(int(tx_stats[gene].get('tx_len', 'NA'))),
                      np.log10(int(tx_stats[gene].get('nearest_gene', 'NA'))),
                      int(tx_stats[gene].get('genes_within_1mb', 'NA')),
                      np.log10(exon_stats[gene]['stats']['summed_size']),
                      exon_stats[gene]['stats']['n'],
                      np.log10(exon_stats[gene]['stats']['min_size']),
                      np.log10(exon_stats[gene]['stats']['median_size']),
                      np.log10(exon_stats[gene]['stats']['max_size']),
                      np.log10(intron_stats[gene]['stats']['min_size']),
                      np.log10(intron_stats[gene]['stats']['median_size']),
                      np.log10(intron_stats[gene]['stats']['max_size'])]
        else:
            gfeats = [int(tx_stats[gene].get('tx_len', 'NA')),
                      tx_stats[gene].get('nearest_gene', 'NA'),
                      tx_stats[gene].get('genes_within_1mb', 'NA'),
                      exon_stats[gene]['stats']['summed_size'],
                      exon_stats[gene]['stats']['n'],
                      exon_stats[gene]['stats']['min_size'],
                      exon_stats[gene]['stats']['median_size'],
                      exon_stats[gene]['stats']['max_size'],
                      intron_stats[gene]['stats']['min_size'],
                      intron_stats[gene]['stats']['median_size'],
                      intron_stats[gene]['stats']['max_size']]
        gfeats_tmp[gene] = gfeats

    # Add chromosome positioning stats, if optioned
    if chrom_stats is not None:
        add_header_cols = 'cen_dist cen_dist_norm tel_dist tel_dist_norm'
        header_cols = header_cols + add_header_cols.split()
        chrom_pos_stats = get_chrom_pos_stats(genes, tx_stats, chrom_stats)
        for gene in genes:
            gfeats = gfeats_tmp[gene]
            gfeats_tmp[gene] = gfeats + list(chrom_pos_stats[gene].values())

    # Add GC content, if optioned
    if ref_fasta is not None:
        header_cols.append('gc_pct')
        for x in txbt:
            gene = x.attrs['gene_name']
            gfeats_tmp[gene].append(calc_gc(x.chrom, x.start, x.end, ref_fasta))

    # Add annotations from specified tracks, if optioned
    if athena_tracks is not None:
        add_header_cols = []
        tx_bed_str = ['\t'.join([x.chrom, str(x.start), str(x.end), x.attrs['gene_name']]) \
                      for x in txbt]
        tx_bed = pbt.BedTool('\n'.join(tx_bed_str), from_string=True)
        with open(athena_tracks) as tsv:
            reader = csv.reader(tsv, delimiter='\t')
            for track, action, tname in reader:
                add_header_cols.append(tname)
                newanno = {x[3]: float(x[4]) for x in \
                           add_local_track(tx_bed, track, action, 8, True)}
                for gene, val in newanno.items():
                    gfeats_tmp[gene].append(val)
        header_cols = header_cols + add_header_cols

    # Format output string of all genomic features per gene
    header = '\t'.join(header_cols)
    genomic_features = {}
    for gene in genes:
        gfeats_str = '\t'.join([str(x) for x in gfeats_tmp[gene]])
        genomic_features[gene] = gfeats_str

    return header, genomic_features


def calc_gtex_stats(gtex_matrix):
    """
    Compute summary stats for every gene from a GTEx expression matrix
    """

    def _summstats(vals):
        """
        Return array of summary statistics for a single gene
        """
        return np.array([np.nanmin(vals),
                         np.nanquantile(vals, q=0.25),
                         np.nanmean(vals),
                         np.nanquantile(vals, q=0.75),
                         np.nanmax(vals),
                         np.nanstd(vals),
                         ((10 ** vals) - 1 > 1).sum()])

    xstats = gtex_matrix.iloc[:, 1:].apply(_summstats, axis=1)
    xstats_df = pd.DataFrame(np.vstack(np.array(xstats)),
                             columns='min q1 mean q3 max sd gt0'.split())

    return pd.concat([gtex_matrix.loc[:, 'gene'], xstats_df], axis=1)


def load_gtex(gtex_matrix, expression_matrix=True):
    """
    Read & clean GTEx expression matrix of gene X covariate, and (optionally) compute summary stats
    """

    gtex = pd.read_csv(gtex_matrix, sep='\t')
    gtex.rename(columns={'#gene' : 'gene'}, inplace=True)

    if expression_matrix:    
        # Drop tissues with NaN values for all genes
        nonnan_cols = gtex.iloc[:, 1:].apply(lambda vals: not all(np.isnan(vals)))
        gtex = gtex.loc[:, ['gene'] + gtex.columns[1:][nonnan_cols].tolist()]

    # Handle duplicate gene symbols
    dups = [g for g, c in gtex.gene.value_counts().to_dict().items() if c > 1]
    if len(dups) > 0:
        if expression_matrix:
            # Sum untransformed values
            # Note: assumes all expression values have been log10(TPM + 1) transformed
            # (this is done by default in preprocess_GTEx.py)
            for gene in dups:
                expr_sum = gtex.loc[gtex.gene == gene, gtex.columns[1:]].\
                                apply(lambda vals: np.log10(np.nansum((10 ** vals) - 1) + 1))
                newrow = pd.Series([gene] + expr_sum.tolist(), index=gtex.columns)
                gtex = gtex.loc[gtex.gene != gene, :]
                gtex = gtex.append(newrow, ignore_index=True)
        else:
            # Otherwise, compute mean of values
            for gene in dups:
                gmean = gtex.loc[gtex.gene == gene, gtex.columns[1:]].\
                             apply(lambda vals: np.nanmean(vals))
                newrow = pd.Series([gene] + gmean.tolist(), index=gtex.columns)
                gtex = gtex.loc[gtex.gene != gene, :]
                gtex = gtex.append(newrow, ignore_index=True)

    # Compute summary stats per gene, if necessary
    if expression_matrix:
        return calc_gtex_stats(gtex)
    else:
        return gtex


def get_expression_features(genes, ensg_ids, gtex_medians, gtex_mads, gtex_pca):
    """
    Collect various expression features per gene
    """

    xfeats_tmp = {g : [] for g in genes}
    header_cols = []

    # Load GTEx medians
    if gtex_medians is not None:
        xmed_df = load_gtex(gtex_medians)
        xmed_cols = 'n_tissues_expressed median_expression_min median_expression_q1 ' + \
                    'median_expression_mean median_expression_q3 median_expression_max ' + \
                    'median_expression_sd'
        xmed_cols = xmed_cols.split()
        for gene in genes:
            if any(xmed_df.gene == gene):
                for v in 'gt0 min q1 mean q3 max sd'.split():
                    xfeats_tmp[gene].append(xmed_df[xmed_df.gene == gene].iloc[:, 1:][v].iloc[0])
            else:
                for col in xmed_cols:
                    xfeats_tmp[gene].append(0)
        header_cols += xmed_cols

    # Load GTEx MADs
    if gtex_mads is not None:
        xmad_df = load_gtex(gtex_mads)
        xmad_cols = 'expression_mad_min expression_mad_q1 expression_mad_mean ' + \
                    'expression_mad_q3 expression_mad_max expression_mad_sd'
        xmad_cols = xmad_cols.split()
        for gene in genes:
            if any(xmad_df.gene == gene):
                for v in 'min q1 mean q3 max sd'.split():
                    xfeats_tmp[gene].append(xmad_df[xmad_df.gene == gene].iloc[:, 1:][v].iloc[0])
            else:
                for col in xmad_cols:
                    xfeats_tmp[gene].append(0)
        header_cols += xmad_cols

    # Load GTEx principal components
    if gtex_pca is not None:
        pca_df = load_gtex(gtex_pca, expression_matrix=False)
        pca_cols = pca_df.columns.tolist()[1:]
        for gene in genes:
            if any(pca_df.gene == gene):
                for v in pca_cols:
                    xfeats_tmp[gene].append(pca_df.loc[pca_df.gene == gene, v].iloc[0])
            else:
                for col in pca_cols:
                    xfeats_tmp[gene].append(0)
        header_cols += xmad_cols

    # Format output string of all expression features per gene
    header = '\t'.join(header_cols)
    expression_features = {}
    for gene in genes:
        xfeats_str = '\t'.join([str(x) for x in xfeats_tmp[gene]])
        expression_features[gene] = xfeats_str

    return header, expression_features


def get_constraint_features(genes, ensg_ids, tx_stats, txbt, exonbt, gene_to_ensg,
                            gnomad_tsv, exac_cnv_tsv, rvis_tsv, eds_tsv, hi_tsv, 
                            ref_fasta, phastcons_url, promoter_size=1000):
    """
    Collect various evolutionary constraint features per gene
    """

    cfeats_tmp = {g : [] for g in genes}

    # Compile feature headers for output file
    header_cols = []

    # Parse gnomAD constraint stats
    if gnomad_tsv is not None:
        # Load gnomAD data
        gnomad = pd.read_csv(gnomad_tsv, delimiter='\t', compression='gzip')
        keep_gnomad_cols = 'gene pLI pNull pRec oe_mis oe_lof oe_mis_upper ' + \
                           'oe_lof_upper mis_z lof_z'
        gnomad = gnomad.loc[gnomad.gene.isin(genes), keep_gnomad_cols.split()]
        # Fill in missing genes and values with overall means
        gnomad_means = gnomad.iloc[:, 1:].apply(np.nanmean).to_dict()
        gnomad.fillna(gnomad_means, axis=0, inplace=True)
        for gene in genes:
            if not any(gnomad.gene == gene):
                newrow = pd.Series([gene] + list(gnomad_means.values()), 
                                   index=gnomad.columns)
                gnomad = gnomad.append(newrow, ignore_index=True)
        # Add values to cfeats per gene
        for gene in genes:
            gvals = gnomad.loc[gnomad.gene == gene, :].values.tolist()[0][1:]
            cfeats_tmp[gene] += gvals
        header_cols += ['gnomad_' + x for x in list(gnomad.columns)[1:]]

    # Add ExAC CNV Z-score
    if exac_cnv_tsv is not None:
        # Load ExAC CNV data
        exac = pd.read_csv(exac_cnv_tsv, delimiter='\t')
        exac_vals = exac.set_index(exac.gene).loc[:, 'cnv_z'].fillna(0).to_dict()
        # Add values to cfeats per gene
        for gene in genes:
            cfeats_tmp[gene].append(exac_vals.get(gene, 0))
        header_cols.append('exac_cnv_z')

    # Add RVIS, if optioned. Assumes RVIS March 2017 release corresponding to gnomAD v2.0
    if rvis_tsv is not None:
        rvis = pd.read_csv(rvis_tsv, delimiter='\t', usecols=[0, 2, 3], skiprows=1,
                           names='gene rvis rvis_pct'.split())
        rvis = rvis.loc[rvis.gene.isin(genes), :]
        rvis_means = rvis.iloc[:, 1:].apply(np.nanmean).values.tolist()
        for gene in genes:
            if any(rvis.gene == gene):
                gvals = rvis.loc[rvis.gene == gene, :].values.tolist()[0][1:]
                cfeats_tmp[gene] += gvals
            else:
                cfeats_tmp[gene] += rvis_means
        header_cols += rvis.columns.tolist()[1:]

    # Add EDS, if optioned
    if eds_tsv is not None:
        eds = pd.read_csv(eds_tsv, delimiter='\t', names='ensg eds'.split(), skiprows=1)
        eds = eds.loc[eds.ensg.isin([x.split('.')[0] for x in gene_to_ensg.values()]), :]
        eds_mean = np.nanmean(eds.eds)
        for gene in genes:
            ensg = gene_to_ensg[gene].split('.')[0]
            if any(eds.ensg == ensg):
                cfeats_tmp[gene].append(eds.eds[eds.ensg == ensg].values[0])
            else:
                cfeats_tmp[gene].append(eds_mean)
        header_cols.append('eds')

    # Add Hurles HI scores, if optioned
    if hi_tsv is not None:
        hi_dat = {f.name.split('|')[0] : float(f.score) for f in pbt.BedTool(hi_tsv)}
        hi_avg = np.nanmean(list(hi_dat.values()))
        for gene in genes:
            cfeats_tmp[gene].append(hi_dat.get(gene, hi_avg))
        header_cols.append('hurles_hi')

    # Make dictionary of promoter coordinates
    promoters = {}
    for gene in genes:
        chrom = tx_stats[gene]['tx_coords'][0].chrom
        if tx_stats[gene]['tx_strand'] == '+':
            prom_end = tx_stats[gene]['tx_coords'][0].start
            prom_start = int(np.nanmax([prom_end - promoter_size, 0]))
        else:
            prom_start = tx_stats[gene]['tx_coords'][0].stop
            prom_end = prom_start + promoter_size
        prom_str = '{}\t{}\t{}\t{}'.format(chrom, prom_start, prom_end, gene)
        prom_bt = pbt.BedTool(prom_str, from_string=True)
        promoters[gene] = {'chrom' : chrom,
                           'start' : prom_start,
                           'end' : prom_end,
                           'prom_bt' : prom_bt}

    # Get CG pct and CpG count for gene promoters if ref_fasta specified
    if ref_fasta is not None:
        for gene in genes:
            prom_gc, prom_cpg = \
                calc_gc(promoters[gene]['chrom'], promoters[gene]['start'], 
                        promoters[gene]['end'], ref_fasta, cpg=True)
            cfeats_tmp[gene].append(prom_gc)
            cfeats_tmp[gene].append(prom_cpg)
        header_cols += 'promoter_gc_pct promoter_cpg_count'.split()

    # Get promoter, exon, and gene body conservation if phastcons_url is provided
    if phastcons_url is not None:
        # Make master pbt.BedTool of all promoters for one-shot conservation calculation
        all_prom_str = ''.join([str(v['prom_bt'][0]) for v in promoters.values()])
        all_prom_bt = pbt.BedTool(all_prom_str, from_string=True)
        all_prom_cons = add_local_track(all_prom_bt, phastcons_url, 'map-mean', 8, True)
        all_prom_cons_df = \
            all_prom_cons.to_dataframe(names='chr start end gene cons'.split())
        # Calculate conservation for all exons
        all_ex_cons = add_local_track(exonbt, phastcons_url, 'map-mean', 8, True)
        all_ex_cons_df = all_ex_cons.to_dataframe().iloc[:, np.r_[0, 3:5, 8:10]]
        all_ex_cons_df.columns = 'chr start end info cons'.split()
        all_ex_cons_df['size'] = all_ex_cons_df['end'] - all_ex_cons_df['start']
        # Calculate conservation for all gene bodies
        all_tx_cons = add_local_track(txbt, phastcons_url, 'map-mean', 8, True)
        all_tx_cons_df = all_tx_cons.\
                             to_dataframe(names='chr x y start end sc st z info cons'.split()).\
                             iloc[:, np.r_[0, 3:5, 8:10]]
        # Annotate for each gene
        for gene in genes:
            prom_cons = all_prom_cons_df.loc[all_prom_cons_df.gene == gene]['cons'].iloc[0]
            ex_keep = all_ex_cons_df['info'].str.contains('gene_name "{}"'.format(gene))
            if ex_keep.sum() > 0:
                ex_cons_df = all_ex_cons_df.loc[ex_keep, 'size cons'.split()]
                ex_cons_df = ex_cons_df.loc[ex_cons_df.cons.astype(str) != '.', :]
                if len(ex_cons_df) > 0:
                    try:
                        ex_cons = np.ma.average(ex_cons_df['cons'].astype(float).to_numpy(), 
                                                weights=ex_cons_df['size'].to_numpy())
                    except:
                        # Debug message
                        print('Failed on exon phastCons calculation for ' + gene)
                        exit(1)
                else:
                    ex_cons = 0
            else:
                ex_cons = 0
            tx_keep = all_tx_cons_df['info'].str.contains('gene_name "{}"'.format(gene))
            if tx_keep.sum() > 0:
                tx_cons = all_tx_cons_df.loc[tx_keep, 'cons'].values[0]
                if str(tx_cons) == '.':
                    tx_cons = 0
            else:
                tx_cons = 0
            cfeats_tmp[gene] += [prom_cons, ex_cons, tx_cons]
        header_cols += 'promoter_phastcons exon_phastcons gene_body_phastcons'.split()

    # Format output string of all constraint features per gene
    header = '\t'.join(header_cols)
    constraint_features = {}
    for gene in genes:
        cfeats_str = '\t'.join([str(x) for x in cfeats_tmp[gene]])
        constraint_features[gene] = cfeats_str

    return header, constraint_features


def write_outbed(outbed, header, genes, txbt, tx_stats, genomic_features,
                 expression_features, constraint_features):
    """
    Format output table of features and write to output BED file
    """

    # Write header
    outbed.write(header + '\n')

    # Write one line per gene with all features optioned at runtime
    for gene in genes:
        gcoords_df = tx_stats[gene]['tx_coords'].to_dataframe()
        gcoords = '\t'.join([tx_stats[gene]['tx_coords'][0].chrom,
                             str(int(np.floor(np.nanmedian(gcoords_df['start'])))),
                             str(int(np.ceil(np.nanmedian(gcoords_df['end']))))])
        outstr = gcoords + '\t' + gene

        if genomic_features is not None:
            outstr = outstr + '\t' + genomic_features[gene]

        if expression_features is not None:
            outstr = outstr + '\t' + expression_features[gene]

        if constraint_features is not None:
            outstr = outstr + '\t' + constraint_features[gene]            

        outbed.write(outstr + '\n')

    # Close connection to flush buffer
    outbed.close()


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gtf', help='GTF of genes to consider.')
    parser.add_argument('--get-genomic', action='store_true', help='Collect genomic ' +
                        'features. [default: False]')
    parser.add_argument('--get-expression', action='store_true', help='Collect ' +
                        'gene expression features. [default: False]')
    parser.add_argument('--get-constraint', action='store_true', help='Collect ' +
                        'evolutionary constraint features. [default: False]')
    parser.add_argument('--centro-telo-bed', help='BED file indicating ' + 
                        'centromere & telomere positions per chromosome.')
    parser.add_argument('--ref-fasta', help='Reference fasta file. Only necessary ' +
                        'to compute GC content if --get-genomic is specified, and ' +
                        'promoter stats if --get-constraint is specified.')
    parser.add_argument('--athena-tracks', help='Annotate transcripts with track(s). ' +
                        'Must be formatted as for athena annotate-bins.')
    parser.add_argument('--gtex-medians', help='GTEx gene X tissue expression medians. ' +
                        'Only used if --get-expression is specified.')
    parser.add_argument('--gtex-mads', help='GTEx gene X tissue expression MADs. ' +
                        'Only used if --get-expression is specified.')
    parser.add_argument('--gtex-pca', help='GTEx gene X tissue principal components. ' +
                        'Only used if --get-expression is specified.')
    parser.add_argument('--gnomad-constraint', help='gnomAD constraint tsv. Only ' +
                        'used if --get-constraint is specified.')
    parser.add_argument('--exac-cnv', help='ExAC CNV constraint tsv. Only used ' +
                        'if --get-constraint is specified.')
    parser.add_argument('--rvis-tsv', help='RVIS tsv. Only used if --get-constraint ' +
                        'is specified.')
    parser.add_argument('--eds-tsv', help='EDS tsv. Only used if --get-constraint ' +
                        'is specified.')
    parser.add_argument('--hi-tsv', help='Hurles HI scores tsv. Only used if ' +
                        '--get-constraint is specified.')
    parser.add_argument('--phastcons-bw-url', help='URL to phastCons bigWig.',
                        default='http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.bw')
    parser.add_argument('--min-intron-size', type=int, default=4, help='Minimum ' +
                        'size of intron to retain (bp). [default: 4]')
    parser.add_argument('--no-scaling', action='store_true', help='Do not perform ' +
                        'scaling & transformations to features.')
    parser.add_argument('-o', '--outbed', help='Path to output BED file. ' +
                        '[default: stdout]')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output BED with bgzip.')
    
    args = parser.parse_args()

    # Open connection to output file
    if args.outbed is None \
    or args.outbed in 'stdout -'.split():
        outbed = stdout
    else:
        if path.splitext(args.outbed)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outbed_path = path.splitext(args.outbed)[0]
        else:
            outbed_path = args.outbed
        outbed = open(outbed_path, 'w')
    outbed_header = '\t'.join('#chr start end gene'.split())

    # Load & filter input GTF
    gtfbt, txbt, exonbt, genes, ensg_ids, transcripts, ensg_to_gene, gene_to_ensg \
        = process_gtf(args.gtf)

    # Get transcript stats
    tx_stats = get_tx_stats(genes, txbt)

    # Get genomic features, if optioned
    if args.get_genomic:
        if args.centro_telo_bed is not None:
            chrom_stats = load_cens_tels(args.centro_telo_bed)
        else:
            chrom_stats = None
        header_add, genomic_features = get_genomic_features(genes, txbt, exonbt, 
                                                            tx_stats,
                                                            args.min_intron_size,
                                                            chrom_stats,
                                                            args.ref_fasta,
                                                            args.athena_tracks,
                                                            args.no_scaling)
        outbed_header = outbed_header + '\t' + header_add
    else:
        genomic_features = None

    # Get expression stats, if optioned
    if args.get_expression:
        header_add, expression_features = \
            get_expression_features(genes, ensg_ids, args.gtex_medians, args.gtex_mads, 
                                    args.gtex_pca)
        outbed_header = outbed_header + '\t' + header_add
    else:
        expression_features = None

    # Get constraint stats, if optioned
    if args.get_constraint:
        header_add, constraint_features = \
            get_constraint_features(genes, ensg_ids, tx_stats, txbt, exonbt, 
                                    gene_to_ensg, args.gnomad_constraint, 
                                    args.exac_cnv, args.rvis_tsv, args.eds_tsv, 
                                    args.hi_tsv, args.ref_fasta, 
                                    args.phastcons_bw_url)
        outbed_header = outbed_header + '\t' + header_add
    else:
        constraint_features = None

    # Format output table of features
    write_outbed(outbed, outbed_header, genes, txbt, tx_stats, genomic_features,
                 expression_features, constraint_features)
    if args.outbed is not None \
    and args.outbed not in 'stdout -'.split() \
    and args.bgzip:
        subprocess.run(['bgzip', '-f', outbed_path])


if __name__ == '__main__':
    main()

