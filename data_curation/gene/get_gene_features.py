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

    return gtfbt, txbt, exonbt, genes, ensg_ids, transcripts


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


def load_gtex(gtex_matrix):
    """
    Read & clean GTEx expression matrix of gene X tissue, and compute summary stats
    """

    gtex = pd.read_csv(gtex_matrix, sep='\t')
    gtex.rename(columns={'#gene' : 'gene'}, inplace=True)
    
    # Drop tissues with NaN values for all genes
    nonnan_cols = gtex.iloc[:, 1:].apply(lambda vals: not all(np.isnan(vals)))
    gtex = gtex.loc[:, ['gene'] + gtex.columns[1:][nonnan_cols].tolist()]

    # Handle duplicate gene symbols by summing their untransformed values
    # Note: assumes all expression values have been log10(TPM + 1) transformed
    # (this is done by default in preprocess_GTEx.py)
    dups = [g for g, c in gtex.gene.value_counts().to_dict().items() if c > 1]
    if len(dups) > 0:
        for gene in dups:
            expr_sum = gtex.loc[gtex.gene == gene, gtex.columns[1:]].\
                           apply(lambda vals: np.log10(np.nansum((10 ** vals) - 1) + 1))
            newrow = pd.Series([gene] + expr_sum.tolist(), index=gtex.columns)
            gtex = gtex.loc[gtex.gene != gene, :]
            gtex = gtex.append(newrow, ignore_index=True)

    # Compute summary stats per gene
    return calc_gtex_stats(gtex)


def get_expression_features(genes, ensg_ids, gtex_medians, gtex_mads):
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

    # Format output string of all expression features per gene
    header = '\t'.join(header_cols)
    expression_features = {}
    for gene in genes:
        xfeats_str = '\t'.join([str(x) for x in xfeats_tmp[gene]])
        expression_features[gene] = xfeats_str

    return header, expression_features


def get_constraint_features(genes, ensg_ids, tx_stats, ref_fasta, 
                            promoter_size=1000):
    """
    Collect various evolutionary constraint features per gene
    """

    cfeats_tmp = {g : [] for g in genes}

    # Compile feature headers for output file
    header_cols = 'promoter_gc_pct promoter_cpg_count'
    header_cols = header_cols.split()

    # Get CG pct and CpG count for gene promoters
    for gene in genes:
        chrom = tx_stats[gene]['tx_coords'][0].chrom
        if tx_stats[gene]['tx_strand'] == '+':
            prom_end = tx_stats[gene]['tx_coords'][0].start
            prom_start = int(np.nanmax([prom_end - 1000, 0]))
        else:
            prom_start = tx_stats[gene]['tx_coords'][0].stop
            prom_end = prom_start + 1000
        prom_gc, prom_cpg = calc_gc(chrom, prom_start, prom_end, ref_fasta, cpg=True)
        cfeats_tmp.append(prom_gc)
        cfeats_tmp.append(prom_cpg)

        import pdb; pdb.set_trace()





def write_outbed(outbed, header, genes, txbt, tx_stats, genomic_features,
                 expression_features):
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
    parser.add_argument('--gtex-medians', help='GTEx gene X tissue expression medians. ' +
                        'Only used if --get-expression is specified.')
    parser.add_argument('--gtex-mads', help='GTEx gene X tissue expression MADs. ' +
                        'Only used if --get-expression is specified.')
    parser.add_argument('--athena-tracks', help='Annotate transcripts with track(s). ' +
                        'Must be formatted as for athena annotate-bins.')
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
    gtfbt, txbt, exonbt, genes, ensg_ids, transcripts = process_gtf(args.gtf)

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
            get_expression_features(genes, ensg_ids, args.gtex_medians, args.gtex_mads)
        outbed_header = outbed_header + '\t' + header_add
    else:
        expression_features = None

    # Get constraint stats, if optioned
    if args.get_constraint:
        header_add, constraint_features = \
            get_constraint_features(genes, ensg_ids, tx_stats, args.ref_fasta)
    else:
        constraint_features = None

    # Format output table of features
    write_outbed(outbed, outbed_header, genes, txbt, tx_stats, genomic_features,
                 expression_features)
    if args.outbed is not None \
    and args.outbed not in 'stdout -'.split() \
    and args.bgzip:
        subprocess.run(['bgzip', '-f', outbed_path])


if __name__ == '__main__':
    main()

