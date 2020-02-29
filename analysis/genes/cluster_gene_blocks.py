#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Refine significant gene associations
"""


from os import path
import gzip
import pandas as pd
import pybedtools as pbt
import networkx as nx
import argparse
import subprocess


def get_bed_header(bedpath):
    """
    Get header from BED file
    """

    if path.splitext(bedpath)[1] in '.gz .bgz .bgzip'.split():
        header = gzip.GzipFile(bedpath).readline().decode('utf-8').rstrip('\n')
    else:
        header = open(bedpath).readline().rstrip('\n')

    return header


def load_sig_df(sig_matrix):
    """
    Load a matrix of True/False significance indications per gene per phenotype
    """

    # Read sig_matrix header to get phenotype keys
    sig_header = get_bed_header(sig_matrix)
    sig_cols = sig_header.rstrip().split('\t')
    sig_cols = [x.split('.')[0].replace('#', '').replace('HP', 'HP:') for x in sig_cols]

    # Read sig matrix
    sig_df = pd.read_table(sig_matrix, names=sig_cols, comment='#')
    sig_df = sig_df[sig_df.iloc[:, 4:].apply(any, axis=1)]
    coords = sig_df.iloc[:, 0:4]
    genes = sig_df['gene'].tolist()
    sig_labs = sig_df.iloc[:, 4:]
    
    # Add column with list of significant HPOs per gene
    sig_hpos_l = []
    for gene in genes:
        sig_hpos = sig_df.columns[4:][sig_df[sig_df['gene'] == gene].iloc[:, 4:].\
                       to_numpy().flatten()].tolist()
        sig_hpos_l.append(';'.join(sig_hpos))
    sig_out = coords
    sig_out['hpos'] = sig_hpos_l

    return sig_out


def cluster_genes(genes_df, dist):
    """
    Cluster genes into blocks, matching on HPO
    """

    genes_bt = pbt.BedTool.from_dataframe(genes_df)

    def _pad_gene(gene, dist):
        gene.start = max([0, gene.start - dist])
        gene.end = gene.end + dist
        return gene

    genes_padded = genes_bt.each(_pad_gene, dist=dist)

    hits = genes_bt.intersect(genes_padded, wa=True, wb=True)

    def _filter_hit(hit):
        # Exclude self-hits
        if hit.name == hit[8]:
            return False
        # Require at least one matching HPO
        else:
            hpos_a = set(hit[4].split(';'))
            hpos_b = set(hit[9].split(';'))
            return len(hpos_a.intersection(hpos_b)) > 0

    hits = hits.filter(_filter_hit).saveas()

    # Define blocks
    blocks = nx.Graph()
    for pair in hits:
        gene_a = pair.name
        gene_b = pair[8]
        blocks.add_nodes_from([gene_a, gene_b])
        blocks.add_edge(gene_a, gene_b)

    # Create bedtool of blocks
    blocks_bt_str = ''
    for genes in nx.connected_components(blocks):
        bchrom = str(genes_df['chr'][genes_df['gene'] == list(genes)[0]].to_numpy()[0])
        bstart = str(min(genes_df['start'][genes_df['gene'].isin(genes)]))
        bend = str(max(genes_df['end'][genes_df['gene'].isin(genes)]))
        blocks_bt_str += '\t'.join([bchrom, bstart, bend, ','.join(list(genes))]) + '\n'

    return pbt.BedTool(blocks_bt_str, from_string=True)


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sig_matrix', help='Matrix of significance labels.')
    parser.add_argument('-d', '--dist', help='Distance to use when clustering.',
                        default=1000000, type=int)
    parser.add_argument('-o', '--outfile', help='Output BED of clustered blocks.')
    parser.add_argument('-z', '--bgzip', help='Bgzip output BED.', 
                        action='store_true', default=False)
    args = parser.parse_args()

    # Load significant genes
    genes_df = load_sig_df(args.sig_matrix)

    # Cluster genes into blocks
    blocks = cluster_genes(genes_df, dist=args.dist)

    # Write blocks to output file
    if args.outfile is None \
    or args.outfile in 'stdout -'.split():
        outfile = stdout
    else:
        if path.splitext(args.outfile)[-1] in '.gz .bz .bgz .bgzip .gzip'.split():
            outfile_path = path.splitext(args.outfile)[0]
        else:
            outfile_path = args.outfile
        outfile = open(outfile_path, 'w')
    for block in blocks:
        outfile.write('\t'.join(block.fields) + '\n')
    outfile.close()

    # Bgzip output, if optioned
    if args.bgzip \
    and args.outfile is not None \
    and args.outfile not in 'stdout -'.split():
        subprocess.run(['bgzip', '-f', outfile_path])


if __name__ == '__main__':
    main()
