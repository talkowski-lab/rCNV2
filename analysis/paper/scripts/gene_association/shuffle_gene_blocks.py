#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Randomly shuffles blocks of colinear genes the genome
"""


import pandas as pd
import pybedtools as pbt
import argparse
from os.path import splitext
from sys import stdout
import subprocess


def load_blocks(blocks_in):
    """
    Load & reformat original set of gene blocks (to be shuffled)
    Returns: dict of block_id : {block attributes} pairs
    """

    blocks = pd.read_csv(blocks_in, sep='\t').iloc[:, :3].astype(str)
    blocks.columns = 'block_id cnv genes'.split()
    blocks.index = blocks['block_id']
    blocks.drop(columns='block_id', inplace=True)
    try:
        blocks_dict = blocks.to_dict(orient='index')
    except:
        import pdb; pdb.set_trace()
    for bid in blocks_dict.keys():
        genes = [g for g in blocks_dict[bid]['genes'].split(';') \
                 if g is not None and g not in 'NA nan NAN'.split() and g != '']
        blocks_dict[bid]['genes'] = genes
        blocks_dict[bid]['n_genes'] = len(genes)

    return blocks_dict


def shuffle_blocks(blocks, univ_df, univ_bt, seed, genome=None, quiet=False):
    """
    Randomly shuffle a dict of colinear gene blocks (output by load_blocks)
    Return: pd.DataFrame of shuffled blocks (matches args.blocks format)
    """

    dist_bt_colnames = 'chrom start end gene chromB startB endB geneB dist'.split()
    newgene_dict = {bid : {'cnv' : bvals['cnv']} for bid, bvals in blocks.items()}

    # Defines 4 x len(blocks) candidate starting points for seed
    # (Need to sample more than len(blocks) due to possibility of collisions)
    index_genes = univ_df.sample(n=4*len(blocks), random_state=seed)
    index_genes.reset_index(inplace=True)
    index_genes.drop(columns='index', inplace=True)
    used_genes = {'DEL' : [], 'DUP' : []}
    k=0

    for bid, bvals in blocks.items():
        cnv = bvals['cnv']
        ngenes = bvals['n_genes']
        
        if ngenes == 0:
            newgene_dict[bid]['genes'] = ''
            continue

        keep_sampling = True
        while keep_sampling:
            index_gene_bt = pbt.BedTool.from_dataframe(index_genes.loc[[k]])
            k += 1
            if genome is not None:
                cbt = index_gene_bt.closest(univ_bt, k=ngenes, g=genome)
            else:
                ichrom = index_gene_bt[0].chrom
                cbt = index_gene_bt.closest(univ_bt.filter(lambda x: x.chrom == ichrom), k=ngenes)
            newgenes = [x[7] for x in cbt]

            if len(set(newgenes).intersection(set(used_genes[cnv]))) > 0:
                if not quiet:
                    print('Overlapping gene collision for ' + bid + '; retrying...')
            else:
                keep_sampling = False

        used_genes[cnv] += newgenes
        newgene_dict[bid]['genes'] = ';'.join(sorted(newgenes))

    newgenes_df = pd.DataFrame.from_dict(newgene_dict, orient='index')
    newgenes_df['#region_id'] = newgenes_df.index.map(lambda x: '_'.join(['perm' + str(seed), x]))

    return newgenes_df.reset_index().loc[:, '#region_id cnv genes'.split()]


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('blocks', help='tsv of gene blocks to be shuffled. Requires ' +
                        'three columns: block id, cnv type, and semicolon-delimited ' +
                        'list of the genes to be shuffled.')
    parser.add_argument('genes', help='BED4 of all genes with coordinates.')
    parser.add_argument('-g', '--genome', required=True, help='BEDTools-style ' +
                        'genome file.')
    parser.add_argument('-n', '--n-perms', type=int, default=1, help='Number of ' +
                        'permutations to perform.')
    parser.add_argument('-s', '--first-seed', type=int, default=1, help='Integer ' +
                        'seed passed to the first permutation. Will be incremented ' +
                        'successively for each subsequent permutation.')
    parser.add_argument('-o', '--outfile', help='Path to output tsv file. [default: ' +
                        'stdout]', default='stdout')
    parser.add_argument('-z', '--gzip', action='store_true', help='Compress ' + 
                        'output tsv file with gzip. [Default: do not compress]')
    parser.add_argument('-q', '--quiet', action='store_true', help='Suppress ' +
                        'verbose output.')
    args = parser.parse_args()

    # Open connection to outfile
    if args.outfile in '- stdout'.split():
        outfile = stdout
    else:
        gz_suf = 'gz gzip bgzip bgz bz'.split()
        if splitext(args.outfile)[-1].replace('.' ,'') in gz_suf:
            outfile_path = splitext(args.outfile)[0]
        else:
            outfile_path = args.outfile
        outfile = open(outfile_path, 'w')

    # Load gene blocks and gene universe to prepare for permutation
    blocks = load_blocks(args.blocks)
    univ_bt = pbt.BedTool(args.genes)
    univ_df = pd.read_csv(args.genes, sep='\t').rename(columns = {'#chr' : 'chr'})

    # Shuffle intervals for each of --n-perms permutations
    for seed in range(args.first_seed, args.first_seed + args.n_perms):
        if not args.quiet:
            print('Starting permutation {}'.format(seed))
        shuffled_blocks = shuffle_blocks(blocks, univ_df, univ_bt, seed=seed,
                                         genome=args.genome, quiet=args.quiet)
        if seed == args.first_seed:
            shuffled_blocks.to_csv(outfile, sep='\t', na_rep='NA', header=True, 
                                   index=False, mode='w')
        else:
            shuffled_blocks.to_csv(outfile, sep='\t', na_rep='NA', header=False, 
                                   index=False, mode='a')
    outfile.close()

    # Gzip, if optioned
    if args.gzip:
        subprocess.run(['gzip', '-f', outfile_path])


if __name__ == '__main__':
    main()
