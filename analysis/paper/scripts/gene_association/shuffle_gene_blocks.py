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
from athena.utils import bgzip


def load_blocks(blocks_in):
    """
    Load & reformat original set of gene blocks (to be shuffled)
    Returns: dict of block_id : {block attributes} pairs
    """

    blocks = pd.read_csv(blocks_in, sep='\t').iloc[:, :3].astype(str)
    blocks.columns = 'block_id cnv genes'.split()
    blocks.index = blocks['block_id']
    blocks.drop(columns='block_id', inplace=True)
    blocks_dict = blocks.to_dict(orient='index')
    for bid in blocks_dict.keys():
        genes = [g for g in blocks_dict[bid]['genes'].split(';') \
                 if g is not None and g not in 'NA nan NAN'.split() and g != '']
        blocks_dict[bid]['genes'] = genes
        blocks_dict[bid]['n_genes'] = len(genes)

    return blocks_dict


def shuffle_blocks(blocks, univ_df, univ_bt, seed):
    """
    Randomly shuffle a dict of colinear gene blocks (output by load_blocks)
    Return: pd.DataFrame of shuffled blocks (matches args.blocks format)
    """

    # Defines 2 x len(blocks) candidate starting points for seed
    index_genes = univ_df.sample(n=2*len(blocks), random_state=seed)
    used_del_genes = []
    used_dup_genes = []

    for bid, bvals in blocks.items():
        ngenes = bvals['n_genes']
        # Select new index gene
        blocks.sample

    # def _seed_prefix(x, seed):
    #     x.name = '_'.join(['perm' + str(seed), x.name])
    #     return x

    import pdb; pdb.set_trace()


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('blocks', help='TSV of gene blocks to be shuffled. Requires ' +
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
    parser.add_argument('-o', '--outfile', help='Path to output BED file. [default: ' +
                        'stdout]', default='stdout')
    parser.add_argument('-z', '--bgzip', action='store_true', help='Compress ' + 
                        'output BED files with bgzip. [Default: do not compress]')
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
        shuffled_blocks = shuffle_blocks(blocks, univ_df, univ_bt, seed=seed)
    #     shuffled_bt = custom_shuffle(segs, seed, args.genome, args.whitelist, 
    #                                  coords_colname=coords_colname)
    #     if seed == args.first_seed:
    #         shuffled_bt.to_csv(outfile, sep='\t', na_rep='NA', header=True, 
    #                            index=False, mode='w')
    #     else:
    #         shuffled_bt.to_csv(outfile, sep='\t', na_rep='NA', header=False, 
    #                            index=False, mode='a')
    # outfile.close()

    # # Bgzip, if optioned
    # if args.bgzip:
    #     bgzip(outfile_path)


if __name__ == '__main__':
    main()
