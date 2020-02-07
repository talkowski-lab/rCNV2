#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Perform local defragmentation of a CNV BED file
"""


import pybedtools as pbt
import numpy as np
import argparse
from sys import stdout


def load_cnvs(inbed):
    """
    Prepare CNVs for defragmentation
    """

    cnvs_orig = pbt.BedTool(inbed)

    def _hack_cnv(cnv):
        """
        Apply hack to chromosome field to allow for easy defragmentation
        """
        cnv.chrom = '_'.join([cnv.chrom, cnv.fields[3], cnv.fields[4]])
        return cnv

    cnvs = cnvs_orig.each(_hack_cnv).cut(range(3)).sort().merge()

    cnvs_str = '\n'.join(['\t'.join([x.chrom, str(x.start), str(x.end), 
                                     '_'.join(x.fields[0:3])]) for x in cnvs])

    return pbt.BedTool(cnvs_str, from_string=True)


def defragment_cnvs(cnvs, maxdist=0.25):
    """
    Defragment CNVs based on extending breakpoints by a fraction of CNV size
    """

    def _extend_cnv(cnv, maxdist):
        """
        Extend CNV as a fraction of size
        """
        ext = np.round(len(cnv) * maxdist)
        cnv.start = np.max([0, cnv.start - ext])
        cnv.end = cnv.end + ext
        return cnv

    cnvs_ext = cnvs.each(_extend_cnv, maxdist=maxdist).sort().\
                    merge(c=4, o='distinct', delim='|')

    def _clean_hit(cnv):
        orig = cnv[3].split('|')
        ostarts = [int(x.split('_')[-2]) for x in orig]
        oends = [int(x.split('_')[-1]) for x in orig]

        cnv.start = np.min(ostarts)
        cnv.end = np.max(oends)
        cnv.name = '.'

        return cnv

    return cnvs_ext.each(_clean_hit)


def reformat_cnvs(cnvs):
    """
    Restructure sample & CNV type columns for final output
    """

    def _unhack_cnv(cnv):
        x = cnv.chrom.split('_')
        chrom = x[0]
        cnv_type = x[1]
        sample_id = '_'.join(x[2:])

        newcnv = '\t'.join([chrom, str(cnv.start), str(cnv.end), 
                            cnv_type, sample_id])

        return newcnv

    return pbt.BedTool('\n'.join([_unhack_cnv(cnv) for cnv in cnvs]), 
                       from_string=True)


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inbed', help='BED file of CNVs to defragment. Must ' + 
                        'be BED5, where col4 = cnv type and col5 = sample ID.')
    parser.add_argument('outbed', help='Output BED file for defragmented CNVs.')
    parser.add_argument('--max-dist', dest='maxdist', type=float, default=0.25, 
                        help='Maximum distance to extend each CNV during ' +
                        'defragmentation, specified as a fraction of total CNV ' +
                        'size.')

    args = parser.parse_args()

    # Load & reformat CNVs for defragmentation
    cnvs = load_cnvs(args.inbed)

    # Defragment CNVs
    defragged_cnvs = defragment_cnvs(cnvs, args.maxdist)

    # Reformat defragged CNVs and write to outfile
    reformat_cnvs(defragged_cnvs).saveas(args.outbed)


if __name__ == '__main__':
    main()
