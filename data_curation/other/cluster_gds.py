#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Cluster genomic disorders across sources
"""


cnvtypes = 'DEL DUP'.split()
elig_chroms = [str(contig) for contig in range(1, 23)]


import csv
import pybedtools as pbt
import argparse
from os import path
from athena.utils import bgzip


def load_sources(tsv_in, cnvtype):
    """
    Load a dict of pbt.BedTools (one per source)
    """

    sources = {}

    with open(tsv_in) as fin:
        reader = csv.reader(fin, delimiter='\t')
        for name, path in reader:
            bed = pbt.BedTool(path).\
                      filter(lambda x: x[3] == cnvtype and x.chrom in elig_chroms).\
                      sort().cut(range(3)).merge().saveas()
            sources[name] = bed
    
    return sources


def density_map(sources, resolution=10000):
    """
    Make map of GD density for all input sources for a single CNV type
    """

    pooled = list(sources.values())[0].cat(*list(sources.values())[0:], 
                                             postmerge=False).sort().saveas()
    counts = pbt.BedTool().makewindows(b=pooled.merge(), w=resolution).\
                 coverage(pooled, counts=True, sorted=True)
    count_bts = []
    for k in range(1, len(sources) + 1):
        hits = counts.filter(lambda x: x[3] == str(k)).saveas()
        if len(hits) > 0:
            count_bts.append(hits.merge(c=4, o='distinct'))

    return count_bts[0].cat(*count_bts[0:], postmerge=False).sort().saveas()


def get_cytobands(bt, cyto_bed):
    """
    Return cytoband nomenclature for all cytobands overlapped by bt
    Note: assumes all entries in bt are on the same chromosome
    """

    chrom = bt[0].chrom
    bands = sorted(list(set([x[-2] for x in bt.intersect(cyto_bed, wb=True)])))
    if len(bands) > 1:
        bandrange = '{}{}-{}'.format(chrom, bands[0], bands[-1])
    else:
        bandrange = chrom + bands[0]

    return bandrange


def finalize_regions(maps, cnv, min_score=0, max_score=10e10, min_size=200000, 
                     max_size=10000000, gfile=None, segdups=None, cyto_bed=None,
                     id_prefix='region', annoprep=False):
    """
    Merge final regions based on count and trim flanking segdups
    """

    regions = maps[cnv].filter(lambda x: int(x[3]) >= min_score and int(x[3]) < max_score).\
                        cut(range(3)).merge()
    if gfile is not None:
        regions = regions.sort(g=gfile).saveas()

    # Trim segdups overlapping ends of each GD, and filter by size
    left_strs = '\n'.join(['\t'.join([x.chrom, str(x.start), str(x.start + 1)]) for x in regions])
    right_strs = '\n'.join(['\t'.join([x.chrom, str(x.end), str(x.end + 1)]) for x in regions])
    ends = pbt.BedTool(left_strs + '\n' + right_strs, from_string=True)
    sd_hits = pbt.BedTool(segdups).intersect(ends, u=True, wa=True)
    regions = regions.subtract(sd_hits).filter(lambda x: len(x) >= min_size and len(x) <= max_size)

    # Add cytoband range as region name and CNV type as fifth column
    all_str = ''
    k = 1
    for gd in regions:
        gd_coord_str = '\t'.join([gd.chrom, str(gd.start), str(gd.end)])
        gd_bt = pbt.BedTool(gd_coord_str, from_string=True)
        if cyto_bed is not None:
            bandrange = get_cytobands(gd_bt, cyto_bed)
            gd_id = '_'.join([id_prefix, str(k), cnv, bandrange])
        else:
            gd_id = '_'.join([id_prefix, cnv, str(k)])
        gd_str = '\t'.join([gd_coord_str, gd_id, cnv])
        # Add extra two columns for gene annotation, if optioned
        if annoprep:
            gd_str += '\t{}:{}-{}\t.'.format(gd.chrom, gd.start, gd.end)
        all_str += gd_str + '\n'
        k += 1
        
    return pbt.BedTool(all_str, from_string=True)


def cnv_filter(gd_bt, cnv_path, cnv, max_cov):
    """
    Exclude GDs based on overlap by common CNVs (while matching on CNV type)
    """

    gd_wcov = gd_bt.coverage(cnv_path)

    gd_bt_filt = gd_wcov.filter(lambda x: x[4] != cnv or float(x[-1]) <= max_cov).\
                         cut(range(len(gd_bt[0].fields))).saveas()

    return gd_bt_filt


def main():
    """
    Main block
    """

    # Parse command-line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sources', help='tsv of input GD BED files. Expects two ' +
                        'columns: source name and path to BED.')
    parser.add_argument('--hc-outfile', required=True, help='Output BED file for ' +
                        'high-confidence GDs. Required.')
    parser.add_argument('--mc-outfile', required=True, help='Output BED file for ' +
                        'medium-confidence GDs. Required.')
    parser.add_argument('--lc-outfile', required=True, help='Output BED file for ' +
                        'low-confidence GDs. Required.')
    parser.add_argument('--hc-cutoff', type=int, default=3,
                        help='Minimum number of sources required for ' +
                        'high-confidence GDs')
    parser.add_argument('--mc-cutoff', type=int, default=2,
                        help='Minimum number of sources required for ' +
                        'medium-confidence GDs')
    parser.add_argument('--lc-cutoff', type=int, default=1,
                        help='Minimum number of sources required for ' +
                        'low-confidence GDs')
    parser.add_argument('-g', '--genome', help='Genome file (for sorting outputs).')
    parser.add_argument('--segdups', help='BED file of segdups (for trimming ends of GDs).')
    parser.add_argument('--cytobands', help='BED file of cytobands (for assigning GD names).')
    parser.add_argument('--common-dels', help='BED file of common deletions ' +
                        '(for filtering weak-effect GDs).')
    parser.add_argument('--common-dups', help='BED file of common duplications ' +
                        '(for filtering weak-effect GDs).')
    parser.add_argument('--common-cnv-cov', type=float, default=0.3, help='Minimum ' +
                        'coverage by common CNVs before discarding GDs.')
    parser.add_argument('--minsize', type=int, default=200000, help='Minimum GD ' +
                        'size to report.')
    parser.add_argument('--maxsize', type=int, default=10000000, help='Maximum GD ' +
                        'size to report.')
    parser.add_argument('--prep-for-gene-anno', action='store_true', help='Format ' + 
                        'output BED files with extra columns for gene annotation.')
    parser.add_argument('-z', '--bgzip', action='store_true', help='Compress ' + 
                        'output BED files with bgzip. [Default: do not compress]')
    args = parser.parse_args()

    # Open paths to output files
    gz_suf = 'gz gzip bgzip bgz bz'.split()
    if path.splitext(args.hc_outfile)[-1].replace('.' ,'') in gz_suf:
        hc_outfile = path.splitext(args.hc_outfile)[0]
    else:
        hc_outfile = args.hc_outfile
    if path.splitext(args.mc_outfile)[-1].replace('.' ,'') in gz_suf:
        mc_outfile = path.splitext(args.mc_outfile)[0]
    else:
        mc_outfile = args.mc_outfile
    if path.splitext(args.lc_outfile)[-1].replace('.' ,'') in gz_suf:
        lc_outfile = path.splitext(args.lc_outfile)[0]
    else:
        lc_outfile = args.lc_outfile

    # Read list of sources
    sources = {cnv : load_sources(args.sources, cnv) for cnv in cnvtypes}

    # Build map of GD density per CNV type
    maps = {cnv : density_map(sources[cnv]) for cnv in cnvtypes}

    # Consolidate hc/mc/lc regions, trim flanking segdups, and write to outfiles
    hc_gd_bts = [finalize_regions(maps, cnv, min_score=args.hc_cutoff,
                                  min_size=args.minsize, max_size=args.maxsize,
                                  gfile=args.genome, segdups=args.segdups, 
                                  cyto_bed=args.cytobands, id_prefix='HC_GD',
                                  annoprep=args.prep_for_gene_anno) \
                 for cnv in cnvtypes]
    hc_gds = hc_gd_bts[0].cat(hc_gd_bts[1], postmerge=False).sort(g=args.genome)
    mc_gd_bts = [finalize_regions(maps, cnv, min_score=args.mc_cutoff,
                                  max_score=args.hc_cutoff, min_size=args.minsize, 
                                  max_size=args.maxsize, gfile=args.genome, 
                                  segdups=args.segdups, cyto_bed=args.cytobands,
                                  id_prefix='MC_GD', annoprep=args.prep_for_gene_anno) \
                 for cnv in cnvtypes]
    mc_gds = mc_gd_bts[0].cat(mc_gd_bts[1], postmerge=False).sort(g=args.genome)
    lc_gd_bts = [finalize_regions(maps, cnv, min_score=args.lc_cutoff,
                                  max_score=args.mc_cutoff, min_size=args.minsize, 
                                  max_size=args.maxsize, gfile=args.genome, 
                                  segdups=args.segdups, cyto_bed=args.cytobands,
                                  id_prefix='LC_GD', annoprep=args.prep_for_gene_anno) \
                 for cnv in cnvtypes]
    lc_gds = lc_gd_bts[0].cat(lc_gd_bts[1], postmerge=False).sort(g=args.genome)

    # Filter vs common CNVs, if optioned
    if args.common_dels is not None:
        hc_gds = cnv_filter(hc_gds, args.common_dels, "DEL", args.common_cnv_cov)
        mc_gds = cnv_filter(mc_gds, args.common_dels, "DEL", args.common_cnv_cov)
        lc_gds = cnv_filter(lc_gds, args.common_dels, "DEL", args.common_cnv_cov)
    if args.common_dups is not None:
        hc_gds = cnv_filter(hc_gds, args.common_dups, "DUP", args.common_cnv_cov)
        mc_gds = cnv_filter(mc_gds, args.common_dups, "DUP", args.common_cnv_cov)
        lc_gds = cnv_filter(lc_gds, args.common_dups, "DUP", args.common_cnv_cov)

    # Write to outfiles
    out_header = '\t'.join('#chr start end gd_id cnv'.split())
    if args.prep_for_gene_anno:
        out_header += '\tcoords\tdummy'
    hc_gds.saveas(hc_outfile, trackline=out_header + '\n')
    mc_gds.saveas(mc_outfile, trackline=out_header + '\n')
    lc_gds.saveas(lc_outfile, trackline=out_header + '\n')

    # Bgzip, if optioned
    if args.bgzip:
        bgzip(hc_outfile)
        bgzip(mc_outfile)
        bgzip(lc_outfile)


if __name__ == '__main__':
    main()

