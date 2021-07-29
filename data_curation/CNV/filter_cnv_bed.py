#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Filter an input CNV BED file
"""


import argparse
from sys import stdin, stdout
from os import path
import pybedtools as pbt
import pysam
from math import ceil
from statsmodels.stats.proportion import proportion_confint
import subprocess
import csv


def freq_filter(cnvsA, cnvsB, nsamp, maxFreq=0.01, ro=0.5, dist=50000, 
                use_upper_binom=True):
    """
    Compare two sets of CNVs, and retain those in set A that have fewer than
    a maximum specified number of hits in set B
    """

    # Save a local copy of cnvsA so generator doesn't exhaust itself
    cnvsA = cnvsA.saveas()

    # Intersect
    xbed = cnvsA.intersect(cnvsB, wa=True, wb=True, f=ro, r=True)

    # Throw out self-hits
    def _no_self(feature):
        if len(feature.fields) >= 10:
            if str(feature[3]) != str(feature[9]):
                return feature
    xbed = xbed.filter(_no_self)

    # Throw out CNV type mismatches
    def _cnv_match(feature):
        if len(feature.fields) >= 11:
            if str(feature[4]) == str(feature[10]) or str(feature[10]) in 'MCNV CNV'.split():
                return feature
    xbed = xbed.filter(_cnv_match)

    # Apply breakpoint distance requirement
    def _bp_dist(feature, dist):
        if len(feature.fields) >= 9:
            if abs(int(feature[1]) - int(feature[7])) <= dist \
            and abs(int(feature[2]) - int(feature[8])) <= dist:
                return feature[0:6]
    xbed = xbed.filter(_bp_dist, dist)

    # Get dictionary of hits per ID
    hits = {}
    for feature in xbed:
        if len(feature) >= 4:
            fname = str(feature[3])
            if fname not in hits.keys():
                hits[fname] = 0
            hits[fname] += 1

    # Get list of CNVs to keep based on freq in xbed
    cutoff = ceil(maxFreq * nsamp)
    # If optioned, set a more conservative cutoff equal to the upper bound of the
    # 95% binomial confidence interval corresponding to maxFreq at nsamp samples
    if use_upper_binom:
        try:
            cutoff = ceil(proportion_confint(count=cutoff, nobs=nsamp,  alpha=0.05)[1] * nsamp)
        except:
            import pdb; pdb.set_trace()
    def _remove_fails(feature, hits, cutoff):
        fname = str(feature[3])
        if hits.get(fname, 0) <= cutoff:
            return feature

    return cnvsA.filter(_remove_fails, hits, cutoff)


def read_vcf(vcfin, maxfreq, af_fields='AF', use_lower_binom=True):
    """
    Reads a SV sites vcf and converts it to a BedTool
    """

    vcf = pysam.VariantFile(vcfin)
    header = vcf.header

    af_fields = af_fields.split(',')
    if 'AF' not in af_fields:
        af_fields.append('AF')

    def _filter_vcf(vcf, maxfreq, use_lower_binom=True):
        vcf_str = ''

        for record in vcf:

            if record.info['SVTYPE'] not in 'DEL DUP MCNV CNV'.split():
                continue

            if record.filter.keys() != ['PASS'] \
            and record.filter.keys() != ['MULTIALLELIC']:
                continue

            def _scrape_afs(record, af_fields, use_lower_binom=True):
                afs = []
                for field in af_fields:

                    # Collect raw AF
                    if field not in record.info.keys():
                        continue
                    else:
                        af = record.info[field]

                    if isinstance(af, tuple):
                        if len(af) > 1:
                            af = sum([f for a, f in zip(list(record.alts), list(af)) if a != '<CN=2>'])
                        else:
                            af = af[0]

                    # If optioned, compute lower bound of 95% binomial confidence 
                    # interval and report that as AF (more conservative / handles
                    # small sample sizes better)
                    if use_lower_binom:
                        if 'CN_NONREF_FREQ' in field:
                            n_samps = record.info[field.replace('NONREF_FREQ', 'NUMBER')]
                        else:
                            if field == 'AF':
                                n_samps = record.info['AN']
                            else:
                                n_samps = record.info[field.replace('_AF', '_AN')]
                        af = proportion_confint(count=ceil(af * n_samps), nobs=n_samps,  alpha=0.05)[0]

                    afs.append(af)
                return afs

            afs = _scrape_afs(record, af_fields, use_lower_binom)
            if len(afs) > 0:
                if max(afs) <= maxfreq:
                    continue

            flist = [str(record.chrom), str(record.pos), str(record.stop), 
                     str(record.id), str(record.info['SVTYPE']).split('_')[0],
                     'HEALTHY_CONTROL']

            fstr = ' '.join(flist) + '\n'

            vcf_str = vcf_str + fstr

        return pbt.BedTool(vcf_str, from_string=True)

    gbed = _filter_vcf(vcf, maxfreq, use_lower_binom)

    return gbed


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inbed', help='Input BED (supports "stdin").')
    parser.add_argument('outbed', help='Output BED (supports "stdout").')
    parser.add_argument('--chr', help='Restrict to a subset of chromosomes. ' +
                        'Specify as comma-separated list.')
    parser.add_argument('--minsize', type=int, help='Minimum CNV size. ' +
                        '[50kb]', default=50000)
    parser.add_argument('--maxsize', type=int, help='Maximum CNV size. ' +
                        '[10Mb]', default=10000000)
    parser.add_argument('--maxfreq', type=float, help='Maximum CNV ' + 
                        'frequency. [default: 0.01]', default=0.01)
    parser.add_argument('--recipoverlap', type=float, help='Reciprocal overlap ' + 
                        'cutoff after which to consider two CNVs matching. ' +
                        '[default: 0.5]',
                        default=0.5)
    parser.add_argument('--dist', type=int, help='Maximum distance permitted ' + 
                        'between breakpoints to consider two CNVs matching. ' +
                        '[default: 50kb]',
                        default=50000)
    parser.add_argument('-w', '--whitelist', action='append', help='BED file ' +
                        'containing regions to exclude from blacklisting and ' +
                        'cross-cohort frequency-based filtering. May be specified ' +
                        'multiple times.')
    parser.add_argument('--wrecip', type=float, help='Minimum reciprocal overlap ' + 
                        'vs. any whitelist interval to automatically include a CNV. ' +
                        '[default: 0.75]', default=0.75)
    parser.add_argument('-x', '--blacklist', action='append', help='BED file ' +
                        'containing regions to blacklist based on CNV coverage. ' +
                        'May be specified multiple times.')
    parser.add_argument('--xcov', type=float, help='Maximum coverage ' + 
                        'by any blacklist before excluding a CNV. [default: 0.3]',
                        default=0.3)
    parser.add_argument('--cohorts-list', help='list of all cohorts to be used ' + 
                        'for freq filtering. Will compare each cohort to input ' + 
                        'bed, one at a time. Requires three tab delimited ' +
                        'columns: name, sample size, and path to BED.')
    parser.add_argument('--allcohorts', help='BED file of all unfiltered ' + 
                        'CNVs across all cohorts, to be used for freq filtering.')
    parser.add_argument('--allcohorts_nsamp', type=int, help='Number of samples ' +
                        'combined across all cohorts. Only used if --allcohorts ' +
                        'is specified.')
    parser.add_argument('-v', '--vcf', action='append', help='SV VCF to use for ' + 
                        'frequency filtering. May be specified multiple times.')
    parser.add_argument('--vcf-af-fields', help='Entry or entries in INFO field ' +
                        'of SV VCF(s) to use as frequency. Separate multiple ' +
                        'entries with commas. If multiple entries are specified, ' +
                        'will take the max. Will default to AF if entry not found.',
                         default='AF', dest='af_fields')
    parser.add_argument('-g', '--genome', help='Genome file (used for sorting).')
    parser.add_argument('-z', '--bgzip', dest='bgzip', action='store_true',
                        help='Compress output BED with bgzip.')

    args = parser.parse_args()

    # Read input bed
    if args.inbed in 'stdin -'.split():
        cnvs = pbt.BedTool(stdin)
    else:
        cnvs = pbt.BedTool(args.inbed)

    # Restrict to a subset of chromosomes, or autosomes if not specified
    if args.chr is None:
        chroms = [i for subl in ['{0} chr{1}'.format(c, c).split() for c in range(1, 23)] for i in subl]
    else:
        chroms = args.chr.split(',')
    cnvs = cnvs.filter(lambda x: x.chrom in chroms)

    # Restrict on size
    cnvs = cnvs.filter(lambda x: len(x) >= args.minsize and len(x) <= args.maxsize)

    # Restrict on VCFs 
    if args.vcf is not None:
        for vcfpath in args.vcf:
            vbed = read_vcf(vcfpath, args.maxfreq, args.af_fields)
            cnvs = freq_filter(cnvs, vbed, 0, maxFreq=args.maxfreq, 
                               ro=args.recipoverlap, dist=args.dist,
                               use_upper_binom=False)

    # Extract CNVs overlapping whitelist and hold them out from all subsequent filtering
    if args.whitelist is not None:
        wbt = pbt.BedTool(args.whitelist[0]).cut(range(3))
        if len(args.whitelist) > 1:
            for inbed in args.whitelist[1:]:
                wbt = wbt.cat(pbt.BedTool(inbed).cut(range(3)), postmerge=False)
        cnvs = cnvs.saveas()
        wcnvs = cnvs.intersect(wbt, r=True, f=args.wrecip, wa=True, u=True).saveas()
        cnvs = cnvs.intersect(wbt, r=True, f=args.wrecip, wa=True, v=True)

    # Restrict on blacklist coverage
    if args.blacklist is not None:
        def _bl_filter(feature, xcov):
            if float(feature[-1]) <= xcov:
                return feature
        for blpath in args.blacklist:
            bl = pbt.BedTool(blpath)
            cnvs = cnvs.coverage(bl).filter(_bl_filter, args.xcov).cut(range(0, 6))

    # Restrict on frequency filtering vs each cohort, one at a time
    if args.cohorts_list is not None:
        with open(args.cohorts_list) as fin:
            for name, N, bcnv_path in csv.reader(fin, delimiter='\t'):
                cnvs = freq_filter(cnvs, bcnv_path, int(N),
                           maxFreq=args.maxfreq, ro=args.recipoverlap, 
                           dist=args.dist)

    # Restrict on global intersection
    if args.allcohorts is not None \
    and args.allcohorts_nsamp is not None:
        cnvs = freq_filter(cnvs, args.allcohorts, args.allcohorts_nsamp,
                           maxFreq=args.maxfreq, ro=args.recipoverlap, 
                           dist=args.dist)

    # Add whitelisted CNVs back into final CNV file
    cnvs = cnvs.saveas()
    if args.genome is not None:
        cnvs = cnvs.cat(wcnvs, postmerge=False).sort(g=args.genome)
    else:
        cnvs = cnvs.cat(wcnvs, postmerge=False).sort()

    # Write filtered CNVs out to file
    outheader = '#chr\tstart\tend\tname\tcnv\tpheno'
    if args.outbed in 'stdout -'.split():
        cnvs.saveas(stdout, trackline=outheader)
    else:
        outbed = args.outbed
        if '.gz' in outbed:
            outbed = path.splitext(outbed)[0]
        cnvs.saveas(outbed, trackline=outheader)
        if args.bgzip:
            subprocess.run(['bgzip', '-f', outbed])

if __name__ == '__main__':
    main()
