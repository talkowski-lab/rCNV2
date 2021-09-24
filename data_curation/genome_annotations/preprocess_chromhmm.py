#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Preprocess Roadmap Epigenomics ChromHMM data and write out as individual BED files
"""


import csv
import pandas as pd
import pybedtools as pbt
import subprocess
from pathlib import Path
import argparse


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


def process_chromhmm_bed(path, states, prefix):
    """
    Splits a ChromHMM BED by state and writes to outfile
    Returns a dict of {state : pbt.BedTool}
    """

    # Load BED as pd.DataFrame
    cdf = pd.read_csv(path, sep='\t', names='chrom start end state'.split())
    cdf['chrom'] = cdf.chrom.apply(lambda x: x.replace('chr', ''))

    # Split by state & convert to pbt.BedTool
    for state in states:
        outpath = '.'.join([prefix, state.replace('/', ''), 'bed'])
        print('Writing state {} to {}.gz...'.format(state, outpath))
        sbt = pbt.BedTool.from_dataframe(cdf.loc[cdf.state == state, :].iloc[:, 0:3]).\
                          sort().merge().saveas(outpath)
        subprocess.run(['bgzip', '-f', outpath])


def process_chromhmm_beds(chromdir, samples, states, prefix, bed_suffix):
    """
    Loads all ChromHMM state BEDs able to be found in chromdir
    Splits each BED by state
    Returns a dict of {eid : {state : pbt.BedTool}}
    """

    for eid in samples.keys():

        # Check for ChromHMM BED, and load if found
        bedpath = chromdir + '/' + eid + bed_suffix
        if Path(bedpath).exists():
            print('Processing ChromHMM BED for {} ({})'.\
                  format(samples[eid]['name'], bedpath))
            out_prefix = '.'.join([prefix, eid])
            process_chromhmm_bed(bedpath, states, out_prefix)
        else:
            print('Warning: unable to locate ChromHMM BED for {} ({})'.\
                  format(samples[eid]['name'], bedpath))
            continue


def main():
    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('chromdir', help='directory containing ChromHMM BEDs')
    parser.add_argument('--state-manifest', required=True, help='ChromHMM state manifest tsv')
    parser.add_argument('--sample-manifest', required=True, help='REP sample manifest tsv')
    parser.add_argument('-p', '--prefix', default='REP_summary', help='Prefix ' +
                        'for output BED files.')
    parser.add_argument('--bed-suffix', help='Suffix for all ChromHMM BEDs.',
                        default='_18_core_K27ac_mnemonics.bed.gz')
    args = parser.parse_args()

    # Load Roadmap manifests
    states = load_state_manifest(args.state_manifest)
    samples = load_sample_manifest(args.sample_manifest)

    # Load Roadmap ChromHMM BEDs, split by state, and write to outdir
    process_chromhmm_beds(args.chromdir, samples, states, args.prefix, args.bed_suffix)


if __name__ == '__main__':
    main()
