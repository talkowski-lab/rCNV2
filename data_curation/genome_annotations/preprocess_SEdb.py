#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Curate SEdb super enhancer and typical enhancer BED files
"""


import csv
from urllib.request import urlopen
import io
import argparse
import subprocess


def load_manifest(manifest, outdir, enh_type):
    """
    Load SEdb sample manifest and open output file handles 
    """

    mfst = {}

    bed_header = '#chrom\tstart\tend\n'

    with open(manifest) as fin:
        for sid, source, biosample, tissue, name in csv.reader(fin, delimiter='\t'):
            name = name.replace('/', '_').replace(' ', '_').replace('(', '').replace(')', '')
            source = source.replace('/', '_').replace(' ', '_').replace('(', '').replace(')', '')
            outname = '{}/SEdb_{}.{}.{}.bed'.format(outdir, enh_type, name, source)
            ofile = open(outname, 'w')
            ofile.write(bed_header)
            mfst[sid] = ofile

    return mfst


def process_stream(url, enh_type, manifest):
    """
    Stream & process SEdb file
    """

    # Stream & parse file
    with urlopen(url) as fin:
        for line in fin:
            fields = line.decode().rstrip().split('\t')
            
            if enh_type == 'SE':
                sid = fields[0]
                if sid in manifest.keys():
                    manifest[sid].write('\t'.join(fields[2:5]) + '\n')

            else:
                sid = fields[3]
                if sid in manifest.keys():
                    manifest[sid].write('\t'.join(fields[0:3]) + '\n')


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('url', help='URL to SEdb file.')
    parser.add_argument('manifest', help='SEdb sample manifest.')
    parser.add_argument('--format', default='SE', choices=('SE', 'TE'), required=True,
    					help='Specify input data type [SE or TE].')
    parser.add_argument('-o', '--outdir', default='./', help='Path to output directory.')
    args = parser.parse_args()

    # Load sample manifest
    mfst = load_manifest(args.manifest, args.outdir, args.format)

    # Stream & process data for all samples
    process_stream(args.url, args.format, mfst)

    # Close and bgzip all outfiles
    for file in mfst.values():
        opath = file.name
        file.close()
        subprocess.run(['bgzip', '-f', opath])


if __name__ == '__main__':
    main()

