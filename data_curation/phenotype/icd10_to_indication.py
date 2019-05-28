#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Convert sample + ICD10 pairs into sample + string indication pairs
"""


import csv
from re import sub
import argparse
from sys import stdout


def read_icd_dict(icd10map):
    """
    Read UKBB official ICD-10 map into dict
    """

    icd_dict = {}
    junk_words = ['OTHER', 'UNSPECIFIED', 'UNKNOWN']
    split_words = ['AND', 'WITH']
    truncate_words = ['WITHOUT']
    junk_phrases = ['PERSONAL_HISTORY_OF_', 'SITE_NOT_SPECIFIED', 'NOT_SPECIFIED', \
                    'NOT_ELSEWHERE_CLASSIFIED']

    with open(icd10map) as infile:

        reader = csv.reader(infile, delimiter='\t')

        for icd, descrip, node_idx, parent_idx, selectable in reader:

            if icd == 'coding' \
            or icd.startswith('Chapter'):
                continue

            # Clean description
            descrip = descrip.upper().split(' ', 1)[1]
            descrip = sub('[^A-Za-z0-9\ ]+', '', descrip)
            descrip = ' '.join([w for w in descrip.split(' ') if w not in junk_words])
            for word in split_words:
                descrip = descrip.replace(' ' + word + ' ', ';')
            for word in truncate_words:
                descrip = descrip.split(word)[0].rstrip()
            descrip = sub('[\ ]+', '_', descrip)
            for phrase in junk_phrases:
                descrip = descrip.replace(phrase, '')
            
            if icd not in icd_dict.keys():
                icd_dict[icd] = descrip
            else:
                icd_dict = ';'.join(icd_dict[icd], descrip)

    return icd_dict


def process_samples(samples, icd_dict, outfile, default='NA', report=False):
    """
    Iterate over a file of sample + ICD-10 pairs and convert to indications
    """

    skip_keywords = ['FAMILY_HISTORY_OF']

    with open(samples) as infile:

        reader = csv.reader(infile, delimiter='\t')

        for sample, icds in reader:

            # Parse ICD-10 codes
            if icds == '':
                pheno = default
            else:
                icds = icds.split(';')
                pheno = ''
                for icd in icds:

                    descrip = icd_dict.get(icd, None)

                    if icd is None:
                        if report:
                            err = 'Failed to find match for {0} for sample {1}'
                            print(err.format(icd, sample))
                    else:
                        if any(keyword in descrip for keyword in skip_keywords):
                            continue
                        else:
                            if pheno == '':
                                pheno = descrip
                            else:
                                pheno = ';'.join([pheno, icd_dict.get(icd, '')])

            if pheno == '':
                pheno = default

            outfile.write('{0}\t{1}\n'.format(sample, pheno))


def main():
    """
    Main block
    """

    # Parse arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('icd_10_map', help='Mapping of ICD-10 ' + 
                        'codes to descriptions. One line per ICD-10. Two ' + \
                        'tab-delimited columns: ICD-10 and description.', 
                        metavar='file')
    parser.add_argument('samples', help='Test file of samples & ICD10 codes. ' + 
                        'One line per sample. Two tab-delimited columns: ' + 
                        'sample ID and string of semicolon-delimited ICD-10s.', 
                        metavar='file')
    parser.add_argument('-d', '--default', help='String to use if no ICD-10 codes ' +
                        'match. [default: HEALTHY_CONTROL]', 
                        default='HEALTHY_CONTROL', metavar='string')
    parser.add_argument('--report-fails', action='store_true', 
                        help='Print warnings for each failed ICD-10 mapping. ' +
                        '[default: silent]')
    parser.add_argument('-o', '--outfile', help='Path to outfile. ' +
                        '[default: stdout]', metavar='file')
    args = parser.parse_args()

    # Read ICD-10 terms
    icd_dict = read_icd_dict(args.icd_10_map)
    
    # Open connection to outfile
    if args.outfile is None:
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Process sample + ICD-10 pairings
    process_samples(args.samples, icd_dict, outfile, args.default,
                    args.report_fails)
    

if __name__ == '__main__':
    main()
