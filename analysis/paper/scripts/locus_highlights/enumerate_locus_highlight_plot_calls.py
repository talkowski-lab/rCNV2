#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# and the Talkowski Laboratory
# Distributed under terms of the MIT license.

"""
Format calls to plot_locus_highlight.R for all
"""


import pandas as pd
import numpy as np
import subprocess
import argparse
from sys import stdout


get_fdr_script = '/opt/rCNV2/analysis/other/estimate_p_for_fdr.R'


def format_segment_calls(segs_in, phenos, outfile, target_dir='./', 
                         bonf_cutoff=10e-6, wex=2/3):
    """
    Generate one call to plot_locus_highlight.R for each row in segs_in 
    and write to outfile
    """

    ss_template = 'meta_stats/{}.rCNV.{}.sliding_window.meta_analysis.stats.bed.gz'

    seg_df = pd.read_csv(segs_in, sep='\t').rename(columns={'#chr' : 'chrom'})

    for sdat in seg_df.itertuples():
        sdat = sdat._asdict()

        # Get basic segment info
        chrom = sdat.get('chrom')
        cnv = sdat.get('cnv')
        o_start = sdat.get('start_min')
        o_end = sdat.get('end_max')
        size = o_end - o_start
        rid = sdat.get('region_id')
        highlights = sdat.get('cred_interval_coords')

        # Scale view window according to overall size
        start = int(np.max([0, o_start - np.floor(wex * size)]))
        end = int(o_end + np.floor(wex * size))
        region = '{}:{}-{}'.format(chrom, start, end)

        # Get associated HPO with largest sample size & its direct parent term
        all_hpos = sdat.get('hpos').split(';')
        largest_hpo = phenos.HPO[phenos.HPO.isin(all_hpos)].values[0]
        hpo_tier = phenos.HPO_tier[phenos.HPO == largest_hpo].values[0]
        if hpo_tier > 1:
            is_parent = phenos.child_terms.astype(str).\
                               apply(lambda x: largest_hpo in x.split(';'))
            main_hpo = None
            k = 0
            while main_hpo is None:
                k += 1
                try:
                    main_hpo = phenos.HPO[(phenos.HPO_tier == hpo_tier - k) & (is_parent)].values[0]
                except:
                    main_hpo = None
                if hpo_tier - k == 1:
                    main_hpo = 'HP:0000118'
            highlight_hpo = largest_hpo
        else:
            main_hpo = largest_hpo
            if len(all_hpos) > 1:
                highlight_hpo = phenos.HPO[phenos.HPO.isin(all_hpos)].values[1]
            else:
                highlight_hpo = None
        if highlight_hpo is not None:
            sumstat_hpo = highlight_hpo
        else:
            sumstat_hpo = main_hpo
        ss_path = ss_template.format(sumstat_hpo.replace(':', ''), cnv)

        # Check FDR-equivalent P-value if necessary
        sig = sdat.get('best_sig_level')
        if sig == 'genome_wide':
            sig_label = 'Genome-wide significance'
            p_cutoff = bonf_cutoff
        else:
            sig_label = 'FDR < 1%'
            fdr_res = subprocess.run([get_fdr_script, ss_path], 
                                      stdout=subprocess.PIPE, 
                                      stderr=subprocess.DEVNULL)
            p_cutoff = float(fdr_res.stdout.decode().rstrip())

        # Write script call to outfile
        outfile.write('echo "Plotting {}..."\n'.format(rid))
        outline = '/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_locus_highlight.R'
        outline += ' --case-hpos "{}"'.format(main_hpo)
        if highlight_hpo is not None:
            outline += ' --highlight-hpo "{}"'.format(highlight_hpo)
        outline += ' --highlights "{}"'.format(highlights)
        outline += ' --sumstats {}'.format(ss_path)
        outline += ' --cytobands refs/GRCh37.cytobands.bed.gz'
        outline += ' --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz'
        outline += ' --gw-sig "{}"'.format(p_cutoff)
        outline += ' --gw-sig-label "{}"'.format(sig_label)
        outline += ' --subgroup-cohorts'
        outline += ' --pdf-height 4'
        outline += ' --cnv-panel-height 1.1'
        outline += ' --cnv-panel-space 0.2'
        outline += ' --standardize-frequencies {} cnvs.input.tsv {}'.format(region, cnv)
        outline += ' refs/HPOs_by_metacohort.table.tsv refs/GRCh37.genome'
        outline += ' {}/{} 1> /dev/null 2>&1\n'.format(target_dir.strip('/'), rid)
        outfile.write(outline)


def format_credset_calls(credsets_in, phenos, outfile, target_dir='./', 
                         bonf_cutoff=10e-6, wex=2/3):
    """
    Generate one call to plot_locus_highlight.R for each row in credsets_in 
    and write to outfile
    """

    ss_template = 'meta_stats/{}.rCNV.{}.gene_burden.meta_analysis.stats.bed.gz'
    pip_template = 'rCNV.{}.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv'

    cs_df = pd.read_csv(credsets_in, sep='\t').rename(columns={'#chr' : 'chrom'})

    for csdat in cs_df.itertuples():
        csdat = csdat._asdict()

        # Get basic credset info
        chrom = csdat.get('chrom')
        cnv = csdat.get('cnv')
        o_start = csdat.get('start')
        o_end = csdat.get('end')
        size = o_end - o_start
        csid = csdat.get('credible_set_id')
        pip_path = pip_template.format(cnv)

        # Scale view window according to overall size
        start = int(np.max([0, o_start - np.floor(wex * size)]))
        end = int(o_end + np.floor(wex * size))
        region = '{}:{}-{}'.format(chrom, start, end)
        highlights = '{}:{}-{}'.format(chrom, o_start, o_end)

        # Get gene information
        cs_genes = csdat.get('all_genes')
        if str(csdat.get('vconf_genes')) != 'nan':
            vconf_genes = csdat.get('vconf_genes').split(';')
        else:
            vconf_genes = []
        if str(csdat.get('conf_genes')) != 'nan':
            conf_genes = csdat.get('conf_genes').split(';')
        else:
            conf_genes = []
        label_genes = list(set(vconf_genes).union(set(conf_genes)))

        # Get associated HPO with largest sample size & its direct parent term
        all_hpos = csdat.get('hpos').split(';')
        largest_hpo = phenos.HPO[phenos.HPO.isin(all_hpos)].values[0]
        hpo_tier = phenos.HPO_tier[phenos.HPO == largest_hpo].values[0]
        if hpo_tier > 1:
            is_parent = phenos.child_terms.astype(str).\
                               apply(lambda x: largest_hpo in x.split(';'))
            main_hpo = None
            k = 0
            while main_hpo is None:
                k += 1
                try:
                    main_hpo = phenos.HPO[(phenos.HPO_tier == hpo_tier - k) & (is_parent)].values[0]
                except:
                    main_hpo = None
                if hpo_tier - k == 1:
                    main_hpo = 'HP:0000118'
            highlight_hpo = largest_hpo
        else:
            main_hpo = largest_hpo
            if len(all_hpos) > 1:
                highlight_hpo = phenos.HPO[phenos.HPO.isin(all_hpos)].values[1]
            else:
                highlight_hpo = None
        if highlight_hpo is not None:
            sumstat_hpo = highlight_hpo
        else:
            sumstat_hpo = main_hpo
        ss_path = ss_template.format(sumstat_hpo.replace(':', ''), cnv)

        # Check FDR-equivalent P-value if necessary
        sig = csdat.get('best_sig_level')
        if sig == 'exome_wide':
            sig_label = 'Exome-wide significance'
            p_cutoff = bonf_cutoff
        else:
            sig_label = 'FDR < 1%'
            fdr_res = subprocess.run([get_fdr_script, ss_path], 
                                      stdout=subprocess.PIPE, 
                                      stderr=subprocess.DEVNULL)
            p_cutoff = float(fdr_res.stdout.decode().rstrip())

        # Write script call to outfile
        outfile.write('echo "Plotting {}..."\n'.format(csid))
        outline = '/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_locus_highlight.R'
        outline += ' --case-hpos "{}"'.format(main_hpo)
        if highlight_hpo is not None:
            outline += ' --highlight-hpo "{}"'.format(highlight_hpo)
            sumstat_hpo = highlight_hpo
        else:
            sumstat_hpo = main_hpo
        outline += ' --highlights "{}"'.format(highlights)
        outline += ' --sumstats ' + ss_template.format(sumstat_hpo.replace(':', ''), cnv)
        outline += ' --cytobands refs/GRCh37.cytobands.bed.gz'
        outline += ' --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz'
        outline += ' --label-genes "{}"'.format(';'.join(label_genes))
        outline += ' --constraint refs/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz'
        outline += ' --pips {}'.format(pip_path)
        outline += ' --gw-sig "{}"'.format(p_cutoff)
        outline += ' --gw-sig-label "{}"'.format(sig_label)
        outline += ' --subgroup-cohorts'
        outline += ' --pdf-height 4.5'
        outline += ' --cnv-panel-height 1.1'
        outline += ' --cnv-panel-space 0.2'
        outline += ' --standardize-frequencies {} cnvs.input.tsv {}'.format(region, cnv)
        outline += ' refs/HPOs_by_metacohort.table.tsv refs/GRCh37.genome'
        outline += ' {}/{} 1> /dev/null 2>&1\n'.format(target_dir.strip('/'), csid)
        outfile.write(outline)


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--segments', help='.bed of large segments')
    parser.add_argument('--credsets', help='.bed of genic credible sets')
    parser.add_argument('--phenotable', required=True, help='HPO table with ' +
                        'parent & child tiers for each term. Required.')
    parser.add_argument('--bonf-cutoff', type=float, default=10e-6, help='P-value ' +
                        'cutoff for genome- or exome-wide significance.')
    parser.add_argument('-o', '--outfile', default='stdout', help='Output .bed of ' +
                        'segments with formatted counts. [default: stdout]')
    parser.add_argument('--target-directory', default='./', help='Output ' +
                        'directory for plots. [default: execution directory]')
    args = parser.parse_args()

    # Require at least one of --segments or --credsets
    if args.segments is None and args.credsets is None:
        exit('Must supply at least one of --segments and/or --credsets')

    # Open connections to output files
    if args.outfile in 'stdout - /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Read phenotable
    phenos = pd.read_csv(args.phenotable, sep='\t').\
                rename(columns={'#HPO_term' : 'HPO'}).\
                sort_values(by='HPO', ascending=True)

    # Format plot script calls
    if args.segments is not None:
        format_segment_calls(args.segments, phenos, outfile, args.target_directory,
                             args.bonf_cutoff)
    if args.credsets is not None:
        format_credset_calls(args.credsets, phenos, outfile, args.target_directory,
                             args.bonf_cutoff)


if __name__ == '__main__':
    main()

