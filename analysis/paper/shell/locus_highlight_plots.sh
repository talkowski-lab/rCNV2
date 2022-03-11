#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block to generate locus-level higlight plots for rCNV2 manuscript


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d2"


# Download necessary data (note: requires permissions)
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/cnv ./
mkdir meta_stats/
gsutil -m cp \
  ${rCNV_bucket}/analysis/sliding_windows/**.rCNV.**.sliding_window.meta_analysis.stats.bed.gz \
  ${rCNV_bucket}/analysis/gene_burden/**.rCNV.**.gene_burden.meta_analysis.stats.bed.gz \
  ${rCNV_bucket}/analysis/crb_burden/**.rCNV.**.crb_burden.meta_analysis.stats.bed.gz \
  meta_stats/
find meta_stats/ -name "*bed.gz" | xargs -I {} tabix -p bed -f {}
gsutil -m cp -r \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/rCNV.*.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv \
  ${rCNV_bucket}/results/* \
  ./


# Download reference files (note: requires permissions)
mkdir refs/ 
gsutil -m cp \
  ${rCNV_bucket}/analysis/analysis_refs/** \
  ${rCNV_bucket}/refs/GRCh37.cytobands.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/hpo/${prefix}.reordered_hpos.txt \
  ${rCNV_bucket}/cleaned_data/phenotypes/hpo_logs_metadata/phenotype_groups.HPO_metadata.txt \
  ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  ${rCNV_bucket}/refs/REP_state_manifest.tsv \
  refs/


# Prepare generic input files
while read cohort; do
  path="cnv/$cohort.rCNV.bed.gz"
  echo -e "$cohort\t$path"
done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt | cut -f1 ) \
> cnvs.input.tsv


# Get genome-wide significance thresholds
example_hpo="HP0012759"
gw_cutoff=$( awk -v FS="\t" -v hpo=${example_hpo} '{ if ($1==hpo) print $2 }' \
             refs/sliding_window.rCNV.DEL.bonferroni_pval.hpo_cutoffs.tsv )
ew_cutoff=$( awk -v FS="\t" -v hpo=${example_hpo} '{ if ($1==hpo) print $2 }' \
             refs/gene_burden.rCNV.DEL.bonferroni_pval.hpo_cutoffs.tsv )


# Plot one standardized locus highlight plot for each large segment and credible set
for dir in all_large_segments all_credsets; do
  if ! [ -e $dir ]; then
    mkdir $dir
  fi
done
/opt/rCNV2/analysis/paper/scripts/locus_highlights/enumerate_locus_highlight_plot_calls.py \
  --segments segment_association/rCNV.final_segments.loci.bed.gz \
  --phenotable refs/phenotype_groups.HPO_metadata.txt \
  --bonf-cutoff "$gw_cutoff" \
  --target-directory all_large_segments/ \
> plot_all_locus_highlights.sh
/opt/rCNV2/analysis/paper/scripts/locus_highlights/enumerate_locus_highlight_plot_calls.py \
  --credsets gene_association/rCNV.final_genes.credible_sets.bed.gz \
  --phenotable refs/phenotype_groups.HPO_metadata.txt \
  --bonf-cutoff "$ew_cutoff" \
  --target-directory all_credsets/ \
>> plot_all_locus_highlights.sh
bash plot_all_locus_highlights.sh


# Declare shared parameters formain highlights
export idio_space=0.2
export pval_space=0.15
export or_space=0.15
export constraint_space=0.15
export pip_space=0.15
export cnv_space=0.1
export cnv_height=1.2
export cnv_dx=300
export pdf_height=3.25


# Highlight panel for GMEB2 duplications
if ! [ -e main_highlights ]; then
  mkdir main_highlights
fi
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_locus_highlight.R \
  --case-hpos "HP:0012639" \
  --highlight-hpo "HP:0002011" \
  --highlights "20:62152076-62251229;20:62218954-62251229" \
  --sumstats meta_stats/HP0002011.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz \
  --genic-sumstats \
  --cytobands refs/GRCh37.cytobands.bed.gz \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --label-genes "GMEB2" \
  --genes-before-sumstats \
  --constraint refs/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  --pips rCNV.DUP.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv \
  --gw-sig "$ew_cutoff" \
  --gw-sig-label "Exome-Wide significance" \
  --standardize-frequencies \
  --collapse-cohorts \
  --idio-panel-space "$idio_space" \
  --pval-panel-space "$pval_space" \
  --or-panel-space "$or_space" \
  --constraint-panel-space "$constraint_space" \
  --pip-panel-space "$pip_space" \
  --cnv-panel-space "$cnv_space" \
  --cnv-panel-height "$cnv_height" \
  --dx "$cnv_dx" \
  --pdf-height "$pdf_height" \
  20:61975000-62425000 \
  cnvs.input.tsv \
  DUP \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  main_highlights/GMEB2


# Highlight panel for KIF13A duplications
if ! [ -e main_highlights ]; then
  mkdir main_highlights
fi
fdr_cutoff=$( /opt/rCNV2/analysis/other/estimate_p_for_fdr.R \
                meta_stats/HP0031466.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz \
                2>/dev/null )
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_locus_highlight.R \
  --case-hpos "HP:0000707" \
  --highlight-hpo "HP:0031466" \
  --highlights "6:17763923-17987800;6:17763923-17987800" \
  --sumstats meta_stats/HP0031466.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz \
  --genic-sumstats \
  --cytobands refs/GRCh37.cytobands.bed.gz \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --label-genes "KIF13A" \
  --genes-before-sumstats \
  --constraint refs/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  --pips rCNV.DUP.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv \
  --gw-sig "$fdr_cutoff" \
  --gw-sig-label "FDR < 1%" \
  --standardize-frequencies \
  --collapse-cohorts \
  --idio-panel-space "$idio_space" \
  --pval-panel-space "$pval_space" \
  --or-panel-space "$or_space" \
  --constraint-panel-space "$constraint_space" \
  --pip-panel-space "$pip_space" \
  --cnv-panel-space "$cnv_space" \
  --cnv-panel-height "$cnv_height" \
  --dx "$cnv_dx" \
  --pdf-height "$pdf_height" \
  6:17250000-18250000 \
  cnvs.input.tsv \
  DUP \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  main_highlights/KIF13A


# Highlight panel for ANKRD11 duplications
if ! [ -e main_highlights ]; then
  mkdir main_highlights
fi
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_locus_highlight.R \
  --case-hpos "HP:0033127" \
  --highlight-hpo "HP:0001507" \
  --highlights "16:89160216-89556969;16:89334037-89556969" \
  --sumstats meta_stats/HP0001507.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz \
  --genic-sumstats \
  --cytobands refs/GRCh37.cytobands.bed.gz \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --label-genes "ANKRD11" \
  --genes-before-sumstats \
  --constraint refs/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  --pips rCNV.DUP.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv \
  --gw-sig "$ew_cutoff" \
  --gw-sig-label "Exome-Wide Significance" \
  --standardize-frequencies \
  --collapse-cohorts \
  --idio-panel-space "$idio_space" \
  --pval-panel-space "$pval_space" \
  --or-panel-space "$or_space" \
  --constraint-panel-space "$constraint_space" \
  --pip-panel-space "$pip_space" \
  --cnv-panel-space "$cnv_space" \
  --cnv-panel-height "$cnv_height" \
  --dx "$cnv_dx" \
  --pdf-height "$pdf_height" \
  16:88900000-89900000 \
  cnvs.input.tsv \
  DUP \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  main_highlights/ANKRD11


# Highlight panel for IGF1R duplications
if ! [ -e main_highlights ]; then
  mkdir main_highlights
fi
fdr_cutoff=$( /opt/rCNV2/analysis/other/estimate_p_for_fdr.R \
                meta_stats/HP0009121.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz \
                2>/dev/null )
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_locus_highlight.R \
  --case-hpos "HP:0033127" \
  --highlight-hpo "HP:0009121" \
  --highlights "15:99192199-100273626;15:99192199-99507759" \
  --sumstats meta_stats/HP0009121.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz \
  --genic-sumstats \
  --cytobands refs/GRCh37.cytobands.bed.gz \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --label-genes "IGF1R" \
  --genes-before-sumstats \
  --constraint refs/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  --pips rCNV.DUP.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv \
  --gw-sig "$fdr_cutoff" \
  --gw-sig-label "FDR < 1%" \
  --standardize-frequencies \
  --collapse-cohorts \
  --idio-panel-space "$idio_space" \
  --pval-panel-space "$pval_space" \
  --or-panel-space "$or_space" \
  --constraint-panel-space "$constraint_space" \
  --pip-panel-space "$pip_space" \
  --cnv-panel-space "$cnv_space" \
  --cnv-panel-height "$cnv_height" \
  --dx "$cnv_dx" \
  --pdf-height "$pdf_height" \
  15:98471248-100994577 \
  cnvs.input.tsv \
  DUP \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  main_highlights/IGF1R


# Highlight panel for SHANK3 deletions
if ! [ -e main_highlights ]; then
  mkdir main_highlights
fi
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_locus_highlight.R \
  --case-hpos "HP:0000707" \
  --highlight-hpo "HP:0012759" \
  --highlights "22:50166930-51171641;22:51113069-51171641" \
  --sumstats meta_stats/HP0012759.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz \
  --genic-sumstats \
  --cytobands refs/GRCh37.cytobands.bed.gz \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --label-genes "SHANK3" \
  --genes-before-sumstats \
  --constraint refs/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  --pips rCNV.DEL.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv \
  --gw-sig-label "Exome-Wide Significance" \
  --gw-sig-label-side "below" \
  --standardize-frequencies \
  --collapse-cohorts \
  --idio-panel-space "$idio_space" \
  --pval-panel-space "$pval_space" \
  --or-panel-space "$or_space" \
  --constraint-panel-space "$constraint_space" \
  --pip-panel-space "$pip_space" \
  --cnv-panel-space "$cnv_space" \
  --cnv-panel-height "$cnv_height" \
  --dx "$cnv_dx" \
  --pdf-height "$pdf_height" \
  22:49800000-51230000 \
  cnvs.input.tsv \
  DEL \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  main_highlights/SHANK3


# Highlight panel for RAI1 deletions
if ! [ -e main_highlights ]; then
  mkdir main_highlights
fi
fdr_cutoff=$( /opt/rCNV2/analysis/other/estimate_p_for_fdr.R \
                meta_stats/HP0012759.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz \
                2>/dev/null )
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_locus_highlight.R \
  --case-hpos "HP:0000118" \
  --highlight-hpo "HP:0012759" \
  --highlights "17:16946073-19320589;17:17584786-17714767" \
  --sumstats meta_stats/HP0012759.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz \
  --genic-sumstats \
  --cytobands refs/GRCh37.cytobands.bed.gz \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --label-genes "RAI1" \
  --genes-before-sumstats \
  --constraint refs/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  --pips rCNV.DEL.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv \
  --gw-sig "$fdr_cutoff" \
  --gw-sig-label "FDR < 1%" \
  --standardize-frequencies \
  --collapse-cohorts \
  --idio-panel-space "$idio_space" \
  --pval-panel-space "$pval_space" \
  --or-panel-space "$or_space" \
  --constraint-panel-space "$constraint_space" \
  --pip-panel-space "$pip_space" \
  --cnv-panel-space "$cnv_space" \
  --cnv-panel-height "$cnv_height" \
  --dx "$cnv_dx" \
  --pdf-height "$pdf_height" \
  17:15700000-20900000 \
  cnvs.input.tsv \
  DEL \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  main_highlights/RAI1


# Copy all plots to final gs:// directory
gsutil -m cp -r \
  all_large_segments \
  all_credsets \
  main_highlights \
  ${rCNV_bucket}/analysis/paper/plots/locus_highlights/

