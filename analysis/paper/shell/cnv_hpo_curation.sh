#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Miscellaneous preprocessing tasks of HPO and CNV data for rCNV formal analyses


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"
export control_hpo="HEALTHY_CONTROL"


# Localize necessary data (note: requires permissions)
mkdir refs/ 
gsutil -m cp \
  ${rCNV_bucket}/analysis/analysis_refs/* \
  ${rCNV_bucket}/cleaned_data/phenotypes/hpo_logs_metadata/phenotype_groups.HPO_metadata.txt \
  refs/
mkdir phenos/
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/phenotypes/filtered/mega.cleaned_phenos.txt \
  phenos/
mkdir cnvs/
gsutil -m cp \
  ${rCNV_bucket}/raw_data/cnv/*raw.bed.gz* \
  ${rCNV_bucket}/cleaned_data/cnv/*rCNV.bed.gz* \
  cnvs/


# Compute Jaccard similarity index for all pairs of phenotypes
/opt/rCNV2/analysis/paper/scripts/misc_setup/hpo_jaccard.py \
  --jaccardfile ${prefix}.hpo_jaccard_matrix.tsv \
  --countsfile ${prefix}.hpo_sample_overlap_counts_matrix.tsv \
  --asymfile ${prefix}.hpo_sample_overlap_fraction_matrix.tsv \
  refs/test_phenotypes.list \
  phenos/mega.cleaned_phenos.txt


# Reorder HPO terms based on hierarchical clustering (for ordering in plots)
/opt/rCNV2/analysis/paper/scripts/misc_setup/reorder_hpo_terms.py \
  --outfile ${prefix}.reordered_hpos.txt \
  refs/phenotype_groups.HPO_metadata.txt \
  ${prefix}.hpo_sample_overlap_fraction_matrix.tsv


# Copy HPO sample overlap matrices and reordered HPO list to GCP (note: requires permissions)
gsutil -m cp \
  ${prefix}.hpo_jaccard_matrix.tsv \
  ${prefix}.hpo_sample_overlap_counts_matrix.tsv \
  ${prefix}.hpo_sample_overlap_fraction_matrix.tsv \
  ${prefix}.reordered_hpos.txt \
  ${rCNV_bucket}/analysis/paper/data/hpo/


# Plot summary figure of HPOs
/opt/rCNV2/analysis/paper/plot/misc/plot_hpo_summary.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  ${prefix}.reordered_hpos.txt \
  refs/phenotype_groups.HPO_metadata.txt \
  refs/HPOs_by_metacohort.table.tsv \
  ${prefix}.hpo_sample_overlap_fraction_matrix.tsv \
  ${prefix}.hpo_summary
gsutil -m cp \
  ${prefix}.hpo_summary*pdf \
  ${rCNV_bucket}/analysis/paper/plots/misc/


# Gather stats per cohort for raw and filtered CNVs
awk -v OFS="\t" '{ print $0, "cnvs/"$1".raw.bed.gz" }' \
  /opt/rCNV2/refs/rCNV_sample_counts.tsv \
| fgrep -v "#" > raw_cnv.input.tsv
/opt/rCNV2/data_curation/CNV/build_cnv_stats_table.py \
  --tsv raw_cnv.stats.tsv \
  raw_cnv.input.tsv
awk -v OFS="\t" '{ print $0, "cnvs/"$1".rCNV.bed.gz" }' \
  /opt/rCNV2/refs/rCNV_sample_counts.tsv \
| fgrep -v "#" > rCNV.input.tsv
/opt/rCNV2/data_curation/CNV/build_cnv_stats_table.py \
  --tsv rCNV.stats.tsv \
  rCNV.input.tsv


# Plot panels for summary figure of CNV filtering pipeline
/opt/rCNV2/analysis/paper/plot/misc/plot_cnv_filtering_summary.R \
  --control-hpo ${control_hpo} \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  raw_cnv.input.tsv \
  raw_cnv.stats.tsv \
  rCNV.input.tsv \
  rCNV.stats.tsv \
  refs/rCNV_metacohort_list.txt \
  ${prefix}.cnv_summary
gsutil -m cp \
  ${prefix}.cnv_summary*pdf \
  ${rCNV_bucket}/analysis/paper/plots/misc/

