#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for secondary analyses of large rCNV-associated segments


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"


# Localize all analysis refs, sliding window meta-analysis stats, and large segment results
mkdir refs/ 
gsutil -m cp \
   ${rCNV_bucket}/analysis/analysis_refs/** \
   refs/
mkdir meta_stats/
gsutil -m cp \
  ${rCNV_bucket}/analysis/sliding_windows/**.rCNV.**.sliding_window.meta_analysis.stats.bed.gz \
  meta_stats/
gsutil -m cp \
  ${rCNV_bucket}/results/segment_association/* \
  ./


# Tabix all meta-analysis stats
find meta_stats/ -name "*meta_analysis.stats.bed.gz" \
| xargs -I {} tabix -f {}


# Compute effect size and max P-value per phenotype per final segment
/opt/rCNV2/analysis/paper/scripts/large_segments/calc_all_seg_stats.py \
  --lnors-out rCNV.final_segments.loci.all_effect_sizes.tsv \
  --pvals-out rCNV.final_segments.loci.all_pvalues.tsv \
  rCNV.final_segments.loci.bed.gz \
  refs/test_phenotypes.list \
  meta_stats
gzip -f rCNV.final_segments.loci.all_effect_sizes.tsv
gzip -f rCNV.final_segments.loci.all_pvalues.tsv


# Plot master grid summarizing segment association across all phenotypes
/opt/rCNV2/analysis/paper/plot/large_segments/plot_association_grid.R
