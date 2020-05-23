#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Helper mini-plots for rCNV2 graphical abstract


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"
example_hpo_nocolon="HP0012759"


# Localize all data
gsutil -m cp \
  ${rCNV_bucket}/results/gene_scoring/rCNV.gene_scores.tsv.gz \
  ${rCNV_bucket}/analysis/sliding_windows/${example_hpo_nocolon}/rCNV/stats/${example_hpo_nocolon}**.meta_analysis.**.bed.gz \
  ./

# Plot mini scatterplot of gene scores
/opt/rCNV2/analysis/paper/plot/misc/plot_mini_gene_score_scatter.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  ${prefix}.mini_genescore_scatter.png

# Plot Miami plot of example phenotype
/opt/rCNV2/analysis/paper/plot/misc/plot_dummy_miami.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  ${example_hpo_nocolon}.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
  ${example_hpo_nocolon}.rCNV.DUP.sliding_window.meta_analysis.stats.bed.gz \
  ${prefix}.example_miami.png

# Copy all plots to output directory
gsutil -m cp \
  ${prefix}.mini_genescore_scatter.png \
  ${prefix}.example_miami.png \
  ${rCNV_bucket}/analysis/paper/plots/misc/

