#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Helper mini-plots for rCNV2 graphical abstract


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d2"
example_hpo_nocolon="HP0012759"


# Localize all data
gsutil -m cp \
  ${rCNV_bucket}/analysis/sliding_windows/${example_hpo_nocolon}/rCNV/stats/${example_hpo_nocolon}**.meta_analysis.**.bed.gz \
  ${rCNV_bucket}/results/gene_scoring/rCNV.gene_scores.tsv.gz \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz \
  ./

# Miami plot of example phenotype
/opt/rCNV2/analysis/paper/plot/misc/plot_dummy_miami.R \
  ${example_hpo_nocolon}.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
  ${example_hpo_nocolon}.rCNV.DUP.sliding_window.meta_analysis.stats.bed.gz \
  ${prefix}.example_miami.png

# Mini scatterplot of gene scores
/opt/rCNV2/analysis/paper/plot/misc/plot_mini_gene_score_scatter.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  ${prefix}.mini_genescore_scatter.png

# Dummy feature matrix for fine mapping
/opt/rCNV2/analysis/paper/plot/misc/plot_mini_feature_matrix.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz \
  ${prefix}.mini_feature_matrix.pdf


# Copy all plots to output directory
gsutil -m cp \
  ${prefix}.example_miami.png \
  ${prefix}.mini_genescore_scatter.png \
  ${prefix}.mini_feature_matrix.pdf \
  ${rCNV_bucket}/analysis/paper/plots/misc/

