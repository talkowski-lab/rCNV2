#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for secondary analyses of gene dosage sensitivity scores


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"


# Download necessary data (note: requires permissions)
gsutil -m cp \
  ${rCNV_bucket}/results/gene_scoring/rCNV.gene_scores.tsv.gz \
  ./
mkdir refs/
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
  refs/


# Plot simple scatterplot distribution of scores
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_scores_scatter.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  ${prefix}


# Plot distributions of features between subgroups of genes


# Perform feature regressions between subgroups of genes


# Copy all plots to final gs:// directory
gsutil -m cp \
  ${prefix}.gene_scores_scatterplot.png \
  ${prefix}.gene_scores_scatterplot.legend.pdf \
  ${rCNV_bucket}/analysis/paper/plots/gene_scores/
