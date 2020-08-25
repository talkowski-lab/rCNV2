#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for secondary analyses of gene association & fine mapping


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"


# Download necessary data (note: requires permissions)
gsutil -m cp \
  ${rCNV_bucket}/results/gene_association/* \
  ./
mkdir refs/
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features*bed.gz \
  gs://rcnv_project/analysis/paper/data/misc/gene_feature_metadata.tsv \
  refs/


# Plot correlation heatmap of gene features
# DEV NOTE: this script is incomplete and does not yet generate a final plot
# TODO: complete this script
/opt/rCNV2/analysis/paper/plot/gene_association/plot_gene_feature_cor_matrix.R \
	ref/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz \
  refs/gencode.v19.canonical.pext_filtered.all_features.no_variation.eigenfeatures.bed.gz \
  refs/gene_feature_metadata.tsv \
  ./gene_features_corplot
