#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Miscellaneous setup tasks for rCNV formal analyses


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"


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


# Compute Jaccard similarity index for all pairs of phenotypes
/opt/rCNV2/analysis/paper/scripts/misc_setup/hpo_jaccard.py \
  --jaccardfile rCNV.hpo_jaccard_matrix.tsv \
  --countsfile rCNV.hpo_sample_overlap_counts_matrix.tsv \
  --asymfile rCNV.hpo_sample_overlap_fraction_matrix.tsv \
  refs/test_phenotypes.list \
  phenos/mega.cleaned_phenos.txt


# Reorder HPO terms based on hierarchical clustering (for ordering in plots)
/opt/rCNV2/analysis/paper/scripts/misc_setup/reorder_hpo_terms.py \
  --outfile ${prefix}.reordered_hpos.txt \
  refs/phenotype_groups.HPO_metadata.txt \
  rCNV.hpo_sample_overlap_fraction_matrix.tsv

