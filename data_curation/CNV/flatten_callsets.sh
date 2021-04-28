#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Flatten aggregated CNV callsets as needed


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv:latest


# Download data (note: permissions required)
gcloud auth login
gsutil -m cp \
  gs://rcnv_project/raw_data/other/cnv_callsets_preflatten/* \
  gs://rcnv_project/cleaned_data/phenotypes/filtered/EstBB.final_cooccurrence_table.tsv.gz \
  gs://rcnv_project/cleaned_data/phenotypes/filtered/BioVU.final_cooccurrence_table.tsv.gz \
  ./


# Flatten EstBB data
/opt/rCNV2/data_curation/CNV/flatten_callset.py \
  --outfile EstBB.raw.bed.gz \
  --hpo-pairs EstBB.final_cooccurrence_table.tsv.gz \
  --cohort EstBB \
  --bgzip \
  EstBB_CNV_QC_woRelatives_defragmented.countsPerHPO.tsv.gz
tabix -f EstBB.raw.bed.gz


# Copy to Google Cloud (note: permissions required)
gsutil -m cp \
  EstBB.raw.bed.gz* \
  gs://rcnv_project/raw_data/cnv/
