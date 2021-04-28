#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Flatten aggregated CNV callsets as needed


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv


# Download data
gcloud auth login
gsutil -m cp \
  gs://rcnv_project/raw_data/other/cnv_callsets_preflatten/* \
  gs://rcnv_project/raw_data/phenotypes/EstBB.HPO_terms_full_cooccurrence_table.tsv.gz \
  ./


# Flatten EstBB data
# NOTE: NEED TO RESTRICT CO-OCCURRENCE TABLE TO ONLY HPO TERMS RETAINED IN ANALYSIS
/opt/rCNV2/data_curation/CNV/flatten_callset.py \
  --outfile EstBB.raw.bed.gz \
  --hpo-pairs EstBB.HPO_terms_full_cooccurrence_table.tsv.gz \
  --cohort EstBB \
  --bgzip \
  EstBB_CNV_QC_woRelatives_defragmented.countsPerHPO.tsv.gz
