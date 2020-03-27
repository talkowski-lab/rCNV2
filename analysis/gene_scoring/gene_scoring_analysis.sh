#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Prediction of dosage sensitivity scores for all genes


# Launch docker image
docker run --rm -it talkowski/rcnv
gcloud auth login


# Test/dev parameters
hpo="HP:0000118"
prefix="HP0000118"
freq_code="rCNV"
# CNV="DEL"
# phenotype_list="refs/test_phenotypes.list"
# metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
# gtf="genes/gencode.v19.canonical.gtf.gz"
rCNV_bucket="gs://rcnv_project"


# Copy all gene lists, metadata, and rCNV association stats
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/gene_lists ./
mkdir gene_metadata && \
  gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/metadata ./gene_metadata
mkdir stats && \
  gsutil -m cp -r ${rCNV_bucket}/analysis/gene_burden/ ./stats


# Define list of high-confidence dosage-sensitive genes
