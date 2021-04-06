#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Extract noncoding CNV subsets from cleaned rCNV data
# Note: assumes CNVs have been filtered using filter_CNV_data.wdl


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Copy necessary data
mkdir cleaned_cnvs/
gsutil -m cp gs://rcnv_project/cleaned_data/cnv/*rCNV* cleaned_cnvs/
gsutil -m cp gs://rcnv_project/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* ./
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/rCNV_metacohort_list.txt ./
gsutil -m cp gs://rcnv_project/cleaned_data/genes/gene_lists/gnomad.v2.1.1.likely_unconstrained.genes.list ./
gsutil -m cp gs://rcnv_project/cleaned_data/genes/gene_lists/HP0000118.HPOdb.genes.list ./


# Subtract OMIM genes from likely unconstrained genes as whitelist
fgrep -wvf \
  HP0000118.HPOdb.genes.list \
  gnomad.v2.1.1.likely_unconstrained.genes.list \
| sort -V | uniq \
> loose_noncoding_whitelist.genes.list


# Extract noncoding CNVs for each metacohort
mkdir noncoding
while read cohort; do
  echo $cohort
  /opt/rCNV2/data_curation/CNV/extract_noncoding_cnvs.py \
    --strict-outbed noncoding/$cohort.rCNV.strict_noncoding.bed.gz \
    --loose-outbed noncoding/$cohort.rCNV.loose_noncoding.bed.gz \
    --whitelist loose_noncoding_whitelist.genes.list \
    --pad-exons 50000 \
    --bgzip \
    cleaned_cnvs/$cohort.rCNV.bed.gz \
    gencode.v19.canonical.pext_filtered.gtf.gz
done < <( cut -f1 rCNV_metacohort_list.txt )


# Copy all cleaned noncoding CNV data to Google Bucket (note: requires permissions)
gsutil -m cp -r noncoding gs://rcnv_project/cleaned_data/cnv/
