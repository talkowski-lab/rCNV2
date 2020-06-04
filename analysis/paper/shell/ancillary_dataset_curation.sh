#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Miscellaneous preprocessing tasks of ancillary datasets for rCNV2 formal analysis


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"


# Download necessary reference files and metadata
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
  ${rCNV_bucket}/raw_data/other/satterstrom_asc_dnms.raw.tsv.gz \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
  ./


# Download & reformat DDD de novo mutation table
wget https://www.biorxiv.org/content/biorxiv/early/2020/04/01/797787/DC4/embed/media-4.txt
/opt/rCNV2/data_curation/other/curate_ddd_dnms.py \
  --dnm-tsv media-4.txt \
  --genes gencode.v19.canonical.pext_filtered.genes.list \
  -o ddd_dnm_counts.tsv.gz \
  -z


# Reformat ASC de novo mutation table
/opt/rCNV2/data_curation/other/curate_asc_dnms.py \
  --dnm-tsv satterstrom_asc_dnms.raw.tsv.gz \
  --genes gencode.v19.canonical.pext_filtered.genes.list \
  -o asc_dnm_counts.tsv.gz \
  -z
/opt/rCNV2/data_curation/other/curate_asc_dnms.py \
  --dnm-tsv satterstrom_asc_dnms.raw.tsv.gz \
  --genes gencode.v19.canonical.pext_filtered.genes.list \
  --controls \
  -o asc_dnm_counts.unaffecteds.tsv.gz \
  -z


# Download & reformat gnomAD mutation rate table
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
/opt/rCNV2/data_curation/other/clean_gnomad_mutation_rates.R \
  gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
  gencode.v19.canonical.pext_filtered.genes.list \
  gene_mutation_rates.tsv
gzip -f gene_mutation_rates.tsv


# Copy curated DNMs and mutation rates to gs:// bucket (note: requires permissions)
gsutil -m cp \
  ddd_dnm_counts.tsv.gz \
  asc_dnm_counts*tsv.gz \
  gene_mutation_rates.tsv.gz \
  ${rCNV_bucket}/analysis/paper/data/misc/
