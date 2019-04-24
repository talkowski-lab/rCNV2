#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens in sliding windows, genome-wide


# Launch docker image
docker run --rm -it talkowski/rcnv


# Copy all filtered CNV data, sliding windows, and other references 
# from the project Google Bucket (note: requires permissions)
gcloud auth login
mkdir cleaned_cnv/
gsutil cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
mkdir windows/
gsutil cp -r gs://rcnv_project/cleaned_data/binned_genome/* windows/


# Gather list of CNV datasets to intersect
while read cohort; do
  for pheno in CASE CTRL; do
    if [ -e cleaned_cnv/$cohort.rCNV.$pheno.bed.gz ] && \
       [ $( zcat cleaned_cnv/$cohort.rCNV.$pheno.bed.gz \
            | fgrep -v "#" | wc -l ) -gt 0 ]; then
      echo -e "./cleaned_cnv/$cohort.rCNV.$pheno.bed.gz"
    fi
  done
done < <( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f1 ) \
> rCNV_bed_paths.list


# Count CNVs in cases and controls per cohort, split by CNV type
# Bins with no annotations
/opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
  -z \
  -o GRCh37.100kb_bins_10kb_steps.raw.rCNV_counts.bed.gz \
  rCNV_bed_paths.list \
  windows/GRCh37.100kb_bins_10kb_steps.raw.bed.gz
for CNV in DEL DUP; do
  /opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
    -z \
    --cnv ${CNV} \
    -o GRCh37.100kb_bins_10kb_steps.raw.rCNV_counts.${CNV}.bed.gz \
    rCNV_bed_paths.list \
    windows/GRCh37.100kb_bins_10kb_steps.raw.bed.gz
done  
# Bins with raw annotations
/opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
  -z \
  -o GRCh37.100kb_bins_10kb_steps.annotated.rCNV_counts.bed.gz \
  rCNV_bed_paths.list \
  windows/GRCh37.100kb_bins_10kb_steps.annotated.bed.gz
for CNV in DEL DUP; do
  /opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
    -z \
    --cnv ${CNV} \
    -o GRCh37.100kb_bins_10kb_steps.annotated.rCNV_counts.${CNV}.bed.gz \
    rCNV_bed_paths.list \
    windows/GRCh37.100kb_bins_10kb_steps.annotated.bed.gz
done  
# Bins with eigenannotations
/opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
  -z \
  -o GRCh37.100kb_bins_10kb_steps.annotated.eigen.rCNV_counts.${CNV}.bed.gz \
  rCNV_bed_paths.list \
  windows/GRCh37.100kb_bins_10kb_steps.annotated.eigen.bed.gz
for CNV in DEL DUP; do
  /opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
    -z \
    --cnv ${CNV} \
    -o GRCh37.100kb_bins_10kb_steps.annotated.eigen.rCNV_counts.${CNV}.bed.gz \
    rCNV_bed_paths.list \
    windows/GRCh37.100kb_bins_10kb_steps.annotated.eigen.bed.gz
done  

