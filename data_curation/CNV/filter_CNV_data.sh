#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Commands executed to filter raw CNV data

# Note: these commands are provided as a reference for what was executed on
# the raw CNV data. In practice, these tasks were parallelized in FireCloud.
# See filter_CNV_data.wdl for more details


# Launch docker image
docker run --rm -it talkowski/rcnv


# Copy all raw CNV data and refs from Google Bucket (note: requires permissions)
gcloud auth login
gsutil cp -r gs://rcnv_project/raw_data/cnv ./
gsutil cp -r gs://rcnv_project/refs ./


# Make master BED file of all raw CNV data
zcat cnv/*bed.gz \
| fgrep -v "#" \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| bgzip -c \
> all_raw_cnvs.bed.gz
allcohorts_nsamp=$( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt \
                    | awk '{ sum+=$2 }END{ print sum }' )


# Filter each cohort for rare CNV callset
mkdir rare_cnv_curated/
while read cohort N; do
  /opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
    --minsize 50000 \
    --maxsize 10000000 \
    --nsamp $N \
    --maxfreq 0.01 \
    --recipoverlap 0.5 \
    --dist 50000 \
    --blacklist refs/GRCh37.segDups_plus_simpleRepeats.bed.gz \
    --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
    --blacklist refs/GRCh37.Nmask.bed.gz \
    --xcov 0.3 \
    --allcohorts all_raw_cnvs.bed.gz \
    --allcohorts_nsamp $allcohorts_nsamp \
    --gnomad refs/gnomAD_v2_SV_MASTER.sites.vcf.gz \
    --gnomad-af-field POPMAX_AF \
    --bgzip \
    cnv/$cohort.raw.bed.gz \
    rare_cnv_curated/$cohort.rCNV.bed.gz
done < <( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f1-2 )


# Filter each cohort for ultrarare CNV callset
mkdir ultrarare_cnv_curated/
for cohort in PGC SSC; do
  /opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
    --minsize 50000 \
    --maxsize 10000000 \
    --nsamp $N \
    --maxfreq 0.0001 \
    --recipoverlap 0.5 \
    --dist 50000 \
    --blacklist refs/GRCh37.segDups_plus_simpleRepeats.bed.gz \
    --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
    --blacklist refs/GRCh37.Nmask.bed.gz \
    --xcov 0.3 \
    --allcohorts all_raw_cnvs.bed.gz \
    --allcohorts_nsamp $allcohorts_nsamp \
    --gnomad refs/gnomAD_v2_SV_MASTER.sites.vcf.gz \
    --gnomad-af-field POPMAX_AF \
    --bgzip \
    cnv/$cohort.raw.bed.gz \
    ultrarare_cnv_curated/$cohort.urCNV.bed.gz
done

