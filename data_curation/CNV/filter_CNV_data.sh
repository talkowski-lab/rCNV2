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



# Sanity check raw CNV data (to ensure proper formatting)
# Check headers
for bed in cnv/*raw.bed.gz; do
  echo $bed
  zcat $bed | head -n1
  echo -e "\n\n"
done
# Check HPO code formatting
for bed in cnv/*raw.bed.gz; do
  echo $bed
  zcat $bed | fgrep -v "#" | cut -f6 | sed 's/\;/\n/g' | sort | uniq -c
  echo -e "\n\n"
done


# # Make master BED file of all raw CNV data
# zcat cnv/*bed.gz \
# | fgrep -v "#" \
# | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
# | bgzip -c \
# > all_raw_cnvs.bed.gz
# allcohorts_nsamp=$( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt \
#                     | awk '{ sum+=$2 }END{ print sum }' )

# Make master list of all raw CNV data per cohort
while read cohort N; do
  if [ -s cnv/$cohort.raw.bed.gz ]; then
    echo "$cohort"
    echo "$N"
    echo "cnv/$cohort.raw.bed.gz"
  fi \
  | paste -s
done < <( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f1-2 ) \
| sort -nk2,2 \
> raw_CNVs.per_cohort.txt


# Filter each cohort for rare CNV callset
mkdir rare_cnv_curated/
while read cohort N; do
  /opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
    --minsize 100000 \
    --maxsize 10000000 \
    --nsamp $N \
    --maxfreq 0.01 \
    --recipoverlap 0.5 \
    --dist 50000 \
    --blacklist refs/GRCh37.segDups_plus_simpleRepeats.bed.gz \
    --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
    --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
    --xcov 0.3 \
    --cohorts-list raw_CNVs.per_cohort.txt \
    --vcf refs/gnomAD_v2_SV_MASTER.sites.vcf.gz \
    --vcf refs/1000Genomes_phase3.sites.vcf.gz \
    --vcf refs/CCDG_Abel_bioRxiv.sites.vcf.gz \
    --vcf-af-fields AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,OTH_AF,POPMAX_AF \
    --bgzip \
    cnv/$cohort.raw.bed.gz \
    rare_cnv_curated/$cohort.rCNV.bed.gz
done < <( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f1-2 )
    


# Filter each cohort for ultrarare CNV callset
mkdir ultrarare_cnv_curated/
while read cohort N; do
  /opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
    --minsize 100000 \
    --maxsize 10000000 \
    --nsamp $N \
    --maxfreq 0.0001 \
    --recipoverlap 0.5 \
    --dist 50000 \
    --blacklist refs/GRCh37.segDups_plus_simpleRepeats.bed.gz \
    --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
    --blacklist refs/GRCh37.Nmask.bed.gz \
    --xcov 0.3 \
    --cohorts-list raw_CNVs.per_cohort.txt \
    --vcf refs/gnomAD_v2_SV_MASTER.sites.vcf.gz \
    --vcf refs/1000Genomes_phase3.sites.vcf.gz \
    --vcf refs/CCDG_Abel_bioRxiv.sites.vcf.gz \
    --vcf-af-fields AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,OTH_AF,POPMAX_AF \
    --bgzip \
    cnv/$cohort.raw.bed.gz \
    ultrarare_cnv_curated/$cohort.urCNV.bed.gz
done < <( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f1-2 )


# Debugging code chunk: chr22 rCNVs for BCH cohort
mkdir debug
/opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
  --chr 22 \
  --minsize 100000 \
  --maxsize 10000000 \
  --nsamp 3591 \
  --maxfreq 0.01 \
  --recipoverlap 0.5 \
  --dist 50000 \
  --blacklist refs/GRCh37.segDups_plus_simpleRepeats.bed.gz \
  --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
  --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
  --xcov 0.3 \
  --cohorts-list raw_CNVs.per_cohort.txt \
  --vcf refs/gnomAD_v2_SV_MASTER.sites.vcf.gz \
  --vcf refs/1000Genomes_phase3.sites.vcf.gz \
  --vcf refs/CCDG_Abel_bioRxiv.sites.vcf.gz \
  --vcf-af-fields AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,OTH_AF,POPMAX_AF \
  --bgzip \
  cnv/BCH.raw.bed.gz \
  debug/BCH.rCNV.bed.gz
# /opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
#   --chr 22 \
#   --minsize 100000 \
#   --maxsize 10000000 \
#   --nsamp 3591 \
#   --maxfreq 0.01 \
#   --recipoverlap 0.5 \
#   --dist 50000 \
#   --blacklist refs/GRCh37.segDups_plus_simpleRepeats.bed.gz \
#   --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
#   --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
#   --xcov 0.3 \
#   --cohorts-list <( fgrep BCH raw_CNVs.per_cohort.txt ) \
#   --vcf-af-fields AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,OTH_AF,POPMAX_AF \
#   --bgzip \
#   cnv/BCH.raw.bed.gz \
#   debug/BCH.rCNV.bed.gz

/opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
  --chr 22 \
  --minsize 100000 \
  --maxsize 10000000 \
  --nsamp 3591 \
  --maxfreq 0.0001 \
  --recipoverlap 0.5 \
  --dist 50000 \
  --blacklist refs/GRCh37.segDups_plus_simpleRepeats.bed.gz \
  --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
  --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
  --xcov 0.3 \
  --cohorts-list raw_CNVs.per_cohort.txt \
  --vcf refs/gnomAD_v2_SV_MASTER.sites.vcf.gz \
  --vcf refs/1000Genomes_phase3.sites.vcf.gz \
  --vcf refs/CCDG_Abel_bioRxiv.sites.vcf.gz \
  --vcf-af-fields AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,OTH_AF,POPMAX_AF \
  --bgzip \
  cnv/BCH.raw.bed.gz \
  debug/BCH.uCNV.bed.gz


# Debugging code chunk: chr22 rCNVs for PGC cohort
/opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
  --chr 22 \
  --minsize 100000 \
  --maxsize 10000000 \
  --nsamp 41371 \
  --maxfreq 0.01 \
  --recipoverlap 0.5 \
  --dist 50000 \
  --blacklist refs/GRCh37.segDups_plus_simpleRepeats.bed.gz \
  --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
  --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
  --xcov 0.3 \
  --cohorts-list raw_CNVs.per_cohort.txt \
  --vcf refs/gnomAD_v2_SV_MASTER.sites.vcf.gz \
  --vcf refs/1000Genomes_phase3.sites.vcf.gz \
  --vcf refs/CCDG_Abel_bioRxiv.sites.vcf.gz \
  --vcf-af-fields AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,OTH_AF,POPMAX_AF \
  --bgzip \
  cnv/PGC.raw.bed.gz \
  rare_cnv_curated/PGC.rCNV.bed.gz

