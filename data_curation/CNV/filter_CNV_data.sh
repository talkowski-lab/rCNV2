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
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv


# Copy all raw CNV data and refs from Google Bucket (note: requires permissions)
gcloud auth login
gsutil -m cp -r gs://rcnv_project/raw_data/cnv ./
gsutil -m cp -r gs://rcnv_project/refs ./
gsutil -m cp gs://analysis/paper/data/large_segments/lit_GDs.*.bed.gz ./refs/



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


# Make list of all populations to filter for allele frequency
for pop in AFR AMR CSA EAS EUR MID OCN OTH SAS; do
  echo $pop | awk -v OFS="\n" '{ print $1"_AF", $1"_NONREF_FREQ" }'
done > all_pop_af_fields.txt


# Filter each cohort for rare CNV callset
mkdir rare_cnv_curated/
while read cohort N; do
  /opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
    --minsize 100000 \
    --maxsize 20000000 \
    --maxfreq 0.01 \
    --recipoverlap 0.5 \
    --dist 100000 \
    --whitelist refs/lit_GDs.hc.bed.gz \
    --whitelist refs/lit_GDs.mc.bed.gz \
    --wrecip 0.75 \
    --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
    --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
    --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
    --xcov 0.5 \
    --cohorts-list raw_CNVs.per_cohort.txt \
    --vcf refs/gnomad_v2.1_sv.nonneuro.sites.vcf.gz \
    --vcf refs/CCDG_Abel_bioRxiv.sites.vcf.gz \
    --vcf refs/1000Genomes_HGSV_highCov.sites.vcf.gz \
    --vcf refs/HGDP.hg19.sites.vcf.gz \
    --vcf-af-fields "AF,$( paste -s -d, all_pop_af_fields.txt )" \
    --genome refs/GRCh37.autosomes.genome \
    --bgzip \
    cnv/$cohort.raw.bed.gz \
    rare_cnv_curated/$cohort.rCNV.bed.gz
done < <( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f1-2 )


# Filter each cohort for very rare CNV callset
mkdir veryrare_cnv_curated/
while read cohort N; do
  /opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
    --minsize 100000 \
    --maxsize 20000000 \
    --maxfreq 0.001 \
    --recipoverlap 0.5 \
    --dist 100000 \
    --whitelist refs/lit_GDs.hc.bed.gz \
    --whitelist refs/lit_GDs.mc.bed.gz \
    --wrecip 0.75 \
    --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
    --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
    --blacklist refs/GRCh37.Nmask.bed.gz \
    --xcov 0.5 \
    --cohorts-list raw_CNVs.per_cohort.txt \
    --vcf refs/gnomad_v2.1_sv.nonneuro.sites.vcf.gz \
    --vcf refs/CCDG_Abel_bioRxiv.sites.vcf.gz \
    --vcf refs/1000Genomes_HGSV_highCov.sites.vcf.gz \
    --vcf refs/HGDP.hg19.sites.vcf.gz \
    --vcf-af-fields "AF,$( paste -s -d, all_pop_af_fields.txt )" \
    --genome refs/GRCh37.autosomes.genome \
    --bgzip \
    cnv/$cohort.raw.bed.gz \
    ultrarare_cnv_curated/$cohort.vCNV.bed.gz
done < <( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f1-2 )


# Filter each cohort for ultrarare CNV callset
mkdir ultrarare_cnv_curated/
while read cohort N; do
  /opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
    --minsize 100000 \
    --maxsize 20000000 \
    --maxfreq 0.0001 \
    --recipoverlap 0.5 \
    --dist 100000 \
    --whitelist refs/lit_GDs.hc.bed.gz \
    --whitelist refs/lit_GDs.mc.bed.gz \
    --wrecip 0.75 \
    --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
    --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
    --blacklist refs/GRCh37.Nmask.bed.gz \
    --xcov 0.5 \
    --cohorts-list raw_CNVs.per_cohort.txt \
    --vcf refs/gnomad_v2.1_sv.nonneuro.sites.vcf.gz \
    --vcf refs/CCDG_Abel_bioRxiv.sites.vcf.gz \
    --vcf refs/1000Genomes_HGSV_highCov.sites.vcf.gz \
    --vcf refs/HGDP.hg19.sites.vcf.gz \
    --vcf-af-fields "AF,$( paste -s -d, all_pop_af_fields.txt )" \
    --genome refs/GRCh37.autosomes.genome \
    --bgzip \
    cnv/$cohort.raw.bed.gz \
    ultrarare_cnv_curated/$cohort.uCNV.bed.gz
done < <( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f1-2 )

