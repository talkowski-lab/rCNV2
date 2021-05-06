#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Curate the HGDP SV callset for use in rCNV2 analyses
# Also converts HGSVC 1000 Genomes callset to hg19


# Launch GATK-SV docker image & authenticate GCP credentials
docker run --rm -it us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline-base:rlc_posthoc_filtering_cnv_mcnv_compatability_9a8561
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
alias "less=zless"
alias "l=ls -ltrha"


# Clone rCNV repo (since these steps are run in GATK-SV docker)
conda install git
git clone https://github.com/talkowski-lab/rCNV2.git /opt/rCNV2 && \
  cd /opt/rCNV2 && \
  git checkout rev1_staging && \
  cd -


# Download HGDP data from Sanger FTP site
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_structural_variation/SV_CALLSET/del.manta.vcf.gz \
     ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_structural_variation/SV_CALLSET/dup.manta.vcf.gz \
     ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_structural_variation/SV_CALLSET/GS.vcf.gz \
     ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_structural_variation/HGDP.metadata


# Recompress raw files with bgzip, not gzip (for pysam compatibility)
for prefix in del.manta dup.manta GS; do
  zcat $prefix.vcf.gz | bgzip -c > $prefix.vcf.bgz
done


# Reformat original HGDP VCFs
/opt/rCNV2/data_curation/other/reformat_hgdp_vcf.py \
  -a manta \
  -c DEL \
  del.manta.vcf.bgz \
  HGDP.manta.DEL.hg38.vcf.gz
/opt/rCNV2/data_curation/other/reformat_hgdp_vcf.py \
  -a manta \
  -c DUP \
  dup.manta.vcf.bgz \
  HGDP.manta.DUP.hg38.vcf.gz
/opt/rCNV2/data_curation/other/reformat_hgdp_vcf.py \
  -a GS \
  GS.vcf.bgz \
  HGDP.GS.hg38.vcf.gz


# Annotate frequencies by continental population
sed '1d' HGDP.metadata | cut -f1,8 \
| sed -e 's/AFRICA/AFR/g' \
      -e 's/AMERICA/AMR/g' \
      -e 's/CENTRAL_SOUTH_ASIA/CSA/g' \
      -e 's/EAST_ASIA/EAS/g' \
      -e 's/EUROPE/EUR/g' \
      -e 's/MIDDLE_EAST/MID/g' \
      -e 's/OCEANIA/OCN/g' \
> HGDP_pop_assignments.tsv
for CNV in DEL DUP; do
  /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py \
    -p HGDP_pop_assignments.tsv \
    HGDP.manta.$CNV.hg38.vcf.gz \
    stdout \
  | cut -f1-8 | bgzip -c \
  > HGDP.manta.$CNV.hg38.wAFs.sites.vcf.gz
done
/opt/sv-pipeline/05_annotation/scripts/compute_AFs.py \
  -p HGDP_pop_assignments.tsv \
  HGDP.GS.hg38.vcf.gz \
  stdout \
| cut -f1-8 | bgzip -c \
> HGDP.GS.hg38.wAFs.sites.vcf.gz


# Push to BED and reverse liftOver to hg19
cd /opt && \
  curl -o liftOver http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && \
  chmod a+x liftOver && \
  cd - && \
  curl -o hg38ToHg19.over.chain.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
for prefix in manta.DEL manta.DUP GS; do
  svtk vcf2bed --no-samples \
    HGDP.$prefix.hg38.wAFs.sites.vcf.gz \
    HGDP.$prefix.hg38.wAFs.sites.bed
  /opt/liftOver \
    -minMatch=0.5 \
    -bedPlus=4 \
    HGDP.$prefix.hg38.wAFs.sites.bed \
    hg38ToHg19.over.chain.gz \
    /dev/stdout \
    HGDP.$prefix.hg38.wAFs.sites.liftFails.bed \
  | sed 's/^chr//g' \
  > HGDP.$prefix.hg19.wAFs.sites.bed
done


# Replace coordinates of all hg38 records with hg19
gsutil -m cp gs://rcnv_project/refs/gnomad_v2.1_sv.nonneuro.sites.vcf.gz* ./
for prefix in manta.DEL manta.DUP GS; do
  cat <( zcat HGDP.$prefix.hg38.wAFs.sites.vcf.gz | grep '^##' | fgrep -v "##contig" ) \
      <( tabix -H gnomad_v2.1_sv.nonneuro.sites.vcf.gz | fgrep "##contig" ) \
      <( zcat HGDP.$prefix.hg38.wAFs.sites.vcf.gz | grep -e '^#' | grep -ve '^##' ) \
  > hg19.header.vcf
  /opt/rCNV2/data_curation/other/rewrite_vcf_coords.py \
    HGDP.$prefix.hg38.wAFs.sites.vcf.gz \
    HGDP.$prefix.hg19.wAFs.sites.bed \
    hg19.header.vcf \
    HGDP.$prefix.hg19.wAFs.unsorted.sites.vcf.gz
  bcftools sort \
    -o HGDP.$prefix.hg19.wAFs.sites.vcf.gz \
    -O z \
    HGDP.$prefix.hg19.wAFs.unsorted.sites.vcf.gz
  tabix -f HGDP.$prefix.hg19.wAFs.sites.vcf.gz
done


# Merge all records from HGDP into single VCF
bcftools concat \
  -a \
  -o HGDP.hg19.sites.vcf.gz \
  -O z \
  HGDP.manta.DEL.hg19.wAFs.sites.vcf.gz \
  HGDP.manta.DUP.hg19.wAFs.sites.vcf.gz \
  HGDP.GS.hg19.wAFs.sites.vcf.gz
tabix -f HGDP.hg19.sites.vcf.gz


# For convenience, also liftOver HGSV 1000 Genomes SV callset 
gsutil -m cp gs://dsmap/data/WGS/HGSV/unfiltered_sites_vcf/HGSV.WGS.wAFs.sites.vcf.gz* ./
svtk vcf2bed --no-samples \
  HGSV.WGS.wAFs.sites.vcf.gz \
  HGSV.WGS.wAFs.sites.bed
/opt/liftOver \
  -minMatch=0.5 \
  -bedPlus=4 \
  HGSV.WGS.wAFs.sites.bed \
  hg38ToHg19.over.chain.gz \
  /dev/stdout \
  HGSV.WGS.wAFs.sites.liftFails.bed \
| sed 's/^chr//g' \
> HGSV.WGS.wAFs.hg19.sites.bed
cat <( zcat HGSV.WGS.wAFs.sites.vcf.gz | grep '^##' | fgrep -v "##contig" ) \
    <( tabix -H gnomad_v2.1_sv.nonneuro.sites.vcf.gz | fgrep "##contig" ) \
    <( zcat HGSV.WGS.wAFs.sites.vcf.gz | grep -e '^#' | grep -ve '^##' ) \
> HGSV.hg19.header.vcf
/opt/rCNV2/data_curation/other/rewrite_vcf_coords.py \
  HGSV.WGS.wAFs.sites.vcf.gz \
  HGSV.WGS.wAFs.hg19.sites.bed \
  HGSV.hg19.header.vcf \
  HGSV.WGS.wAFs.hg19.unsorted.sites.vcf.gz
bcftools sort \
  -o 1000Genomes_HGSV_highCov.sites.vcf.gz \
  -O z \
  HGSV.WGS.wAFs.hg19.unsorted.sites.vcf.gz
tabix -f 1000Genomes_HGSV_highCov.sites.vcf.gz


# Copy all data to rCNV bucket (note: requires permissions)
gsutil -m cp \
  HGDP.hg19.sites.vcf.gz \
  HGDP.hg19.sites.vcf.gz.tbi \
  1000Genomes_HGSV_highCov.sites.vcf.gz \
  1000Genomes_HGSV_highCov.sites.vcf.gz.tbi \
  ${rCNV_bucket}/refs/
