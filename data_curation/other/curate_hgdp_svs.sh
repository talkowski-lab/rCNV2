#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Curate the HGDP SV callset for use in rCNV2 analyses


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
# TODO: ADD GENOME STRIP CURATION


# Annotate frequencies by continental population
sed '1d' HGDP.metadata | cut -f1,8 > HGDP_pop_assignments.tsv
for CNV in DEL DUP; do
  /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py \
    -p HGDP_pop_assignments.tsv \
    HGDP.manta.$CNV.hg38.vcf.gz \
    stdout \
  | cut -f1-8 | bgzip -c \
  > HGDP.manta.$CNV.hg38.wAFs.sites.vcf.gz
done




# Push to BED and reverse liftOver to hg19


# Replace coordinates of all records


