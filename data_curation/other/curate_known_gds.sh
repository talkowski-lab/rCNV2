#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Curate list of previously reported genomic disorder regions for rCNV2 analyses


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"


# Download necessary reference files
mkdir refs
gsutil -m cp \
  ${rCNV_bucket}/refs/*bed.gz \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
  refs/
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz


# Reformat segdup file to collapse overlapping intervals
zcat genomicSuperDups.txt.gz \
| cut -f2-4 \
| sed 's/^chr//g' \
| grep -e '^[0-9]' \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
| bgzip -c \
> segdups.merged.bed.gz
bedtools merge -d 10000 -i segdups.merged.bed.gz \
| bgzip -c \
> segdups.merged.10kb_slop.bed.gz 


# Download & parse ClinGen CNV regions
wget ftp://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh37.tsv
/opt/rCNV2/data_curation/other/parse_clingen_regions.R \
  ClinGen_region_curation_list_GRCh37.tsv \
  ./clingen_regions
for conf in hmc all; do
  cat <( grep -e '^#' ./clingen_regions.DEL_GD_${conf}.bed ) \
      <( cat ./clingen_regions.*_GD_${conf}.bed | fgrep -v "#" | grep -e '^[0-9]' \
         | sort -Vk1,1 -k2,2n -k3,3n -k5,5V ) \
  | bgzip -c \
  > ClinGen_GD.${conf}.bed.gz
done


# Integrate DECIPHER, ClinGen, Owen, Girirajan, and Dittwald
TAB=$( printf '\t' )
cat << EOF > cluster_gds.input.tsv
DECIPHER${TAB}/opt/rCNV2/refs/Decipher_GD.bed.gz
ClinGen${TAB}/ClinGen_GD.hmc.bed.gz
Owen${TAB}/opt/rCNV2/refs/UKBB_GD.Owen_2018.bed.gz
Girirajan${TAB}/opt/rCNV2/refs/Girirajan_2012_GD.bed.gz
Dittwald${TAB}/opt/rCNV2/refs/Dittwald_2013_GD.bed.gz
EOF
# Clusters GDs (and formats them in preparation for gene annotation, below)
/opt/rCNV2/data_curation/other/cluster_gds.py \
  --hc-outfile lit_GDs.hc.no_genes.bed.gz \
  --lc-outfile lit_GDs.lc.no_genes.bed.gz \
  --hc-cutoff 4 \
  --lc-cutoff 2 \
  --genome refs/GRCh37.autosomes.genome \
  --segdups segdups.merged.10kb_slop.bed.gz \
  --cytobands refs/GRCh37.cytobands.bed.gz \
  --minsize 200000 \
  --maxsize 10000000 \
  --prep-for-gene-anno \
  --bgzip \
  cluster_gds.input.tsv


# Annotate GDs with genes
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  ./
for conf in hc lc; do
  /opt/rCNV2/analysis/sliding_windows/get_genes_per_region.py \
      -o lit_GDs.${conf}.w_genes.bed \
      lit_GDs.${conf}.no_genes.bed.gz \
      gencode.v19.canonical.pext_filtered.gtf.gz
  cut -f1-5,8-9 lit_GDs.${conf}.w_genes.bed \
  | bgzip -c > lit_GDs.${conf}.bed.gz
done


# Build comparison table of predicted NAHR-mediated CNVs
/opt/rCNV2/data_curation/other/predict_nahr_cnvs.sh


# Copy all files to analysis data bucket
gsutil -m cp \
  lit_GDs.hc.bed.gz \
  lit_GDs.lc.bed.gz \
  clustered_nahr_regions.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/

