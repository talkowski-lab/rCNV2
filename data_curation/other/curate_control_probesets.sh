#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Curate SNP array probesets used for all control samples included in rCNV2 analyses


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"


# Prep output directory
mkdir control_probesets


# Localize all raw probe data stored in rCNV gs:// bucket (note: requires permissions)
gsutil -m cp \
  gs://rcnv_project/raw_data/other/200K_GSAMD.probes.txt.gz \
  gs://rcnv_project/raw_data/other/Axiom_BioBank1-na35-bed.zip \
  gs://rcnv_project/raw_data/other/CytoScanHD_Array-na32-3-bed.zip \
  gs://rcnv_project/raw_data/other/GenomeDx_v5.1_086081_NoSNP.txt \
  gs://rcnv_project/raw_data/other/031780_PrenatalSNPv1_no_SNP.txt \
  ./


# Curate probesets stored in rCNV gs:// bucket
# GSA
zcat 200K_GSAMD.probes.txt.gz | sed '1d' \
| awk -v OFS="\t" '{ print $2, $3-1, $3 }' \
| sort -Vk1,1 -k2,2n -k3,3n | bgzip -c \
> control_probesets/illumina_gsa.bed.gz
# Axiom
unzip Axiom_BioBank1-na35-bed.zip
fgrep -v "#" Axiom_BioBank1.na35.bed | sed '1d' | cut -f1-3 \
| sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | uniq | bgzip -c \
> control_probesets/ukbbaxiom.bed.gz
# CytoScanHD
unzip CytoScanHD_Array-na32-3-bed.zip
sed '1d' CytoScanHD_Array.na32.3.bed | cut -f1-3 \
| sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | uniq | bgzip -c \
> control_probesets/affy_cyto_hd.bed.gz
# GDX_Agilent_60k
sed 's/^chr//g' 031780_PrenatalSNPv1_no_SNP.txt | cut -f1-3 \
| sort -Vk1,1 -k2,2n -k3,3n | uniq | bgzip -c \
> control_probesets/agilent_gdx_60k.bed.gz
# GDX_Agilent_180k
sed 's/^chr//g' GenomeDx_v5.1_086081_NoSNP.txt | cut -f1-3 \
| sort -Vk1,1 -k2,2n -k3,3n | uniq | bgzip -c \
> control_probesets/agilent_gdx_180k.bed.gz


# Download & clean probesets hosted by Illumina
# 610k Quad
wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/human610/human610-quadv1_h.zip
unzip human610-quadv1_h.zip
sed '1d' Human610-Quadv1_H.bed \
| cut -f1-3 | sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | uniq | bgzip -c \
> control_probesets/illumina_610k_quad.bed.gz
# Omni Express
wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanomniexpress-24/v1-3/infinium-omniexpress-24-v1-3-a1-bed.zip
unzip infinium-omniexpress-24-v1-3-a1-bed.zip
sed '1d' InfiniumOmniExpress-24v1-3_A1.bed \
| cut -f1-3 | sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | uniq | bgzip -c \
> control_probesets/omniexpress.bed.gz
# Omni 2.5
wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanomni25/v1-5/infinium-omni2-5-8v1-5-a1-bed.zip
unzip infinium-omni2-5-8v1-5-a1-bed.zip
sed '1d' InfiniumOmni2-5-8v1-5_A1.bed \
| cut -f1-3 | sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | uniq | bgzip -c \
> control_probesets/omni_2.5.bed.gz
# MEGA
wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/multiethnic-global/multi-ethnic-global-8-d1-bed.zip
unzip multi-ethnic-global-8-d1-bed.zip
sed '1d' Multi-EthnicGlobal_D1.bed \
| cut -f1-3 | sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | uniq | bgzip -c \
> control_probesets/illumina_mega.bed.gz


# Download all array probesets stored as UCSC tables
for table in snpArrayAffy250Nsp snpArrayAffy250Sty snpArrayAffy5 snpArrayAffy6 \
             snpArrayIllumina300 snpArrayIllumina550 snpArrayIllumina650 \
             snpArrayIllumina1M; do
  wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/$table.txt.gz
done


# Clean all array probesets downloaded from UCSC
# Affy 500k
zcat snpArrayAffy250Nsp.txt.gz snpArrayAffy250Sty.txt.gz \
| cut -f2-4 | sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | uniq | bgzip -c \
> control_probesets/affy_500k.bed.gz
# All others
while read table prefix; do
  zcat $table.txt.gz \
  | cut -f2-4 | sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | uniq | bgzip -c \
  > control_probesets/$prefix.bed.gz
done < <( echo -e "snpArrayAffy5\taffy_5.0\nsnpArrayAffy6\taffy_6.0\nsnpArrayIllumina300\tillumina_300k\nsnpArrayIllumina550\tillumina_550k\nsnpArrayIllumina650\tillumina_650k\nsnpArrayIllumina1M\tillumina_1m_duo" )


# Tabix all curated probesets
for file in control_probesets/*bed.gz; do
  tabix -p bed -f $file
done


# Copy curated probesets to rCNV gs:// bucket for storage
gsutil -m cp -r control_probesets/* gs://rcnv_project/cleaned_data/control_probesets/

