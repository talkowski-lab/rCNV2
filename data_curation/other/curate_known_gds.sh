#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Curate files related to genomic disorder analyses for rCNV2 paper, including:
#   1. List of previously reported genomic disorder regions
#   2. List of predicted NAHR-mediated CNVs
#   3. List of regions with common (freq ≥ 1%) control CNVs


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"


# Download necessary reference files
mkdir refs
gsutil -m cp \
  ${rCNV_bucket}/refs/*bed.gz \
  ${rCNV_bucket}/refs/*vcf.gz* \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
  ${rCNV_bucket}/analysis/analysis_refs/rCNV_metacohort_list.txt \
  ${rCNV_bucket}/analysis/analysis_refs/HPOs_by_metacohort.table.tsv \
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


# Create common CNV blacklist from 1000Genomes, gnomAD, and CCDG
while read af suffix; do
  for CNV in DEL DUP; do
    athena vcf-filter \
      --minAF $af \
      --minAC 1 \
      --include-chroms $( seq 1 22 | paste -s -d, ) \
      --svtypes ${CNV},CNV,MCNV \
      --vcf-filters PASS,MULTIALLELIC \
      --af-field AF \
      --bgzip \
      refs/1000Genomes_phase3.sites.vcf.gz \
      1000Genomes_phase3.common_cnvs.${CNV}.$suffix.vcf.gz
    athena vcf-filter \
      --minAF $af \
      --minAC 1 \
      --include-chroms $( seq 1 22 | paste -s -d, ) \
      --svtypes ${CNV},CNV,MCNV \
      --vcf-filters PASS,MULTIALLELIC \
      --af-field AF \
      --bgzip \
      refs/gnomad_v2.1_sv.nonneuro.sites.vcf.gz \
      gnomAD.common_cnvs.${CNV}.$suffix.vcf.gz
    athena vcf-filter \
      --minAF $af \
      --minAC 1 \
      --include-chroms $( seq 1 22 | paste -s -d, ) \
      --svtypes ${CNV},CNV,MCNV \
      --vcf-filters PASS \
      --af-field AF \
      --bgzip \
      refs/CCDG_Abel_bioRxiv.sites.vcf.gz \
      CCDG.common_cnvs.${CNV}.$suffix.vcf.gz
    # Merge filtered VCFs and convert to BED
    /opt/rCNV2/data_curation/other/vcf2bed_merge.py \
      --genome refs/GRCh37.autosomes.genome \
      --outfile wgs_common_cnvs.${CNV}.$suffix.bed.gz \
      *.common_cnvs.${CNV}.$suffix.vcf.gz
  done
done < <( echo -e "0.01\t1pct\n0.001\t01pct" )

# Create common CNV blacklist from raw rCNV2 controls 
# (raw CNVs are necessary because curated rCNV2 controls are frequency-filtered)
mkdir raw_cnvs/
gsutil -m cp ${rCNV_bucket}/raw_data/cnv/*bed.gz raw_cnvs/
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/binned_genome/GRCh37.200kb_bins_10kb_steps.raw.bed.gz \
  refs/
# Step 1: count raw control CNVs per 200kb window per metacohort
for CNV in DEL DUP; do
  while read meta cohorts; do
    echo -e "$meta\t$CNV"
    # Pool control CNVs
    while read cohort; do
      zcat raw_cnvs/$cohort.raw.bed.gz \
      | fgrep -w "HEALTHY_CONTROL" \
      | fgrep -w $CNV
    done < <( echo -e "$cohorts" | sed 's/;/\n/g' ) \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bgzip -c \
    > $meta.raw_control_cnvs.$CNV.bed.gz
    # Count CNVs per window
    /opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
      --fraction 0.5 \
      -t $CNV \
      -o $meta.raw_control_cnv_counts.$CNV.bed.gz \
      --bgzip \
      $meta.raw_control_cnvs.$CNV.bed.gz \
      refs/GRCh37.200kb_bins_10kb_steps.raw.bed.gz
  done < <( fgrep meta refs/rCNV_metacohort_list.txt )
done
# Step 2: normalize CNV counts by control sample size and reduce to nonredundant intervals
while read af suffix; do
  for CNV in DEL DUP; do
    while read cohort; do
      cohort_idx=$( head -n1 refs/HPOs_by_metacohort.table.tsv | sed 's/\t/\n/g' \
                    | awk -v cohort=$cohort '{ if ($1==cohort) print NR }' )
      n_controls=$( awk -v FS="\t" -v idx=$cohort_idx \
                    '{ if ($1=="HEALTHY_CONTROL") print $idx }' \
                    refs/HPOs_by_metacohort.table.tsv )
      /opt/rCNV2/data_curation/other/get_common_control_cnv_regions.py \
        --n-controls $n_controls \
        --min-freq $af \
        --genome refs/GRCh37.autosomes.genome \
        $cohort.raw_control_cnv_counts.$CNV.bed.gz
    done < <( fgrep meta refs/rCNV_metacohort_list.txt | cut -f1 ) \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | grep -ve '^#' \
    | bedtools merge -i - \
    | cat <( echo -e "#chr\tstart\tend" ) - \
    | bgzip -c \
    > rCNV2_common_cnvs.${CNV}.$suffix.bed.gz
  done
done < <( echo -e "0.01\t1pct\n0.001\t01pct" )


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
# Prep input files
TAB=$( printf '\t' )
cat << EOF > cluster_gds.input.tsv
DECIPHER${TAB}/opt/rCNV2/refs/Decipher_GD.bed.gz
ClinGen${TAB}/ClinGen_GD.hmc.bed.gz
Owen${TAB}/opt/rCNV2/refs/UKBB_GD.Owen_2018.bed.gz
Girirajan${TAB}/opt/rCNV2/refs/Girirajan_2012_GD.bed.gz
Dittwald${TAB}/opt/rCNV2/refs/Dittwald_2013_GD.bed.gz
EOF
# af_suffix="01pct"
# for CNV in DEL DUP; do
#   zcat \
#     wgs_common_cnvs.$CNV.$af_suffix.bed.gz \
#     rCNV2_common_cnvs.$CNV.$af_suffix.bed.gz \
#   | grep -ve '^#' \
#   | cut -f1-3 \
#   | sort -Vk1,1 -k2,2n -k3,3n \
#   | bedtools merge -i -\
#   | bgzip -c \
#   > combined_common_cnvs.$CNV.$af_suffix.bed.gz
# done
# Clusters GDs (and formats them in preparation for gene annotation, below)
# Note: no longer apply control frequency filter (this can be handled with the 
# benign annotation in the final segments table later in analysis)
/opt/rCNV2/data_curation/other/cluster_gds.py \
  --hc-outfile lit_GDs.hc.no_genes.bed.gz \
  --mc-outfile lit_GDs.mc.no_genes.bed.gz \
  --lc-outfile lit_GDs.lc.no_genes.bed.gz \
  --hc-cutoff 4 \
  --mc-cutoff 2 \
  --lc-cutoff 1 \
  --genome refs/GRCh37.autosomes.genome \
  --segdups segdups.merged.10kb_slop.bed.gz \
  --cytobands refs/GRCh37.cytobands.bed.gz \
  --minsize 200000 \
  --maxsize 10000000 \
  --prep-for-gene-anno \
  --bgzip \
  cluster_gds.input.tsv


# Annotate GDs with genes
for conf in hc mc lc; do
  /opt/rCNV2/analysis/sliding_windows/get_genes_per_region.py \
      -o lit_GDs.${conf}.w_genes.bed \
      lit_GDs.${conf}.no_genes.bed.gz \
      refs/gencode.v19.canonical.pext_filtered.gtf.gz
  cut -f1-5,8-9 lit_GDs.${conf}.w_genes.bed \
  | bgzip -c > lit_GDs.${conf}.bed.gz
done


# Build comparison table of predicted NAHR-mediated CNVs
/opt/rCNV2/data_curation/other/predict_nahr_cnvs.sh


# Copy all files to analysis data bucket
gsutil -m cp \
  wgs_common_cnvs.*.bed.gz \
  rCNV2_common_cnvs.*.bed.gz \
  lit_GDs.hc.bed.gz \
  lit_GDs.mc.bed.gz \
  lit_GDs.lc.bed.gz \
  clustered_nahr_regions.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/

