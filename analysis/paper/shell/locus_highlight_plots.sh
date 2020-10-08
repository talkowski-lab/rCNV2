#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block to generate locus-level higlight plots for rCNV2 manuscript


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"


# Download necessary data (note: requires permissions)
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/cnv ./
mkdir meta_stats/
gsutil -m cp \
  ${rCNV_bucket}/analysis/sliding_windows/**.rCNV.**.sliding_window.meta_analysis.stats.bed.gz \
  ${rCNV_bucket}/analysis/gene_burden/**.rCNV.**.gene_burden.meta_analysis.stats.bed.gz \
  ${rCNV_bucket}/analysis/crb_burden/**.rCNV.**.crb_burden.meta_analysis.stats.bed.gz \
  meta_stats/
find meta_stats/ -name "*bed.gz" | xargs -I {} tabix -p bed -f {}
gsutil -m cp -r \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/rCNV.*.gene_fine_mapping.gene_stats.merged_no_variation_features.all_genes_from_blocks.tsv \
  ./


# Download reference files (note: requires permissions)
mkdir refs/ 
gsutil -m cp \
  ${rCNV_bucket}/analysis/analysis_refs/** \
  ${rCNV_bucket}/analysis/paper/data/hpo/${prefix}.reordered_hpos.txt \
  ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  refs/


# Prepare generic input files
while read cohort; do
  path="cnv/$cohort.rCNV.bed.gz"
  echo -e "$cohort\t$path"
done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt | cut -f1 ) \
> cnvs.input.tsv


#####################
#  CBLN2 Deletions  #
#####################
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_CBLN2_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  cnvs.input.tsv \
  meta_stats/HP0012759.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}



####################
#  IER5 Deletions  #
####################
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_IER5_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  cnvs.input.tsv \
  meta_stats/HP0012759.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}


#######################
#  RAF1 Duplications  #
#######################
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_RAF1_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  cnvs.input.tsv \
  meta_stats/UNKNOWN.rCNV.DUP.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}


#######################
#  1q44 Duplications  #
#######################
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_1q44_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  cnvs.input.tsv \
  meta_stats/HP0001626.rCNV.DUP.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}


########################
#  GMEB2 Duplications  #
########################
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_GMEB2_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --pips rCNV.DUP.gene_fine_mapping.gene_stats.merged_no_variation_features.all_genes_from_blocks.tsv \
  cnvs.input.tsv \
  meta_stats/HP0012639.rCNV.DUP.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}


##########################
#  ANKRD11 Duplications  #
##########################
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_ANKRD11_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --pips rCNV.DUP.gene_fine_mapping.gene_stats.merged_no_variation_features.all_genes_from_blocks.tsv \
  cnvs.input.tsv \
  meta_stats/HP0001507.rCNV.DUP.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}


######################
#  SHANK3 Deletions  #
######################
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_SHANK3_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --pips rCNV.DEL.gene_fine_mapping.gene_stats.merged_no_variation_features.all_genes_from_blocks.tsv \
  cnvs.input.tsv \
  meta_stats/HP0012759.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}


#######################
#  SMARCA2 Deletions  #
#######################
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_SMARCA2_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  cnvs.input.tsv \
  meta_stats/HP0012759.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}


###################
#  QKI Deletions  #
###################
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_QKI_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --pips rCNV.DEL.gene_fine_mapping.gene_stats.merged_no_variation_features.all_genes_from_blocks.tsv \
  cnvs.input.tsv \
  meta_stats/UNKNOWN.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}


######################
#  SLC2A3 Deletions  #
######################
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_SLC2A3_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --pips rCNV.DEL.gene_fine_mapping.gene_stats.merged_no_variation_features.all_genes_from_blocks.tsv \
  cnvs.input.tsv \
  meta_stats/HP0001250.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}


#####################
#  CADM2 Deletions  #
#####################
# Generate 3kb bins around CADM2
athena make-bins \
  --exclude-chroms $( cut -f1 refs/GRCh37.genome | fgrep -wv 3 | paste -s -d, ) \
  --bgzip \
  refs/GRCh37.genome \
  3000 \
  chr3.3kb_bins.bed.gz
bedtools intersect -header -wa \
  -a chr3.3kb_bins.bed.gz \
  -b <( echo -e "3\t83892686\t87239026" ) \
| awk '{ print "chr"$0 }' \
| bgzip -c \
> CADM2.3kb_bins.bed.gz
# Download & process fetal cortex RNAseq data from ENCODE
wget https://www.encodeproject.org/files/ENCFF862KEW/@@download/ENCFF862KEW.bigWig
bigWigToBedGraph \
  -chrom=chr3 -start=83892686 -end=87239026 \
  ENCFF862KEW.bigWig \
  ENCFF862KEW.CADM2.bg
wget https://www.encodeproject.org/files/ENCFF572HVL/@@download/ENCFF572HVL.bigWig
bigWigToBedGraph \
  -chrom=chr3 -start=83892686 -end=87239026 \
  ENCFF572HVL.bigWig \
  ENCFF572HVL.CADM2.bg
bedtools map -c 4 -o mean \
  -a CADM2.3kb_bins.bed.gz \
  -b ENCFF862KEW.CADM2.bg \
| bedtools map -c 4 -o mean \
  -a - \
  -b ENCFF572HVL.CADM2.bg \
| awk -v OFS="\t" '{ if ($4==".") $4=0; print }' \
| awk -v OFS="\t" '{ if ($5==".") $5=0; print }' \
| awk -v OFS="\t" '{ print $1, $2, $3, $4, $5, ($4+$5)/2 }' \
| sort -Vk1,1 -k2,2n -k3,3n \
| sed 's/^chr//g' \
| cat <( echo -e "#chr\tstart\tend\tfemale_signal\tmale_signal\tavg_signal" ) - \
| bgzip -c \
> CADM2.fetal_cortex_RNAseq.3kb_bins.bed.gz
# Extract all CADM2 transcripts
curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz \
| gunzip -c | grep -e '^chr3' | fgrep -w CADM2 | sed 's/^chr//g' | sort -Vk1,1 -k4,4n -k5,5n | bgzip -c \
> gencode.v19.annotation.CADM2.gtf.gz
tabix -p gff -f gencode.v19.annotation.CADM2.gtf.gz
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_CADM2_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --rnaseq CADM2.fetal_cortex_RNAseq.3kb_bins.bed.gz \
  cnvs.input.tsv \
  meta_stats/HP0001250.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}


# Copy all plots to final gs:// directory
gsutil -m cp -r \
  ${prefix}.locus_highlight.*.pdf \
  ${rCNV_bucket}/analysis/paper/plots/locus_highlights/

