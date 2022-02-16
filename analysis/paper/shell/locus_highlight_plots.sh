#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block to generate locus-level higlight plots for rCNV2 manuscript


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d2"


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
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/rCNV.*.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv \
  ${rCNV_bucket}/results/* \
  ./


# Download reference files (note: requires permissions)
mkdir refs/ 
gsutil -m cp \
  ${rCNV_bucket}/analysis/analysis_refs/** \
  ${rCNV_bucket}/analysis/paper/data/hpo/${prefix}.reordered_hpos.txt \
  ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  ${rCNV_bucket}/refs/REP_state_manifest.tsv \
  refs/


# Prepare generic input files
while read cohort; do
  path="cnv/$cohort.rCNV.bed.gz"
  echo -e "$cohort\t$path"
done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt | cut -f1 ) \
> cnvs.input.tsv


# Get genome-wide significance thresholds
example_hpo="HP0012759"
gw_cutoff=$( awk -v FS="\t" -v hpo=${example_hpo} '{ if ($1==hpo) print $2 }' \
             refs/sliding_window.rCNV.DEL.bonferroni_pval.hpo_cutoffs.tsv )


# Simple dev case for new script - SHANK3 deletions
Rscript -e "install.packages('opt/rCNV2/source/rCNV2_0.1.0.tar.gz', repos=NULL)" && \
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_locus_highlight.R \
  --case-hpos "HP:0000707" \
  --highlights "22:50166930-51171641" \
  --sumstats meta_stats/HP0012759.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --pips rCNV.DEL.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv \
  --gw-sig-label-side bottom \
  --standardize-frequencies \
  --collapse-cohorts \
  --cnv-panel-height 1.2 \
  --pdf-height 3.5 \
  22:49310000-51270000 \
  cnvs.input.tsv \
  "DEL" \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  locus_highlight_test



#####################
#  CBLN2 Deletions  #
#####################
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_CBLN2_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --gw-sig "$gw_cutoff" \
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
  --gw-sig "$gw_cutoff" \
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
  --gw-sig "$gw_cutoff" \
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
  --gw-sig "$gw_cutoff" \
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
  --gw-sig "$gw_cutoff" \
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
  --gw-sig "$gw_cutoff" \
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
  --gw-sig "$gw_cutoff" \
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
  --gw-sig "$gw_cutoff" \
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
  --gw-sig "$gw_cutoff" \
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
  --gw-sig "$gw_cutoff" \
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
tabix -p bed -f CADM2.fetal_cortex_RNAseq.3kb_bins.bed.gz
# Extract all CADM2 transcripts
curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz \
| gunzip -c | grep -e '^chr3' | fgrep -w CADM2 | sed 's/^chr//g' | sort -Vk1,1 -k4,4n -k5,5n | bgzip -c \
> gencode.v19.annotation.CADM2.gtf.gz
tabix -p gff -f gencode.v19.annotation.CADM2.gtf.gz
# Prep ChromHMM tracks
for eid in E067 E069 E072 E073; do
  while read number abbrev junk; do
    sid="${number}_$( echo $abbrev | sed 's/\///g' )"
    gsutil -m cat \
      ${rCNV_bucket}/cleaned_data/genome_annotations/chromhmm_beds/roadmap_chromhmm.$eid.$sid.bed.gz \
    | gunzip -c | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - \
    | awk -v number=$number -v OFS="\t" '{ print $1, $2, $3, number }'
  done < <( fgrep -v "#" refs/REP_state_manifest.tsv ) \
  | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
  | cat <( echo -e "#chrom\tstart\tend\tstate" ) - \
  | bgzip -c \
  > $eid.chromhmm.bed.gz
  tabix -p bed -f $eid.chromhmm.bed.gz
  echo -e "$eid.chromhmm.bed.gz"
done > CADM2.chromhmm_paths.tsv
# Plot locus
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_CADM2_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf gencode.v19.annotation.CADM2.gtf.gz \
  --rnaseq CADM2.fetal_cortex_RNAseq.3kb_bins.bed.gz \
  --chromhmm-tracks CADM2.chromhmm_paths.tsv \
  --chromhmm-manifest refs/REP_state_manifest.tsv \
  cnvs.input.tsv \
  meta_stats/HP0000752.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}


########################
#  10p14 Duplications  #
########################
# Prep lincRNAs
curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz \
| gunzip -c | grep -e '^chr10' | sed 's/^chr//g' | fgrep -v protein_coding | sort -Vk1,1 -k4,4n -k5,5n | bgzip -c \
> gencode.v19.annotation.chr10_ncRNAs.gtf.gz
tabix -p gff -f gencode.v19.annotation.chr10_ncRNAs.gtf.gz
# Generate 1kb bins near 10p14 locus
athena make-bins \
  --exclude-chroms $( cut -f1 refs/GRCh37.genome | fgrep -wv 10 | paste -s -d, ) \
  --bgzip \
  refs/GRCh37.genome \
  1000 \
  chr10.1kb_bins.bed.gz
bedtools intersect -header -wa \
  -a chr10.1kb_bins.bed.gz \
  -b <( echo -e "10\t5520000\t8280000" ) \
| awk '{ print "chr"$0 }' \
| bgzip -c \
> 10p14.1kb_bins.bed.gz
# Prep ENCODE H3K27ac ChIP-seq data from ESC lines
wget https://www.encodeproject.org/files/ENCFF661NZB/@@download/ENCFF661NZB.bam #H1
samtools index ENCFF661NZB.bam
wget https://www.encodeproject.org/files/ENCFF930EEK/@@download/ENCFF930EEK.bam #H9
samtools index ENCFF930EEK.bam
# wget https://www.encodeproject.org/files/ENCFF156LTI/@@download/ENCFF156LTI.bam #HUES48
# samtools index ENCFF156LTI.bam
# wget https://www.encodeproject.org/files/ENCFF326JSN/@@download/ENCFF326JSN.bam #HUES6
# samtools index ENCFF326JSN.bam
# wget https://www.encodeproject.org/files/ENCFF599FVC/@@download/ENCFF599FVC.bam #HUES64
# samtools index ENCFF599FVC.bam
# wget https://www.encodeproject.org/files/ENCFF465XYF/@@download/ENCFF465XYF.bam #female
# samtools index ENCFF465XYF.bam
# wget https://www.encodeproject.org/files/ENCFF042WIR/@@download/ENCFF042WIR.bam #male
# samtools index ENCFF042WIR.bam
bedtools coverage -counts \
  -a 10p14.1kb_bins.bed.gz \
  -b <( samtools view -b ENCFF661NZB.bam chr10 ) \
| bedtools coverage -counts \
  -a - \
  -b <( samtools view -b ENCFF930EEK.bam chr10 ) \
| awk -v OFS="\t" '{ if ($4==".") $4=0; print }' \
| awk -v OFS="\t" '{ if ($5==".") $5=0; print }' \
| awk -v OFS="\t" '{ print $1, $2, $3, $4, $5, ($4+$5)/2 }' \
| sort -Vk1,1 -k2,2n -k3,3n \
| sed 's/^chr//g' \
| cat <( echo -e "#chr\tstart\tend\tH1_cov\tH9_cov\tavg_cov" ) - \
| bgzip -c \
> 10p14.ESC_H3k27ac.1kb_bins.bed.gz
tabix -p bed -f 10p14.ESC_H3k27ac.1kb_bins.bed.gz
# Prep ESC ChromHMM tracks
for eid in E008 E015 E014 E016 E003; do
  while read number abbrev junk; do
    sid="${number}_$( echo $abbrev | sed 's/\///g' )"
    gsutil -m cat \
      ${rCNV_bucket}/cleaned_data/genome_annotations/chromhmm_beds/roadmap_chromhmm.$eid.$sid.bed.gz \
    | gunzip -c | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - \
    | awk -v number=$number -v OFS="\t" '{ print $1, $2, $3, number }'
  done < <( fgrep -v "#" refs/REP_state_manifest.tsv ) \
  | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
  | cat <( echo -e "#chrom\tstart\tend\tstate" ) - \
  | bgzip -c \
  > $eid.chromhmm.bed.gz
  tabix -p bed -f $eid.chromhmm.bed.gz
  echo -e "$eid.chromhmm.bed.gz"
done > 10p14.chromhmm_paths.tsv
# Plot locus
/opt/rCNV2/analysis/paper/plot/locus_highlights/plot_10p14_locus.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --gtf refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  --ncrna-gtf gencode.v19.annotation.chr10_ncRNAs.gtf.gz \
  --chipseq 10p14.ESC_H3k27ac.1kb_bins.bed.gz \
  --esc-chromhmm-tracks 10p14.chromhmm_paths.tsv \
  --adult-chromhmm-tracks CADM2.chromhmm_paths.tsv \
  --chromhmm-manifest refs/REP_state_manifest.tsv \
  cnvs.input.tsv \
  meta_stats/HP0012759.rCNV.DUP.sliding_window.meta_analysis.stats.bed.gz \
  refs/HPOs_by_metacohort.table.tsv \
  refs/GRCh37.genome \
  ${prefix}


# Copy all plots to final gs:// directory
gsutil -m cp -r \
  ${prefix}.locus_highlight.*.pdf \
  ${rCNV_bucket}/analysis/paper/plots/locus_highlights/

