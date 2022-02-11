#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Miscellaneous helper analyses for rCNV2 manuscript


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d2"
export finemap_dist=1000000


# Download necessary data (note: requires permissions)
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/cnv \
  ${rCNV_bucket}/analysis/analysis_refs/HPOs_by_metacohort.table.tsv \
  ${rCNV_bucket}/results/gene_association/* \
  ${rCNV_bucket}/results/segment_association/* \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/rCNV.*.gene_fine_mapping.gene_stats.merged_no_variation_features.tsv \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  ${rCNV_bucket}/analysis/gene_scoring/gene_lists \
  ${rCNV_bucket}/analysis/paper/data/global/${prefix}.global_burden_stats.tsv.gz \
  ./


# Count total number of CNVs
for CNV in DEL DUP; do
  zcat cnv/mega.rCNV.bed.gz | fgrep -w $CNV | wc -l
done


# Get total number of samples, cases. and controls
n_control=$( fgrep -w HEALTHY_CONTROL HPOs_by_metacohort.table.tsv | cut -f3 )
n_case=$( fgrep -w "HP:0000118" HPOs_by_metacohort.table.tsv | cut -f3 )
echo -e "$n_case\t$n_control" | awk -v OFS="\n" '{ print $1+$2, $1, $2 }'


# Compute total number of significant loci
for CNV in DEL DUP; do
  zcat rCNV.final_segments.loci.bed.gz \
       rCNV.final_genes.credible_sets.bed.gz \
  | fgrep -w $CNV \
  | cut -f1-3 \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  | wc -l
done
# (As above, but split by significance)
for CNV in DEL DUP; do
  # Genome- or exome-wide
  zcat rCNV.final_segments.loci.bed.gz \
       rCNV.final_genes.credible_sets.bed.gz \
  | fgrep -w $CNV | grep -e '[genome|exome]_wide' \
  | cut -f1-3 \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  > gw_ew.$CNV.bed
  cat gw_ew.$CNV.bed | wc -l
  # FDR
  zcat rCNV.final_segments.loci.bed.gz \
       rCNV.final_genes.credible_sets.bed.gz \
  | fgrep -w $CNV | fgrep -w FDR \
  | cut -f1-3 \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  | bedtools intersect -v -a - -b gw_ew.$CNV.bed \
  | wc -l
  rm gw_ew.$CNV.bed
done


# Compute total number of significant associations
for CNV in DEL DUP; do
  zcat rCNV.final_segments.associations.bed.gz \
       rCNV.final_genes.credible_sets.bed.gz \
  | fgrep -w $CNV \
  | awk -v FS="\t" -v OFS="\t" '{ print $6"_"$1, $2, $3 }' \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  | wc -l
done
# (As above, but split by significance)
for CNV in DEL DUP; do
  # Genome- or exome-wide
  zcat rCNV.final_segments.associations.bed.gz \
       rCNV.final_genes.credible_sets.bed.gz \
  | fgrep -w $CNV | grep -e '[genome|exome]_wide' \
  | awk -v FS="\t" -v OFS="\t" '{ print $6"_"$1, $2, $3 }' \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  > gw_ew.$CNV.bed
  cat gw_ew.$CNV.bed | wc -l
  # FDR
  zcat rCNV.final_segments.associations.bed.gz \
       rCNV.final_genes.credible_sets.bed.gz \
  | fgrep -w $CNV | fgrep -w FDR \
  | awk -v FS="\t" -v OFS="\t" '{ print $6"_"$1, $2, $3 }' \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  | bedtools intersect -v -a - -b gw_ew.$CNV.bed \
  | wc -l
  rm gw_ew.$CNV.bed
done


# Get number of candidate novel TS disease genes
zcat rCNV.final_genes.genes.bed.gz | fgrep -w DUP | cut -f4 \
| fgrep -wvf gene_lists/DDG2P.all_gof.genes.list \
| fgrep -wvf gene_lists/ClinGen.all_triplosensitive.genes.list \
| sort | uniq | wc -l



# Get total number of fine-mapped candidate driver genes by CNV type
zcat rCNV.final_genes.genes.bed.gz \
| fgrep -v "#" | cut -f4 | sort | uniq | wc -l 


# Get gold-standard dosage-sensitive genes & PIPs to higlight in bottom of graphical abstract
cat \
  <( zcat rCNV.final_genes.genes.bed.gz | fgrep -w DEL \
     | fgrep -wf  gene_lists/gold_standard.haploinsufficient.genes.list ) \
  <( zcat rCNV.final_genes.genes.bed.gz | fgrep -w DUP \
     | fgrep -wf  gene_lists/gold_standard.triplosensitive.genes.list ) \
| cut -f4,5,12 | sort -nrk3,3 | head -n10


# Get constrained genes & PIPs to higlight in bottom of graphical abstract
cat gene_lists/gnomad.v2.1.1.*_constrained.genes.list | sort | uniq \
| fgrep -wf - <( zcat rCNV.final_genes.genes.bed.gz ) | cut -f4,5,12 \
| sort -nrk3,3 | head -n10


# Get median & IQR of GD effect sizes across phenotypes
# (Note: category 3 = GDs & 7 = constrained genes outside of GDs)
for categ in 3 7; do
  cat << EOF > get_OR_distrib.cat$categ.R
x <- read.table("${prefix}.global_burden_stats.tsv.gz", header=T, comment.char="")
round(summary(exp(x[which(x\$category==$categ & x\$CNV=="CNV"), "meta_lnOR"])), 2)
EOF
Rscript get_OR_distrib.cat$categ.R
rm get_OR_distrib.cat$categ.R
done

# Get interquartile range of CNV sizes
Rscript -e "summary(abs(apply(read.table('cnv/mega.rCNV.bed.gz', \
            header=T, sep=\"\\\t\", comment.char=\"\")[, 2:3], 1, diff)))"

