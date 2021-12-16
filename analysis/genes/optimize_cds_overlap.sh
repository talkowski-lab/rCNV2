#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Optimization of CDS overlap fractions for gene-based rCNV association testing


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv:latest
gcloud auth login


# Set parameters
export rCNV_bucket="gs://rcnv_project"


# Copy all filtered CNV data and other references from the project Google Bucket (note: requires permissions)
mkdir cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
mkdir refs/
gsutil -m cp -r \
  gs://rcnv_project/refs/GRCh37.*.bed.gz \
  gs://rcnv_project/analysis/analysis_refs/* \
  gs://rcnv_project/cleaned_data/genes/gene_lists \
  gs://rcnv_project/cleaned_data/genes/*gtf* \
  gs://rcnv_project/analysis/paper/data/large_segments/clustered_nahr_regions.bed.gz \
  refs/
gsutil -m cp \
  gs://rcnv_project/analysis/gene_scoring/gene_lists/* \
  refs/gene_lists/


# Note: must also localize probe-based conditional exclusion BED, either from
# FireCloud or regenerated locally (see gene_burden_analysis.sh)
export exclusion_bed=gencode.v19.canonical.cohort_exclusion.bed.gz


# Subset GTF to genes of interest for deletions & duplications
opt/rCNV2/utils/filter_gtf_by_genelist.py \
  -o gencode.v19.canonical.pext_filtered.DEL.gtf.gz \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  refs/gene_lists/gold_standard.haploinsufficient.genes.list
opt/rCNV2/utils/filter_gtf_by_genelist.py \
  -o gencode.v19.canonical.pext_filtered.DUP.gtf.gz \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  refs/gene_lists/gold_standard.triplosensitive.genes.list


# Prepare CNV type-specific exclusion lists of NAHR-mediated GDs
for CNV in DEL DUP; do
  bedtools intersect -u -r -f 0.25 \
    -a /opt/rCNV2/refs/lit_GDs.all.$CNV.bed.gz \
    -b refs/clustered_nahr_regions.bed.gz \
  | bgzip -c \
  > refs/lit_GDs.all.$CNV.NAHR_only.bed.gz
done


# Annotate gene & CNV overlap for developmental HPOs while excluding known NAHR-mediated GDs
# Also apply cohort-specific probe-based conditional exclusion at this step (rather than meta step)
mkdir cnv_counts
while read meta cohorts; do
  # Extract genes to be excluded from this cohort
  zcat ${exclusion_bed} | fgrep -w $meta | cut -f1-3 \
  | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | bgzip -c \
  > conditional_exclusion.$meta.bed.gz
  # Count deletions & duplications separately
  for CNV in DEL DUP; do
    /opt/rCNV2/analysis/genes/count_cnvs_per_gene.py \
      --min-cds-ovr 0 \
      -t $CNV \
      --hpo "$( paste -d\; -s refs/rCNV2.hpos_by_severity.developmental.list )" \
      --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
      --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
      --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
      --blacklist refs/lit_GDs.all.$CNV.NAHR_only.bed.gz \
      --blacklist conditional_exclusion.$meta.bed.gz \
      -z \
      --verbose \
      -o /dev/null \
      --cnvs-out cnv_counts/${meta}.${CNV}.genes_per_cnv.bed.gz \
      --annotate-cds-per-gene \
      cleaned_cnv/${meta}.rCNV.bed.gz \
      gencode.v19.canonical.pext_filtered.${CNV}.gtf.gz
  done
done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt )


# Summarize counts & effect size per cohort at various CDS cutoffs
mkdir cds_optimization_data
while read meta cohorts; do
  ncase=$( fgrep -w $meta refs/rCNV2.hpos_by_severity.developmental.counts.tsv | cut -f2 )
  cidx=$( head -n1 refs/HPOs_by_metacohort.table.tsv | sed 's/\t/\n/g' | awk -v meta="$meta" '{ if ($1==meta) print NR }' )
  nctrl=$( awk -v FS="\t" -v cidx=$cidx '{ if ($1=="HEALTHY_CONTROL") print $(cidx) }' refs/HPOs_by_metacohort.table.tsv )
  for CNV in DEL DUP; do
    /opt/rCNV2/analysis/other/count_cnvs_by_cds.R \
      --n-cases $ncase \
      --n-controls $nctrl \
      cnv_counts/${meta}.${CNV}.genes_per_cnv.bed.gz \
      cds_optimization_data/${meta}.${CNV}.counts_per_cds.tsv
  done
done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt )


# Meta-analyze CNV counts across cohorts at each CDS cutoff
for CNV in DEL DUP; do
  while read meta cohorts; do
    echo -e "$meta\tcds_optimization_data/${meta}.${CNV}.counts_per_cds.tsv"
  done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt ) \
  > ${CNV}.meta_analysis.input.txt
  /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
    --model "fe" \
    --keep-n-columns 4 \
    --p-is-neg-log10 \
    --min-cases 300 \
    ${CNV}.meta_analysis.input.txt \
    /dev/stdout \
  | cut -f4- | sed -e 's/^mincds_//g' -e 's/^min_cds/#min_cds/g' \
  | awk -v FS="\t" '{ if ($10!="") print $0 }' \
  > cds_optimization_data/${CNV}.cds_optimization.meta_analysis.stats.tsv
done


# Optimize CDS cutoffs
/opt/rCNV2/analysis/other/optimize_min_cds.R \
  --optimize-power \
  cds_optimization_data/DEL.cds_optimization.meta_analysis.stats.tsv \
  cds_optimization_data/DEL
/opt/rCNV2/analysis/other/optimize_min_cds.R \
  cds_optimization_data/DUP.cds_optimization.meta_analysis.stats.tsv \
  cds_optimization_data/DUP


# Copy optimization data to Google bucket
gsutil -m cp -r \
  cds_optimization_data \
  gs://rcnv_project/analysis/gene_burden/
