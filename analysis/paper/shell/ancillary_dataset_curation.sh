#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Miscellaneous preprocessing tasks of ancillary datasets for rCNV2 formal analysis


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"


# Download necessary reference files and metadata
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
  ${rCNV_bucket}/raw_data/other/satterstrom_asc_dnms.raw.tsv.gz \
  ${rCNV_bucket}/raw_data/other/redin_2017.bca_breakpoints.all.tsv.gz \
  ${rCNV_bucket}/analysis/analysis_refs/GRCh37.genome \
  ${rCNV_bucket}/refs/gnomad_v2.1_sv.nonneuro.sites.vcf* \
  ${rCNV_bucket}/raw_data/other/asc_spark_denovo_cnvs \
  ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs.*.bed.gz \
  ./


# Download & reformat DDD de novo mutation table
wget https://www.biorxiv.org/content/biorxiv/early/2020/04/01/797787/DC4/embed/media-4.txt
/opt/rCNV2/data_curation/other/curate_ddd_dnms.py \
  --dnm-tsv media-4.txt \
  --genes gencode.v19.canonical.pext_filtered.genes.list \
  -o ddd_dnm_counts.tsv.gz \
  -z


# Reformat ASC de novo mutation table
/opt/rCNV2/data_curation/other/curate_asc_dnms.py \
  --dnm-tsv satterstrom_asc_dnms.raw.tsv.gz \
  --genes gencode.v19.canonical.pext_filtered.genes.list \
  -o asc_dnm_counts.tsv.gz \
  -z
/opt/rCNV2/data_curation/other/curate_asc_dnms.py \
  --dnm-tsv satterstrom_asc_dnms.raw.tsv.gz \
  --genes gencode.v19.canonical.pext_filtered.genes.list \
  --controls \
  -o asc_dnm_counts.unaffecteds.tsv.gz \
  -z


# Reformat Redin 2017 translocation & inversion breakpoints
/opt/rCNV2/data_curation/other/curate_redin_bcas.py \
  --redin-tsv redin_2017.bca_breakpoints.all.tsv.gz \
  --gtf gencode.v19.canonical.pext_filtered.gtf.gz \
  --genes gencode.v19.canonical.pext_filtered.genes.list \
  --genome GRCh37.genome \
  -o redin_bca_counts.tsv.gz \
  --outbed redin_bca_breakpoints.bed.gz \
  -z


# Count functional SVs per gene from non-neuro gnomAD-SV subset
/opt/rCNV2/data_curation/other/curate_gnomad_sv_gene_counts.py \
  --gnomad-sv-vcf gnomad_v2.1_sv.nonneuro.sites.vcf.gz \
  --genes gencode.v19.canonical.pext_filtered.genes.list \
  -o gnomad_sv_nonneuro_counts.tsv.gz \
  -z


# Download & reformat gnomAD mutation rate table
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
/opt/rCNV2/data_curation/other/clean_gnomad_mutation_rates.R \
  gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
  gencode.v19.canonical.pext_filtered.genes.list \
  gene_mutation_rates.tsv
gzip -f gene_mutation_rates.tsv


# Parse ASC/SPARK pedigree for de novo CNV analyses
/opt/rCNV2/data_curation/other/parse_asc_spark_cnv_pedigree.R \
  asc_spark_denovo_cnvs/pedigree_cnv_07_27_2020.csv \
  asc_spark_child_phenotypes.list


# Process ASC/SPARK de novo CNVs, including:
# 1. reverse liftover from hg38 to GRCh37
# 2. restrict to autosomes
# 3. exclude known genomic disorders
# 4. re-annotate vs genes used in rCNV (further requiring copy-gain for dups)
# 5. restrict to samples in cleaned list of child phenotypes (see above)
/opt/rCNV2/data_curation/other/curate_asc_spark_denovo_cnvs.py \
  --bgzip \
  --outbed asc_spark_denovo_cnvs.cleaned.hg38.bed.gz \
  asc_spark_denovo_cnvs/dnv_cnv_07_27_2020.txt \
  asc_spark_child_phenotypes.list
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
liftOver -minMatch=0.5 -bedPlus=3 \
  asc_spark_denovo_cnvs.cleaned.hg38.bed.gz \
  hg38ToHg19.over.chain.gz \
  asc_spark_denovo_cnvs.cleaned.hg19.bed \
  asc_spark_denovo_cnvs.cleaned.hg38.liftFail.txt
for CNV in DEL DUP; do
  case $CNV in
    DEL)
      cds_ovr=0.05
      ;;
    DUP)
      cds_ovr=1.0
      ;;
  esac
  sed 's/^chr//g' asc_spark_denovo_cnvs.cleaned.hg19.bed \
  | sort -Vk1,1 -k2,2n -k3,3n \
  > asc_spark_denovo_cnvs.cleaned.b37.bed
  /opt/rCNV2/analysis/genes/count_cnvs_per_gene.py \
    --min-cds-ovr $cds_ovr \
    --type $CNV \
    --control-hpo Control \
    --outbed /dev/null \
    --cnvs-out asc_spark_denovo_cnvs.cleaned.b37.${CNV}.annotated.bed.gz \
    --bgzip \
    asc_spark_denovo_cnvs.cleaned.b37.bed \
    gencode.v19.canonical.pext_filtered.gtf.gz
done
for CNV in DEL DUP; do
  bedtools intersect -f 0.5 -r -v -wa \
    -a asc_spark_denovo_cnvs.cleaned.b37.${CNV}.annotated.bed.gz \
    -b <( zcat lit_GDs.*.bed.gz | fgrep ${CNV} ) \
  | bgzip -c \
  > asc_spark_denovo_cnvs.cleaned.b37.${CNV}.annotated.noGDs.bed.gz
done
zcat asc_spark_denovo_cnvs.cleaned.b37.DEL.annotated.noGDs.bed.gz \
     asc_spark_denovo_cnvs.cleaned.b37.DUP.annotated.noGDs.bed.gz \
| fgrep -v "#" \
| awk -v FS="\t" -v OFS="\t" '{ if ($7>0) print $1, $2, $3, $5, $4, $6, $7, $9 }' \
| sort -Vk1,1 -k2,2 -k3,3n \
| cat <( echo -e "#chr\tstart\tend\tcnv\tchild_id\tpheno\tn_genes\tgenes" ) - \
| bgzip -c \
> asc_spark_denovo_cnvs.cleaned.b37.annotated.bed.gz


# Copy curated DNMs, BCAs, mutation rates, and de novo CNVs to gs:// bucket (note: requires permissions)
gsutil -m cp \
  ddd_dnm_counts.tsv.gz \
  asc_dnm_counts*tsv.gz \
  redin_bca_counts.tsv.gz \
  redin_bca_breakpoints.bed.gz \
  gnomad_sv_nonneuro_counts.tsv.gz \
  gene_mutation_rates.tsv.gz \
  asc_spark_child_phenotypes.list \
  asc_spark_denovo_cnvs.cleaned.b37.annotated.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/misc/

