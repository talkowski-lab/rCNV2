#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Partition HPO terms into developmental and adult-onset subsets based on empirical CNV effect sizes


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv:latest
gcloud auth login


# Copy all filtered CNV data and other references from the project Google Bucket (note: requires permissions)
mkdir cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
mkdir phenos/
gsutil -m cp -r \
  gs://rcnv_project/cleaned_data/phenotypes/filtered/* \
  phenos/
mkdir refs/
gsutil -m cp -r \
  gs://rcnv_project/refs/GRCh37.*.bed.gz \
  gs://rcnv_project/analysis/analysis_refs/* \
  gs://rcnv_project/cleaned_data/genes/gene_lists \
  gs://rcnv_project/cleaned_data/phenotypes/hpo_logs_metadata/phenotype_groups.HPO_metadata.txt \
  refs/


# Extract constrained genes from GTF
gsutil -m cat \
  gs://rcnv_project/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
| gunzip -c | cut -f1-4 | fgrep -v "#" | sort -k4,4 \
| join -1 4 -2 1 -t $'\t' - refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
| awk -v FS="\t" -v OFS="\t" '{ print $2, $3, $4, $1 }' \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| cat <( echo -e "#chr\tstart\tend\tgene" ) - \
> gencode.v19.canonical.constrained.coords.bed


# Collect counts of whole-gene deletions of constrained genes per phenotype
# Note: for an unknown reason, the below bedtools command does not work 
while read meta cohorts; do
  cnv_bed=cleaned_cnv/$meta.rCNV.bed.gz
  if [ -e $cnv_bed ]; then
    bedtools intersect -wa -u -F 1.0 \
      -a <( zcat $cnv_bed | fgrep -w DEL ) \
      -b gencode.v19.canonical.constrained.coords.bed \
    | cut -f6 | sed 's/;/\n/g' | sort -Vk1,1 | uniq -c \
    | awk -v OFS="\t" '{ print $2, $1 }' \
    > $meta.constrained_del_counts.per_hpo.tsv
  fi
done < <( fgrep -v mega refs/rCNV_metacohort_list.txt )


# Compute single-cohort effect sizes for meta-analysis
while read prefix hpo; do
  while read meta cohorts; do
    if [ -e $meta.constrained_del_counts.per_hpo.tsv ]; then
      # Make synthetic counts
      n_ctrl=$( awk -v hpo="HEALTHY_CONTROL" '{ if ($1==hpo) print $2; else print "0"}' \
                  $meta.constrained_del_counts.per_hpo.tsv \
                | sort -nrk1,1 | uniq | head -n1 )
      n_case=$( awk -v hpo="$hpo" '{ if ($1==hpo) print $2; else print "0"}' \
                  $meta.constrained_del_counts.per_hpo.tsv \
                | sort -nrk1,1 | uniq | head -n1 )
      echo -e "1\t1\t2\tall_constrained_genes\t$( cat gencode.v19.canonical.constrained.coords.bed | wc -l )" \
      | awk -v n_ctrl="$n_ctrl" -v n_case="$n_case" -v OFS="\t" \
        '{ print $0, n_ctrl, n_ctrl, n_case, n_case }' \
      | cat meta_tmp_files/counts_header.tsv - \
      | bgzip -c > meta_tmp_files/$meta.$prefix.constrained_del_counts.bed.gz

      # Compute single-cohort stats
      /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
        --pheno-table refs/HPOs_by_metacohort.table.tsv \
        --cohort-name $meta \
        --case-hpo "$hpo" \
        --keep-n-columns 4 \
        --bgzip \
        "meta_tmp_files/$meta.$prefix.constrained_del_counts.bed.gz" \
        "meta_tmp_files/$meta.$prefix.constrained_del_stats.bed.gz"
    fi
  done < <( fgrep -v mega refs/rCNV_metacohort_list.txt )
done < refs/test_phenotypes.list


# Compute meta-analysis effect size per HPO
mkdir meta_tmp_files
echo -e "#chr\tstart\tend\tgene\tcds\tcontrol_cnvs\tcontrol_cnvs_weighted\tcase_cnvs\tcase_cnvs_weighted" \
> meta_tmp_files/counts_header.tsv
while read prefix hpo; do
  while read meta cohorts; do
    counts_file="meta_tmp_files/$meta.$prefix.constrained_del_stats.bed.gz"
    if [ -e $counts_file ]; then
      echo -e "${meta}\t${counts_file}"
    fi
  done < <( fgrep -v mega refs/rCNV_metacohort_list.txt ) \
  > meta_tmp_files/$prefix.meta_input.tsv
  if [ $( cat meta_tmp_files/$prefix.meta_input.tsv | wc -l ) -gt 1 ]; then
    /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
      --model "fe" \
      --p-is-neg-log10 \
      --min-cases 300 \
      --keep-n-columns 4 \
      meta_tmp_files/$prefix.meta_input.tsv \
      meta_tmp_files/$prefix.stats.tsv
  fi
done < refs/test_phenotypes.list


# Collapse effect size estimates into single table for analysis
while read prefix hpo; do
  if [ -e meta_tmp_files/$prefix.stats.tsv ]; then
    sed '1d' meta_tmp_files/$prefix.stats.tsv | cut -f8-12 \
    | awk -v OFS="\t" -v hpo="$hpo" '{ print hpo, $0 }'
  fi
done < refs/test_phenotypes.list | sort -nrk4,4 \
| cat <( echo -e "hpo\tcase_freq\tcontrol_freq\tlnOR\tlnOR_lower\tlnOR_upper" ) - \
> constrained_gene_del_stats.all_hpos.tsv


# Partition phenotypes by effect size of whole-gene deletions of constrained genes
/opt/rCNV2/analysis/other/split_phenos_by_severity.R \
  --out-prefix rCNV2.hpos_by_severity \
  --threshold 2 \
  --use-lower \
  constrained_gene_del_stats.all_hpos.tsv


# Get counts of cases per metacohort matching at least one vs. zero developmental terms
for file in pheno_inputs.tsv precomp_pairs.tsv; do
  if [ -e $file ]; then rm $file; fi
done
while read meta cohorts; do
  if [ -e "phenos/${cohorts}.final_cooccurrence_table.tsv.gz" ]; then
    echo -e "${meta}\tphenos/${cohorts}.final_cooccurrence_table.tsv.gz" >> precomp_pairs.tsv
  else
    echo -e "${meta}\tphenos/${meta}.cleaned_phenos.txt" >> pheno_inputs.tsv
  fi
done < <( fgrep -v mega refs/rCNV_metacohort_list.txt )
/opt/rCNV2/data_curation/phenotype/count_samples_by_pheno_list.py \
  --hpo-pair-cohorts precomp_pairs.tsv \
  --hpo-tier-metadata refs/phenotype_groups.HPO_metadata.txt \
  --outfile rCNV2.hpos_by_severity.developmental.counts.tsv \
  pheno_inputs.tsv \
  rCNV2.hpos_by_severity.developmental.list
/opt/rCNV2/data_curation/phenotype/count_samples_by_pheno_list.py \
  --hpo-pair-cohorts precomp_pairs.tsv \
  --invert \
  --hpo-tier-metadata refs/phenotype_groups.HPO_metadata.txt \
  --outfile rCNV2.hpos_by_severity.adult.counts.tsv \
  pheno_inputs.tsv \
  rCNV2.hpos_by_severity.developmental.list


# Copy phenotype classifications to Google bucket
gsutil -m cp \
  rCNV2.hpos_by_severity.*list \
  rCNV2.hpos_by_severity.*tsv \
  gs://rcnv_project/analysis/analysis_refs/

