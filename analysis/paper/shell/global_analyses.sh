#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for global secondary analyses for rCNV manuscript


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export study_prefix="rCNV2_analysis_d2"
export control_hpo="HEALTHY_CONTROL"


# Localize all CNV data and analysis references
mkdir cnvs/
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/cnv/* cnvs/
mkdir refs/
gsutil -m cp -r \
  ${rCNV_bucket}/refs/GRCh37.*.bed.gz \
  ${rCNV_bucket}/analysis/analysis_refs/* \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  ${rCNV_bucket}/cleaned_data/phenotypes/hpo_logs_metadata/phenotype_groups.HPO_metadata.txt \
  ${rCNV_bucket}/analysis/paper/data/hpo/${study_prefix}.reordered_hpos.txt \
  refs/
gsutil -m cp \
  ${rCNV_bucket}/analysis/gene_scoring/gene_lists/* \
  refs/gene_lists/


# Prep directory tree and helper files
for dir in counts stats meta cnv_ids filtered_cnv_beds; do
  if [ -e $dir ]; then
    rm -rf $dir
  fi
  mkdir $dir
done
echo -e "#chr\tstart\tend\tgene\tcds\tcontrol_cnvs\tcontrol_cnvs_weighted\tcase_cnvs\tcase_cnvs_weighted" \
> counts/counts_header.tsv


# Helper function to format synthetic counts
format_counts () {
  # Positional args: hpo, k, n_ctrl, n_case, outfile
  echo -e "$1\t1\t2\t${1}_category_${2}\t1" \
  | awk -v n_ctrl="$3" -v n_case="$4" -v OFS="\t" \
    '{ print $0, n_ctrl, n_ctrl, n_case, n_case }' \
  | cat counts/counts_header.tsv - | bgzip -c > $5
}


# Collect counts per HPO for non-GTF-based comparisons
i=0
while read prefix hpo; do

  i=$((i+1))
  echo -e "$i\t$hpo"

  while read meta cohorts; do

    echo  "$meta"
    cnvbed=cnvs/$meta.rCNV.bed.gz

    # Collect counts for deletions and duplications separately
    for CNV in DEL DUP; do
      
      echo "$CNV"
      case $CNV in
        DEL)
          cds_frac=0.02
          ds_genes=refs/gene_lists/gold_standard.haploinsufficient.genes.list
          ;;
        DUP)
          cds_frac=0.84
          ds_genes=refs/gene_lists/gold_standard.triplosensitive.genes.list
          ;;
      esac

      # Category 1: Any rCNV
      n_ctrl=$( zcat $cnvbed | fgrep -w $CNV | fgrep -w $control_hpo | wc -l )
      n_case=$( zcat $cnvbed | fgrep -w $CNV | fgrep -w $hpo | wc -l )
      format_counts $hpo 1 $n_ctrl $n_case counts/$meta.$prefix.$CNV.1.bed.gz

      # Category 3: Known GDs
      zcat $cnvbed | fgrep -w $CNV \
      | bedtools intersect -u -F 0.5 -a - -b /opt/rCNV2/refs/lit_GDs.all.$CNV.bed.gz \
      | cut -f4 > cnv_ids/$meta.$CNV.3.ids.txt
      n_ctrl=$( zcat $cnvbed | fgrep -wf cnv_ids/$meta.$CNV.3.ids.txt \
                | fgrep -w $control_hpo | wc -l )
      n_case=$( zcat $cnvbed | fgrep -wf cnv_ids/$meta.$CNV.3.ids.txt \
                | fgrep -w $hpo | wc -l )
      format_counts $hpo 3 $n_ctrl $n_case counts/$meta.$prefix.$CNV.3.bed.gz

      # Category 10: All noncoding rCNVs
      n_ctrl=$( zcat cnvs/noncoding/$meta.rCNV.strict_noncoding.bed.gz \
                | fgrep -w $CNV | fgrep -w $control_hpo | wc -l )
      n_case=$( zcat cnvs/noncoding/$meta.rCNV.strict_noncoding.bed.gz \
                | fgrep -w $CNV | fgrep -w $hpo | wc -l )
      format_counts $hpo 10 $n_ctrl $n_case counts/$meta.$prefix.$CNV.10.bed.gz

    done

    # Combine counts for deletions and duplications
    for k in 1 3 10; do
      n_ctrl=$( zcat counts/$meta.$prefix.D*.$k.bed.gz | fgrep -v "#" \
                | awk -v FS="\t" '{ sum+=$6 }END{ print sum }' )
      n_case=$( zcat counts/$meta.$prefix.D*.$k.bed.gz | fgrep -v "#" \
                | awk -v FS="\t" '{ sum+=$8 }END{ print sum }' )
      format_counts $hpo $k $n_ctrl $n_case counts/$meta.$prefix.CNV.$k.bed.gz
    done

  done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt )

done < refs/test_phenotypes.list


# Filter GTF to various gene subsets
# Category 4: known DS genes
/opt/rCNV2/utils/filter_gtf_by_genelist.py \
  --o gencode.v19.canonical.pext_filtered.known_HI.gtf.gz \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  refs/gene_lists/gold_standard.haploinsufficient.genes.list
/opt/rCNV2/utils/filter_gtf_by_genelist.py \
  --o gencode.v19.canonical.pext_filtered.known_TS.gtf.gz \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  refs/gene_lists/gold_standard.triplosensitive.genes.list
# Category 5: Known DS genes outside of known GDs
gsutil -m cat \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
| gunzip -c | cut -f1-4 | fgrep -v "#" | sort -k4,4 \
| join -1 4 -2 1 -t $'\t' - <( sort refs/gene_lists/gold_standard.haploinsufficient.genes.list ) \
| awk -v FS="\t" -v OFS="\t" '{ print $2, $3, $4, $1 }' \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| bedtools intersect -v -a - -b /opt/rCNV2/refs/lit_GDs.all.DEL.bed.gz \
| cut -f4 | sort -Vk1,1 | uniq \
> refs/gene_lists/gold_standard.haploinsufficient.no_GDs.genes.list
/opt/rCNV2/utils/filter_gtf_by_genelist.py \
  --o gencode.v19.canonical.pext_filtered.known_HI.no_GDs.gtf.gz \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  refs/gene_lists/gold_standard.haploinsufficient.no_GDs.genes.list
gsutil -m cat \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
| gunzip -c | cut -f1-4 | fgrep -v "#" | sort -k4,4 \
| join -1 4 -2 1 -t $'\t' - <( sort refs/gene_lists/gold_standard.triplosensitive.genes.list ) \
| awk -v FS="\t" -v OFS="\t" '{ print $2, $3, $4, $1 }' \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| bedtools intersect -v -a - -b /opt/rCNV2/refs/lit_GDs.all.DUP.bed.gz \
| cut -f4 | sort -Vk1,1 | uniq \
> refs/gene_lists/gold_standard.triplosensitive.no_GDs.genes.list
/opt/rCNV2/utils/filter_gtf_by_genelist.py \
  --o gencode.v19.canonical.pext_filtered.known_TS.no_GDs.gtf.gz \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  refs/gene_lists/gold_standard.triplosensitive.no_GDs.genes.list
# Category 6: LoF-constrained genes
/opt/rCNV2/utils/filter_gtf_by_genelist.py \
  --o gencode.v19.canonical.pext_filtered.constrained.gtf.gz \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list
# Category 7: Constrained genes outside of known GDs
gsutil -m cat \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
| gunzip -c | cut -f1-4 | fgrep -v "#" | sort -k4,4 \
| join -1 4 -2 1 -t $'\t' - <( sort refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list ) \
| awk -v FS="\t" -v OFS="\t" '{ print $2, $3, $4, $1 }' \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| bedtools intersect -v -a - -b <( zcat /opt/rCNV2/refs/lit_GDs.all.DEL.bed.gz | fgrep -v "#" ) \
| cut -f4 | sort -Vk1,1 | uniq \
> refs/gene_lists/gold_standard.constrained.no_GDs.genes.list
/opt/rCNV2/utils/filter_gtf_by_genelist.py \
  --o gencode.v19.canonical.pext_filtered.constrained.no_GDs.gtf.gz \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  refs/gene_lists/gold_standard.constrained.no_GDs.genes.list
gsutil -m cat \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
| gunzip -c | cut -f1-4 | fgrep -v "#" | sort -k4,4 \
| join -1 4 -2 1 -t $'\t' - <( sort refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list ) \
| awk -v FS="\t" -v OFS="\t" '{ print $2, $3, $4, $1 }' \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| bedtools intersect -v -a - -b <( zcat /opt/rCNV2/refs/lit_GDs.all.DUP.bed.gz | fgrep -v "#" ) \
| cut -f4 | sort -Vk1,1 | uniq \
> refs/gene_lists/gold_standard.constrained.no_GDs.genes.list
# Category 8: Constrained genes outside of known GDs and known DS genes
gsutil -m cat \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
| gunzip -c | cut -f1-4 | fgrep -v "#" | sort -k4,4 \
| join -1 4 -2 1 -t $'\t' - <( sort refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
                               | fgrep -wvf refs/gene_lists/gold_standard.haploinsufficient.genes.list ) \
| awk -v FS="\t" -v OFS="\t" '{ print $2, $3, $4, $1 }' \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| bedtools intersect -v -a - -b <( zcat /opt/rCNV2/refs/lit_GDs.all.DEL.bed.gz | fgrep -v "#" ) \
| cut -f4 | sort -Vk1,1 | uniq \
> refs/gene_lists/gold_standard.constrained.no_GDs_no_HI.genes.list
/opt/rCNV2/utils/filter_gtf_by_genelist.py \
  --o gencode.v19.canonical.pext_filtered.constrained.no_GDs_no_HI.gtf.gz \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  refs/gene_lists/gold_standard.constrained.no_GDs_no_HI.genes.list
gsutil -m cat \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
| gunzip -c | cut -f1-4 | fgrep -v "#" | sort -k4,4 \
| join -1 4 -2 1 -t $'\t' - <( sort refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
                               | fgrep -wvf refs/gene_lists/gold_standard.triplosensitive.genes.list ) \
| awk -v FS="\t" -v OFS="\t" '{ print $2, $3, $4, $1 }' \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| bedtools intersect -v -a - -b <( zcat /opt/rCNV2/refs/lit_GDs.all.DUP.bed.gz | fgrep -v "#" ) \
| cut -f4 | sort -Vk1,1 | uniq \
> refs/gene_lists/gold_standard.constrained.no_GDs_no_TS.genes.list
/opt/rCNV2/utils/filter_gtf_by_genelist.py \
  --o gencode.v19.canonical.pext_filtered.constrained.no_GDs_no_TS.gtf.gz \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  refs/gene_lists/gold_standard.constrained.no_GDs_no_TS.genes.list
# Category 9: All genes outside of known GDs, known DS genes, and constrained genes
gsutil -m cat \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
| gunzip -c | cut -f1-4 | fgrep -v "#" \
| bedtools intersect -u -wa -a - -b <( zcat /opt/rCNV2/refs/lit_GDs.all.DEL.bed.gz | fgrep -v "#" ) \
| cut -f4 | cat - refs/gene_lists/gold_standard.haploinsufficient.genes.list \
                  refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
| sort -V | uniq \
> refs/gene_lists/all_GDs_or_HI_or_constrained.genes.list
/opt/rCNV2/utils/filter_gtf_by_genelist.py \
  --invert \
  --o gencode.v19.canonical.pext_filtered.other_genes.no_GDs_no_HI_no_constrained.gtf.gz \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  refs/gene_lists/all_GDs_or_HI_or_constrained.genes.list
gsutil -m cat \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
| gunzip -c | cut -f1-4 | fgrep -v "#" \
| bedtools intersect -u -wa -a - -b <( zcat /opt/rCNV2/refs/lit_GDs.all.DUP.bed.gz | fgrep -v "#" ) \
| cut -f4 | cat - refs/gene_lists/gold_standard.triplosensitive.genes.list \
                  refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
| sort -V | uniq \
> refs/gene_lists/all_GDs_or_TS_or_constrained.genes.list
/opt/rCNV2/utils/filter_gtf_by_genelist.py \
  --invert \
  --o gencode.v19.canonical.pext_filtered.other_genes.no_GDs_no_TS_no_constrained.gtf.gz \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.gtf.gz \
  refs/gene_lists/all_GDs_or_TS_or_constrained.genes.list


# Helper function to call count_cnvs_per_gene.py
count_genic_hits () {
  # Positional args: cds_frac, CNV, cnvbed, gtf, cnv_ids_outfile, outfile
      /opt/rCNV2/analysis/genes/count_cnvs_per_gene.py \
        --min-cds-ovr $1 \
        -t $2 \
        --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
        --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
        --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
        --verbose \
        --cnvs-out /dev/stdout \
        -o /dev/null \
        $3 \
        $4 \
      | fgrep -v "#" | fgrep -v "input" \
      | awk -v OFS="\t" '{ if ($7>0) print $4 }' > $5
      zcat $3 | fgrep -wf $5 | cut -f6 | sed 's/;/\n/g' | sort -Vk1,1 \
      | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' > $6
}


# Pre-compute lookup tables of CNVs overlapping gene sets per HPO per metacohort
while read meta cohorts; do

    echo  "$meta"
    cnvbed=cnvs/$meta.rCNV.bed.gz

    for CNV in DEL DUP; do
      
      echo "$CNV"
      case $CNV in
        DEL)
          cds_frac=0.02
          ds_label=HI
          ;;
        DUP)
          cds_frac=0.84
          ds_label=TS
          ;;
      esac

      # Category 2: Genic rCNVs
      count_genic_hits \
        $cds_frac $CNV $cnvbed \
        refs/gencode.v19.canonical.pext_filtered.gtf.gz \
        cnv_ids/$meta.$CNV.2.ids.txt \
        counts/$meta.genic_CNVs_per_hpo.$CNV.2.tsv

      # Category 4: Known DS genes
      count_genic_hits \
        $cds_frac $CNV $cnvbed \
        gencode.v19.canonical.pext_filtered.known_$ds_label.gtf.gz \
        cnv_ids/$meta.$CNV.4.ids.txt \
        counts/$meta.genic_CNVs_per_hpo.$CNV.4.tsv

      # Category 5: Known DS genes outside of known GDs
      zcat $cnvbed | fgrep -wvf cnv_ids/$meta.$CNV.3.ids.txt | bgzip -c \
      > filtered_cnv_beds/$meta.rCNV.no_3.bed.gz
      count_genic_hits \
        $cds_frac $CNV filtered_cnv_beds/$meta.rCNV.no_3.bed.gz \
        gencode.v19.canonical.pext_filtered.known_$ds_label.no_GDs.gtf.gz \
        cnv_ids/$meta.$CNV.5.ids.txt \
        counts/$meta.genic_CNVs_per_hpo.$CNV.5.tsv

      # Category 6: Constrained genes
      count_genic_hits \
        $cds_frac $CNV $cnvbed \
        gencode.v19.canonical.pext_filtered.constrained.gtf.gz \
        cnv_ids/$meta.$CNV.6.ids.txt \
        counts/$meta.genic_CNVs_per_hpo.$CNV.6.tsv

      # Category 7: Constrained genes outside of known GDs
      count_genic_hits \
        $cds_frac $CNV filtered_cnv_beds/$meta.rCNV.no_3.bed.gz \
        gencode.v19.canonical.pext_filtered.constrained.no_GDs_no_${ds_label}.gtf.gz \
        cnv_ids/$meta.$CNV.7.ids.txt \
        counts/$meta.genic_CNVs_per_hpo.$CNV.7.tsv

      # Category 8: Constrained genes outside of known GDs and known DS genes
      zcat $cnvbed | fgrep -wvf cnv_ids/$meta.$CNV.3.ids.txt \
      | fgrep -wvf cnv_ids/$meta.$CNV.5.ids.txt | bgzip -c \
      > filtered_cnv_beds/$meta.rCNV.no_3_5.bed.gz
      count_genic_hits \
        $cds_frac $CNV filtered_cnv_beds/$meta.rCNV.no_3_5.bed.gz \
        gencode.v19.canonical.pext_filtered.constrained.no_GDs_no_${ds_label}.gtf.gz \
        cnv_ids/$meta.$CNV.8.ids.txt \
        counts/$meta.genic_CNVs_per_hpo.$CNV.8.tsv

      # Category 9: All genes outside of known GDs, known DS genes, and constrained genes
      zcat $cnvbed | fgrep -wvf cnv_ids/$meta.$CNV.3.ids.txt \
      | fgrep -wvf cnv_ids/$meta.$CNV.5.ids.txt \
      | fgrep -wvf cnv_ids/$meta.$CNV.8.ids.txt | bgzip -c \
      > filtered_cnv_beds/$meta.rCNV.no_3_5_8.bed.gz
      count_genic_hits \
        $cds_frac $CNV filtered_cnv_beds/$meta.rCNV.no_3_5_8.bed.gz \
        gencode.v19.canonical.pext_filtered.other_genes.no_GDs_no_${ds_label}_no_constrained.gtf.gz \
        cnv_ids/$meta.$CNV.9.ids.txt \
        counts/$meta.genic_CNVs_per_hpo.$CNV.9.tsv

    done

done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt )


# Collect counts per HPO for all GTF-based comparisons
i=0
while read prefix hpo; do

  i=$((i+1))
  echo -e "$i\t$hpo"

  while read meta cohorts; do

    echo  "$meta"
    cnvbed=cnvs/$meta.rCNV.bed.gz

    # Collect counts for deletions and duplications separately
    for CNV in DEL DUP; do
      
      echo "$CNV"
      case $CNV in
        DEL)
          cds_frac=0.02
          ds_genes=refs/gene_lists/gold_standard.haploinsufficient.genes.list
          ;;
        DUP)
          cds_frac=0.84
          ds_genes=refs/gene_lists/gold_standard.triplosensitive.genes.list
          ;;
      esac

      for k in 2 $( seq 4 9 ); do
        n_ctrl=$( fgrep -w $control_hpo counts/$meta.genic_CNVs_per_hpo.$CNV.$k.tsv | cut -f2 )
        if [ -z $n_ctrl ]; then n_ctrl=0; fi
        n_case=$( fgrep -w $hpo counts/$meta.genic_CNVs_per_hpo.$CNV.$k.tsv | cut -f2 )
        if [ -z $n_case ]; then n_case=0; fi
        format_counts $hpo $k $n_ctrl $n_case counts/$meta.$prefix.$CNV.$k.bed.gz
      done

    done

    # Combine counts for deletions and duplications
    for k in 2 $( seq 4 9 ); do
      n_ctrl=$( zcat counts/$meta.$prefix.D*.$k.bed.gz | fgrep -v "#" \
                | awk -v FS="\t" '{ sum+=$6 }END{ print sum }' )
      n_case=$( zcat counts/$meta.$prefix.D*.$k.bed.gz | fgrep -v "#" \
                | awk -v FS="\t" '{ sum+=$8 }END{ print sum }' )
      format_counts $hpo $k $n_ctrl $n_case counts/$meta.$prefix.CNV.$k.bed.gz
    done

  done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt )

done < refs/test_phenotypes.list


# Compute single-cohort effect size estimates for each HPO & metacohort
i=0
while read prefix hpo; do

  i=$((i+1))
  echo -e "$i\t$hpo"

  while read meta cohorts; do

    echo "$meta"

    for CNV in DEL DUP CNV; do

      echo "$CNV"

      for k in $( seq 1 10 ); do
        /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
          --pheno-table refs/HPOs_by_metacohort.table.tsv \
          --cohort-name $meta \
          --case-hpo "$hpo" \
          --keep-n-columns 4 \
          --bgzip \
          "counts/$meta.$prefix.$CNV.$k.bed.gz" \
          "stats/$meta.$prefix.$CNV.$k.stats.bed.gz" \
        &> /dev/null
      done

    done

  done < <( fgrep -v mega refs/rCNV_metacohort_list.txt )

done < refs/test_phenotypes.list


# Conduct meta-analysis for each HPO & category
i=0
while read prefix hpo; do

  i=$((i+1))
  echo -e "$i\t$hpo"

  for CNV in DEL DUP CNV; do

    echo "$CNV"

    for k in $( seq 1 10 ); do

      while read meta cohorts; do
        counts_file="stats/$meta.$prefix.$CNV.$k.stats.bed.gz"
        if [ -e $counts_file ]; then
          echo -e "${meta}\t${counts_file}"
        fi
      done < <( fgrep -v mega refs/rCNV_metacohort_list.txt ) \
      > meta/$prefix.$CNV.$k.meta_input.tsv

      /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
        --model "fe" \
        --p-is-neg-log10 \
        --min-cases 300 \
        --keep-n-columns 4 \
        meta/$prefix.$CNV.$k.meta_input.tsv \
        meta/$prefix.$CNV.$k.stats.tsv \
      &> /dev/null

    done

  done

done < refs/test_phenotypes.list


# Pool all meta-analysis results into single files for analysis
head -n1 meta/$( head -n1 refs/test_phenotypes.list | cut -f1 ).DEL.1.stats.tsv \
| cut -f5- | paste <( echo -e "#HPO\tCNV\tcategory" ) - \
> meta/meta_stats.header.tsv
for CNV in DEL DUP CNV; do

  while read prefix hpo; do

    for k in $( seq 1 10 ); do
      fgrep -v "#" meta/$prefix.$CNV.$k.stats.tsv | cut -f5- \
      | paste <( echo -e "${hpo}\t${CNV}\t${k}" ) -
    done

  done < refs/test_phenotypes.list

done | sort -Vk1,1 -k2,2V -k3,3n | cat meta/meta_stats.header.tsv - | gzip -c \
> ${study_prefix}.global_burden_stats.tsv.gz


# Copy global stats to Google Bucket (note: requires permissions)
gsutil -m cp \
  ${study_prefix}.global_burden_stats.tsv.gz \
  ${rCNV_bucket}/analysis/paper/data/global/


# Plot large grid summarizing results
/opt/rCNV2/analysis/paper/plot/misc/plot_global_summary.R \
  refs/${study_prefix}.reordered_hpos.txt \
  refs/phenotype_groups.HPO_metadata.txt \
  refs/HPOs_by_metacohort.table.tsv \
  ${study_prefix}.global_burden_stats.tsv.gz \
  ${study_prefix}.global_summary


# Copy plot to Google Bucket (note: requires permissions)
gsutil -m cp \
  ${study_prefix}.global_summary.pdf \
  ${rCNV_bucket}/analysis/paper/plots/misc/

