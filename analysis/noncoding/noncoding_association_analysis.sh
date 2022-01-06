#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control noncoding CNV burdens for cis-regulatory blocks (CRBs)


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Copy all filtered CNV data, gene coordinates, and other references 
# from the project Google Bucket (note: requires permissions)
mkdir cleaned_cnv/
gsutil -m cp -r \
  gs://rcnv_project/cleaned_data/cnv/noncoding/* \
  gs://rcnv_project/cleaned_data/cnv/*bed.gz* \
  cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/genome_annotations/*bed.gz* ./
mkdir refs/
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/


# Test/dev parameters
# Params for seizures
hpo="HP:0001250"
prefix="HP0001250"
meta="meta2"
# General params
freq_code="rCNV"
noncoding_filter="loose"
CNV="DEL"
phenotype_list="refs/test_phenotypes.list"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
crbs="rCNV.crbs.bed.gz"
crb_elements="rCNV.crb_elements.bed.gz"
rCNV_bucket="gs://rcnv_project"
pad_controls=0
min_element_ovr=1.0
min_frac_all_elements=0.05
p_cutoff=0.000005076658
meta_p_cutoff=0.000005076658
max_manhattan_neg_log10_p=30
n_pheno_perms=50
meta_model_prefix="fe"
i=1
winsorize_meta_z=0.99
meta_min_cases=300



### Determine probe density-based conditional exclusion list
# NOTE: This code must be run using a DIFFERENT DOCKER: us.gcr.io/broad-dsmap/athena-cloud
# Test/dev parameters
crbs_prefix="rCNV.crbs" #Note: this can be inferred in WDL as basename(crbs, ".bed.gz")
min_probes_per_crb=10
min_frac_controls_probe_exclusion=0.9
metacohort_list="refs/rCNV_metacohort_list.txt"
rCNV_bucket="gs://rcnv_project"
freq_code="rCNV"

# Download probesets (note: requires permissions)
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/control_probesets \
  ./

# Subset crb coordinates to minimal BED4
zcat ${crbs} | cut -f1-4 | bgzip -c > crb_coords.bed.gz

# Compute conditional cohort exclusion mask
for file in control_probesets/*bed.gz; do
  echo -e "$file\t$( basename $file | sed 's/\.bed\.gz//g' )"
done > probeset_tracks.tsv
/opt/rCNV2/data_curation/other/probe_based_exclusion.py \
  --outfile ${crbs_prefix}.cohort_exclusion.bed.gz \
  --probecounts-outfile ${crbs_prefix}.probe_counts.bed.gz \
  --control-mean-counts-outfile ${crbs_prefix}.mean_probe_counts_per_cohort.bed.gz \
  --frac-pass-outfile ${crbs_prefix}.frac_passing.bed.gz \
  --min-probes ${min_probes_per_crb} \
  --min-frac-samples ${min_frac_controls_probe_exclusion} \
  --keep-n-columns 4 \
  --bgzip \
  crb_coords.bed.gz \
  probeset_tracks.tsv \
  control_probesets/rCNV.control_counts_by_array.tsv \
  <( fgrep -v mega ${metacohort_list} )



# Count CNVs in cases and controls per phenotype, split by metacohort and CNV type
# Iterate over phenotypes
while read pheno hpo; do
  # Set HPO-specific parameters
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )

  # Iterate over metacohorts
  while read meta cohorts; do
    echo $meta

    # Set metacohort-specific parameters
    cnv_bed="cleaned_cnv/$meta.${freq_code}.${noncoding_filter}_noncoding.bed.gz"
    meta_idx=$( head -n1 "${metacohort_sample_table}" \
                | sed 's/\t/\n/g' \
                | awk -v meta="$meta" '{ if ($1==meta) print NR }' )
    ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
             | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
    nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
             | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
             | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
    title="$descrip (${hpo})\n$ncase cases vs $nctrl controls in '$meta' cohort"

    # Iterate over CNV types
    for CNV in DEL DUP; do
      echo $CNV

      # Count CNVs
      /opt/rCNV2/analysis/noncoding/count_cnvs_per_crb.py \
        --cnvs $cnv_bed \
        --crbs ${crbs} \
        --elements ${crb_elements} \
        --pad-controls ${pad_controls} \
        --min-element-ovr ${min_element_ovr} \
        --min-frac-all-elements ${min_frac_all_elements} \
        -t $CNV \
        --hpo ${hpo} \
        -z \
        -o "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz" 
      tabix -f "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz"

      # Perform burden test
      /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
        --pheno-table ${metacohort_sample_table} \
        --cohort-name $meta \
        --case-hpo ${hpo} \
        --keep-n-columns 4 \
        --bgzip \
        "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz" \
        "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"
      tabix -f "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "fisher_neg_log10_p" \
        --p-is-neg-log10 \
        --max-neg-log10-p ${max_manhattan_neg_log10_p} \
        --cutoff ${p_cutoff} \
        --label-prefix "$CNV" \
        --title "$title" \
        "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden"
    done

    # Generate Miami & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --miami \
      --p-col-name "fisher_neg_log10_p" \
      --p-is-neg-log10 \
      --max-neg-log10-p ${max_manhattan_neg_log10_p} \
      --cutoff ${p_cutoff} \
      --label-prefix "DUP" \
      --label-prefix-2 "DEL" \
      --title "$title" \
      "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.DUP.crb_burden.stats.bed.gz" \
      "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.DEL.crb_burden.stats.bed.gz" \
      "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.crb_burden"
  done < ${metacohort_list}
done < ${phenotype_list}




# Count *all* CNVs in cases and controls per phenotype, split by metacohort and CNV type
# NOTE: without restricting on noncoding CNVs
# Iterate over phenotypes
while read pheno hpo; do
  # Iterate over metacohorts
  while read meta cohorts; do
    echo $meta
    cnv_bed="cleaned_cnv/$meta.${freq_code}.bed.gz"

    # Iterate over CNV types
    for CNV in DEL DUP; do
      echo $CNV

      # Count CNVs
      /opt/rCNV2/analysis/noncoding/count_cnvs_per_crb.py \
        --cnvs $cnv_bed \
        --crbs ${crbs} \
        --elements ${crb_elements} \
        --pad-controls ${pad_controls} \
        --min-element-ovr ${min_element_ovr} \
        --min-frac-all-elements ${min_frac_all_elements} \
        -t $CNV \
        --hpo ${hpo} \
        -z \
        -o "$meta.${prefix}.${freq_code}.$CNV.crb_burden.counts.bed.gz" 
      tabix -f "$meta.${prefix}.${freq_code}.$CNV.crb_burden.counts.bed.gz"

      # Perform burden test
      /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
        --pheno-table ${metacohort_sample_table} \
        --cohort-name $meta \
        --case-hpo ${hpo} \
        --keep-n-columns 4 \
        --bgzip \
        "$meta.${prefix}.${freq_code}.$CNV.crb_burden.counts.bed.gz" \
        "$meta.${prefix}.${freq_code}.$CNV.crb_burden.stats.bed.gz"
      tabix -f "$meta.${prefix}.${freq_code}.$CNV.crb_burden.stats.bed.gz"
    done
  done < ${metacohort_list}
done < ${phenotype_list}




# Permute CNVs in cases and controls per phenotype, split by metacohort and CNV type
# Iterate over phenotypes
while read prefix hpo; do
  for i in $( seq 1 ${n_pheno_perms} ); do

    # Shuffle phenotypes for each metacohort CNV dataset, and restrict CNVs to
    # phenotype of relevance
    if ! [ -e shuffled_cnv/ ]; then
      mkdir shuffled_cnv/
    fi
    yes $i | head -n1000000 > seed_$i.txt
    while read meta cohorts; do

      cnvbed="cleaned_cnv/$meta.${freq_code}.${noncoding_filter}_noncoding.bed.gz"

      # Shuffle phenotypes, matching by CNV type
      for CNV in DEL DUP; do
        zcat $cnvbed | fgrep -w $CNV \
        > cnv_subset.bed
        paste <( cut -f1-5 cnv_subset.bed ) \
              <( cut -f6 cnv_subset.bed | shuf --random-source seed_$i.txt ) \
        | awk -v hpo=${hpo} '{ if ($NF ~ "HEALTHY_CONTROL" || $NF ~ hpo) print $0 }'
        rm cnv_subset.bed
      done \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat <( tabix -H $cnvbed ) - \
      | bgzip -c \
      > shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz
      tabix -f shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz

    done < <( fgrep -v "mega" ${metacohort_list} )

    # Iterate over CNV types
    for CNV in DEL DUP; do

      # Perform association test for each metacohort
      while read meta cohorts; do

        echo -e "[$( date )] Starting permutation $i for $CNV in ${hpo} from $meta...\n"

        cnv_bed="shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz"

        # Count CNVs
        /opt/rCNV2/analysis/noncoding/count_cnvs_per_crb.py \
          --cnvs $cnv_bed \
          --crbs ${crbs} \
          --elements ${crb_elements} \
          --pad-controls ${pad_controls} \
          --min-element-ovr ${min_element_ovr} \
          --min-frac-all-elements ${min_frac_all_elements} \
          -t $CNV \
          --hpo ${hpo} \
          -z \
          -o "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz"

        # Perform burden test
        /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --case-hpo ${hpo} \
          --keep-n-columns 4 \
          --bgzip \
          "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz" \
          "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"

      done < <( fgrep -v "mega" ${metacohort_list} )

      # Perform meta-analysis
      while read meta cohorts; do
        echo -e "$meta\t$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"
      done < <( fgrep -v mega ${metacohort_list} ) \
      > ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.input.txt
      /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
        --model ${meta_model_prefix} \
        --conditional-exclusion ${exclusion_bed} \
        --p-is-neg-log10 \
        --spa \
        --spa-exclude /opt/rCNV2/refs/lit_GDs.all.$CNV.bed.gz \
        --winsorize ${winsorize_meta_z} \
        --min-cases ${meta_min_cases} \
        --keep-n-columns 4 \
        ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.input.txt \
        ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.perm_$i.bed
      bgzip -f ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.perm_$i.bed
      tabix -f ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.perm_$i.bed.gz

  done
done < ${phenotype_list}

# Gather all permutation results and compute FDR CDFs
mkdir perm_res/
while read prefix hpo; do
  gsutil -m cp \
    "${rCNV_bucket}/analysis/crb_burden/$prefix/${freq_code}/permutations/$prefix.${freq_code}.${noncoding_filter}_noncoding.${CNV}.crb_burden.meta_analysis.stats.perm_*.bed.gz" \
    perm_res/
  for i in $( seq 1 ${n_pheno_perms} ); do
    p_idx=$( zcat perm_res/$prefix.${freq_code}.${noncoding_filter}_noncoding.${CNV}.crb_burden.meta_analysis.stats.perm_$i.bed.gz \
             | sed -n '1p' | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' \
             | fgrep -w ${p_val_column_name} | cut -f2 )
    zcat perm_res/$prefix.${freq_code}.${noncoding_filter}_noncoding.${CNV}.crb_burden.meta_analysis.stats.perm_$i.bed.gz \
    | grep -ve '^#' \
    | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
    | cat <( echo "$prefix.${CNV}.$i" ) - \
    > perm_res/$prefix.${freq_code}.${noncoding_filter}_noncoding.${CNV}.crb_burden.meta_analysis.permuted_p_values.$i.txt
  done
  rm perm_res/$prefix.${freq_code}.${noncoding_filter}_noncoding.${CNV}.crb_burden.meta_analysis.stats.perm_*.bed.gz
done < ${phenotype_list}
paste perm_res/*.crb_burden.meta_analysis.permuted_p_values.*.txt \
| gzip -c \
> ${freq_code}.${noncoding_filter}_noncoding.${CNV}.permuted_pval_matrix.txt.gz

# Analyze p-values and compute FDR
/opt/rCNV2/analysis/sliding_windows/calc_empirical_fdr.R \
  --cnv ${CNV} \
  --fdr-target ${meta_p_cutoff} \
  --linear-fit \
  --flat-ladder \
  --plot crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}_permutation_results.png \
  ${freq_code}.${noncoding_filter}_noncoding.${CNV}.permuted_pval_matrix.txt.gz \
  ${metacohort_sample_table} \
  crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}

# Also produce an optional table of flat Bonferroni P-value cutoffs
awk -v pval=${meta_p_cutoff} -v FS="\t" -v OFS="\t" \
  '{ print $1, pval }' ${phenotype_list} \
> crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.bonferroni_pval.hpo_cutoffs.tsv



# Run meta-analysis for each phenotype
# Dev: copy precomputed counts (so the entire upper block doesn't have to be run)
gsutil -m cp ${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/stats/meta*.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz ./
exclusion_bed="rCNV.crbs.cohort_exclusion.bed.gz" #Note: this file must be generated above
while read prefix hpo; do

  # Get metadata for meta-analysis (while accounting for cohorts below inclusion criteria)
  last_cohort_col=$( head -n1 "${metacohort_sample_table}" | awk '{ print NF-1 }' )
  keep_cols=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | cut -f4-$last_cohort_col \
               | sed 's/\t/\n/g' \
               | awk -v min_n=${meta_min_cases} '{ if ($1>=min_n) print NR+3 }' \
               | paste -s -d, )
  ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
           | cut -f$keep_cols \
           | sed 's/\t/\n/g' \
           | awk '{ sum+=$1 }END{ print sum }' )
  nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
           | cut -f$keep_cols \
           | sed 's/\t/\n/g' \
           | awk '{ sum+=$1 }END{ print sum }' )
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )
  title="$descrip (${hpo})\nMeta-analysis of $ncase cases and $nctrl controls"
  DEL_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                  crb_burden.${freq_code}.${noncoding_filter}_noncoding.DEL.bonferroni_pval.hpo_cutoffs.tsv )
  DUP_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                  crb_burden.${freq_code}.${noncoding_filter}_noncoding.DUP.bonferroni_pval.hpo_cutoffs.tsv )

  # Set HPO-specific parameters
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )

  # Run meta-analysis for each CNV type
  for CNV in DEL DUP; do
    # Set CNV-specific parameters
    case "$CNV" in
      DEL)
        meta_p_cutoff=$DEL_p_cutoff
        ;;
      DUP)
        meta_p_cutoff=$DUP_p_cutoff
        ;;
    esac

    # Perform meta-analysis of CNV counts
    while read meta cohorts; do
      echo -e "$meta\t$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"
    done < <( fgrep -v mega ${metacohort_list} ) \
    > ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.input.txt
    /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
      --model ${meta_model_prefix} \
      --conditional-exclusion ${exclusion_bed} \
      --p-is-neg-log10 \
      --spa \
      --spa-exclude /opt/rCNV2/refs/lit_GDs.all.$CNV.bed.gz \
      --winsorize ${winsorize_meta_z} \
      --min-cases ${meta_min_cases} \
      --keep-n-columns 4 \
      ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.input.txt \
      ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed
    bgzip -f ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed
    tabix -f ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed.gz

    # Generate Manhattan & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --p-col-name "meta_neg_log10_p" \
      --p-is-neg-log10 \
      --max-neg-log10-p ${max_manhattan_neg_log10_p} \
      --cutoff $meta_p_cutoff \
      --label-prefix "$CNV" \
      --title "$title" \
      "${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed.gz" \
      "${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis"
  done

  # Generate Miami & QQ plots
  /opt/rCNV2/utils/plot_manhattan_qq.R \
    --miami \
    --p-col-name "meta_neg_log10_p" \
    --p-is-neg-log10 \
    --max-neg-log10-p ${max_manhattan_neg_log10_p} \
    --cutoff $DUP_p_cutoff \
    --label-prefix "DUP" \
    --cutoff-2 $DEL_p_cutoff \
    --label-prefix-2 "DEL" \
    --title "$title" \
    "${prefix}.${freq_code}.${noncoding_filter}_noncoding.DUP.crb_burden.meta_analysis.stats.bed.gz" \
    "${prefix}.${freq_code}.${noncoding_filter}_noncoding.DEL.crb_burden.meta_analysis.stats.bed.gz" \
    "${prefix}.${freq_code}.${noncoding_filter}_noncoding.crb_burden.meta_analysis"

done < refs/test_phenotypes.list

# Collapse all meta-analysis p-values into single matrix for visualizing calibration
mkdir meta_res/
# Download data
while read prefix hpo; do
  for CNV in DEL DUP; do
    echo -e "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/stats/${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed.gz"
  done
done < ${phenotype_list} \
| gsutil -m cp -I meta_res/
# Primary p-values
p_val_column_name="meta_neg_log10_p"
while read prefix hpo; do
  echo -e "$prefix\n\n"
  for CNV in DEL DUP; do
      stats=meta_res/${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed.gz
      if [ -e $stats ]; then
        p_idx=$( zcat $stats | sed -n '1p' | sed 's/\t/\n/g' \
                 | awk -v OFS="\t" '{ print $1, NR }' \
                 | fgrep -w ${p_val_column_name} | cut -f2 )
        zcat $stats | grep -ve '^#' \
        | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
        | cat <( echo "$prefix.$CNV" ) - \
        > meta_res/$prefix.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.p_values.txt
      fi
  done
done < ${phenotype_list}
paste meta_res/*.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.p_values.txt \
| gzip -c \
> ${freq_code}.observed_pval_matrix.txt.gz
# Secondary p-values
p_val_column_name="meta_neg_log10_p_secondary"
while read prefix hpo; do
  echo -e "$prefix\n\n"
  for CNV in DEL DUP; do
      stats=meta_res/${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed.gz
      if [ -e $stats ]; then
        p_idx=$( zcat $stats | sed -n '1p' | sed 's/\t/\n/g' \
                 | awk -v OFS="\t" '{ print $1, NR }' \
                 | fgrep -w ${p_val_column_name} | cut -f2 )
        zcat $stats | grep -ve '^#' \
        | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
        | cat <( echo "$prefix.$CNV" ) - \
        > meta_res/$prefix.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.secondary_p_values.txt
      fi
  done
done < ${phenotype_list}
paste meta_res/*.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.secondary_p_values.txt \
| gzip -c \
> ${freq_code}.observed_pval_matrix_secondary.txt.gz
# DEV NOTE: these p-values can be visualized with plot_gene_burden_meta_analysis_p_values.R,
#           which is currently just a code snippet referencing local filepaths




# Run unfiltered meta-analysis (including coding CNVs) for each phenotype
while read prefix hpo; do

  # Get metadata for meta-analysis
  DEL_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                  crb_burden.${freq_code}.${noncoding_filter}_noncoding.DEL.bonferroni_pval.hpo_cutoffs.tsv )
  DUP_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                  crb_burden.${freq_code}.${noncoding_filter}_noncoding.DUP.bonferroni_pval.hpo_cutoffs.tsv )

  # Run meta-analysis for each CNV type
  for CNV in DEL DUP; do

    # Perform meta-analysis of CNV counts
    while read meta cohorts; do
      echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.crb_burden.stats.bed.gz"
    done < <( fgrep -v mega ${metacohort_list} ) \
    > ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.input.txt
    
    /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
      --model ${meta_model_prefix} \
      --conditional-exclusion ${exclusion_bed} \
      --p-is-neg-log10 \
      --spa \
      --spa-exclude /opt/rCNV2/refs/lit_GDs.all.$CNV.bed.gz \
      --winsorize ${winsorize_meta_z} \
      --min-cases ${meta_min_cases} \
      --keep-n-columns 4 \
      ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.input.txt \
      ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.stats.bed
    bgzip -f ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.stats.bed
    tabix -f ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.stats.bed.gz

  done

done < refs/test_phenotypes.list




# Extract CRBs meeting all significance criteria
# TODO: NEED TO ADD FDR SUBSET TO THIS
mkdir meta_res/
# Download data
while read prefix hpo; do
  for CNV in DEL DUP; do
    echo -e "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/stats/${prefix}.${freq_code}.**$CNV.crb_burden.meta_analysis.stats.bed.gz"
  done
done < ${phenotype_list} \
| gsutil -m cp -I meta_res/
for CNV in DEL DUP; do
  # Build input files
  for x in coding noncoding; do
    if [ -e ${freq_code}.${x}.sig_crb_input.${CNV}.tsv ]; then
      rm ${freq_code}.${x}.sig_crb_input.${CNV}.tsv
    fi
  done
  while read prefix hpo; do
    echo -e "${hpo}\tmeta_res/${prefix}.${freq_code}.${noncoding_filter}_noncoding.${CNV}.crb_burden.meta_analysis.stats.bed.gz" \
    >> ${freq_code}.noncoding.sig_crb_input.${CNV}.tsv
    echo -e "${hpo}\tmeta_res/${prefix}.${freq_code}.${CNV}.crb_burden.meta_analysis.stats.bed.gz" \
    >> ${freq_code}.coding.sig_crb_input.${CNV}.tsv
  done < ${phenotype_list}
  # Extract significant CRBs
  /opt/rCNV2/analysis/noncoding/get_sig_crbs.py \
    --sumstats ${freq_code}.noncoding.sig_crb_input.${CNV}.tsv \
    --primary-p ${meta_p_cutoff} \
    --secondary-p 0.05 \
    --n-nominal 2 \
    --secondary-or-nominal \
    --coding-sumstats ${freq_code}.coding.sig_crb_input.${CNV}.tsv \
    --coding-p ${meta_p_cutoff} \
    --cnv ${CNV} \
    --outfile ${freq_code}.sig_CRBs.${CNV}.bed.gz \
    --bgzip
done
# DEV: get significant CRBs
# for CNV in DEL DUP; do
#   echo -e "\n\n\n${CNV}"
#   while read prefix hpo; do
#     zcat meta_res/${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed.gz \
#     | fgrep -v "#" \
#     | awk -v FS="\t" -v OFS="\t" -v prefix=$prefix -v CNV=$CNV \
#       '{ if ($13>=5.491278 && $13!="NA" && ($5>1 || ($18>=1.30103 && $18!="NA"))) print $1, $2, $3, $4, CNV, prefix, $9 }'
#   done < ${phenotype_list} \
#   | sort -k4,4V -Vk1,1 -k2,2n -k3,3n -k6,6V -k5,5V 
# done
