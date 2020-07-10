#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control noncoding CNV burdens for cis-regulatory blocks (CRBs)


# Launch docker image
docker run --rm -it talkowski/rcnv
gcloud auth login


# Copy all filtered CNV data, gene coordinates, and other references 
# from the project Google Bucket (note: requires permissions)
mkdir cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/noncoding/* cleaned_cnv/
gsutil -m cp gs://rcnv_project/cleaned_data/cnv/*bed.gz* cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/genome_annotations/*bed.gz* ./
mkdir refs/
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/


# Test/dev parameters
# Params for seizures
hpo="HP:0001250"
prefix="HP0001250"
meta="meta2"
# Params for NDDs
hpo="HP:0012759"
prefix="HP0012759"
meta="meta1"
# Params for UNKNOWN
hpo="UNKNOWN"
prefix="UNKNOWN"
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
p_cutoff=0.000003226431
meta_p_cutoff=0.000003226431
max_manhattan_phred_p=30
n_pheno_perms=50
meta_model_prefix="fe"
i=1


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
      /opt/rCNV2/analysis/noncoding/crb_burden_test.R \
        --pheno-table ${metacohort_sample_table} \
        --cohort-name $meta \
        --cnv $CNV \
        --case-hpo ${hpo} \
        --bgzip \
        "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz" \
        "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"
      tabix -f "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "fisher_phred_p" \
        --p-is-phred \
        --max-phred-p ${max_manhattan_phred_p} \
        --cutoff ${p_cutoff} \
        --label-prefix "$CNV" \
        --title "$title" \
        "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden"
    done

    # Generate Miami & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --miami \
      --p-col-name "fisher_phred_p" \
      --p-is-phred \
      --max-phred-p ${max_manhattan_phred_p} \
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
      /opt/rCNV2/analysis/noncoding/crb_burden_test.R \
        --pheno-table ${metacohort_sample_table} \
        --cohort-name $meta \
        --cnv $CNV \
        --case-hpo ${hpo} \
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

    # Shuffle phenotypes for each metacohort CNV dataset, and restrict CNVs from
    # phenotype of relevance
    if ! [ -e shuffled_cnv/ ]; then
      mkdir shuffled_cnv/
    fi
    yes $i | head -n1000000 > seed_$i.txt
    while read meta cohorts; do

      cnvbed="cleaned_cnv/$meta.${freq_code}.${noncoding_filter}_noncoding.bed.gz"

      # Determine CNV size quintiles
      for CNV in DEL DUP; do
        zcat $cnvbed \
        | awk -v CNV="$CNV" '{ if ($1 !~ "#" && $5==CNV) print $3-$2 }' \
        | /opt/rCNV2/utils/quantiles.py \
          --quantiles "0,0.2,0.4,0.6,0.8,1.0" \
          --no-header \
        > $meta.$CNV.quantiles.tsv
      done

      # Shuffle phenotypes, matching by CNV type and size quintile
      for CNV in DEL DUP; do
        for qr in $( seq 1 5 ); do
          smin=$( awk -v qr="$qr" '{ if (NR==qr) print $2 }' $meta.$CNV.quantiles.tsv )
          smax=$( awk -v qr="$qr" '{ if (NR==(qr+1)) print $2 + 1 }' $meta.$CNV.quantiles.tsv )
          zcat $cnvbed | sed '1d' \
          | awk -v smin="$smin" -v smax="$smax" -v CNV="$CNV" \
            '{ if ($3-$2>=smin && $3-$2<smax && $5==CNV) print $0 }' \
          > cnv_subset.bed
          paste <( cut -f1-5 cnv_subset.bed ) \
                <( cut -f6 cnv_subset.bed | shuf --random-source seed_$i.txt ) \
          | awk -v hpo=${hpo} '{ if ($NF ~ "HEALTHY_CONTROL" || $NF ~ hpo) print $0 }'
          rm cnv_subset.bed
        done
      done \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat <( tabix -H $cnvbed ) - \
      | bgzip -c \
      > shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz
      tabix -f shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz

    done < ${metacohort_list}

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
        /opt/rCNV2/analysis/noncoding/crb_burden_test.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --cnv $CNV \
          --case-hpo ${hpo} \
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
      /opt/rCNV2/analysis/genes/gene_meta_analysis.R \
        --model ${meta_model_prefix} \
        --p-is-phred \
        --spa \
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
  --fdr-target ${fdr_target} \
  --linear-fit \
  --flat-ladder \
  --plot crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}_permutation_results.png \
  ${freq_code}.${noncoding_filter}_noncoding.${CNV}.permuted_pval_matrix.txt.gz \
  ${metacohort_sample_table} \
  crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}




# Run meta-analysis for each phenotype
while read prefix hpo; do

  # Get metadata for meta-analysis
  mega_idx=$( head -n1 "${metacohort_sample_table}" \
              | sed 's/\t/\n/g' \
              | awk '{ if ($1=="mega") print NR }' )
  ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
           | awk -v FS="\t" -v mega_idx="$mega_idx" '{ print $mega_idx }' \
           | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
  nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
           | awk -v FS="\t" -v mega_idx="$mega_idx" '{ print $mega_idx }' \
           | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
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
    /opt/rCNV2/analysis/noncoding/crb_meta_analysis.R \
      --or-corplot ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.or_corplot_grid.jpg \
      --model ${meta_model_prefix} \
      --p-is-phred \
      --spa \
      ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.input.txt \
      ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed
    bgzip -f ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed
    tabix -f ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed.gz

    # Generate Manhattan & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --p-col-name "meta_phred_p" \
      --p-is-phred \
      --max-phred-p ${max_manhattan_phred_p} \
      --cutoff $meta_p_cutoff \
      --label-prefix "$CNV" \
      --title "$title" \
      "${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed.gz" \
      "${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis"
  done

  # Generate Miami & QQ plots
  /opt/rCNV2/utils/plot_manhattan_qq.R \
    --miami \
    --p-col-name "meta_phred_p" \
    --p-is-phred \
    --max-phred-p ${max_manhattan_phred_p} \
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
p_val_column_name="meta_phred_p"
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
p_val_column_name="meta_phred_p_secondary"
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
    
    /opt/rCNV2/analysis/noncoding/crb_meta_analysis.R \
      --or-corplot ${prefix}.${freq_code}.$CNV.crb_burden.or_corplot_grid.jpg \
      --model ${meta_model_prefix} \
      --p-is-phred \
      --spa \
      ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.input.txt \
      ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.stats.bed
    bgzip -f ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.stats.bed
    tabix -f ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.stats.bed.gz

  done

done < refs/test_phenotypes.list




# Extract CRBs meeting all significance criteria
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
    --cnv ${CNV} \
    --outfile ${freq_code}.sig_CRBs.${CNV}.bed.gz \
    --bgzip
    # --coding-sumstats ${freq_code}.coding.sig_crb_input.${CNV}.tsv \
    # --coding-p ${meta_p_cutoff} \
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
