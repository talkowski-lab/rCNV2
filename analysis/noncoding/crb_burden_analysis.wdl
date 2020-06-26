######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens per cis-regulatory block (CRB)


import "https://api.firecloud.org/ga4gh/v1/tools/rCNV:scattered_crb_burden_perm_test/versions/2/plain-WDL/descriptor" as scattered_perm


workflow crb_burden_analysis {
  File phenotype_list
  File metacohort_list
  File metacohort_sample_table
  File crbs
  File crb_elements
  String noncoding_filter
  Int n_pheno_perms
  Int pad_controls
  String meta_model_prefix
  Float min_element_ovr
  Float min_frac_all_elements
  Float p_cutoff
  Int max_manhattan_phred_p
  Float meta_secondary_p_cutoff
  Int meta_nominal_cohorts_cutoff
  String rCNV_bucket
  String fisher_cache_string
  String perm_cache_string
  String meta_cache_string

  Array[Array[String]] phenotypes = read_tsv(phenotype_list)

  Array[String] cnv_types = ["DEL", "DUP"]

  # Scatter over phenotypes
  scatter ( pheno in phenotypes ) {
    # Run assocation tests per phenotype for rCNV
    call burden_test as rCNV_burden_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        noncoding_filter=noncoding_filter,
        crbs=crbs,
        crb_elements=crb_elements,
        pad_controls=pad_controls,
        min_element_ovr=min_element_ovr,
        min_frac_all_elements=min_frac_all_elements,
        p_cutoff=p_cutoff,
        max_manhattan_phred_p=max_manhattan_phred_p,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0],
        cache_string=fisher_cache_string
    }

    # Permute phenotypes to estimate empirical FDR
    call scattered_perm.scattered_crb_burden_perm_test as rCNV_perm_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        noncoding_filter=noncoding_filter,
        crbs=crbs,
        crb_elements=crb_elements,
        pad_controls=pad_controls,
        min_element_ovr=min_element_ovr,
        min_frac_all_elements=min_frac_all_elements,
        p_cutoff=p_cutoff,
        n_pheno_perms=n_pheno_perms,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0],
        cache_string=perm_cache_string
    }
  }


  # Determine appropriate primary P-value thresholds for meta-analysis
  scatter ( cnv in cnv_types ) {
    # Genome-wide, primary
    call calc_meta_p_cutoff as calc_genome_wide_cutoffs {
      input:
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        noncoding_filter=noncoding_filter,
        CNV=cnv,
        n_pheno_perms=n_pheno_perms,
        fdr_target=p_cutoff,
        rCNV_bucket=rCNV_bucket,
        dummy_completion_markers=rCNV_perm_test.completion_marker,
        fdr_table_suffix="empirical_genome_wide_pval",
        p_val_column_name="meta_phred_p"
    }
  }

  # Perform meta-analysis of rCNV association statistics
  scatter ( pheno in phenotypes ) {
    call meta_analysis as rCNV_meta_analysis {
      input:
        stats_beds=rCNV_burden_test.stats_beds,
        stats_bed_idxs=rCNV_burden_test.stats_bed_idxs,
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        noncoding_filter=noncoding_filter,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        max_manhattan_phred_p=max_manhattan_phred_p,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0],
        cache_string=meta_cache_string
    }
  }

  output {}
}


# Run burden test for a single phenotype for all metacohorts
task burden_test {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  File crbs
  File crb_elements
  String noncoding_filter
  Int pad_controls
  Float min_element_ovr
  Float min_frac_all_elements
  Float p_cutoff
  Int max_manhattan_phred_p
  String rCNV_bucket
  String prefix
  String cache_string

  command <<<
    set -e

    # Copy CNV data and other references
    mkdir cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/noncoding/* cleaned_cnv/
    mkdir refs/
    gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

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


    # Copy results to output bucket
    gsutil -m cp *.crb_burden.stats.bed.gz* \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/counts/"
    gsutil -m cp *.crb_burden.stats.bed.gz* \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.crb_burden.*.png \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/plots/"

    echo "${cache_string}" > completion.txt
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:cbeb85322bbd6b8e249c1e9135fc860d9fc34cc7758e45cc9790ee2d83a30f18"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] stats_beds = glob("*.crb_burden.stats.bed.gz")
    Array[File] stats_bed_idxs = glob("*.crb_burden.stats.bed.gz.tbi")
    File completion_token = "completion.txt"
  }
}


# Aggregate all permutation results to determine empirical P-value cutoff
task calc_meta_p_cutoff {
  File phenotype_list
  File metacohort_sample_table
  String freq_code
  String CNV
  String noncoding_filter
  Int n_pheno_perms
  Float fdr_target
  String rCNV_bucket
  Array[File] dummy_completion_markers #Must delocalize something or Cromwell will bypass permutation test
  String fdr_table_suffix
  String p_val_column_name

  command <<<
    set -e 

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

    # Also produce an optional table of flat Bonferroni P-value cutoffs
    awk -v pval=${fdr_target} -v FS="\t" -v OFS="\t" \
      '{ print $1, pval }' ${phenotype_list} \
    > crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.bonferroni_pval.hpo_cutoffs.tsv

    # Copy cutoff tables to output bucket
    gsutil -m cp crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}*tsv \
      "${rCNV_bucket}/analysis/analysis_refs/"
    gsutil -m cp crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.bonferroni_pval.hpo_cutoffs.tsv \
      "${rCNV_bucket}/analysis/analysis_refs/"
    gsutil -m cp crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}_permutation_results.png \
      "${rCNV_bucket}/analysis/crb_burden/empirical_fdr_plots/"
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:cbeb85322bbd6b8e249c1e9135fc860d9fc34cc7758e45cc9790ee2d83a30f18"
    preemptible: 1
    memory: "32 GB"
    disks: "local-disk 275 HDD"
    bootDiskSizeGb: "40"
  }

  output {
    File perm_results_plot = "crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}_permutation_results.png"
    File p_cutoff_table = "crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}.hpo_cutoffs.tsv"
    File p_cutoff_ladder = "crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}.ncase_cutoff_ladder.tsv"
    File bonferroni_cutoff_table = "crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.bonferroni_pval.hpo_cutoffs.tsv "
  }  
}

# Run meta-analysis (both weighted and raw) across metacohorts for a single phenotype
task meta_analysis {
  Array[Array[File]] stats_beds
  Array[Array[File]] stats_bed_idxs
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  String noncoding_filter
  Array[File] meta_p_cutoff_tables
  Int max_manhattan_phred_p
  String meta_model_prefix
  String rCNV_bucket
  String prefix
  String cache_string

  command <<<
    set -e

    # Copy burden counts & gene coordinates
    find / -name "*${prefix}.${freq_code}.${noncoding_filter}_noncoding.*.crb_burden.stats.bed.gz*" \
    | xargs -I {} mv {} ./
    find / -name "*crb_burden.${freq_code}.${noncoding_filter}_noncoding.*.bonferroni_pval.hpo_cutoffs.tsv*" \
    | xargs -I {} mv {} ./
    mkdir refs/
    gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

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

    # Copy results to output bucket
    gsutil -m cp *.crb_burden.*meta_analysis.stats.bed.gz* \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.crb_burden.*or_corplot_grid.jpg \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/plots/"
    gsutil -m cp *.crb_burden.*meta_analysis.*.png \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/plots/"

    # Must delocalize completion marker to prevent caching of final step
    echo "${cache_string}" > completion.txt
  >>>

  output {
    File completion_token = "completion.txt"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:cbeb85322bbd6b8e249c1e9135fc860d9fc34cc7758e45cc9790ee2d83a30f18"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }
}
