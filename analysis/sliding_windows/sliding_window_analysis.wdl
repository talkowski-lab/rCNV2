######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens in sliding windows, genome-wide


import "https://api.firecloud.org/ga4gh/v1/tools/rCNV:scattered_sliding_window_perm_test/versions/14/plain-WDL/descriptor" as scattered_perm


workflow sliding_window_analysis {
  File phenotype_list
  File metacohort_list
  File metacohort_sample_table
  File binned_genome
  Float bin_overlap
  Int pad_controls
  Float p_cutoff
  Int n_pheno_perms
  String meta_model_prefix
  Float meta_secondary_p_cutoff
  Float meta_or_cutoff
  Int meta_nominal_cohorts_cutoff
  Float credible_interval
  Int sig_window_pad
  Int refine_max_cnv_size
  File contigfile
  String rCNV_bucket

  Array[Array[String]] phenotypes = read_tsv(phenotype_list)

  Array[Array[String]] contigs = read_tsv(contigfile)

  Array[String] cnv_types = ["DEL", "DUP"]

  # Scatter over phenotypes
  scatter ( pheno in phenotypes ) {
    # Run rCNV assocation tests per phenotype
    call burden_test as rCNV_burden_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        binned_genome=binned_genome,
        bin_overlap=bin_overlap,
        pad_controls=pad_controls,
        p_cutoff=p_cutoff,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0]
    }

    # Permute phenotypes to estimate empirical FDR
    call scattered_perm.scattered_sliding_window_perm_test as rCNV_perm_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        binned_genome=binned_genome,
        bin_overlap=bin_overlap,
        pad_controls=pad_controls,
        p_cutoff=p_cutoff,
        n_pheno_perms=n_pheno_perms,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0]
    }
  }
  
  # Determine appropriate P-value thresholds for primary and secondary meta-analysis
  scatter ( cnv in cnv_types ) {
    # Genome-wide, primary
    call calc_meta_p_cutoff as calc_genome_wide_cutoffs {
      input:
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        n_pheno_perms=n_pheno_perms,
        fdr_target=p_cutoff,
        rCNV_bucket=rCNV_bucket,
        dummy_completion_markers=rCNV_perm_test.completion_marker,
        fdr_table_suffix="empirical_genome_wide_pval",
        p_val_column_name="meta_phred_p"
    }

    # 1% FDR, primary
    call calc_meta_p_cutoff as calc_fdr_1pct_cutoffs {
      input:
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        n_pheno_perms=n_pheno_perms,
        fdr_target=0.01,
        rCNV_bucket=rCNV_bucket,
        dummy_completion_markers=rCNV_perm_test.completion_marker,
        fdr_table_suffix="empirical_fdr_1pct_pval",
        p_val_column_name="meta_phred_p"
    }

    # 5% FDR, primary
    call calc_meta_p_cutoff as calc_fdr_5pct_cutoffs {
      input:
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        n_pheno_perms=n_pheno_perms,
        fdr_target=0.05,
        rCNV_bucket=rCNV_bucket,
        dummy_completion_markers=rCNV_perm_test.completion_marker,
        fdr_table_suffix="empirical_fdr_5pct_pval",
        p_val_column_name="meta_phred_p"
    }

    # Genome-wide, secondary
    call calc_meta_p_cutoff as calc_genome_wide_cutoffs_secondary {
      input:
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        n_pheno_perms=n_pheno_perms,
        fdr_target=p_cutoff,
        rCNV_bucket=rCNV_bucket,
        dummy_completion_markers=rCNV_perm_test.completion_marker,
        fdr_table_suffix="empirical_genome_wide_pval_secondary",
        p_val_column_name="meta_phred_p_secondary"
    }

    # 1% FDR, secondary
    call calc_meta_p_cutoff as calc_fdr_1pct_cutoffs_secondary {
      input:
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        n_pheno_perms=n_pheno_perms,
        fdr_target=0.01,
        rCNV_bucket=rCNV_bucket,
        dummy_completion_markers=rCNV_perm_test.completion_marker,
        fdr_table_suffix="empirical_fdr_1pct_pval_secondary",
        p_val_column_name="meta_phred_p_secondary"
    }

    # 5% FDR, secondary
    call calc_meta_p_cutoff as calc_fdr_5pct_cutoffs_secondary {
      input:
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        n_pheno_perms=n_pheno_perms,
        fdr_target=0.05,
        rCNV_bucket=rCNV_bucket,
        dummy_completion_markers=rCNV_perm_test.completion_marker,
        fdr_table_suffix="empirical_fdr_5pct_pval_secondary",
        p_val_column_name="meta_phred_p_secondary"
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
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.p_cutoff_table,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0]
    }
  }

  # Refine minimal credible regions
  scatter ( cnv in cnv_types ) {
    call prep_refinement {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_list=metacohort_list,
        binned_genome=binned_genome,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.p_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_or_cutoff=meta_or_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        sig_window_pad=sig_window_pad,
        rCNV_bucket=rCNV_bucket
    }
  }
  scatter ( contig in contigs ) {
    # DEL
    call refine_regions as refine_rCNV_regions_DEL {
      input:
        contig=contig[0],
        regions_to_refine=prep_refinement.regions_to_refine[0],
        metacohort_info_tsv=prep_refinement.metacohort_info_tsv[0],
        pval_matrix=prep_refinement.pval_matrix[0],
        labeled_windows=prep_refinement.labeled_windows[0],
        freq_code="rCNV",
        CNV="DEL",
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.p_cutoff_table,
        meta_p_ladder_cutoff_tables=calc_genome_wide_cutoffs.p_cutoff_ladder,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_or_cutoff=meta_or_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        meta_model_prefix=meta_model_prefix,
        credible_interval=credible_interval,
        refine_max_cnv_size=refine_max_cnv_size,
        rCNV_bucket=rCNV_bucket
    }
    # DUP
    call refine_regions as refine_rCNV_regions_DUP {
      input:
        contig=contig[0],
        regions_to_refine=prep_refinement.regions_to_refine[1],
        metacohort_info_tsv=prep_refinement.metacohort_info_tsv[1],
        pval_matrix=prep_refinement.pval_matrix[1],
        labeled_windows=prep_refinement.labeled_windows[1],
        freq_code="rCNV",
        CNV="DUP",
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.p_cutoff_table,
        meta_p_ladder_cutoff_tables=calc_genome_wide_cutoffs.p_cutoff_ladder,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_or_cutoff=meta_or_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        meta_model_prefix=meta_model_prefix,
        credible_interval=credible_interval,
        refine_max_cnv_size=refine_max_cnv_size,
        rCNV_bucket=rCNV_bucket
    }
  }
  call merge_refinements as merge_refinements_DEL {
    input:
      loci=refine_rCNV_regions_DEL.loci,
      associations=refine_rCNV_regions_DEL.associations,
      logfiles=refine_rCNV_regions_DEL.logfile,
      freq_code="rCNV",
      CNV="DEL",
      rCNV_bucket=rCNV_bucket
  }
  call merge_refinements as merge_refinements_DUP {
    input:
      loci=refine_rCNV_regions_DUP.loci,
      associations=refine_rCNV_regions_DUP.associations,
      logfiles=refine_rCNV_regions_DUP.logfile,
      freq_code="rCNV",
      CNV="DUP",
      rCNV_bucket=rCNV_bucket
  }

  # Plot summary metrics for final credible regions
  call plot_region_summary as plot_rCNV_regions {
    input:
      freq_code="rCNV",
      DEL_regions=merge_refinements_DEL.final_loci,
      DUP_regions=merge_refinements_DUP.final_loci,
      rCNV_bucket=rCNV_bucket
  }

  output {
    Array[File] final_sig_regions = [merge_refinements_DEL.final_loci, merge_refinements_DUP.final_loci]
    Array[File] final_sig_associations = [merge_refinements_DEL.final_associations, merge_refinements_DUP.final_associations]
  }
}


# Run burden test for a single phenotype for all metacohorts
task burden_test {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  File binned_genome
  Float bin_overlap
  Int pad_controls
  Float p_cutoff
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Copy CNV data
    mkdir cleaned_cnv/
    gsutil -m cp ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/

    # Iterate over metacohorts
    while read meta cohorts; do

      # Set metacohort parameters
      cnv_bed="cleaned_cnv/$meta.${freq_code}.bed.gz"
      descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
                 | awk -v FS="\t" '{ print $2 }' )
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
        # Set CNV-specific parameters
        case "$CNV" in
          DEL)
            highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz
            highlight_title="Known DEL GDs (Owen 2018)"
            ;;
          DUP)
            highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.DUP.bed.gz
            highlight_title="Known DUP GDs (Owen 2018)"
            ;;
        esac

        # Count CNVs
        /opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
          --fraction ${bin_overlap} \
          --pad-controls ${pad_controls} \
          -t $CNV \
          --hpo ${hpo} \
          -z \
          -o "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz" \
          $cnv_bed \
          ${binned_genome}
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz"

        # Perform burden test
        /opt/rCNV2/analysis/sliding_windows/window_burden_test.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --case-hpo ${hpo} \
          --bgzip \
          "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"

        # Generate Manhattan & QQ plots
        /opt/rCNV2/utils/plot_manhattan_qq.R \
          --p-col-name "fisher_phred_p" \
          --p-is-phred \
          --max-phred-p 100 \
          --cutoff ${p_cutoff} \
          --highlight-bed "$highlight_bed" \
          --highlight-name "$highlight_title" \
          --label-prefix "$CNV" \
          --title "$title" \
          "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.sliding_window"
      done

      # Generate Miami & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --miami \
        --p-col-name "fisher_phred_p" \
        --p-is-phred \
        --max-phred-p 100 \
        --cutoff ${p_cutoff} \
        --highlight-bed /opt/rCNV2/refs/UKBB_GD.Owen_2018.DUP.bed.gz \
        --highlight-name "Known DUP GDs (Owen 2018)" \
        --label-prefix "DUP" \
        --highlight-bed-2 /opt/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz \
        --highlight-name-2 "Known DEL GDs (Owen 2018)" \
        --label-prefix-2 "DEL" \
        --title "$title" \
        "$meta.${prefix}.${freq_code}.DUP.sliding_window.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.DEL.sliding_window.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.sliding_window"
    done < ${metacohort_list}

    # Copy results to output bucket
    gsutil -m cp *.sliding_window.stats.bed.gz* \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.sliding_window.*.png \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/plots/"
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:3f404f25d6821167744a6041b1159ece926d3720d0f33c75abbad1178775b242"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] stats_beds = glob("*.sliding_window.stats.bed.gz")
    Array[File] stats_bed_idxs = glob("*.sliding_window.stats.bed.gz.tbi")
    # Array[File] plots = glob("*.sliding_window.*.png")
  }
}


# Aggregate all permutation results to determine empirical P-value cutoff
task calc_meta_p_cutoff {
  File phenotype_list
  File metacohort_sample_table
  String freq_code
  String CNV
  Int n_pheno_perms
  Float fdr_target
  String rCNV_bucket
  Array[File] dummy_completion_markers #Must delocalize something or Cromwell will bypass permutation test
  String fdr_table_suffix
  String p_val_column_name

  command <<<
    # Gather all permutation results and compute FDR CDFs
    mkdir perm_res/
    while read prefix hpo; do
      echo -e "$prefix\n\n"
      gsutil -m cp \
        "${rCNV_bucket}/analysis/sliding_windows/$prefix/${freq_code}/permutations/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.perm_*.bed.gz" \
        perm_res/
      for i in $( seq 1 ${n_pheno_perms} ); do
        p_idx=$( zcat perm_res/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.perm_$i.bed.gz \
                 | sed -n '1p' | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' \
                 | fgrep -w ${p_val_column_name} | cut -f2 )
        zcat perm_res/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.perm_$i.bed.gz \
        | grep -ve '^#' \
        | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
        | cat <( echo "$prefix.${CNV}.$i" ) - \
        > perm_res/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.permuted_p_values.$i.txt
      done
    done < ${phenotype_list}
    paste perm_res/*.sliding_window.meta_analysis.permuted_p_values.*.txt \
    | gzip -c \
    > ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz

    # Analyze p-values and compute FDR
    /opt/rCNV2/analysis/sliding_windows/calc_empirical_fdr.R \
      --cnv ${CNV} \
      --fdr-target ${fdr_target} \
      --plot sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
      ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz \
      ${metacohort_sample_table} \
      sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}

    # Copy cutoff tables to output bucket
    gsutil -m cp sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}*tsv \
      "${rCNV_bucket}/analysis/analysis_refs/"
    gsutil -m cp sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
      "${rCNV_bucket}/analysis/sliding_windows/empirical_fdr_plots/"
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:01175c8d4c0d1ca32ae57b35c85902e0eaee9427eaeca5b9f9b57c733267453c"
    preemptible: 1
    memory: "16 GB"
    disks: "local-disk 100 HDD"
    bootDiskSizeGb: "20"
  }

  output {
    File perm_results_plot = "sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png"
    File p_cutoff_table = "sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}.hpo_cutoffs.tsv"
    File p_cutoff_ladder = "sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}.ncase_cutoff_ladder.tsv"
  }  
}


# Run meta-analysis across metacohorts for a single phenotype
task meta_analysis {
  Array[Array[File]] stats_beds
  Array[Array[File]] stats_bed_idxs
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  Array[File] meta_p_cutoff_tables
  String meta_model_prefix
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Copy burden stats & p-value cutoff tables
    find / -name "*${prefix}.${freq_code}.*.sliding_window.stats.bed.gz*" \
    | xargs -I {} mv {} ./
    find / -name "*sliding_window.${freq_code}.*.empirical_genome_wide_pval.hpo_cutoffs.tsv*" \
    | xargs -I {} mv {} ./
    # gsutil -m cp \
    #   ${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/** \
    #   ./

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
                    sliding_window.${freq_code}.DEL.empirical_genome_wide_pval.hpo_cutoffs.tsv )
    DUP_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                    sliding_window.${freq_code}.DUP.empirical_genome_wide_pval.hpo_cutoffs.tsv )

    # Run meta-analysis for each CNV type
    for CNV in DEL DUP; do
      # Set CNV-specific parameters
      case "$CNV" in
        DEL)
          highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz
          highlight_title="Known DEL GDs (Owen 2018)"
          meta_p_cutoff=$DEL_p_cutoff
          ;;
        DUP)
          highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.DUP.bed.gz
          highlight_title="Known DUP GDs (Owen 2018)"
          meta_p_cutoff=$DUP_p_cutoff
          ;;
      esac

      # Perform meta-analysis
      while read meta cohorts; do
        echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"
      done < <( fgrep -v mega ${metacohort_list} ) \
      > ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt
      /opt/rCNV2/analysis/sliding_windows/window_meta_analysis.R \
        --or-corplot ${prefix}.${freq_code}.$CNV.sliding_window.or_corplot_grid.jpg \
        --model ${meta_model_prefix} \
        --p-is-phred \
        ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt \
        ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed
      bgzip -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed
      tabix -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "meta_phred_p" \
        --p-is-phred \
        --cutoff $meta_p_cutoff \
        --highlight-bed "$highlight_bed" \
        --highlight-name "$highlight_title" \
        --label-prefix "$CNV" \
        --title "$title" \
        "${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz" \
        "${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis"
    done

    # Generate Miami & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --miami \
      --p-col-name "meta_phred_p" \
      --p-is-phred \
      --cutoff $DUP_p_cutoff \
      --highlight-bed /opt/rCNV2/refs/UKBB_GD.Owen_2018.DUP.bed.gz \
      --highlight-name "Known DUP GDs (Owen 2018)" \
      --label-prefix "DUP" \
      --cutoff-2 $DEL_p_cutoff \
      --highlight-bed-2 /opt/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz \
      --highlight-name-2 "Known DEL GDs (Owen 2018)" \
      --label-prefix-2 "DEL" \
      --title "$title" \
      "${prefix}.${freq_code}.DUP.sliding_window.meta_analysis.stats.bed.gz" \
      "${prefix}.${freq_code}.DEL.sliding_window.meta_analysis.stats.bed.gz" \
      "${prefix}.${freq_code}.sliding_window.meta_analysis"

    # Copy results to output bucket
    gsutil -m cp *.sliding_window.meta_analysis.stats.bed.gz* \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.sliding_window.or_corplot_grid.jpg \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/plots/"
    gsutil -m cp *.sliding_window.meta_analysis.*.png \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/plots/"

    # Must delocalize completion marker to prevent caching of final step
    echo "DONE" > completion.txt
  >>>

  output {
    File completion_token = "completion.txt"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:01175c8d4c0d1ca32ae57b35c85902e0eaee9427eaeca5b9f9b57c733267453c"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }
}


# Prepare data for refinement
task prep_refinement {
  Array[File] completion_tokens # Must delocalize something from meta-analysis step to prevent caching
  File phenotype_list
  File metacohort_list
  File binned_genome
  String freq_code
  String CNV
  Array[File] meta_p_cutoff_tables
  Float meta_secondary_p_cutoff
  Float meta_or_cutoff
  Int meta_nominal_cohorts_cutoff
  Int sig_window_pad
  String rCNV_bucket

  command <<<
    set -e

    # Copy p-value cutoff tables
    find / -name "*sliding_window.${freq_code}.*.empirical_genome_wide_pval.hpo_cutoffs.tsv*" \
    | xargs -I {} mv {} ./

    # Download all meta-analysis stats files and necessary data
    mkdir cleaned_cnv/
    gsutil -m cp -r ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/
    mkdir stats/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/sliding_windows/**.${freq_code}.**.sliding_window.meta_analysis.stats.bed.gz \
      stats/
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/* refs/
    mkdir phenos/
    gsutil -m cp ${rCNV_bucket}/cleaned_data/phenotypes/filtered/* phenos/

    # Iterate over phenotypes and make matrix of p-values, odds ratios (lower 95% CI), and nominal sig cohorts
    mkdir pvals/
    mkdir secondary_pvals/
    mkdir ors/
    mkdir nomsig/
    while read pheno hpo; do
      stats=stats/$pheno.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.bed.gz
      p_idx=$( zcat $stats | sed -n '1p' | sed 's/\t/\n/g' \
             | awk -v OFS="\t" '{ if ($1=="meta_phred_p") print NR }' )
      secondary_idx=$( zcat $stats | sed -n '1p' | sed 's/\t/\n/g' \
                       | awk -v OFS="\t" '{ if ($1=="meta_phred_p_secondary") print NR }' )
      lnor_lower_idx=$( zcat $stats | sed -n '1p' | sed 's/\t/\n/g' \
                        | awk -v OFS="\t" '{ if ($1=="meta_lnOR_lower") print NR }' )
      nom_idx=$( zcat $stats | sed -n '1p' | sed 's/\t/\n/g' \
                 | awk -v OFS="\t" '{ if ($1=="n_nominal_cohorts") print NR }' )
      zcat $stats | awk -v FS="\t" -v idx=$p_idx '{ if ($1 !~ "#") print $(idx) }' \
      | cat <( echo "$pheno.${CNV}" ) - \
      > pvals/$pheno.${CNV}.pvals.txt
      zcat $stats | awk -v FS="\t" -v idx=$secondary_idx '{ if ($1 !~ "#") print $(idx) }' \
      | cat <( echo "$pheno.${CNV}" ) - \
      > pvals/$pheno.${CNV}.secondary_pvals.txt
      zcat $stats | awk -v FS="\t" -v idx=$lnor_lower_idx '{ if ($1 !~ "#") print $(idx) }' \
      | cat <( echo "$pheno.${CNV}" ) - \
      > ors/$pheno.${CNV}.lnOR_lower.txt
      zcat $stats | awk -v FS="\t" -v idx=$nom_idx '{ if ($1 !~ "#") print $(idx) }' \
      | cat <( echo "$pheno.${CNV}" ) - \
      > nomsig/$pheno.${CNV}.nomsig_counts.txt
    done < ${phenotype_list}
    paste <( zcat ${binned_genome} | cut -f1-3 ) \
          pvals/*.${CNV}.pvals.txt \
    | bgzip -c \
    > ${CNV}.pval_matrix.bed.gz
    paste <( zcat ${binned_genome} | cut -f1-3 ) \
          pvals/*.${CNV}.secondary_pvals.txt \
    | bgzip -c \
    > ${CNV}.secondary_pval_matrix.bed.gz
    paste <( zcat ${binned_genome} | cut -f1-3 ) \
          ors/*.${CNV}.lnOR_lower.txt \
    | bgzip -c \
    > ${CNV}.lnOR_lower_matrix.bed.gz
    paste <( zcat ${binned_genome} | cut -f1-3 ) \
          nomsig/*.${CNV}.nomsig_counts.txt \
    | bgzip -c \
    > ${CNV}.nominal_cohort_counts.bed.gz

    # Get matrix of window significance labels
    /opt/rCNV2/analysis/sliding_windows/get_significant_windows.R \
      --pvalues ${CNV}.pval_matrix.bed.gz \
      --secondary-pvalues ${CNV}.secondary_pval_matrix.bed.gz \
      --p-is-phred \
      --p-cutoffs sliding_window.${freq_code}.${CNV}.empirical_genome_wide_pval.hpo_cutoffs.tsv \
      --odds-ratios ${CNV}.lnOR_lower_matrix.bed.gz \
      --or-is-ln \
      --min-secondary-p ${meta_secondary_p_cutoff} \
      --min-or ${meta_or_cutoff} \
      --nominal-counts ${CNV}.nominal_cohort_counts.bed.gz \
      --min-nominal ${meta_nominal_cohorts_cutoff} \
      --out-prefix ${freq_code}.${CNV}. \
      ${binned_genome}
    bgzip -f ${freq_code}.${CNV}.all_windows_labeled.bed
    bgzip -f ${freq_code}.${CNV}.significant_windows.bed

    # Define regions to be refined (sig windows padded by $sig_window_pad and merged)
    zcat ${freq_code}.${CNV}.significant_windows.bed.gz \
    | fgrep -v "#" \
    | awk -v buf=${sig_window_pad} -v OFS="\t" '{ print $1, $2-buf, $3+buf }' \
    | awk -v OFS="\t" '{ if ($2<0) $2=0; print $1, $2, $3 }' \
    | sort -Vk1,1 -k2,2V -k3,3V \
    | bedtools merge -i - \
    | bgzip -c \
    > ${freq_code}.${CNV}.sig_regions_to_refine.bed.gz

    # Prep input file
    while read meta; do
      echo -e "$meta\tcleaned_cnv/$meta.${freq_code}.bed.gz\tphenos/$meta.cleaned_phenos.txt"
    done < <( cut -f1 ${metacohort_list} | fgrep -v "mega" )\
    > window_refinement.${freq_code}_metacohort_info.tsv
  >>>

  output {
    File regions_to_refine = "${freq_code}.${CNV}.sig_regions_to_refine.bed.gz"
    File metacohort_info_tsv = "window_refinement.${freq_code}_metacohort_info.tsv"
    File pval_matrix = "${CNV}.pval_matrix.bed.gz"
    File labeled_windows = "${freq_code}.${CNV}.all_windows_labeled.bed.gz"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:e0c7c3d285a39bb285f655b02747336d6497acd1bd5783d83ec93a53a1e9ea89"
    preemptible: 1
    memory: "8 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }
}


# # Refine associated regions to minimal credible regions
task refine_regions {
  String contig
  File regions_to_refine
  File metacohort_info_tsv
  File pval_matrix
  File labeled_windows
  String freq_code
  String CNV
  Array[File] meta_p_cutoff_tables
  Array[File] meta_p_ladder_cutoff_tables
  Float meta_secondary_p_cutoff
  Float meta_or_cutoff
  Int meta_nominal_cohorts_cutoff
  String meta_model_prefix
  Float credible_interval
  Int refine_max_cnv_size
  String rCNV_bucket

  command <<<
    set -e

    # Copy p-value cutoff tables
    find / -name "*sliding_window.${freq_code}.*.empirical_genome_wide_pval.hpo_cutoffs.tsv" \
    | xargs -I {} mv {} ./
    find / -name "*sliding_window.${freq_code}.*.empirical_genome_wide_pval.ncase_cutoff_ladder.tsv" \
    | xargs -I {} mv {} ./

    # Download all meta-analysis stats files and necessary data
    mkdir cleaned_cnv/
    gsutil -m cp -r ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/
    mkdir stats/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/sliding_windows/**.${freq_code}.**.sliding_window.meta_analysis.stats.bed.gz \
      stats/
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/* refs/
    mkdir phenos/
    gsutil -m cp ${rCNV_bucket}/cleaned_data/phenotypes/filtered/* phenos/

    # Tabix input to single chromosome
    tabix -f ${regions_to_refine}
    tabix -h ${regions_to_refine} ${contig} | bgzip -c > regions_to_refine.bed.gz
    tabix -f ${pval_matrix}
    tabix -h ${pval_matrix} ${contig} | bgzip -c > pval_matrix.bed.gz
    tabix -f ${labeled_windows}
    tabix -h ${labeled_windows} ${contig} | bgzip -c > labeled_windows.bed.gz

    # Perform refinement
    if [ $( zcat regions_to_refine.bed.gz | fgrep -v "#" | wc -l ) -gt 0 ]; then
      /opt/rCNV2/analysis/sliding_windows/refine_significant_regions.py \
        --cnv-type ${CNV} \
        --model ${meta_model_prefix} \
        --hpo-p-cutoffs sliding_window.${freq_code}.${CNV}.empirical_genome_wide_pval.hpo_cutoffs.tsv \
        --p-cutoff-ladder sliding_window.${freq_code}.${CNV}.empirical_genome_wide_pval.ncase_cutoff_ladder.tsv \
        --p-is-phred \
        --secondary-p-cutoff ${meta_secondary_p_cutoff} \
        --min-or-lower ${meta_or_cutoff} \
        --retest-min-or-lower ${meta_or_cutoff} \
        --max-cnv-size ${refine_max_cnv_size} \
        --min-nominal ${meta_nominal_cohorts_cutoff} \
        --credible-interval ${credible_interval} \
        --prefix "${freq_code}_${CNV}" \
        --log ${freq_code}.${CNV}.region_refinement.${contig}.log \
        regions_to_refine.bed.gz \
        ${metacohort_info_tsv} \
        pval_matrix.bed.gz \
        labeled_windows.bed.gz \
        ${freq_code}.${CNV}.final_regions.associations.${contig}.bed \
        ${freq_code}.${CNV}.final_regions.loci.${contig}.bed
    else
      touch ${freq_code}.${CNV}.final_regions.associations.${contig}.bed
      touch ${freq_code}.${CNV}.final_regions.loci.${contig}.bed
      touch ${freq_code}.${CNV}.region_refinement.${contig}.log
    fi
    bgzip -f ${freq_code}.${CNV}.final_regions.associations.${contig}.bed
    bgzip -f ${freq_code}.${CNV}.final_regions.loci.${contig}.bed
  >>>

  output {
    File loci = "${freq_code}.${CNV}.final_regions.loci.${contig}.bed.gz"
    File associations = "${freq_code}.${CNV}.final_regions.associations.${contig}.bed.gz"
    File logfile = "${freq_code}.${CNV}.region_refinement.${contig}.log"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:caaf44768dc13701026c761972685ca317788ccf2e749f4beffdb3823cd39d7c"
    preemptible: 1
    memory: "8 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }
}


# Merge refined loci across all chromosomes
task merge_refinements {
  Array[File] loci
  Array[File] associations
  Array[File] logfiles
  String freq_code
  String CNV
  String rCNV_bucket

  command <<<
    set -e

    # Merge loci
    zcat ${loci[0]} | sed -n '1p' > header.tsv
    zcat ${sep=" " loci} | fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n \
    | cat header.tsv - | bgzip -c \
    > ${freq_code}.${CNV}.final_regions.loci.bed.gz

    # Merge associations
    zcat ${associations[0]} | sed -n '1p' > header.tsv
    zcat ${sep=" " associations} | fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n -k7,7V \
    | cat header.tsv - | bgzip -c \
    > ${freq_code}.${CNV}.final_regions.associations.bed.gz

    # Merge logfiles
    cat ${sep=" " logfiles} > ${freq_code}.${CNV}.region_refinement.log

    # Annotate final regions with genes
    gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes ./
    /opt/rCNV2/analysis/sliding_windows/get_genes_per_region.py \
      -o ${freq_code}.${CNV}.final_regions.loci.bed \
      ${freq_code}.${CNV}.final_regions.loci.bed.gz \
      genes/gencode.v19.canonical.gtf.gz
    bgzip -f ${freq_code}.${CNV}.final_regions.loci.bed

    # Copy results to output bucket
    gsutil -m cp ${freq_code}.${CNV}.final_regions.loci.bed.gz \
      "${rCNV_bucket}/results/sliding_windows/"
    gsutil -m cp ${freq_code}.${CNV}.final_regions.associations.bed.gz \
      "${rCNV_bucket}/results/sliding_windows/"
    gsutil -m cp ${freq_code}.${CNV}.region_refinement.log \
      "${rCNV_bucket}/results/sliding_windows/"
  >>>

  output {
    File final_loci = "${freq_code}.${CNV}.final_regions.loci.bed.gz"
    File final_associations = "${freq_code}.${CNV}.final_regions.associations.bed.gz"
    File merged_logfile = "${freq_code}.${CNV}.region_refinement.log"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:caaf44768dc13701026c761972685ca317788ccf2e749f4beffdb3823cd39d7c"
    preemptible: 1
    memory: "8 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }
}


task plot_region_summary {
  String freq_code
  File DEL_regions
  File DUP_regions
  String rCNV_bucket

  command <<<
    /opt/rCNV2/analysis/sliding_windows/regions_summary.plot.R \
      -o "${freq_code}.final_regions." \
      ${DEL_regions} \
      ${DUP_regions}

    gsutil -m cp \
      ./*.jpg \
      "${rCNV_bucket}/results/sliding_windows/plots/"
    gsutil -m cp \
      "${freq_code}.final_regions.multipanel_summary.jpg" \
      ${rCNV_bucket}/public/
    gsutil acl ch -u AllUsers:R ${rCNV_bucket}/public/*.jpg
  >>>

  output {
    File summary_plot = "${freq_code}.final_regions.multipanel_summary.jpg"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:caaf44768dc13701026c761972685ca317788ccf2e749f4beffdb3823cd39d7c"
    preemptible: 1
  }
}

