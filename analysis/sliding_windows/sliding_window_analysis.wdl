######################
#    rCNV Project    #
######################

# Copyright (c) 2019-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens in sliding windows, genome-wide


import "https://api.firecloud.org/ga4gh/v1/tools/rCNV:scattered_sliding_window_perm_test/versions/23/plain-WDL/descriptor" as scattered_perm


workflow sliding_window_analysis {
  File phenotype_list
  File metacohort_list
  File metacohort_sample_table
  File binned_genome
  Float bin_overlap
  Int pad_controls
  Float p_cutoff
  Int max_manhattan_phred_p
  Int n_pheno_perms
  Int min_probes_per_window
  Float min_frac_controls_probe_exclusion
  String meta_model_prefix
  Float meta_secondary_p_cutoff
  Float meta_or_cutoff
  Int meta_nominal_cohorts_cutoff
  Float credible_interval
  Int sig_window_pad
  Float FDR_cutoff
  File gtf
  File contigfile
  String rCNV_bucket
  String rCNV_docker
  String rCNV_docker_refine # This is specified separately for convenience, but could be converted to the same docker for final run
  String athena_cloud_docker
  String fisher_cache_string
  String perm_cache_string
  String meta_cache_string

  Array[Array[String]] phenotypes = read_tsv(phenotype_list)

  Array[Array[String]] contigs = read_tsv(contigfile)

  Array[String] cnv_types = ["DEL", "DUP"]

  # Determine conditional cohort exclusion list based on probe density
  call build_exclusion_list {
    input:
      binned_genome=binned_genome,
      binned_genome_prefix=basename(binned_genome, '.bed.gz'),
      min_probes_per_window=min_probes_per_window,
      min_frac_controls_probe_exclusion=min_frac_controls_probe_exclusion,
      metacohort_list=metacohort_list,
      rCNV_bucket=rCNV_bucket,
      athena_cloud_docker=athena_cloud_docker,
      freq_code="rCNV"
  }

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
        max_manhattan_phred_p=max_manhattan_phred_p,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker,
        prefix=pheno[0],
        cache_string=fisher_cache_string
    }

    # Permute phenotypes to estimate empirical FDR
    call scattered_perm.scattered_sliding_window_perm_test as rCNV_perm_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        exclusion_bed=build_exclusion_list.exclusion_bed,
        freq_code="rCNV",
        binned_genome=binned_genome,
        bin_overlap=bin_overlap,
        pad_controls=pad_controls,
        p_cutoff=p_cutoff,
        n_pheno_perms=n_pheno_perms,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker,
        prefix=pheno[0],
        cache_string=perm_cache_string
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
        rCNV_docker=rCNV_docker,
        dummy_completion_markers=rCNV_perm_test.completion_marker,
        fdr_table_suffix="empirical_genome_wide_pval",
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
        rCNV_docker=rCNV_docker,
        dummy_completion_markers=rCNV_perm_test.completion_marker,
        fdr_table_suffix="empirical_genome_wide_pval_secondary",
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
        exclusion_bed=build_exclusion_list.exclusion_bed,
        freq_code="rCNV",
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        max_manhattan_phred_p=max_manhattan_phred_p,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker,
        prefix=pheno[0],
        cache_string=meta_cache_string
    }
  }

  # Refine minimal credible regions
  call refine_regions as refine_DEL {
    input:
      completion_tokens=rCNV_meta_analysis.completion_token,
      phenotype_list=phenotype_list,
      metacohort_list=metacohort_list,
      metacohort_sample_table=metacohort_sample_table,
      freq_code="rCNV",
      CNV="DEL",
      meta_p_cutoffs_tsv=calc_genome_wide_cutoffs.bonferroni_cutoff_table[0],
      meta_secondary_p_cutoff=meta_secondary_p_cutoff,
      meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
      sig_window_pad=sig_window_pad,
      credset=credible_interval,
      FDR_cutoff=FDR_cutoff,
      gtf=gtf,
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker_refine
  }
  call refine_regions as refine_DUP {
    input:
      completion_tokens=rCNV_meta_analysis.completion_token,
      phenotype_list=phenotype_list,
      metacohort_list=metacohort_list,
      metacohort_sample_table=metacohort_sample_table,
      freq_code="rCNV",
      CNV="DUP",
      meta_p_cutoffs_tsv=calc_genome_wide_cutoffs.bonferroni_cutoff_table[1],
      meta_secondary_p_cutoff=meta_secondary_p_cutoff,
      meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
      sig_window_pad=sig_window_pad,
      credset=credible_interval,
      FDR_cutoff=FDR_cutoff,
      gtf=gtf,
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker_refine
  }

  # Merge refined associations & regions
  call merge_refined_regions {
    input:
      assoc_beds=[refine_DEL.associations, refine_DUP.associations],
      loci_beds=[refine_DEL.loci, refine_DUP.loci],
      freq_code="rCNV",
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker_refine
  }

  # Plot summary metrics for final credible regions
  call plot_region_summary as plot_rCNV_regions {
    input:
      freq_code="rCNV",
      DEL_regions=refine_DEL.loci,
      DUP_regions=refine_DUP.loci,
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker_refine
  }

  output {
    File final_sig_loci = merge_refined_regions.final_loci
    File final_sig_associations = merge_refined_regions.final_associations
  }
}


# Build list of cohorts to exclude per locus based on inadequate probe density
task build_exclusion_list {
  File binned_genome
  String binned_genome_prefix
  Int min_probes_per_window
  Float min_frac_controls_probe_exclusion
  File metacohort_list
  String rCNV_bucket
  String athena_cloud_docker
  String freq_code

  command <<<
    set -euo pipefail

    # Download control probesets
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/control_probesets \
      ./

    # Make input for conditional exclusion script
    for file in control_probesets/*bed.gz; do
      echo -e "$file\t$( basename $file | sed 's/\.bed\.gz//g' )"
    done > probeset_tracks.tsv

    # Clone rCNV2 repo (not present in athena-cloud Docker)
    git clone https://github.com/talkowski-lab/rCNV2.git

    # Build conditional exclusion list
    rCNV2/data_curation/other/probe_based_exclusion.py \
      --outfile ${binned_genome_prefix}.cohort_exclusion.bed.gz \
      --probecounts-outfile ${binned_genome_prefix}.probe_counts.bed.gz \
      --control-mean-counts-outfile ${binned_genome_prefix}.mean_probe_counts_per_cohort.bed.gz \
      --frac-pass-outfile ${binned_genome_prefix}.frac_passing.bed.gz \
      --min-probes ${min_probes_per_window} \
      --min-frac-samples ${min_frac_controls_probe_exclusion} \
      --keep-n-columns 3 \
      --bgzip \
      ${binned_genome} \
      probeset_tracks.tsv \
      control_probesets/rCNV.control_counts_by_array.tsv \
      <( fgrep -v mega ${metacohort_list} )

    # Copy to Google bucket for storage
    gsutil -m cp \
      ${binned_genome_prefix}.cohort_exclusion.bed.gz \
      ${binned_genome_prefix}.probe_counts.bed.gz \
      ${binned_genome_prefix}.frac_passing.bed.gz \
      ${rCNV_bucket}/analysis/analysis_refs/
  >>>

  runtime {
    docker: "${athena_cloud_docker}"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    File exclusion_bed = "${binned_genome_prefix}.cohort_exclusion.bed.gz"
    File probe_counts = "${binned_genome_prefix}.probe_counts.bed.gz"
    File control_mean_counts = "${binned_genome_prefix}.mean_probe_counts_per_cohort.bed.gz"
    File frac_passing = "${binned_genome_prefix}.frac_passing.bed.gz"
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
  Int max_manhattan_phred_p
  String rCNV_bucket
  String rCNV_docker
  String prefix
  String cache_string

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
        /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
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
          --max-phred-p ${max_manhattan_phred_p} \
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
        --max-phred-p ${max_manhattan_phred_p} \
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

    echo "${cache_string}" > completion.txt
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] stats_beds = glob("*.sliding_window.stats.bed.gz")
    Array[File] stats_bed_idxs = glob("*.sliding_window.stats.bed.gz.tbi")
    File completion_token = "completion.txt"
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
  String rCNV_docker
  Array[File] dummy_completion_markers #Must delocalize something or Cromwell will bypass permutation test
  String fdr_table_suffix
  String p_val_column_name

  command <<<
    set -e 

    # Gather all permutation results and compute FDR CDFs
    mkdir perm_res/
    while read prefix hpo; do
      echo -e "\nSTARTING $prefix\n"
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
      rm perm_res/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.perm_*.bed.gz
      echo -e "\nFINISHED $prefix\n"
    done < ${phenotype_list}
    echo -e "\nMAKING P-VALUE MATRIX\n"
    paste perm_res/*.sliding_window.meta_analysis.permuted_p_values.*.txt \
    | gzip -c \
    > ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz

    # Analyze p-values and compute FDR
    echo -e "\nANALYZING P-VALUE MATRIX\n"
    /opt/rCNV2/analysis/sliding_windows/calc_empirical_fdr.R \
      --cnv ${CNV} \
      --fdr-target ${fdr_target} \
      --linear-fit \
      --flat-ladder \
      --plot sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
      ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz \
      ${metacohort_sample_table} \
      sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}

    # Also produce an optional table of flat Bonferroni P-value cutoffs
    awk -v pval=${fdr_target} -v FS="\t" -v OFS="\t" \
      '{ print $1, pval }' ${phenotype_list} \
    > sliding_window.${freq_code}.${CNV}.bonferroni_pval.hpo_cutoffs.tsv

    # Copy cutoff tables to output bucket
    gsutil -m cp sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}*tsv \
      "${rCNV_bucket}/analysis/analysis_refs/"
    gsutil -m cp sliding_window.${freq_code}.${CNV}.bonferroni_pval.hpo_cutoffs.tsv \
      "${rCNV_bucket}/analysis/analysis_refs/"
    gsutil -m cp sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
      "${rCNV_bucket}/analysis/sliding_windows/empirical_fdr_plots/"
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "32 GB"
    disks: "local-disk 100 HDD"
    bootDiskSizeGb: "20"
  }

  output {
    File perm_results_plot = "sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png"
    File p_cutoff_table = "sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}.hpo_cutoffs.tsv"
    File p_cutoff_ladder = "sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}.ncase_cutoff_ladder.tsv"
    File bonferroni_cutoff_table = "sliding_window.${freq_code}.${CNV}.bonferroni_pval.hpo_cutoffs.tsv"
  }  
}


# Run meta-analysis across metacohorts for a single phenotype
task meta_analysis {
  Array[Array[File]] stats_beds
  Array[Array[File]] stats_bed_idxs
  String hpo
  File metacohort_list
  File metacohort_sample_table
  File exclusion_bed
  String freq_code
  Array[File] meta_p_cutoff_tables
  Int max_manhattan_phred_p
  String meta_model_prefix
  String rCNV_bucket
  String rCNV_docker
  String prefix
  String cache_string

  command <<<
    set -e

    # Copy burden stats & p-value cutoff tables
    find / -name "*${prefix}.${freq_code}.*.sliding_window.stats.bed.gz*" \
    | xargs -I {} mv {} ./
    find / -name "*sliding_window.${freq_code}.*.bonferroni_pval.hpo_cutoffs.tsv*" \
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
                    sliding_window.${freq_code}.DEL.bonferroni_pval.hpo_cutoffs.tsv )
    DUP_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                    sliding_window.${freq_code}.DUP.bonferroni_pval.hpo_cutoffs.tsv )

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
      /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
        --or-corplot ${prefix}.${freq_code}.$CNV.sliding_window.or_corplot_grid.jpg \
        --model ${meta_model_prefix} \
        --conditional-exclusion ${exclusion_bed} \
        --p-is-phred \
        --spa \
        ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt \
        ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed
      bgzip -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed
      tabix -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "meta_phred_p" \
        --p-is-phred \
        --max-phred-p ${max_manhattan_phred_p} \
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
      --max-phred-p ${max_manhattan_phred_p} \
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
    echo "${cache_string}" > completion.txt
  >>>

  output {
    File completion_token = "completion.txt"
  }

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }
}


# Refine associated regions to minimal credible regions
task refine_regions {
  Array[File] completion_tokens
  String freq_code
  String CNV
  File phenotype_list
  File metacohort_list
  File metacohort_sample_table
  String rCNV_bucket
  String rCNV_docker
  File meta_p_cutoffs_tsv
  Float meta_secondary_p_cutoff
  Int meta_nominal_cohorts_cutoff
  Float FDR_cutoff
  Int sig_window_pad
  Float credset
  File gtf

  command <<<
    set -e

    # Download all meta-analysis stats files and necessary data
    mkdir stats/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/sliding_windows/**.${freq_code}.**.sliding_window.meta_analysis.stats.bed.gz \
      stats/
    mkdir refs/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/analysis_refs/* \
      ${rCNV_bucket}/refs/GRCh37.cytobands.bed.gz \
      ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs*.bed.gz \
      refs/

    # Write tsv inputs
    while read prefix hpo; do
      for wrapper in 1; do
        echo "$hpo"
        echo "stats/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.subset.bed.gz"
        awk -v x=$prefix -v FS="\t" '{ if ($1==x) print $2 }' \
          ${meta_p_cutoffs_tsv}
      done | paste -s
    done < ${phenotype_list} \
    > ${freq_code}.${CNV}.segment_refinement.stats_input.tsv
    zcat refs/lit_GDs.*.bed.gz | fgrep -w ${CNV} | cut -f1-5 | \
    sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V | bedtools merge -i - \
    > all_GDs.${CNV}.bed
    echo "all_GDs.${CNV}.bed" > known_causal_loci_lists.${CNV}.tsv

    # Apply an initial loose mask per HPO to P<0.1 regions ±sig_window_pad
    # to reduce I/O time reading sumstats files in refinement
    # Also add all known GD regions for null variance estimation
    while read prefix hpo; do
      echo -e "Subsetting $prefix summary stats prior to refinement"
      allstats="stats/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.bed.gz"
      primary_p_idx=$( zcat $allstats | sed -n '1p' | sed 's/\t/\n/g' \
                       | awk -v FS="\t" '{ if ($1=="meta_phred_p") print NR }' )
      zcat $allstats \
      | grep -ve '^#' \
      | awk -v FS="\t" -v OFS="\t" -v idx=$primary_p_idx \
        '{ if ($(idx) >= 1 && $(idx) != "NA") print $1, $2, $3 }' \
      | cat - all_GDs.${CNV}.bed \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | bedtools merge -i - \
      | awk -v OFS="\t" -v dist=${sig_window_pad} \
        '{ if ($2-dist < 1) print $1, "1", $3+dist; else print $1, $2-dist, $3+dist }' \
      | bedtools intersect -wa -u -header \
        -a $allstats \
        -b - \
      | bgzip -c \
      > stats/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.subset.bed.gz
      tabix -f stats/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.subset.bed.gz
    done < ${phenotype_list}

    # Refine significant segments
    /opt/rCNV2/analysis/sliding_windows/refine_significant_regions.py \
      --cnv ${CNV} \
      --secondary-p-cutoff ${meta_secondary_p_cutoff} \
      --min-nominal ${meta_nominal_cohorts_cutoff} \
      --secondary-or-nominal \
      --fdr-q-cutoff ${FDR_cutoff} \
      --secondary-for-fdr \
      --credible-sets ${credset} \
      --distance ${sig_window_pad} \
      --known-causal-loci-list known_causal_loci_lists.${CNV}.tsv \
      --cytobands refs/GRCh37.cytobands.bed.gz \
      --sig-loci-bed ${freq_code}.${CNV}.final_segments.loci.pregenes.bed \
      --sig-assoc-bed ${freq_code}.${CNV}.final_segments.associations.pregenes.bed \
      ${freq_code}.${CNV}.segment_refinement.stats_input.tsv \
      ${metacohort_sample_table}

    # Annotate final regions with genes & sort by coordinates
    for entity in loci associations; do
      /opt/rCNV2/analysis/sliding_windows/get_genes_per_region.py \
        -o ${freq_code}.${CNV}.final_segments.$entity.bed \
        ${freq_code}.${CNV}.final_segments.$entity.pregenes.bed \
        ${gtf}
    done
  >>>

  output {
    File loci = "${freq_code}.${CNV}.final_segments.loci.bed"
    File associations = "${freq_code}.${CNV}.final_segments.associations.bed"
  }

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }
}


# Merge final refined regions results for public export
task merge_refined_regions {
  Array[File] assoc_beds
  Array[File] loci_beds
  String freq_code
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -e 

    # Get headers
    sed -n '1p' ${assoc_beds[0]} > assoc_header.tsv
    sed -n '1p' ${loci_beds[0]} > loci_header.tsv

    # Merge associations
    cat ${sep=" " assoc_beds} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n -k5,5V -k6,6V \
    | cat assoc_header.tsv - \
    | bgzip -c \
    > ${freq_code}.final_segments.associations.bed.gz

    # Merge segments
    cat ${sep=" " loci_beds} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n -k5,5V \
    | cat loci_header.tsv - \
    | bgzip -c \
    > ${freq_code}.final_segments.loci.bed.gz

    # Copy final files to results bucket (note: requires permissions)
    gsutil -m cp \
      ${freq_code}.final_segments.*.bed.gz \
      ${rCNV_bucket}/results/segment_association/
  >>>

  output {
    File final_associations = "${freq_code}.final_segments.associations.bed.gz"
    File final_loci = "${freq_code}.final_segments.loci.bed.gz"
  }

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    bootDiskSizeGb: "20"
  }
}


task plot_region_summary {
  String freq_code
  File DEL_regions
  File DUP_regions
  String rCNV_bucket
  String rCNV_docker

  command <<<
    /opt/rCNV2/analysis/sliding_windows/regions_summary.plot.R \
      -o "${freq_code}.final_segments." \
      ${DEL_regions} \
      ${DUP_regions}

    gsutil -m cp \
      "${freq_code}.final_segments.multipanel_summary.jpg" \
      ${rCNV_bucket}/public/
    gsutil acl ch -u AllUsers:R ${rCNV_bucket}/public/*.jpg
  >>>

  output {
    File summary_plot = "${freq_code}.final_segments.multipanel_summary.jpg"
  }

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
  }
}

