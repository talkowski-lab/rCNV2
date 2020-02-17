######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens per gene


workflow gene_burden_analysis {
  File phenotype_list
  File metacohort_list
  File metacohort_sample_table
  File gtf
  Int n_pheno_perms
  Int pad_controls
  String meta_model_prefix
  String weight_mode
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  Float p_cutoff
  String rCNV_bucket

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
        gtf=gtf,
        pad_controls=pad_controls,
        weight_mode=weight_mode,
        min_cds_ovr_del=min_cds_ovr_del,
        min_cds_ovr_dup=min_cds_ovr_dup,
        max_genes_per_cnv=max_genes_per_cnv,
        p_cutoff=p_cutoff,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0]
    }

    # Permute phenotypes to estimate empirical FDR
    call scattered_perm.scattered_gene_burden_perm_test as rCNV_perm_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        gtf=gtf,
        pad_controls=pad_controls,
        weight_mode=weight_mode,
        min_cds_ovr_del=min_cds_ovr_del,
        min_cds_ovr_dup=min_cds_ovr_dup,
        max_genes_per_cnv=max_genes_per_cnv,
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

    # # Genome-wide, secondary
    # call calc_meta_p_cutoff as calc_genome_wide_cutoffs_secondary {
    #   input:
    #     phenotype_list=phenotype_list,
    #     metacohort_sample_table=metacohort_sample_table,
    #     freq_code="rCNV",
    #     CNV=cnv,
    #     n_pheno_perms=n_pheno_perms,
    #     fdr_target=p_cutoff,
    #     rCNV_bucket=rCNV_bucket,
    #     dummy_completion_markers=rCNV_perm_test.completion_marker,
    #     fdr_table_suffix="empirical_genome_wide_pval_secondary",
    #     p_val_column_name="meta_phred_p_secondary"
    # }
  }

  # Perform meta-analysis of rCNV association statistics
  scatter ( pheno in phenotypes ) {
    call meta_analysis as rCNV_meta_analysis {
      input:
        count_beds=rCNV_burden_test.count_beds,
        count_bed_idxs=rCNV_burden_test.count_bed_idxs,
        stats_beds=rCNV_burden_test.stats_beds,
        stats_bed_idxs=rCNV_burden_test.stats_bed_idxs,
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        meta_p_cutoff=p_cutoff,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0]
    }
  }
}


# Run burden test for a single phenotype for all metacohorts
task burden_test {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  File gtf
  Int pad_controls
  String weight_mode
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  Float p_cutoff
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Copy CNV data and constrained gene coordinates
    mkdir cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/genes ./
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/refs/GRCh37.*.bed.gz refs/
    gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

    # Set HPO-specific parameters
    descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | awk -v FS="\t" '{ print $2 }' )
    zcat refs/gencode.v19.canonical.constrained.bed.gz \
    | fgrep -wf genes/gene_lists/${prefix}.HPOdb.constrained.genes.list \
    | cat <( echo -e "#chr\tstart\tend\tgene" ) - \
    > ${prefix}.highlight_regions.bed

    # Iterate over metacohorts
    while read meta cohorts; do
      echo $meta

      # Set metacohort-specific parameters
      cnv_bed="cleaned_cnv/$meta.${freq_code}.bed.gz"
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

        # Set CNV-specific parameters
        case "$CNV" in
          DEL)
            min_cds_ovr=${min_cds_ovr_del}
            ;;
          DUP)
            min_cds_ovr=${min_cds_ovr_dup}
            ;;
        esac

        # Count CNVs
        /opt/rCNV2/analysis/genes/count_cnvs_per_gene.py \
          --pad-controls ${pad_controls} \
          --weight-mode ${weight_mode} \
          --min-cds-ovr $min_cds_ovr \
          --max-genes ${max_genes_per_cnv} \
          -t $CNV \
          --hpo ${hpo} \
          --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
          --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
          --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
          -z \
          --verbose \
          -o "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz" \
          "$cnv_bed" \
          ${gtf}
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz"

        # Perform burden test
        /opt/rCNV2/analysis/genes/gene_burden_test.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --cnv $CNV \
          --case-hpo ${hpo} \
          --bgzip \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"

        # Generate Manhattan & QQ plots
        /opt/rCNV2/utils/plot_manhattan_qq.R \
          --p-col-name "fisher_phred_p" \
          --p-is-phred \
          --max-phred-p 100 \
          --cutoff ${p_cutoff} \
          --highlight-bed "${prefix}.highlight_regions.bed" \
          --highlight-name "Constrained genes associated with this phenotype" \
          --title "$title" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden"
      done

      # Generate Miami & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --miami \
        --p-col-name "fisher_phred_p" \
        --p-is-phred \
        --max-phred-p 100 \
        --cutoff ${p_cutoff} \
        --highlight-bed "${prefix}.highlight_regions.bed" \
        --highlight-name "Constrained genes associated with this phenotype" \
        --label-prefix "DUP" \
        --highlight-bed-2 "${prefix}.highlight_regions.bed" \
        --highlight-name-2 "Constrained genes associated with this phenotype" \
        --label-prefix-2 "DEL" \
        --title "$title" \
        "$meta.${prefix}.${freq_code}.DUP.gene_burden.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.DEL.gene_burden.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.gene_burden"
    done < ${metacohort_list}

    # Copy results to output bucket
    gsutil -m cp *.gene_burden.counts.bed.gz* \
      "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/counts/"
    gsutil -m cp *.gene_burden.stats.bed.gz* \
      "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.gene_burden.*.png \
      "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/plots/"
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:90840bcba5f6f5ef8ab183a4879d6eb25971845c79959386fca917cc0b155e11"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] stats_beds = glob("*.gene_burden.stats.bed.gz")
    Array[File] stats_bed_idxs = glob("*.gene_burden.stats.bed.gz.tbi")
    Array[File] count_beds = glob("*.gene_burden.counts.bed.gz")
    Array[File] count_bed_idxs = glob("*.gene_burden.counts.bed.gz.tbi")
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
    set -e 

    # Gather all permutation results and compute FDR CDFs
    mkdir perm_res/
    while read prefix hpo; do
      echo -e "\nSTARTING $prefix\n"
      gsutil -m cp \
        "${rCNV_bucket}/analysis/gene_burden/$prefix/${freq_code}/permutations/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.perm_*.bed.gz" \
        perm_res/
      for i in $( seq 1 ${n_pheno_perms} ); do
        p_idx=$( zcat perm_res/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.perm_$i.bed.gz \
                 | sed -n '1p' | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' \
                 | fgrep -w ${p_val_column_name} | cut -f2 )
        zcat perm_res/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.perm_$i.bed.gz \
        | grep -ve '^#' \
        | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
        | cat <( echo "$prefix.${CNV}.$i" ) - \
        > perm_res/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.permuted_p_values.$i.txt
      done
      rm perm_res/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.perm_*.bed.gz
      echo -e "\nFINISHED $prefix\n"
    done < ${phenotype_list}
    echo -e "\nMAKING P-VALUE MATRIX\n"
    paste perm_res/*.gene_burden.meta_analysis.permuted_p_values.*.txt \
    | gzip -c \
    > ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz

    # Analyze p-values and compute FDR
    echo -e "\nANALYZING P-VALUE MATRIX\n"
    /opt/rCNV2/analysis/sliding_windows/calc_empirical_fdr.R \
      --cnv ${CNV} \
      --fdr-target ${fdr_target} \
      --plot gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
      ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz \
      ${metacohort_sample_table} \
      gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}

    # Copy cutoff tables to output bucket
    gsutil -m cp gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}*tsv \
      "${rCNV_bucket}/analysis/analysis_refs/"
    gsutil -m cp gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
      "${rCNV_bucket}/analysis/gene_burden/empirical_fdr_plots/"
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:a96cc89624ac64474737d020ccba7484a406899e46268d5103a240f618d02b20"
    preemptible: 1
    memory: "32 GB"
    disks: "local-disk 275 HDD"
    bootDiskSizeGb: "40"
  }

  output {
    File perm_results_plot = "gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png"
    File p_cutoff_table = "gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}.hpo_cutoffs.tsv"
    File p_cutoff_ladder = "gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}.ncase_cutoff_ladder.tsv"
  }  
}

# Run meta-analysis (both weighted and raw) across metacohorts for a single phenotype
task meta_analysis {
  Array[File] count_beds
  Array[File] count_bed_idxs
  Array[File] stats_beds
  Array[File] stats_bed_idxs
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  Float meta_p_cutoff
  String meta_model_prefix
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Copy burden counts & gene coordinates
    find / -name "*${prefix}.${freq_code}.*.gene_burden.counts.bed.gz*" \
    | xargs -I {} mv {} ./
    find / -name "*${prefix}.${freq_code}.*.gene_burden.stats.bed.gz*" \
    | xargs -I {} mv {} ./
    gsutil -m cp -r gs://rcnv_project/cleaned_data/genes ./
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/refs/GRCh37.*.bed.gz refs/
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

    # Set HPO-specific parameters
    descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | awk -v FS="\t" '{ print $2 }' )
    zcat refs/gencode.v19.canonical.constrained.bed.gz \
    | fgrep -wf genes/gene_lists/${prefix}.HPOdb.constrained.genes.list \
    | cat <( echo -e "#chr\tstart\tend\tgene" ) - \
    > ${prefix}.highlight_regions.bed

    # Run meta-analysis for each CNV type
    for CNV in DEL DUP; do
      # Perform meta-analysis for unweighted CNVs
      while read meta cohorts; do
        echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"
      done < <( fgrep -v mega ${metacohort_list} ) \
      > ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt
      /opt/rCNV2/analysis/genes/gene_meta_analysis.R \
        --or-corplot ${prefix}.${freq_code}.$CNV.gene_burden.or_corplot_grid.jpg \
        --model ${meta_model_prefix} \
        --p-is-phred \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed
      bgzip -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed
      tabix -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "meta_phred_p" \
        --p-is-phred \
        --cutoff ${meta_p_cutoff} \
        --highlight-bed "${prefix}.highlight_regions.bed" \
        --highlight-name "Constrained genes associated with this phenotype" \
        --label-prefix "$CNV" \
        --title "$title" \
        "${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz" \
        "${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis"
    done

    # Generate Miami & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --miami \
      --p-col-name "meta_phred_p" \
      --p-is-phred \
      --cutoff ${meta_p_cutoff} \
      --highlight-bed "${prefix}.highlight_regions.bed" \
      --highlight-name "Constrained genes associated with this phenotype" \
      --label-prefix "DUP" \
      --highlight-bed-2 "${prefix}.highlight_regions.bed" \
      --highlight-name-2 "Constrained genes associated with this phenotype" \
      --label-prefix-2 "DEL" \
      --title "$title" \
      "${prefix}.${freq_code}.DUP.gene_burden.meta_analysis.stats.bed.gz" \
      "${prefix}.${freq_code}.DEL.gene_burden.meta_analysis.stats.bed.gz" \
      "${prefix}.${freq_code}.gene_burden.meta_analysis"

    # Copy results to output bucket
    gsutil -m cp *.gene_burden.*meta_analysis.stats.bed.gz* \
      "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.gene_burden.*or_corplot_grid.jpg \
      "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/plots/"
    gsutil -m cp *.gene_burden.*meta_analysis.*.png \
      "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/plots/"

    # Must delocalize completion marker to prevent caching of final step
    echo "DONE" > completion.txt
  >>>

  output {
    File completion_token = "completion.txt"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:0d25626d340a68933f598c0a5f3aa9078cf5f1c6b51ec42943e654c38d14e0c9"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }
}

