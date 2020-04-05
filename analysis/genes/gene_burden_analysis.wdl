######################
#    rCNV Project    #
######################

# Copyright (c) 2019-2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens per gene


import "https://api.firecloud.org/ga4gh/v1/tools/rCNV:scattered_gene_burden_perm_test/versions/8/plain-WDL/descriptor" as scattered_perm


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
  Int max_manhattan_phred_p
  Float meta_secondary_p_cutoff
  Int meta_nominal_cohorts_cutoff
  Float finemap_elnet_alpha
  Float finemap_elnet_l1_l2_mix
  Int finemap_distance
  File finemap_genomic_features
  File finemap_expression_features
  File finemap_constraint_features
  File finemap_merged_features
  File raw_finemap_merged_features
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
        gtf=gtf,
        pad_controls=pad_controls,
        weight_mode=weight_mode,
        min_cds_ovr_del=min_cds_ovr_del,
        min_cds_ovr_dup=min_cds_ovr_dup,
        max_genes_per_cnv=max_genes_per_cnv,
        p_cutoff=p_cutoff,
        max_manhattan_phred_p=max_manhattan_phred_p,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0],
        cache_string=fisher_cache_string
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
        freq_code="rCNV",
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.p_cutoff_table,
        max_manhattan_phred_p=max_manhattan_phred_p,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0],
        cache_string=meta_cache_string
    }
  }

  # Fine-map significant genes
  scatter ( cnv in cnv_types ) {
    # Genomic features
    call finemap_genes as finemap_genomic {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.p_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_output_label="genomic_features",
        gene_features=finemap_genomic_features,
        rCNV_bucket=rCNV_bucket
    }
    
    # Expression features
    call finemap_genes as finemap_expression {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.p_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_output_label="expression_features",
        gene_features=finemap_expression_features,
        rCNV_bucket=rCNV_bucket
    }
    
    # Constraint features
    call finemap_genes as finemap_constraint {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.p_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_output_label="constraint_features",
        gene_features=finemap_constraint_features,
        rCNV_bucket=rCNV_bucket
    }
    
    # Merged features
    call finemap_genes as finemap_merged {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.p_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_output_label="merged_features",
        gene_features=finemap_merged_features,
        rCNV_bucket=rCNV_bucket
    }
  }

  # Once complete, plot finemap results
  scatter( cnv in cnv_types ) {
    call plot_finemap_res {
      input:
        completion_tokens=[finemap_genomic.completion_token, finemap_expression.completion_token, finemap_constraint.completion_token, finemap_merged.completion_token],
        freq_code="rCNV",
        CNV=cnv,
        raw_features_merged=raw_finemap_merged_features,
        phenotype_list=phenotype_list,
        rCNV_bucket=rCNV_bucket
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
  File gtf
  Int pad_controls
  String weight_mode
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  Float p_cutoff
  Int max_manhattan_phred_p
  String rCNV_bucket
  String prefix
  String cache_string

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
          --max-phred-p ${max_manhattan_phred_p} \
          --cutoff ${p_cutoff} \
          --highlight-bed "${prefix}.highlight_regions.bed" \
          --highlight-name "Constrained genes associated with this phenotype" \
          --label-prefix "$CNV" \
          --title "$title" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden"
      done

      # Generate Miami & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --miami \
        --p-col-name "fisher_phred_p" \
        --p-is-phred \
        --max-phred-p ${max_manhattan_phred_p} \
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
    # gsutil -m cp *.gene_burden.counts.bed.gz* \
    #   "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/counts/"
    gsutil -m cp *.gene_burden.stats.bed.gz* \
      "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.gene_burden.*.png \
      "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/plots/"

    echo "${cache_string}" > completion.txt
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:fde16e782393b0e6319dc904e59137dcede0ec0d16e634ccb91a4b8de7a1565f"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] stats_beds = glob("*.gene_burden.stats.bed.gz")
    Array[File] stats_bed_idxs = glob("*.gene_burden.stats.bed.gz.tbi")
    File completion_token = "completion.txt"
    # Array[File] count_beds = glob("*.gene_burden.counts.bed.gz")
    # Array[File] count_bed_idxs = glob("*.gene_burden.counts.bed.gz.tbi")
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
    done < ${phenotype_list}
    paste perm_res/*.gene_burden.meta_analysis.permuted_p_values.*.txt \
    | gzip -c \
    > ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz

    # Analyze p-values and compute FDR
    /opt/rCNV2/analysis/sliding_windows/calc_empirical_fdr.R \
      --cnv ${CNV} \
      --fdr-target ${fdr_target} \
      --linear-fit \
      --flat-ladder \
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
    docker: "talkowski/rcnv@sha256:6519871cadff358c95c52f0b443ce4197089ec1d881178c0ada97f1744e4cc2b"
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
  Array[Array[File]] stats_beds
  Array[Array[File]] stats_bed_idxs
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  Array[File] meta_p_cutoff_tables
  Int max_manhattan_phred_p
  String meta_model_prefix
  String rCNV_bucket
  String prefix
  String cache_string

  command <<<
    set -e

    # Copy burden counts & gene coordinates
    find / -name "*${prefix}.${freq_code}.*.gene_burden.stats.bed.gz*" \
    | xargs -I {} mv {} ./
    find / -name "*gene_burden.${freq_code}.*.empirical_genome_wide_pval.hpo_cutoffs.tsv*" \
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
    DEL_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                    gene_burden.${freq_code}.DEL.empirical_genome_wide_pval.hpo_cutoffs.tsv )
    DUP_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                    gene_burden.${freq_code}.DUP.empirical_genome_wide_pval.hpo_cutoffs.tsv )

    # Set HPO-specific parameters
    descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | awk -v FS="\t" '{ print $2 }' )
    zcat refs/gencode.v19.canonical.constrained.bed.gz \
    | fgrep -wf genes/gene_lists/${prefix}.HPOdb.constrained.genes.list \
    | cat <( echo -e "#chr\tstart\tend\tgene" ) - \
    > ${prefix}.highlight_regions.bed


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

      # Perform meta-analysis for unweighted CNVs
      while read meta cohorts; do
        echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"
      done < <( fgrep -v mega ${metacohort_list} ) \
      > ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt
      /opt/rCNV2/analysis/genes/gene_meta_analysis.R \
        --or-corplot ${prefix}.${freq_code}.$CNV.gene_burden.or_corplot_grid.jpg \
        --model ${meta_model_prefix} \
        --p-is-phred \
        --spa \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed
      bgzip -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed
      tabix -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "meta_phred_p" \
        --p-is-phred \
        --max-phred-p ${max_manhattan_phred_p} \
        --cutoff $meta_p_cutoff \
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
      --max-phred-p ${max_manhattan_phred_p} \
      --cutoff $DUP_p_cutoff \
      --highlight-bed "${prefix}.highlight_regions.bed" \
      --highlight-name "Constrained genes associated with this phenotype" \
      --label-prefix "DUP" \
      --cutoff-2 $DEL_p_cutoff \
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
    echo "${cache_string}" > completion.txt
  >>>

  output {
    File completion_token = "completion.txt"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:6519871cadff358c95c52f0b443ce4197089ec1d881178c0ada97f1744e4cc2b"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }
}


# Fine-map genes
task finemap_genes {
  Array[File] completion_tokens # Must delocalize something from meta-analysis step to prevent caching
  File phenotype_list
  File metacohort_sample_table
  String freq_code
  String CNV
  Array[File] meta_p_cutoff_tables
  Float meta_secondary_p_cutoff
  Int meta_nominal_cohorts_cutoff
  Float finemap_elnet_alpha
  Float finemap_elnet_l1_l2_mix
  Int finemap_distance
  String finemap_output_label
  File gene_features
  String rCNV_bucket

  command <<<
    set -e

    # Copy p-value cutoff tables
    find / -name "*gene_burden.${freq_code}.*.empirical_genome_wide_pval.hpo_cutoffs.tsv*" \
    | xargs -I {} mv {} ./

    # Copy association stats & gene lists from the project Google Bucket (note: requires permissions)
    mkdir stats
    gsutil -m cp \
      ${rCNV_bucket}/analysis/gene_burden/**.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz \
      stats/
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/genes/gene_lists \
      ./

    # Write tsv input
    while read prefix hpo; do
      for wrapper in 1; do
        echo "$hpo"
        echo "stats/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz"
        awk -v x=$prefix -v FS="\t" '{ if ($1==x) print $2 }' \
          gene_burden.${freq_code}.${CNV}.empirical_genome_wide_pval.hpo_cutoffs.tsv
      done | paste -s
    done < ${phenotype_list} \
    > ${freq_code}.${CNV}.gene_fine_mapping.stats_input.tsv

    # Make lists of a priori "true causal" genes for null variance estimation
    case ${CNV} in
      "DEL")
        for wrapper in 1; do
          echo "gene_lists/DDG2P.hmc_lof.genes.list"
          echo "gene_lists/DDG2P.hmc_other.genes.list"
          echo "gene_lists/ClinGen.hmc_haploinsufficient.genes.list"
          echo "gene_lists/HP0000118.HPOdb.constrained.genes.list"
          echo "gene_lists/gnomad.v2.1.1.lof_constrained.genes.list"
        done > known_causal_gene_lists.tsv
        ;;
      "DUP")
        for wrapper in 1; do
          echo "gene_lists/DDG2P.all_gof.genes.list"
          echo "gene_lists/DDG2P.hmc_other.genes.list"
          echo "gene_lists/ClinGen.all_triplosensitive.genes.list"
          echo "gene_lists/HP0000118.HPOdb.constrained.genes.list"
          echo "gene_lists/gnomad.v2.1.1.lof_constrained.genes.list"
        done > known_causal_gene_lists.tsv
        ;;
    esac

    # Run functional fine-mapping procedure
    /opt/rCNV2/analysis/genes/finemap_genes.py \
      --secondary-p-cutoff ${meta_secondary_p_cutoff} \
      --min-nominal ${meta_nominal_cohorts_cutoff} \
      --secondary-or-nominal \
      --regularization-alpha ${finemap_elnet_alpha} \
      --regularization-l1-l2-mix ${finemap_elnet_l1_l2_mix} \
      --distance ${finemap_distance} \
      --known-causal-gene-lists known_causal_gene_lists.tsv \
      --outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.tsv \
      --all-genes-outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.all_genes_from_blocks.tsv \
      --naive-outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.${finemap_output_label}.tsv \
      --genetic-outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.${finemap_output_label}.tsv \
      --coeffs-out ${freq_code}.${CNV}.gene_fine_mapping.logit_coeffs.${finemap_output_label}.tsv \
      ${freq_code}.${CNV}.gene_fine_mapping.stats_input.tsv \
      ${gene_features}

    # Copy results to output bucket
    gsutil -m cp \
      ${freq_code}.${CNV}.gene_fine_mapping.*.${finemap_output_label}*tsv \
      ${rCNV_bucket}/analysis/gene_burden/fine_mapping/

    # Make completion token to track caching
    echo "Done" > completion.txt
  >>>

  output {
    File finemapped_output = "${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.tsv"
    File finemapped_output_all_genes = "${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.all_genes_from_blocks.tsv"
    File naive_output = "${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.${finemap_output_label}.tsv"
    File genetic_output = "${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.${finemap_output_label}.tsv"
    File logit_coeffs = "${freq_code}.${CNV}.gene_fine_mapping.logit_coeffs.${finemap_output_label}.tsv"
    File completion_token = "completion.txt"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:6519871cadff358c95c52f0b443ce4197089ec1d881178c0ada97f1744e4cc2b"
    preemptible: 1
    memory: "8 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }
}


# Plot fine-mapping results
task plot_finemap_res {
  Array[Array[File]] completion_tokens
  String freq_code
  String CNV
  File phenotype_list
  File raw_features_genomic
  File raw_features_expression
  File raw_features_constraint
  File raw_features_merged
  String rCNV_bucket

  command <<<
    set -e

    # Copy all fine-mapped gene lists
    mkdir finemap_stats/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/gene_burden/fine_mapping/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.*.tsv \
      finemap_stats/

    # Make input tsvs
    for wrapper in 1; do
      echo -e "Prior\tgrey70\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.genomic_features.tsv"
      echo -e "Posterior\t'#264653'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.genomic_features.tsv"
      echo -e "Genomic features\t'#E76F51'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genomic_features.tsv"
      echo -e "Gene expression\t'#E9C46A'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.expression_features.tsv"
      echo -e "Gene constraint\t'#F4A261'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.constraint_features.tsv"
      echo -e "Full model\t'#2A9D8F'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.merged_features.tsv"
    done > finemap_roc_input.tsv
    for wrapper in 1; do
      echo -e "Genomic features\t'#E76F51'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genomic_features.all_genes_from_blocks.tsv\t${raw_features_genomic}"
      echo -e "Gene expression\t'#E9C46A'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.expression_features.all_genes_from_blocks.tsv\t${raw_features_expression}"
      echo -e "Gene constraint\t'#F4A261'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.constraint_features.all_genes_from_blocks.tsv\t${raw_features_constraint}"
      echo -e "Full model\t'#2A9D8F'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.merged_features.all_genes_from_blocks.tsv\t${raw_features_merged}"
    done > finemap_feature_cor_input.tsv

    # Make all gene truth sets
    gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/gene_lists ./

    # HPO-associated
    while read prefix hpo; do
      awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }' \
        gene_lists/$prefix.HPOdb.genes.list
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > hpo_truth_set.tsv

    # Make CNV type-dependent truth sets
    case ${CNV} in
      "DEL")
        # Union (ClinGen HI + DDG2P dominant LoF)
        while read prefix hpo; do
          cat gene_lists/ClinGen.hmc_haploinsufficient.genes.list \
              gene_lists/DDG2P.hmc_lof.genes.list \
          | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
        done < ${phenotype_list} \
        | sort -Vk1,1 -k2,2V | uniq \
        | cat <( echo -e "#HPO\tgene" ) - \
        > union_truth_set.lof.tsv

        # Intersection (ClinGen HI + DDG2P dominant LoF)
        while read prefix hpo; do
          fgrep -wf \
            gene_lists/ClinGen.hmc_haploinsufficient.genes.list \
            gene_lists/DDG2P.hmc_lof.genes.list \
          | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
        done < ${phenotype_list} \
        | sort -Vk1,1 -k2,2V | uniq \
        | cat <( echo -e "#HPO\tgene" ) - \
        > intersection_truth_set.lof.tsv

        # ClinGen HI alone
        while read prefix hpo; do
          awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }' \
            gene_lists/ClinGen.all_haploinsufficient.genes.list
        done < ${phenotype_list} \
        | sort -Vk1,1 -k2,2V | uniq \
        | cat <( echo -e "#HPO\tgene" ) - \
        > clingen_truth_set.lof.tsv

        # DDG2P lof + other alone
        while read prefix hpo; do
          cat gene_lists/DDG2P.all_lof.genes.list \
              gene_lists/DDG2P.all_other.genes.list \
          | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
        done < ${phenotype_list} \
        | sort -Vk1,1 -k2,2V | uniq \
        | cat <( echo -e "#HPO\tgene" ) - \
        > ddg2p_truth_set.lof.tsv

        # Combine all truth sets into master union
        cat union_truth_set.lof.tsv \
            clingen_truth_set.lof.tsv \
            ddg2p_truth_set.lof.tsv \
            hpo_truth_set.tsv \
        | fgrep -v "#" | sort -Vk1,1 -k2,2V | uniq \
        | cat <( echo -e "#HPO\tgene" ) - \
        > master_union_truth_set.lof.tsv

        # Write truth set input tsv
        for wrapper in 1; do
          echo -e "Union of all truth sets\tmaster_union_truth_set.lof.tsv"
          echo -e "ClinGen HI & DECIPHER LoF (union)\tunion_truth_set.lof.tsv"
          echo -e "ClinGen HI & DECIPHER LoF (intersection)\tintersection_truth_set.lof.tsv"
          echo -e "ClinGen HI (any confidence)\tclingen_truth_set.lof.tsv"
          echo -e "DECIPHER dominant LoF/unk. (any confidence)\tddg2p_truth_set.lof.tsv"
          echo -e "HPO-matched disease genes\thpo_truth_set.tsv"
        done > finemap_roc_truth_sets.tsv
        ;;

      "DUP")
        # Union (ClinGen HI + DDG2P dominant CG)
        while read prefix hpo; do
          cat gene_lists/ClinGen.hmc_triplosensitive.genes.list \
              gene_lists/DDG2P.hmc_gof.genes.list \
          | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
        done < ${phenotype_list} \
        | sort -Vk1,1 -k2,2V | uniq \
        | cat <( echo -e "#HPO\tgene" ) - \
        > union_truth_set.gof.tsv

        # ClinGen triplo alone
        while read prefix hpo; do
          awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }' \
            gene_lists/ClinGen.all_triplosensitive.genes.list
        done < ${phenotype_list} \
        | sort -Vk1,1 -k2,2V | uniq \
        | cat <( echo -e "#HPO\tgene" ) - \
        > clingen_truth_set.triplo.tsv

        # DDG2P gof + other alone
        while read prefix hpo; do
          cat gene_lists/DDG2P.all_gof.genes.list \
              gene_lists/DDG2P.all_other.genes.list \
          | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
        done < ${phenotype_list} \
        | sort -Vk1,1 -k2,2V | uniq \
        | cat <( echo -e "#HPO\tgene" ) - \
        > ddg2p_truth_set.gof.tsv

        # Combine all truth sets into master union
        cat union_truth_set.gof.tsv \
            clingen_truth_set.gof.tsv \
            ddg2p_truth_set.gof.tsv \
            hpo_truth_set.tsv \
        | fgrep -v "#" | sort -Vk1,1 -k2,2V | uniq \
        | cat <( echo -e "#HPO\tgene" ) - \
        > master_union_truth_set.gof.tsv

        # Write truth set input tsv
        for wrapper in 1; do
          echo -e "Union of all truth sets\tmaster_union_truth_set.gof.tsv"
          echo -e "ClinGen TS & DECIPHER GoF (union)\tunion_truth_set.gof.tsv"
          echo -e "ClinGen TS (any confidence)\tclingen_truth_set.triplo.tsv"
          echo -e "DECIPHER dominant GoF/unk. (any confidence)\tddg2p_truth_set.gof.tsv"
          echo -e "HPO-matched disease genes\thpo_truth_set.tsv"
        done > finemap_roc_truth_sets.tsv
        ;;

    esac

    # Plot finemapping QC & feature correlations
    mkdir ${freq_code}_${CNV}_finemap_plots/
    /opt/rCNV2/analysis/genes/plot_finemap_results.R \
      finemap_roc_input.tsv \
      finemap_roc_truth_sets.tsv \
      ${freq_code}_${CNV}_finemap_plots/${freq_code}.${CNV}.finemap_results
    /opt/rCNV2/analysis/genes/plot_finemap_coefficients.R \
      finemap_feature_cor_input.tsv \
      ${freq_code}_${CNV}_finemap_plots/${freq_code}.${CNV}.finemap_feature_cors

    # Compress results
    tar -czvf ${freq_code}_${CNV}_finemap_plots.tgz ${freq_code}_${CNV}_finemap_plots
  >>>

  output {
    File plots_tarball = "${freq_code}_${CNV}_finemap_plots.tgz"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:6519871cadff358c95c52f0b443ce4197089ec1d881178c0ada97f1744e4cc2b"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }
}
