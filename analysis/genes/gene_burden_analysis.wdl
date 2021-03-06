######################
#    rCNV Project    #
######################

# Copyright (c) 2019-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens per gene


import "https://api.firecloud.org/ga4gh/v1/tools/rCNV:scattered_gene_burden_perm_test/versions/9/plain-WDL/descriptor" as scattered_perm


workflow gene_burden_analysis {
  File phenotype_list
  File metacohort_list
  File metacohort_sample_table
  File gtf
  Int n_pheno_perms
  Int pad_controls
  Int min_probes_per_gene
  Float min_frac_controls_probe_exclusion
  String meta_model_prefix
  String weight_mode
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  Float p_cutoff
  Int max_manhattan_phred_p
  Float meta_secondary_p_cutoff
  Int meta_nominal_cohorts_cutoff
  Float FDR_cutoff
  Float finemap_elnet_alpha
  Float finemap_elnet_l1_l2_mix
  Int finemap_distance
  Float finemap_conf_pip
  Float finemap_vconf_pip
  File finemap_genomic_features
  File finemap_expression_features
  File finemap_chromatin_features
  File finemap_constraint_features
  File finemap_variation_features
  File finemap_merged_no_variation_features
  File finemap_merged_features
  File raw_finemap_genomic_features
  File raw_finemap_expression_features
  File raw_finemap_chromatin_features
  File raw_finemap_constraint_features
  File raw_finemap_variation_features
  File raw_finemap_merged_no_variation_features
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
  }

  # Determine conditional cohort exclusion list based on probe density
  call build_exclusion_list {
    input:
      genes_bed=rCNV_burden_test.stats_beds[0],
      gtf_prefix=basename(gtf, '.gtf.gz'),
      min_probes_per_gene=min_probes_per_gene,
      min_frac_controls_probe_exclusion=min_frac_controls_probe_exclusion,
      metacohort_list=metacohort_list,
      rCNV_bucket=rCNV_bucket,
      freq_code=freq_code
  }

  # Scatter over phenotypes
  scatter ( pheno in phenotypes ) {
    # Permute phenotypes to estimate empirical FDR
    call scattered_perm.scattered_gene_burden_perm_test as rCNV_perm_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        exclusion_bed=build_exclusion_list.exclusion_bed,
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
        prefix=pheno[0],
        cache_string=meta_cache_string
    }
  }

  # Fine-map significant genes
  scatter ( cnv in cnv_types ) {
    # Genomic features
    call finemap_genes as finemap_genomic_gw {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="genomic_features",
        gene_features=finemap_genomic_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="FALSE",
        output_suffix="strict_gw_sig",
        rCNV_bucket=rCNV_bucket
    }
    call finemap_genes as finemap_genomic_fdr {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=1,
        meta_nominal_cohorts_cutoff=0,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="genomic_features",
        gene_features=finemap_genomic_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="TRUE",
        output_suffix="fdr",
        rCNV_bucket=rCNV_bucket
    }
    
    # Expression features
    call finemap_genes as finemap_expression_gw {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="expression_features",
        gene_features=finemap_expression_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="FALSE",
        output_suffix="strict_gw_sig",
        rCNV_bucket=rCNV_bucket
    }
    call finemap_genes as finemap_expression_fdr {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=1,
        meta_nominal_cohorts_cutoff=0,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="expression_features",
        gene_features=finemap_expression_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="TRUE",
        output_suffix="fdr",
        rCNV_bucket=rCNV_bucket
    }
    
    # Chromatin features
    call finemap_genes as finemap_chromatin_gw {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="chromatin_features",
        gene_features=finemap_chromatin_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="FALSE",
        output_suffix="strict_gw_sig",
        rCNV_bucket=rCNV_bucket
    }
    call finemap_genes as finemap_chromatin_fdr {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=1,
        meta_nominal_cohorts_cutoff=0,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="chromatin_features",
        gene_features=finemap_chromatin_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="TRUE",
        output_suffix="fdr",
        rCNV_bucket=rCNV_bucket
    }

    # Constraint features
    call finemap_genes as finemap_constraint_gw {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="constraint_features",
        gene_features=finemap_constraint_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="FALSE",
        output_suffix="strict_gw_sig",
        rCNV_bucket=rCNV_bucket
    }
    call finemap_genes as finemap_constraint_fdr {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=1,
        meta_nominal_cohorts_cutoff=0,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="constraint_features",
        gene_features=finemap_constraint_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="TRUE",
        output_suffix="fdr",
        rCNV_bucket=rCNV_bucket
    }

    # Variation features
    call finemap_genes as finemap_variation_gw {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="variation_features",
        gene_features=finemap_variation_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="FALSE",
        output_suffix="strict_gw_sig",
        rCNV_bucket=rCNV_bucket
    }
    call finemap_genes as finemap_variation_fdr {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=1,
        meta_nominal_cohorts_cutoff=0,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="variation_features",
        gene_features=finemap_variation_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="TRUE",
        output_suffix="fdr",
        rCNV_bucket=rCNV_bucket
    }
    
    # Merged features
    call finemap_genes as finemap_merged_gw {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="merged_features",
        gene_features=finemap_merged_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="FALSE",
        output_suffix="strict_gw_sig",
        rCNV_bucket=rCNV_bucket
    }
    call finemap_genes as finemap_merged_fdr {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=1,
        meta_nominal_cohorts_cutoff=0,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="merged_features",
        gene_features=finemap_merged_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="TRUE",
        output_suffix="fdr",
        rCNV_bucket=rCNV_bucket
    }
    
    # Merged features (no variation)
    call finemap_genes as finemap_merged_no_variation_gw {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="merged_no_variation_features",
        gene_features=finemap_merged_no_variation_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="FALSE",
        output_suffix="strict_gw_sig",
        rCNV_bucket=rCNV_bucket
    }
    call finemap_genes as finemap_merged_no_variation_fdr {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=1,
        meta_nominal_cohorts_cutoff=0,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_distance=finemap_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="merged_no_variation_features",
        gene_features=finemap_merged_no_variation_features,
        FDR_cutoff=FDR_cutoff,
        use_FDR="TRUE",
        output_suffix="fdr",
        rCNV_bucket=rCNV_bucket
    }
  }

  # Combine merged features (without variation) BED files as final output
  call merge_finemap_res as merge_finemap_res_gw {
    input:
      gene_beds=finemap_merged_no_variation_gw.sig_genes_bed,
      assoc_beds=finemap_merged_no_variation_gw.sig_assocs_bed,
      credset_beds=finemap_merged_no_variation_gw.cred_sets_bed,
      freq_code="rCNV",
      output_suffix="strict_gw_sig",
      rCNV_bucket=rCNV_bucket
  }
  call merge_finemap_res as merge_finemap_res_fdr {
    input:
      gene_beds=finemap_merged_no_variation_fdr.sig_genes_bed,
      assoc_beds=finemap_merged_no_variation_fdr.sig_assocs_bed,
      credset_beds=finemap_merged_no_variation_fdr.cred_sets_bed,
      freq_code="rCNV",
      output_suffix="fdr",
      rCNV_bucket=rCNV_bucket
  }

  # Once complete, plot finemap results
  scatter( cnv in cnv_types ) {
    call plot_finemap_res as plot_finemap_res_gw {
      input:
        completion_tokens=[finemap_genomic_gw.completion_token, finemap_expression_gw.completion_token, finemap_chromatin_gw.completion_token, finemap_constraint_gw.completion_token, finemap_variation_gw.completion_token, finemap_merged_no_variation_gw.completion_token, finemap_merged_gw.completion_token],
        freq_code="rCNV",
        CNV=cnv,
        raw_features_genomic=raw_finemap_genomic_features,
        raw_features_expression=raw_finemap_expression_features,
        raw_features_chromatin=raw_finemap_chromatin_features,
        raw_features_constraint=raw_finemap_constraint_features,
        raw_features_variation=raw_finemap_variation_features,
        raw_features_merged_no_variation=raw_finemap_merged_no_variation_features,
        raw_features_merged=raw_finemap_merged_features,
        phenotype_list=phenotype_list,
        output_suffix="strict_gw_sig",
        rCNV_bucket=rCNV_bucket
    }
    call plot_finemap_res as plot_finemap_res_fdr {
      input:
        completion_tokens=[finemap_genomic_fdr.completion_token, finemap_expression_fdr.completion_token, finemap_chromatin_fdr.completion_token, finemap_constraint_fdr.completion_token, finemap_variation_fdr.completion_token, finemap_merged_no_variation_fdr.completion_token, finemap_merged_fdr.completion_token],
        freq_code="rCNV",
        CNV=cnv,
        raw_features_genomic=raw_finemap_genomic_features,
        raw_features_expression=raw_finemap_expression_features,
        raw_features_chromatin=raw_finemap_chromatin_features,
        raw_features_constraint=raw_finemap_constraint_features,
        raw_features_variation=raw_finemap_variation_features,
        raw_features_merged_no_variation=raw_finemap_merged_no_variation_features,
        raw_features_merged=raw_finemap_merged_features,
        phenotype_list=phenotype_list,
        output_suffix="fdr",
        rCNV_bucket=rCNV_bucket
    }
  }

  output {
    File final_sig_genes_gw = merge_finemap_res_gw.final_genes
    File final_sig_associations_gw = merge_finemap_res_gw.final_associations
    File final_credible_sets_gw = merge_finemap_res_gw.final_credsets
    File final_sig_genes_fdr = merge_finemap_res_fdr.final_genes
    File final_sig_associations_fdr = merge_finemap_res_fdr.final_associations
    File final_credible_sets_fdr = merge_finemap_res_fdr.final_credsets
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
        /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --case-hpo ${hpo} \
          --keep-n-columns 4 \
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
    docker: "talkowski/rcnv@sha256:93ec0fee2b0ad415143eda627c2b3c8d2e1ef3c8ff4d3d620767637614fee5f8"
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


# Build list of cohorts to exclude per gene based on inadequate probe density
task build_exclusion_list {
  File genes_bed
  String gtf_prefix
  Int min_probes_per_gene
  Float min_frac_controls_probe_exclusion
  File metacohort_list
  String rCNV_bucket
  String freq_code

  command <<<
    set -euo pipefail

    # Download control probesets
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/control_probesets \
      ./

    # Make inputs for conditional exclusion script
    zcat ${genes_bed} | cut -f1-4 | bgzip -c > gene_coords.bed.gz
    for file in control_probesets/*bed.gz; do
      echo -e "$file\t$( basename $file | sed 's/\.bed\.gz//g' )"
    done > probeset_tracks.tsv

    # Build conditional exclusion list
    /opt/rCNV2/data_curation/other/probe_based_exclusion.py \
      --outfile ${gtf_prefix}.cohort_exclusion.bed.gz \
      --probecounts-outfile ${gtf_prefix}.probe_counts.bed.gz \
      --frac-pass-outfile ${gtf_prefix}.frac_passing.bed.gz \
      --min-probes ${min_probes_per_gene} \
      --min-frac-samples ${min_frac_controls_probe_exclusion} \
      --keep-n-columns 4 \
      --bgzip \
      gene_coords.bed.gz \
      probeset_tracks.tsv \
      control_probesets/rCNV.control_counts_by_array.tsv \
      <( fgrep -v mega ${metacohort_list} )
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:db7a75beada57d8e2649ce132581f675eb47207de489c3f6ac7f3452c51ddb6e"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    File exclusion_bed = "${gtf_prefix}.cohort_exclusion.bed.gz"
    File probe_counts = "${gtf_prefix}.probe_counts.bed.gz"
    File frac_passing = "${gtf_prefix}.frac_passing.bed.gz"
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

    # Also produce an optional table of flat Bonferroni P-value cutoffs
    awk -v pval=${fdr_target} -v FS="\t" -v OFS="\t" \
      '{ print $1, pval }' ${phenotype_list} \
    > gene_burden.${freq_code}.${CNV}.bonferroni_pval.hpo_cutoffs.tsv

    # Copy cutoff tables to output bucket
    gsutil -m cp gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}*tsv \
      "${rCNV_bucket}/analysis/analysis_refs/"
    gsutil -m cp gene_burden.${freq_code}.${CNV}.bonferroni_pval.hpo_cutoffs.tsv \
      "${rCNV_bucket}/analysis/analysis_refs/"
    gsutil -m cp gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
      "${rCNV_bucket}/analysis/gene_burden/empirical_fdr_plots/"
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:2942c7386b43479d02b29506ad1f28fdcff17bdf8b279f2e233be0c4d2cd50fa"
    preemptible: 1
    memory: "32 GB"
    disks: "local-disk 275 HDD"
    bootDiskSizeGb: "40"
  }

  output {
    File perm_results_plot = "gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png"
    File p_cutoff_table = "gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}.hpo_cutoffs.tsv"
    File p_cutoff_ladder = "gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}.ncase_cutoff_ladder.tsv"
    File bonferroni_cutoff_table = "gene_burden.${freq_code}.${CNV}.bonferroni_pval.hpo_cutoffs.tsv "
  }  
}

# Run meta-analysis (both weighted and raw) across metacohorts for a single phenotype
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
  String prefix
  String cache_string

  command <<<
    set -e

    # Copy burden counts & gene coordinates
    find / -name "*${prefix}.${freq_code}.*.gene_burden.stats.bed.gz*" \
    | xargs -I {} mv {} ./
    find / -name "*gene_burden.${freq_code}.*.bonferroni_pval.hpo_cutoffs.tsv*" \
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
                    gene_burden.${freq_code}.DEL.bonferroni_pval.hpo_cutoffs.tsv )
    DUP_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                    gene_burden.${freq_code}.DUP.bonferroni_pval.hpo_cutoffs.tsv )

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
      /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
        --or-corplot ${prefix}.${freq_code}.$CNV.gene_burden.or_corplot_grid.jpg \
        --model ${meta_model_prefix} \
        --conditional-exclusion ${exclusion_bed} \
        --keep-n-columns 4 \
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
    docker: "talkowski/rcnv@sha256:2942c7386b43479d02b29506ad1f28fdcff17bdf8b279f2e233be0c4d2cd50fa"
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
  Float FDR_cutoff
  String use_FDR
  Float finemap_elnet_alpha
  Float finemap_elnet_l1_l2_mix
  Int finemap_distance
  Float finemap_conf_pip
  Float finemap_vconf_pip
  String finemap_output_label
  File gene_features
  String output_suffix
  String rCNV_bucket

  command <<<
    set -e

    # Copy p-value cutoff tables
    find / -name "*gene_burden.${freq_code}.*.bonferroni_pval.hpo_cutoffs.tsv*" \
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
        if [ ${use_FDR} == "TRUE" ]; then
          echo "${FDR_cutoff}"
        else
          awk -v x=$prefix -v FS="\t" '{ if ($1==x) print $2 }' \
            gene_burden.${freq_code}.${CNV}.bonferroni_pval.hpo_cutoffs.tsv
        fi
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
    if [ ${use_FDR} == "TRUE" ]; then
      /opt/rCNV2/analysis/genes/finemap_genes.py \
        --cnv ${CNV} \
        --use-fdr \
        --secondary-p-cutoff 1 \
        --min-nominal 0 \
        --secondary-or-nominal \
        --regularization-alpha ${finemap_elnet_alpha} \
        --regularization-l1-l2-mix ${finemap_elnet_l1_l2_mix} \
        --distance ${finemap_distance} \
        --confident-pip ${finemap_conf_pip} \
        --very-confident-pip ${finemap_vconf_pip} \
        --known-causal-gene-lists known_causal_gene_lists.tsv \
        --outfile ${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.${finemap_output_label}.tsv \
        --sig-genes-bed ${freq_code}.${CNV}.final_genes.${output_suffix}.genes.bed \
        --sig-assoc-bed ${freq_code}.${CNV}.final_genes.${output_suffix}.associations.bed \
        --sig-credsets-bed ${freq_code}.${CNV}.final_genes.${output_suffix}.credible_sets.bed \
        --all-genes-outfile ${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.${finemap_output_label}.all_genes_from_blocks.tsv \
        --naive-outfile ${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.naive_priors.${finemap_output_label}.tsv \
        --genetic-outfile ${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.genetics_only.${finemap_output_label}.tsv \
        --coeffs-out ${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.logit_coeffs.${finemap_output_label}.tsv \
        ${freq_code}.${CNV}.gene_fine_mapping.stats_input.tsv \
        ${gene_features} \
        ${metacohort_sample_table}
    else
      /opt/rCNV2/analysis/genes/finemap_genes.py \
        --cnv ${CNV} \
        --secondary-p-cutoff ${meta_secondary_p_cutoff} \
        --min-nominal ${meta_nominal_cohorts_cutoff} \
        --secondary-or-nominal \
        --regularization-alpha ${finemap_elnet_alpha} \
        --regularization-l1-l2-mix ${finemap_elnet_l1_l2_mix} \
        --distance ${finemap_distance} \
        --confident-pip ${finemap_conf_pip} \
        --very-confident-pip ${finemap_vconf_pip} \
        --known-causal-gene-lists known_causal_gene_lists.tsv \
        --outfile ${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.${finemap_output_label}.tsv \
        --sig-genes-bed ${freq_code}.${CNV}.final_genes.${output_suffix}.genes.bed \
        --sig-assoc-bed ${freq_code}.${CNV}.final_genes.${output_suffix}.associations.bed \
        --sig-credsets-bed ${freq_code}.${CNV}.final_genes.${output_suffix}.credible_sets.bed \
        --all-genes-outfile ${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.${finemap_output_label}.all_genes_from_blocks.tsv \
        --naive-outfile ${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.naive_priors.${finemap_output_label}.tsv \
        --genetic-outfile ${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.genetics_only.${finemap_output_label}.tsv \
        --coeffs-out ${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.logit_coeffs.${finemap_output_label}.tsv \
        ${freq_code}.${CNV}.gene_fine_mapping.stats_input.tsv \
        ${gene_features} \
        ${metacohort_sample_table}
    fi

    # Copy results to output bucket
    gsutil -m cp \
      ${freq_code}.${CNV}.gene_fine_mapping.*.${finemap_output_label}*tsv \
      ${rCNV_bucket}/analysis/gene_burden/fine_mapping/

    # Make completion token to track caching
    echo "Done" > completion.txt
  >>>

  output {
    File finemapped_output = "${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.${finemap_output_label}.tsv"
    File finemapped_output_all_genes = "${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.${finemap_output_label}.all_genes_from_blocks.tsv"
    File naive_output = "${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.naive_priors.${finemap_output_label}.tsv"
    File genetic_output = "${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.genetics_only.${finemap_output_label}.tsv"
    File logit_coeffs = "${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.logit_coeffs.${finemap_output_label}.tsv"
    File sig_genes_bed = "${freq_code}.${CNV}.final_genes.${output_suffix}.genes.bed"
    File sig_assocs_bed = "${freq_code}.${CNV}.final_genes.${output_suffix}.associations.bed"
    File cred_sets_bed = "${freq_code}.${CNV}.final_genes.${output_suffix}.credible_sets.bed"
    File completion_token = "completion.txt"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:5289b227ad6a283da7ffff23e305bb19db8454164cbca9c9a7c5c8aa78e0463e"
    preemptible: 1
    memory: "8 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }
}


# Merge final gene association results for public export
task merge_finemap_res {
  Array[File] gene_beds
  Array[File] assoc_beds
  Array[File] credset_beds
  String freq_code
  String output_suffix
  String rCNV_bucket

  command <<<
    set -e 

    # Get headers
    sed -n '1p' ${gene_beds[0]} > gene_header.tsv
    sed -n '1p' ${assoc_beds[0]} > assoc_header.tsv
    sed -n '1p' ${credset_beds[0]} > credset_header.tsv

    # Merge genes
    cat ${sep=" " gene_beds} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k13,13V \
    | cat gene_header.tsv - \
    | bgzip -c \
    > ${freq_code}.final_genes.${output_suffix}.genes.bed.gz

    # Merge associations
    cat ${sep=" " assoc_beds} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k6,6V \
    | cat assoc_header.tsv - \
    | bgzip -c \
    > ${freq_code}.final_genes.${output_suffix}.associations.bed.gz

    # Merge genes
    cat ${sep=" " credset_beds} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n -k5,5V -k6,6V \
    | cat credset_header.tsv - \
    | bgzip -c \
    > ${freq_code}.final_genes.${output_suffix}.credible_sets.bed.gz

    # Copy final files to results bucket (note: requires permissions)
    gsutil -m cp \
      ${freq_code}.final_genes.${output_suffix}.*.bed.gz \
      ${rCNV_bucket}/results/gene_association/
  >>>

  output {
    File final_genes = "${freq_code}.final_genes.${output_suffix}.genes.bed.gz"
    File final_associations = "${freq_code}.final_genes.${output_suffix}.associations.bed.gz"
    File final_credsets = "${freq_code}.final_genes.${output_suffix}.credible_sets.bed.gz"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:91ee6094cfae27626def40ff4557f30701cb4df99e2221991cc87a5361aaf796"
    preemptible: 1
    bootDiskSizeGb: "20"
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
  File raw_features_chromatin
  File raw_features_constraint
  File raw_features_variation
  File raw_features_merged_no_variation
  File raw_features_merged
  String output_suffix
  String rCNV_bucket

  command <<<
    set -e

    # Copy all fine-mapped gene lists
    mkdir finemap_stats/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/gene_burden/fine_mapping/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.*.tsv \
      finemap_stats/

    # Make input tsvs
    for wrapper in 1; do
      echo -e "Prior\tgrey70\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.naive_priors.genomic_features.tsv"
      echo -e "Posterior\t'#264653'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.genetics_only.genomic_features.tsv"
      echo -e "Genomic features\t'#331C8C'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.genomic_features.tsv"
      echo -e "Gene expression\t'#87C9F2'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.expression_features.tsv"
      echo -e "Chromatin\t'#3DAB9C'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.chromatin_features.tsv"
      echo -e "Gene constraint\t'#027831'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.constraint_features.tsv"
      echo -e "Variation\t'#E0CC70'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.variation_features.tsv"
      echo -e "Full (no var.)\t'#CF6576'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.merged_no_variation_features.tsv"
      echo -e "Full model\t'#AE3F9D'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.merged_features.tsv"
    done > finemap_roc_input.tsv
    for wrapper in 1; do
      echo -e "Genomic features\t'#331C8C'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.genomic_features.all_genes_from_blocks.tsv\t${raw_features_genomic}"
      echo -e "Gene expression\t'#87C9F2'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.expression_features.all_genes_from_blocks.tsv\t${raw_features_expression}"
      echo -e "Chromatin\t'#3DAB9C'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.chromatin_features.all_genes_from_blocks.tsv\t${raw_features_chromatin}"
      echo -e "Gene constraint\t'#027831'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.constraint_features.all_genes_from_blocks.tsv\t${raw_features_constraint}"
      echo -e "Variation\t'#E0CC70'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.variation_features.all_genes_from_blocks.tsv\t${raw_features_variation}"
      echo -e "Full (no var.)\t'#CF6576'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.merged_no_variation_features.all_genes_from_blocks.tsv\t${raw_features_merged_no_variation}"
      echo -e "Full model\t'#AE3F9D'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.${output_suffix}.gene_stats.merged_features.all_genes_from_blocks.tsv\t${raw_features_merged}"
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
    mkdir ${freq_code}_${CNV}_${output_suffix}_finemap_plots/
    /opt/rCNV2/analysis/genes/plot_finemap_results.R \
      finemap_roc_input.tsv \
      finemap_roc_truth_sets.tsv \
      ${freq_code}_${CNV}_${output_suffix}_finemap_plots/${freq_code}.${CNV}.${output_suffix}.finemap_results
    /opt/rCNV2/analysis/genes/plot_finemap_coefficients.R \
      finemap_feature_cor_input.tsv \
      ${freq_code}_${CNV}_${output_suffix}_finemap_plots/${freq_code}.${CNV}.${output_suffix}.finemap_feature_cors

    # Compress results
    tar -czvf \
      ${freq_code}_${CNV}_${output_suffix}_finemap_plots.tgz \
      ${freq_code}_${CNV}_${output_suffix}_finemap_plots
  >>>

  output {
    File plots_tarball = "${freq_code}_${CNV}_${output_suffix}_finemap_plots.tgz"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:91ee6094cfae27626def40ff4557f30701cb4df99e2221991cc87a5361aaf796"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }
}
