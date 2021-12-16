######################
#    rCNV Project    #
######################

# Copyright (c) 2019-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens per gene


import "https://api.firecloud.org/ga4gh/v1/tools/rCNV:scattered_gene_burden_perm_test/versions/10/plain-WDL/descriptor" as scattered_perm


workflow gene_burden_analysis {
  File phenotype_list
  File metacohort_list
  File metacohort_sample_table
  File gtf
  Int n_pheno_perms
  Int pad_controls
  Int min_probes_per_gene
  Float min_frac_controls_probe_exclusion
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  Float p_cutoff
  Int max_manhattan_neg_log10_p
  String meta_model_prefix
  Float winsorize_meta_z
  Int meta_min_cases
  Float meta_secondary_p_cutoff
  Int meta_nominal_cohorts_cutoff
  Float FDR_cutoff
  Float finemap_elnet_alpha
  Float finemap_elnet_l1_l2_mix
  Int finemap_cluster_distance
  Int finemap_nonsig_distance
  Float finemap_conf_pip
  Float finemap_vconf_pip
  File finemap_genomic_features
  File finemap_expression_features
  File finemap_chromatin_features
  File finemap_protein_features
  File finemap_constraint_features
  File finemap_variation_features
  File finemap_merged_no_variation_features
  File finemap_merged_features
  File raw_finemap_genomic_features
  File raw_finemap_expression_features
  File raw_finemap_chromatin_features
  File raw_finemap_protein_features
  File raw_finemap_constraint_features
  File raw_finemap_variation_features
  File raw_finemap_merged_no_variation_features
  File raw_finemap_merged_features
  String rCNV_bucket
  String rCNV_docker_assoc     # Note: strictly for dev/caching convenience; can collapse into a single docker at a later date
  String rCNV_docker_finemap   # Note: strictly for dev/caching convenience; can collapse into a single docker at a later date
  String athena_cloud_docker
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
        min_cds_ovr_del=min_cds_ovr_del,
        min_cds_ovr_dup=min_cds_ovr_dup,
        max_genes_per_cnv=max_genes_per_cnv,
        p_cutoff=p_cutoff,
        max_manhattan_neg_log10_p=max_manhattan_neg_log10_p,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_assoc,
        prefix=pheno[0],
        cache_string=fisher_cache_string
    }
  }

  # Determine conditional cohort exclusion list based on probe density
  call build_exclusion_list {
    input:
      genes_bed=flatten(rCNV_burden_test.stats_beds)[0],
      gtf_prefix=basename(gtf, '.gtf.gz'),
      min_probes_per_gene=min_probes_per_gene,
      min_frac_controls_probe_exclusion=min_frac_controls_probe_exclusion,
      metacohort_list=metacohort_list,
      rCNV_bucket=rCNV_bucket,
      athena_cloud_docker=athena_cloud_docker
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
        min_cds_ovr_del=min_cds_ovr_del,
        min_cds_ovr_dup=min_cds_ovr_dup,
        max_genes_per_cnv=max_genes_per_cnv,
        p_cutoff=p_cutoff,
        n_pheno_perms=n_pheno_perms,
        meta_model_prefix=meta_model_prefix,
        winsorize_meta_z=winsorize_meta_z,
        meta_min_cases=meta_min_cases,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_assoc,
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
        rCNV_docker=rCNV_docker_assoc,
        dummy_completion_markers=rCNV_perm_test.completion_marker,
        fdr_table_suffix="empirical_genome_wide_pval",
        p_val_column_name="meta_neg_log10_p"
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
        max_manhattan_neg_log10_p=max_manhattan_neg_log10_p,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_assoc,
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
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        meta_secondary_p_cutoff=meta_secondary_p_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        finemap_elnet_alpha=finemap_elnet_alpha,
        finemap_elnet_l1_l2_mix=finemap_elnet_l1_l2_mix,
        finemap_cluster_distance=finemap_cluster_distance,
        finemap_nonsig_distance=finemap_nonsig_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="genomic_features",
        gene_features=finemap_genomic_features,
        FDR_cutoff=FDR_cutoff,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_finemap
    }
    
    # Expression features
    call finemap_genes as finemap_expression {
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
        finemap_cluster_distance=finemap_cluster_distance,
        finemap_nonsig_distance=finemap_nonsig_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="expression_features",
        gene_features=finemap_expression_features,
        FDR_cutoff=FDR_cutoff,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_finemap
    }
    
    # Chromatin features
    call finemap_genes as finemap_chromatin {
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
        finemap_cluster_distance=finemap_cluster_distance,
        finemap_nonsig_distance=finemap_nonsig_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="chromatin_features",
        gene_features=finemap_chromatin_features,
        FDR_cutoff=FDR_cutoff,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_finemap
    }
    
    # Protein features
    call finemap_genes as finemap_protein {
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
        finemap_cluster_distance=finemap_cluster_distance,
        finemap_nonsig_distance=finemap_nonsig_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="protein_features",
        gene_features=finemap_protein_features,
        FDR_cutoff=FDR_cutoff,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_finemap
    }

    # Constraint features
    call finemap_genes as finemap_constraint {
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
        finemap_cluster_distance=finemap_cluster_distance,
        finemap_nonsig_distance=finemap_nonsig_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="constraint_features",
        gene_features=finemap_constraint_features,
        FDR_cutoff=FDR_cutoff,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_finemap
    }

    # Variation features
    call finemap_genes as finemap_variation {
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
        finemap_cluster_distance=finemap_cluster_distance,
        finemap_nonsig_distance=finemap_nonsig_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="variation_features",
        gene_features=finemap_variation_features,
        FDR_cutoff=FDR_cutoff,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_finemap
    }
    
    # Merged features
    call finemap_genes as finemap_merged {
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
        finemap_cluster_distance=finemap_cluster_distance,
        finemap_nonsig_distance=finemap_nonsig_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="merged_features",
        gene_features=finemap_merged_features,
        FDR_cutoff=FDR_cutoff,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_finemap
    }
    
    # Merged features (no variation)
    call finemap_genes as finemap_merged_no_variation {
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
        finemap_cluster_distance=finemap_cluster_distance,
        finemap_nonsig_distance=finemap_nonsig_distance,
        finemap_conf_pip=finemap_conf_pip,
        finemap_vconf_pip=finemap_vconf_pip,
        finemap_output_label="merged_no_variation_features",
        gene_features=finemap_merged_no_variation_features,
        FDR_cutoff=FDR_cutoff,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_finemap
    }
  }

  # Combine merged features (without variation) BED files as final output
  call merge_finemap_res as merge_finemap_res {
    input:
      gene_beds=finemap_merged_no_variation.sig_genes_bed,
      assoc_beds=finemap_merged_no_variation.sig_assocs_bed,
      credset_beds=finemap_merged_no_variation.cred_sets_bed,
      freq_code="rCNV",
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker_finemap
  }

  # Once complete, plot finemap results
  scatter( cnv in cnv_types ) {
    call plot_finemap_res as plot_finemap_res {
      input:
        completion_tokens=[finemap_genomic.completion_token, finemap_expression.completion_token, 
                           finemap_chromatin.completion_token, finemap_protein.completion_token, 
                           finemap_constraint.completion_token, finemap_variation.completion_token, 
                           finemap_merged_no_variation.completion_token, finemap_merged.completion_token],
        freq_code="rCNV",
        CNV=cnv,
        raw_features_genomic=raw_finemap_genomic_features,
        raw_features_expression=raw_finemap_expression_features,
        raw_features_chromatin=raw_finemap_chromatin_features,
        raw_features_protein=raw_finemap_protein_features,
        raw_features_constraint=raw_finemap_constraint_features,
        raw_features_variation=raw_finemap_variation_features,
        raw_features_merged_no_variation=raw_finemap_merged_no_variation_features,
        raw_features_merged=raw_finemap_merged_features,
        phenotype_list=phenotype_list,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_finemap
    }
  }

  output {
    File final_sig_genes = merge_finemap_res.final_genes
    File final_sig_associations = merge_finemap_res.final_associations
    File final_credible_sets = merge_finemap_res.final_credsets
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
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  Float p_cutoff
  Int max_manhattan_neg_log10_p
  String rCNV_bucket
  String rCNV_docker
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
          --p-col-name "fisher_neg_log10_p" \
          --p-is-neg-log10 \
          --max-neg-log10-p ${max_manhattan_neg_log10_p} \
          --cutoff ${p_cutoff} \
          --highlight-bed "${prefix}.highlight_regions.bed" \
          --highlight-name "Constrained genes associated with ${hpo}" \
          --label-prefix "$CNV" \
          --title "$title" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden"
      done

      # Generate Miami & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --miami \
        --p-col-name "fisher_neg_log10_p" \
        --p-is-neg-log10 \
        --max-neg-log10-p ${max_manhattan_neg_log10_p} \
        --cutoff ${p_cutoff} \
        --highlight-bed "${prefix}.highlight_regions.bed" \
        --highlight-name "Constrained genes associated with ${hpo}" \
        --label-prefix "DUP" \
        --highlight-bed-2 "${prefix}.highlight_regions.bed" \
        --highlight-name-2 "Constrained genes associated with ${hpo}" \
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
    docker: "${rCNV_docker}"
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
  String athena_cloud_docker

  command <<<
    set -euo pipefail

    # Download probesets (note: requires permissions)
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/control_probesets \
      ./

    # Subset gene coordinates to minimal BED4
    zcat ${genes_bed} | cut -f1-4 | bgzip -c > gene_coords.bed.gz

    # Make input for conditional exclusion script
    for file in control_probesets/*bed.gz; do
      echo -e "$file\t$( basename $file | sed 's/\.bed\.gz//g' )"
    done > probeset_tracks.tsv

    # Clone rCNV2 repo (not present in athena-cloud Docker)
    git clone https://github.com/talkowski-lab/rCNV2.git

    # Build conditional exclusion list
    rCNV2/data_curation/other/probe_based_exclusion.py \
      --outfile ${gtf_prefix}.cohort_exclusion.bed.gz \
      --probecounts-outfile ${gtf_prefix}.probe_counts.bed.gz \
      --control-mean-counts-outfile ${gtf_prefix}.mean_probe_counts_per_cohort.bed.gz \
      --frac-pass-outfile ${gtf_prefix}.frac_passing.bed.gz \
      --min-probes ${min_probes_per_gene} \
      --min-frac-samples ${min_frac_controls_probe_exclusion} \
      --keep-n-columns 4 \
      --bgzip \
      gene_coords.bed.gz \
      probeset_tracks.tsv \
      control_probesets/rCNV.control_counts_by_array.tsv \
      <( fgrep -v mega ${metacohort_list} )

    # Copy to Google bucket for storage
    gsutil -m cp \
      ${gtf_prefix}.cohort_exclusion.bed.gz \
      ${gtf_prefix}.probe_counts.bed.gz \
      ${gtf_prefix}.frac_passing.bed.gz \
      ${rCNV_bucket}/analysis/analysis_refs/
  >>>

  runtime {
    docker: "${athena_cloud_docker}"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    File exclusion_bed = "${gtf_prefix}.cohort_exclusion.bed.gz"
    File probe_counts = "${gtf_prefix}.probe_counts.bed.gz"
    File control_mean_counts = "${gtf_prefix}.mean_probe_counts_per_cohort.bed.gz"
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
  String rCNV_docker
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
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "32 GB"
    disks: "local-disk 100 HDD"
    bootDiskSizeGb: "20"
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
  Int max_manhattan_neg_log10_p
  String meta_model_prefix
  Float winsorize_meta_z
  Int meta_min_cases
  String rCNV_bucket
  String rCNV_docker
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
        --model ${meta_model_prefix} \
        --conditional-exclusion ${exclusion_bed} \
        --p-is-neg-log10 \
        --spa \
        --spa-exclude /opt/rCNV2/refs/lit_GDs.all.$CNV.bed.gz \
        --winsorize ${winsorize_meta_z} \
        --min-cases ${meta_min_cases} \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed
      bgzip -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed
      tabix -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "meta_neg_log10_p" \
        --p-is-neg-log10 \
        --max-neg-log10-p ${max_manhattan_neg_log10_p} \
        --cutoff $meta_p_cutoff \
        --highlight-bed "${prefix}.highlight_regions.bed" \
        --highlight-name "Constrained genes associated with ${hpo}" \
        --label-prefix "$CNV" \
        --title "$title" \
        "${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz" \
        "${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis"
    done

    # Generate Miami & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --miami \
      --p-col-name "meta_neg_log10_p" \
      --p-is-neg-log10 \
      --max-neg-log10-p ${max_manhattan_neg_log10_p} \
      --cutoff $DUP_p_cutoff \
      --highlight-bed "${prefix}.highlight_regions.bed" \
      --highlight-name "Constrained genes associated with ${hpo}" \
      --label-prefix "DUP" \
      --cutoff-2 $DEL_p_cutoff \
      --highlight-bed-2 "${prefix}.highlight_regions.bed" \
      --highlight-name-2 "Constrained genes associated with ${hpo}" \
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
    docker: "${rCNV_docker}"
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
  Float finemap_elnet_alpha
  Float finemap_elnet_l1_l2_mix
  Int finemap_cluster_distance
  Int finemap_nonsig_distance
  Float finemap_conf_pip
  Float finemap_vconf_pip
  String finemap_output_label
  File gene_features
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -e

    # Copy p-value cutoff tables
    find / -name "*gene_burden.${freq_code}.*.bonferroni_pval.hpo_cutoffs.tsv*" \
    | xargs -I {} mv {} ./

    # Copy association stats & other reference data from the project Google Bucket (note: requires permissions)
    mkdir stats
    gsutil -m cp \
      ${rCNV_bucket}/analysis/gene_burden/**.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz \
      stats/
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/genes/gene_lists \
      ./
    gsutil -m cp \
      gs://rcnv_project/analysis/gene_scoring/gene_lists/* \
      ./gene_lists/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/paper/data/large_segments/clustered_nahr_regions.bed.gz \
      ./
    gsutil -m cp \
      ${rCNV_bucket}/analysis/analysis_refs/rCNV2.hpos_by_severity.*list \
      ./

    # Make list of genes from predicted NAHR-mediated CNV regions for training exclusion
    zcat clustered_nahr_regions.bed.gz | fgrep -v "#" \
    | awk -v FS="\t" '{ if ($5>0) print $NF }' \
    | sed 's/;/\n/g' | sort | uniq > nahr.genes.list

    # Write tsv inputs for developmental & adult-onset subgroups (to be finemapped separately)
    for subgroup in developmental adult; do
      while read prefix hpo; do
        for wrapper in 1; do
          echo "$hpo"
          echo "stats/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz"
          awk -v x=$prefix -v FS="\t" '{ if ($1==x) print $2 }' \
            gene_burden.${freq_code}.${CNV}.bonferroni_pval.hpo_cutoffs.tsv
        done | paste -s
      done < <( fgrep -wf rCNV2.hpos_by_severity.$subgroup.list ${phenotype_list} ) \
      > ${freq_code}.${CNV}.gene_fine_mapping.stats_input.$subgroup.tsv
    done

    # Make lists of a priori "true causal" genes for null variance estimation
    case ${CNV} in
      "DEL")
        echo "gene_lists/gold_standard.haploinsufficient.genes.list" > known_causal_gene_lists.developmental.tsv
        ;;
      "DUP")
        echo "gene_lists/gold_standard.triplosensitive.genes.list" > known_causal_gene_lists.developmental.tsv
        ;;
    esac
    echo "gene_lists/HP0000118.HPOdb.genes.list" > known_causal_gene_lists.adult.tsv

    # Run functional fine-mapping procedure
    for subgroup in developmental adult; do
      echo -e "\n...Starting fine-mapping for $subgroup phenotypes...\n"
      /opt/rCNV2/analysis/genes/finemap_genes.py \
        --cnv ${CNV} \
        --secondary-p-cutoff ${meta_secondary_p_cutoff} \
        --min-nominal ${meta_nominal_cohorts_cutoff} \
        --secondary-or-nominal \
        --fdr-q-cutoff ${FDR_cutoff} \
        --secondary-for-fdr \
        --regularization-alpha ${finemap_elnet_alpha} \
        --regularization-l1-l2-mix ${finemap_elnet_l1_l2_mix} \
        --distance ${finemap_cluster_distance} \
        --nonsig-distance ${finemap_nonsig_distance} \
        --training-exclusion nahr.genes.list \
        --use-max-pip-per-gene \
        --confident-pip ${finemap_conf_pip} \
        --very-confident-pip ${finemap_vconf_pip} \
        --known-causal-gene-lists known_causal_gene_lists.$subgroup.tsv \
        --outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.$subgroup.tsv \
        --sig-genes-bed ${freq_code}.${CNV}.final_genes.genes.${finemap_output_label}.$subgroup.bed \
        --sig-assoc-bed ${freq_code}.${CNV}.final_genes.associations.${finemap_output_label}.$subgroup.bed \
        --sig-credsets-bed ${freq_code}.${CNV}.final_genes.credible_sets.${finemap_output_label}.$subgroup.bed \
        --all-genes-outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.${finemap_output_label}.$subgroup.tsv \
        --naive-outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.$subgroup.tsv \
        --genetic-outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.$subgroup.tsv \
        --coeffs-out ${freq_code}.${CNV}.gene_fine_mapping.logit_coeffs.${finemap_output_label}.$subgroup.tsv \
        ${freq_code}.${CNV}.gene_fine_mapping.stats_input.$subgroup.tsv \
        ${gene_features} \
        ${metacohort_sample_table}
    done

    # Merge & sort outputs for developmental & adult-onset subgroups
    mkdir finemapping_merged_outputs/
    cat <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.developmental.tsv ) \
        <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.adult.tsv ) \
    | sort -nrk4,4 \
    | cat <( grep -e '^#' ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.developmental.tsv ) - \
    > finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.tsv
    cat <( fgrep -v "#" ${freq_code}.${CNV}.final_genes.associations.${finemap_output_label}.developmental.bed ) \
        <( fgrep -v "#" ${freq_code}.${CNV}.final_genes.associations.${finemap_output_label}.adult.bed ) \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k6,6V \
    | cat <( grep -e '^#' ${freq_code}.${CNV}.final_genes.associations.${finemap_output_label}.developmental.bed ) - \
    > finemapping_merged_outputs/${freq_code}.${CNV}.final_genes.associations.${finemap_output_label}.bed
    cat <( fgrep -v "#" ${freq_code}.${CNV}.final_genes.credible_sets.${finemap_output_label}.developmental.bed ) \
        <( fgrep -v "#" ${freq_code}.${CNV}.final_genes.credible_sets.${finemap_output_label}.adult.bed ) \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k6,6V \
    | cat <( grep -e '^#' ${freq_code}.${CNV}.final_genes.credible_sets.${finemap_output_label}.developmental.bed ) - \
    > finemapping_merged_outputs/${freq_code}.${CNV}.final_genes.credible_sets.${finemap_output_label}.bed
    cat <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.${finemap_output_label}.developmental.tsv ) \
        <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.${finemap_output_label}.adult.tsv ) \
    | sort -nrk4,4 \
    | cat <( grep -e '^#' ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.${finemap_output_label}.developmental.tsv ) - \
    > finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.${finemap_output_label}.tsv
    cat <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.developmental.tsv ) \
        <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.adult.tsv ) \
    | sort -nrk4,4 \
    | cat <( grep -e '^#' ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.developmental.tsv ) - \
    > finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.tsv
    cat <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.developmental.tsv ) \
        <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.adult.tsv ) \
    | sort -nrk4,4 \
    | cat <( grep -e '^#' ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.developmental.tsv ) - \
    > finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.tsv
    cat <( fgrep -v "#" ${freq_code}.${CNV}.final_genes.genes.${finemap_output_label}.developmental.bed \
           | awk -v OFS="\t" '{ print $0, "developmental" }' ) \
        <( fgrep -v "#" ${freq_code}.${CNV}.final_genes.genes.${finemap_output_label}.adult.bed \
           | awk -v OFS="\t" '{ print $0, "adult" }' ) \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k6,6V \
    | cat <( grep -e "^#" ${freq_code}.${CNV}.final_genes.genes.${finemap_output_label}.developmental.bed \
           | awk -v OFS="\t" '{ print $0, "finemapping_model" }' ) \
    > finemapping_merged_outputs/${freq_code}.${CNV}.final_genes.genes.${finemap_output_label}.bed

    # Copy results to output bucket
    gsutil -m cp \
      finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.*tsv \
      ${rCNV_bucket}/analysis/gene_burden/fine_mapping/

    # Make completion token to track caching
    echo "Done" > completion.txt
  >>>

  output {
    File finemapped_output = "finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.tsv"
    File finemapped_output_all_genes = "finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.${finemap_output_label}.tsv"
    File naive_output = "finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.tsv"
    File genetic_output = "finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.tsv"
    File logit_coeffs_developmental = "${freq_code}.${CNV}.gene_fine_mapping.logit_coeffs.${finemap_output_label}.developmental.tsv"
    File logit_coeffs_adult = "${freq_code}.${CNV}.gene_fine_mapping.logit_coeffs.${finemap_output_label}.adult.tsv"
    File sig_genes_bed = "finemapping_merged_outputs/${freq_code}.${CNV}.final_genes.genes.${finemap_output_label}.bed"
    File sig_assocs_bed = "finemapping_merged_outputs/${freq_code}.${CNV}.final_genes.associations.${finemap_output_label}.bed"
    File cred_sets_bed = "finemapping_merged_outputs/${freq_code}.${CNV}.final_genes.credible_sets.${finemap_output_label}.bed"
    File completion_token = "completion.txt"
  }

  runtime {
    docker: "${rCNV_docker}"
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
  String rCNV_bucket
  String rCNV_docker

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
    > ${freq_code}.final_genes.genes.bed.gz

    # Merge associations
    cat ${sep=" " assoc_beds} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k6,6V \
    | cat assoc_header.tsv - \
    | bgzip -c \
    > ${freq_code}.final_genes.associations.bed.gz

    # Merge genes
    cat ${sep=" " credset_beds} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n -k5,5V -k6,6V \
    | cat credset_header.tsv - \
    | bgzip -c \
    > ${freq_code}.final_genes.credible_sets.bed.gz

    # Copy final files to results bucket (note: requires permissions)
    gsutil -m cp \
      ${freq_code}.final_genes.*.bed.gz \
      ${rCNV_bucket}/results/gene_association/
  >>>

  output {
    File final_genes = "${freq_code}.final_genes.genes.bed.gz"
    File final_associations = "${freq_code}.final_genes.associations.bed.gz"
    File final_credsets = "${freq_code}.final_genes.credible_sets.bed.gz"
  }

  runtime {
    docker: "${rCNV_docker}"
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
  File raw_features_protein
  File raw_features_constraint
  File raw_features_variation
  File raw_features_merged_no_variation
  File raw_features_merged
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -e

    # Copy all fine-mapped gene lists
    mkdir finemap_stats/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/gene_burden/fine_mapping/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.*.tsv \
      finemap_stats/

    # Make input tsvs
    for wrapper in 1; do
      echo -e "Prior\tgrey70\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.tsv"
      echo -e "Posterior\t'#264653'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.tsv"
      echo -e "Genomic features\t'#490C65'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genomic_features.tsv"
      echo -e "Gene expression\t'#BA7FD0'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.expression_features.tsv"
      echo -e "Chromatin\t'#001588'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.chromatin_features.tsv"
      echo -e "Protein\t'#0180C9'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.protein_features.tsv"
      echo -e "Gene constraint\t'#F6313E'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.constraint_features.tsv"
      echo -e "Variation\t'#FFA300'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.variation_features.tsv"
      echo -e "Full (no var.)\t'#46A040'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.merged_no_variation_features.tsv"
      echo -e "Full model\t'#00441B'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.merged_features.tsv"
    done > finemap_roc_input.tsv
    for wrapper in 1; do
      echo -e "Genomic features\t'#490C65'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.genomic_features.tsv\t${raw_features_genomic}"
      echo -e "Gene expression\t'#BA7FD0'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.expression_features.tsv\t${raw_features_expression}"
      echo -e "Chromatin\t'#001588'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.chromatin_features.tsv\t${raw_features_chromatin}"
      echo -e "Protein\t'#0180C9'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.protein_features.tsv\t${raw_features_protein}"
      echo -e "Gene constraint\t'#F6313E'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.constraint_features.tsv\t${raw_features_constraint}"
      echo -e "Variation\t'#FFA300'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.variation_features.tsv\t${raw_features_variation}"
      echo -e "Full (no var.)\t'#46A040'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv\t${raw_features_merged_no_variation}"
      echo -e "Full model\t'#00441B'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_features.tsv\t${raw_features_merged}"
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

    fgrep -wvf \
      gene_lists/HP0000118.HPOdb.genes.list \
      gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list \
    > true_negatives.genes.tsv

    # Plot finemapping QC & feature correlations
    for subgroup in developmental adult; do
      mkdir ${freq_code}_${CNV}_finemap_plots_$subgroup/
      /opt/rCNV2/analysis/genes/plot_finemap_results.R \
        finemap_roc_input.tsv \
        finemap_roc_truth_sets.tsv \
        true_negatives.genes.tsv \
        $subgroup \
        ${freq_code}_${CNV}_finemap_plots_$subgroup/${freq_code}.${CNV}.$subgroup.finemap_results
      /opt/rCNV2/analysis/genes/plot_finemap_coefficients.R \
        finemap_feature_cor_input.tsv \
        $subgroup \
        ${freq_code}_${CNV}_finemap_plots_$subgroup/${freq_code}.${CNV}.$subgroup.finemap_feature_cors

      # Compress results
      tar -czvf \
        ${freq_code}_${CNV}_finemap_plots_$subgroup.tgz \
        ${freq_code}_${CNV}_finemap_plots_$subgroup
    done
  >>>

  output {
    File adult_plots_tarball = "${freq_code}_${CNV}_finemap_plots_adult.tgz"
    File developmental_plots_tarball = "${freq_code}_${CNV}_finemap_plots_developmental.tgz"
  }

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }
}
