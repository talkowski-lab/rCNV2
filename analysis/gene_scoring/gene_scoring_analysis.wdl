######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Gene scoring for haploinsufficiency and triplosensitivity


workflow gene_burden_analysis {
  String prefix
  File metacohort_list
  File metacohort_sample_table
  File gtf
  File gtf_idx
  Int pad_controls
  Int max_cnv_size
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  String meta_model_prefix
  Float winsorize_meta_z
  Int meta_min_cases
  Float max_standard_error
  Float prior_frac
  Float prior_lnor_thresholding_pct
  File training_excludelist
  File gene_features
  File raw_gene_features
  Float max_true_bfdp
  Float min_false_bfdp
  String rCNV_bucket
  String rCNV_docker
  String rCNV_docker_scoring #Note: this is separate just for development purposes; can be collapsed in at a later date
  File contiglist

  Array[Array[String]] contigs = read_tsv(contiglist)

  Array[String] cnv_types = ["DEL", "DUP"]

  Array[String] models = ["logit", "svm", "randomforest", "lda", "naivebayes", "neuralnet", "gbdt", "knn"]

  # Scatter over contigs (for speed)
  scatter ( contig in contigs ) {
    # Gather data required for parameter optimization
    call get_annotated_cnvs as get_optimization_data {
      input:
        hpo="HP:0000118",
        prefix="HP0000118",
        contig=contig[0],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        gtf=gtf,
        gtf_idx=gtf_idx,
        pad_controls=pad_controls,
        min_cds_ovr_del=min_cds_ovr_del,
        min_cds_ovr_dup=min_cds_ovr_dup,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker
    }
    # Run assocation tests per cohort
    call burden_test as rCNV_burden_test {
      input:
        prefix=prefix,
        contig=contig[0],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        gtf=gtf,
        gtf_idx=gtf_idx,
        pad_controls=pad_controls,
        max_cnv_size=max_cnv_size,
        min_cds_ovr_del=min_cds_ovr_del,
        min_cds_ovr_dup=min_cds_ovr_dup,
        max_genes_per_cnv=max_genes_per_cnv,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker
    }
  }

  # Merge annotated CNVs for parameter optimization
  scatter ( cnv in cnv_types ) {
    call merge_annotated_cnvs {
      input:
        annotated_cnv_beds=get_optimization_data.annotated_cnvs_bed,
        annotated_cnv_bed_idxs=get_optimization_data.annotated_cnvs_idx,
        metacohort_list=metacohort_list,
        CNV=cnv,
        freq_code="rCNV",
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker
    }
  }

  # Merge association statistics and run meta-analysis
  scatter ( CNV in cnv_types ) {
    call merge_and_meta_analysis as rCNV_meta_analysis {
      input:
        stats_beds=rCNV_burden_test.stats_beds,
        stats_beds_idxs=rCNV_burden_test.stats_bed_idxs,
        prefix=prefix,
        CNV=CNV,
        metacohort_list=metacohort_list,
        meta_model_prefix=meta_model_prefix,
        winsorize_meta_z=winsorize_meta_z,
        meta_min_cases=meta_min_cases,
        freq_code="rCNV",
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker
    }
    call get_underpowered_genes {
      input:
        meta_stats=rCNV_meta_analysis.meta_stats_bed,
        prefix=prefix,
        CNV=CNV,
        max_standard_error=max_standard_error,
        freq_code="rCNV",
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_scoring
    }
  }

  # Estimate prior effect sizes and compute BFDPs
  call calc_priors_bfdp {
    input:
      del_meta_stats=rCNV_meta_analysis.meta_stats_bed[0],
      dup_meta_stats=rCNV_meta_analysis.meta_stats_bed[1],
      underpowered_genes=get_underpowered_genes.underpowered_genes,
      freq_code="rCNV",
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker_scoring,
      prior_frac=prior_frac,
      prior_lnor_thresholding_pct=prior_lnor_thresholding_pct
  }

  # Score genes for each model
  scatter ( model in models ) {
    call score_genes as score_genes_DEL {
      input:
        CNV="DEL",
        BFDP_stats=calc_priors_bfdp.del_bfdp,
        excludelist=training_excludelist,
        underpowered_genes=get_underpowered_genes.underpowered_genes[0],
        gene_features=gene_features,
        model=model,
        max_true_bfdp=max_true_bfdp,
        min_false_bfdp=min_false_bfdp,
        freq_code="rCNV",
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_scoring
    }
    call score_genes as score_genes_DUP {
      input:
        CNV="DUP",
        BFDP_stats=calc_priors_bfdp.dup_bfdp,
        excludelist=training_excludelist,
        underpowered_genes=get_underpowered_genes.underpowered_genes[1],
        gene_features=gene_features,
        model=model,
        max_true_bfdp=max_true_bfdp,
        min_false_bfdp=min_false_bfdp,
        freq_code="rCNV",
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker_scoring
    }
  }

  # Score with ensemble classifier, and return updated array of all scores + ensemble
  call score_ensemble as score_ensemble_DEL {
    input:
      CNV="DEL",
      scores=score_genes_DEL.scores_tsv,
      freq_code="rCNV",
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker_scoring
  }
  call score_ensemble as score_ensemble_DUP {
    input:
      CNV="DUP",
      scores=score_genes_DUP.scores_tsv,
      freq_code="rCNV",
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker_scoring
  }

  # Determine best model & QC final scores
  call qc_scores {
    input:
      del_scores=score_ensemble_DEL.all_scores,
      dup_scores=score_ensemble_DUP.all_scores,
      models=models,
      raw_gene_features=raw_gene_features,
      freq_code="rCNV",
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker_scoring
  }

  output {
    File gene_scores = qc_scores.final_scores
  }

}


# Collect data necessary for parameter optimization
task get_annotated_cnvs {
  String hpo
  String prefix
  String contig
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  File gtf
  File gtf_idx
  Int pad_controls
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -e

    # Copy CNV data and constrained gene coordinates
    mkdir cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/genes ./
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/refs/GRCh37.*.bed.gz refs/
    gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

    # Extract contig of interest from GTF
    tabix ${gtf} ${contig} | bgzip -c > subset.gtf.gz

    # Iterate over metacohorts
    while read meta cohorts; do
      echo $meta

      # Extract CNVs of interest from BED
      tabix -f cleaned_cnv/$meta.${freq_code}.bed.gz
      tabix -h cleaned_cnv/$meta.${freq_code}.bed.gz ${contig} \
      | bgzip -c > $meta.${contig}.bed.gz
      cnv_bed="$meta.${contig}.bed.gz"

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
          --max-cnv-size 300000000 \
          --min-cds-ovr $min_cds_ovr \
          --max-genes 25000 \
          -t $CNV \
          --hpo ${hpo} \
          --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
          --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
          --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
          -z \
          --verbose \
          --cnvs-out "$meta.${prefix}.${freq_code}.$CNV.genes_per_cnv.${contig}.bed.gz" \
          -o "$meta.${prefix}.${freq_code}.$CNV.genes_per_cnv.${contig}.bed.gz" \
          "$cnv_bed" \
          subset.gtf.gz
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.genes_per_cnv.${contig}.bed.gz"
      done
    done < ${metacohort_list}
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] annotated_cnvs_bed = glob("*.genes_per_cnv.${contig}.bed.gz")
    Array[File] annotated_cnvs_idx = glob("*.genes_per_cnv.${contig}.bed.gz.tbi")
  }
}


# Run burden test for a single chromosome for all cohorts
task burden_test {
  String prefix
  String contig
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  File gtf
  File gtf_idx
  Int pad_controls
  Int max_cnv_size
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -e

    # Copy CNV data and constrained gene coordinates
    mkdir cleaned_cnv/ refs/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/genes ./
    gsutil -m cp \
      ${rCNV_bucket}/refs/GRCh37.*.bed.gz \
      ${rCNV_bucket}/analysis/paper/data/large_segments/loose_unclustered_nahr_regions.bed.gz \
      gs://rcnv_project/analysis/analysis_refs/* \
      refs/

    # Extract contig of interest from GTF
    tabix ${gtf} ${contig} | bgzip -c > ${contig}.gtf.gz

    # Create string of HPOs to keep
    keep_hpos=$( cat refs/rCNV2.hpos_by_severity.developmental.list | paste -s -d\; )

    # Iterate over metacohorts to compute single-cohort stats
    while read meta cohorts; do
      echo $meta

      # Exclude all NAHR CNVs from burden testing
      bedtools intersect -v -r -f 0.5 \
        -a cleaned_cnv/$meta.${freq_code}.bed.gz \
        -b refs/loose_unclustered_nahr_regions.bed.gz \
      | bgzip -c \
      > cleaned_cnv/$meta.${freq_code}.no_NAHR.bed.gz
      tabix -f cleaned_cnv/$meta.${freq_code}.no_NAHR.bed.gz

      # Set metacohort-specific parameters
      cnv_bed="cleaned_cnv/$meta.${freq_code}.no_NAHR.bed.gz"
      effective_case_n=$( fgrep -w $meta refs/rCNV2.hpos_by_severity.developmental.counts.tsv | cut -f2 )

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
          --max-cnv-size ${max_cnv_size} \
          --min-cds-ovr $min_cds_ovr \
          --max-genes ${max_genes_per_cnv} \
          -t $CNV \
          --hpo "$keep_hpos" \
          --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
          --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
          --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
          -z \
          --verbose \
          -o "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.${contig}.bed.gz" \
          "$cnv_bed" \
          ${contig}.gtf.gz
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.${contig}.bed.gz"

        # Perform burden test
        /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --effective-case-n $effective_case_n \
          --keep-n-columns 4 \
          --bgzip \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.${contig}.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.${contig}.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.${contig}.bed.gz"
      done
    done < <( fgrep -v "mega" ${metacohort_list} )
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] stats_beds = glob("*.gene_burden.stats.${contig}.bed.gz")
    Array[File] stats_bed_idxs = glob("*.gene_burden.stats.${contig}.bed.gz.tbi")
  }
}


# Merge annotated CNVs per metacohort for parameter optimization
task merge_annotated_cnvs {
  Array[Array[File]] annotated_cnv_beds
  Array[Array[File]] annotated_cnv_bed_idxs
  File metacohort_list
  String freq_code
  String CNV
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -e

    # Make list of stats files to be considered
    find / -name "*.${freq_code}.*.genes_per_cnv.*.bed.gz" \
    > stats.paths.list
    # Debug: print paths to stdout
    cat stats.paths.list

    # Merge list of genes hit per cohort, per phenotype
    mkdir optimization_data/
    zcat $( sed -n '1p' stats.paths.list ) | sed -n '1p' > header.tsv
    while read meta cohort; do
      zcat $( fgrep $meta stats.paths.list ) \
      | fgrep -w ${CNV} \
      | grep -ve '^#' \
      | sort -Vk4,4 \
      | uniq \
      | cat header.tsv - \
      | awk -v FS="\t" -v OFS="\t" '{ print $4, $6, $7, $9 }' \
      | gzip -c \
      > optimization_data/${freq_code}.${CNV}.$meta.genes_per_cnv.tsv.gz
    done < ${metacohort_list}

    # Copy output to Google bucket
    gsutil -m cp -r \
      optimization_data \
      ${rCNV_bucket}/analysis/gene_scoring/
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
  }

  output {}
}


# Merge burden tests and perform meta-analysis per cohort per CNV type
task merge_and_meta_analysis {
  Array[Array[File]] stats_beds
  Array[Array[File]] stats_beds_idxs
  String prefix
  String CNV
  File metacohort_list
  String meta_model_prefix
  Float winsorize_meta_z
  Int meta_min_cases
  String freq_code
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -euo pipefail

    # Copy necessary reference files
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/* refs/
    export exclusion_bed=refs/gencode.v19.canonical.pext_filtered.cohort_exclusion.bed.gz

    # Make list of stats files to be considered
    find / -name "*.${prefix}.${freq_code}.${CNV}.gene_burden.stats.*.bed.gz" \
    > stats.paths.list

    # Merge burden stats per cohort
    zcat $( sed -n '1p' stats.paths.list ) | sed -n '1p' > header.tsv
    while read meta cohort; do
      zcat $( fgrep $meta stats.paths.list ) \
      | grep -ve '^#' \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat header.tsv - \
      | bgzip -c \
      > "$meta.${prefix}.${freq_code}.${CNV}.gene_burden.stats.bed.gz"
      tabix -f "$meta.${prefix}.${freq_code}.${CNV}.gene_burden.stats.bed.gz"
    done < <( fgrep -v mega ${metacohort_list} )

    # Make input for meta-analysis
    while read meta cohorts; do
      echo -e "$meta\t$meta.${prefix}.${freq_code}.${CNV}.gene_burden.stats.bed.gz"
    done < <( fgrep -v mega ${metacohort_list} ) \
    > ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt
    
    # Run meta-analysis
    /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
      --model ${meta_model_prefix} \
      --conditional-exclusion $exclusion_bed \
      --p-is-neg-log10 \
      --spa \
      --spa-exclude /opt/rCNV2/refs/lit_GDs.all.${CNV}.bed.gz \
      --winsorize ${winsorize_meta_z} \
      --min-cases ${meta_min_cases} \
      --keep-n-columns 4 \
      ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt \
      ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed
    bgzip -f ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed
    tabix -p bed -f ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz

    # Copy meta-analysis results to Google bucket
    gsutil -m cp \
      ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz* \
      ${rCNV_bucket}/analysis/gene_scoring/data/
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "8 GB"
    bootDiskSizeGb: "30"
    disks: "local-disk 50 HDD"
  }

  output {
    File meta_stats_bed = "${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz"
    File meta_stats_bed_idx = "${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz.tbi"
  }
}


# Get list of genes with insufficient CNV data to be used for model training
task get_underpowered_genes {
  File meta_stats
  String prefix
  String CNV
  Float max_standard_error
  String freq_code
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -euo pipefail

    # Extract list of genes with standard error < max_standard_error (to be used as excludelist later for training)
    /opt/rCNV2/analysis/gene_scoring/get_underpowered_genes.R \
      --max-se "${max_standard_error}" \
      ${meta_stats} \
      ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed
    awk -v OFS="\t" '{ print $1, $2, $3, $4 }' \
      ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed \
    | bgzip -c \
    > ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz
    tabix -p bed -f ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz

    # Copy meta-analysis results to Google bucket
    gsutil -m cp \
      ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz \
      ${rCNV_bucket}/analysis/gene_scoring/data/
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "8 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File underpowered_genes = "${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz"
  }
}


# Estimate prior effect sizes and compute BFDPs per gene
task calc_priors_bfdp {
  File del_meta_stats
  File dup_meta_stats
  Array[File] underpowered_genes
  String freq_code
  String rCNV_bucket
  String rCNV_docker
  Float prior_frac
  Float prior_lnor_thresholding_pct

  command <<<
    set -e

    # Localize necessary references
    gsutil -m cp -r \
      ${rCNV_bucket}/analysis/gene_scoring/refs/${freq_code}.gene_scoring.training_gene_excludelist.bed.gz \
      ${rCNV_bucket}/analysis/gene_scoring/optimization_data \
      ${rCNV_bucket}/analysis/gene_scoring/gene_lists/*.genes.list \
      ${rCNV_bucket}/cleaned_data/genes/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
      ./
    mkdir refs/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/analysis_refs/* \
      refs/

    # Locate underpowered genes
    find / -name "*.gene_burden.underpowered_genes.bed.gz" | xargs -I {} mv {} ./

    # Create CNV-type-specific gene excludelists
    for CNV in DEL DUP; do
      zcat *.$CNV.gene_burden.underpowered_genes.bed.gz \
        ${freq_code}.gene_scoring.training_gene_excludelist.bed.gz \
      | fgrep -v "#" | cut -f4 | sort -V | uniq \
      > ${freq_code}.$CNV.training_excludelist.genes.list
    done
    
    # Prepare files for on-the fly meta-analysis of CNV effect sizes
    echo -e "cohort\tn_case\tn_control\tDEL_path\tDUP_path" \
    > prior_estimation.meta_inputs.tsv
    zcat optimization_data/rCNV.DEL.meta1.genes_per_cnv.tsv.gz | head -n1 \
    > optimization_data/opt_data.header.tsv
    while read meta cohorts; do
      # Subset CNVs to developmental cases & controls
      for CNV in DEL DUP; do
        zcat optimization_data/rCNV.$CNV.$meta.genes_per_cnv.tsv.gz \
        | fgrep -wf <( cat refs/rCNV2.hpos_by_severity.developmental.list \
                           <( echo "HEALTHY_CONTROL" ) ) \
        | sort -Vk1,1 | cat optimization_data/opt_data.header.tsv - | gzip -c \
        > optimization_data/rCNV.$meta.annotated_developmental_and_control.$CNV.tsv.gz
      done

      # Write info to meta-analysis input
      for dummy in 1; do
        fgrep -w $meta refs/rCNV2.hpos_by_severity.developmental.counts.tsv
        cidx=$( head -n1 refs/HPOs_by_metacohort.table.tsv | sed 's/\t/\n/g' \
                | awk -v OFS="\t" '{ print NR, $0 }' | fgrep -w $meta | cut -f1 )
        fgrep -w "HEALTHY_CONTROL" refs/HPOs_by_metacohort.table.tsv | cut -f$cidx
        for CNV in DEL DUP; do
          echo optimization_data/rCNV.$meta.annotated_developmental_and_control.$CNV.tsv.gz
        done
      done | paste -s >> prior_estimation.meta_inputs.tsv
    done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt )

    # Compute prior effect sizes
    /opt/rCNV2/analysis/gene_scoring/estimate_prior_effect_sizes.R \
      prior_estimation.meta_inputs.tsv \
      ${freq_code}.gene_scoring.training_gene_excludelist.bed.gz \
      gnomad.v2.1.1.lof_constrained.genes.list \
      gold_standard.haploinsufficient.genes.list \
      gold_standard.haplosufficient.genes.list \
      gold_standard.triplosensitive.genes.list \
      gold_standard.triploinsensitive.genes.list \
      ${freq_code}.prior_estimation
    awk -v FS="\t" '{ if ($1=="theta0" && $2=="DEL") print $3 }' \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
    > theta0_del.tsv
    awk -v FS="\t" '{ if ($1=="theta0" && $2=="DUP") print $3 }' \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
    > theta0_dup.tsv
    awk -v FS="\t" '{ if ($1=="theta1" && $2=="DEL") print $3 }' \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
    > theta1_del.tsv
    awk -v FS="\t" '{ if ($1=="theta1" && $2=="DUP") print $3 }' \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
    > theta1_dup.tsv

    # Compute BFDP per gene
    for CNV in DEL DUP; do
      # Set CNV-specific variables
      case $CNV in
        "DEL")
          theta0=$( cat theta0_del.tsv )
          theta1=$( cat theta1_del.tsv )
          statsbed=${del_meta_stats}
          ;;
        "DUP")
          theta0=$( cat theta0_dup.tsv )
          theta1=$( cat theta1_dup.tsv )
          statsbed=${dup_meta_stats}
          ;;
      esac

      # Compute BF & BFDR for all genes
      /opt/rCNV2/analysis/gene_scoring/calc_gene_bfs.py \
        --theta0 $theta0 \
        --theta1 $theta1 \
        --var0 1 \
        --prior ${prior_frac} \
        --blacklist ${freq_code}.$CNV.training_excludelist.genes.list \
        --outfile ${freq_code}.$CNV.gene_abfs.tsv \
        "$statsbed"
    done

    # Copy results to Google bucket
    gsutil -m cp \
      ${freq_code}.*.gene_abfs.tsv \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
      ${rCNV_bucket}/analysis/gene_scoring/data/
    gsutil -m cp \
      *.pdf \
      ${rCNV_bucket}/analysis/gene_scoring/plots/
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    File priors_tsv = "${freq_code}.prior_estimation.empirical_prior_estimates.tsv"
    Float theta0_del = read_float("theta0_del.tsv")
    Float theta0_dup = read_float("theta0_dup.tsv")
    Float theta1_del = read_float("theta1_del.tsv")
    Float theta1_dup = read_float("theta1_dup.tsv")
    Float var0 = 1.0
    File del_bfdp = "${freq_code}.DEL.gene_abfs.tsv"
    File dup_bfdp = "${freq_code}.DUP.gene_abfs.tsv"
    Array[File] prior_plots = glob("${freq_code}.prior_estimation*pdf")
  }
}


# Score all genes for one model & CNV type
task score_genes {
  String CNV
  File BFDP_stats
  File excludelist
  File underpowered_genes
  File gene_features
  Float max_true_bfdp
  Float min_false_bfdp
  String model
  String freq_code
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -e

    # Copy necessary references
    gsutil -m cp \
      ${rCNV_bucket}/refs/GRCh37.centromeres_telomeres.bed.gz \
      ${rCNV_bucket}/analysis/gene_scoring/gene_lists/*genes.list \
      ./

    # Merge excludelist and list of underpowered genes
    zcat ${excludelist} ${underpowered_genes} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
    | uniq \
    | cat <( zcat ${excludelist} | grep -e '^#' ) - \
    | bgzip -c  \
    > excludelist_plus_underpowered.bed.gz

    # Set CNV type-specific parameters
    case ${CNV} in
      "DEL")
        true_pos="gold_standard.haploinsufficient.genes.list"
        true_neg="gold_standard.haplosufficient.genes.list"
        ;;
      "DUP")
        true_pos="gold_standard.triplosensitive.genes.list"
        true_neg="gold_standard.triploinsensitive.genes.list"
        ;;
    esac

    # Score genes
    /opt/rCNV2/analysis/gene_scoring/score_genes.py \
      --centromeres GRCh37.centromeres_telomeres.bed.gz \
      --blacklist excludelist_plus_underpowered.bed.gz \
      --true-positives $true_pos \
      --true-negatives $true_neg \
      --model ${model} \
      --max-true-bfdp ${max_true_bfdp} \
      --min-false-bfdp ${min_false_bfdp} \
      --chromsplit \
      --no-out-of-sample-prediction \
      --outfile ${freq_code}.${CNV}.gene_scores.${model}.tsv \
      ${BFDP_stats} \
      ${gene_features}

    # Copy scores to Google bucket
    gsutil -m cp \
      ${freq_code}.${CNV}.gene_scores.${model}.tsv \
      ${rCNV_bucket}/analysis/gene_scoring/all_models/
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "16 GB"
    cpu: "8"
    bootDiskSizeGb: "20"
  }

  output {
    File scores_tsv = "${freq_code}.${CNV}.gene_scores.${model}.tsv"
  }
}


# Score genes with ensemble classifier of all individual models
task score_ensemble {
  String CNV
  Array[File] scores
  String freq_code
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -e

    # Localize scores to working directory
    find / -name "${freq_code}.${CNV}.gene_scores.*.tsv" | xargs -I {} mv {} ./

    # Copy gene lists
    gsutil -m cp \
      ${rCNV_bucket}/analysis/gene_scoring/gene_lists/gold_standard.*.genes.list \
      ./

    # Make inputs for ensemble classifier
    find ./ -name "${freq_code}.${CNV}.gene_scores.*.tsv" \
    > ensemble_input.tsv
    case ${CNV} in
      DEL)
        pos_genes="gold_standard.haploinsufficient.genes.list"
        neg_genes="gold_standard.haplosufficient.genes.list"
        ;;
      DUP)
        pos_genes="gold_standard.triplosensitive.genes.list"
        neg_genes="gold_standard.triploinsensitive.genes.list"
        ;;
    esac

    # Run ensemble classifier
    /opt/rCNV2/analysis/gene_scoring/ensemble_classifier.R \
      ensemble_input.tsv \
      $pos_genes \
      $neg_genes \
      ${freq_code}.${CNV}.gene_scores.ensemble.tsv

    # Copy scores to Google bucket
    gsutil -m cp \
      ${freq_code}.${CNV}.gene_scores.ensemble.tsv \
      ${rCNV_bucket}/analysis/gene_scoring/all_models/
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] all_scores = glob("${freq_code}.${CNV}.gene_scores.*.tsv")
  }
}


# Compare scores between models, determine best model, and QC final set of scores
task qc_scores {
  Array[File] del_scores
  Array[File] dup_scores
  Array[String] models
  File raw_gene_features
  String freq_code
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -e

    # Localize scores to working directory
    find / -name "*.gene_scores.*.tsv" | xargs -I {} mv {} ./

    # Copy gene lists
    gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/gene_lists ./
    gsutil -m cp \
      ${rCNV_bucket}/analysis/gene_scoring/gene_lists/gold_standard.*.genes.list \
      gene_lists/

    # Evaluate every model to determine best overall predictor
    compdir=${freq_code}_gene_scoring_model_comparisons
    if ! [ -e $compdir ]; then
      mkdir $compdir
    fi
    for CNV in DEL DUP; do
      for wrapper in 1; do
        for model in ${sep=" " models} ensemble; do
          echo $model
          echo ${freq_code}.$CNV.gene_scores.$model.tsv
        done | paste - -
      done > ${freq_code}.$CNV.model_evaluation.input.tsv
    done
    /opt/rCNV2/analysis/gene_scoring/compare_models.R \
      ${freq_code}.DEL.model_evaluation.input.tsv \
      ${freq_code}.DUP.model_evaluation.input.tsv \
      gene_lists/gold_standard.haploinsufficient.genes.list \
      gene_lists/gold_standard.triplosensitive.genes.list \
      gene_lists/gold_standard.haplosufficient.genes.list \
      gene_lists/gold_standard.triploinsensitive.genes.list \
      $compdir/${freq_code}_gene_scoring_model_comparisons

    # Compute harmonic mean of AUCs and merge scores from ensemble model
    sed -n '2p' $compdir/${freq_code}_gene_scoring_model_comparisons.summary_table.tsv \
    | cut -f1 > best_average_auc.tsv
    best_model="ensemble"
    /opt/rCNV2/analysis/gene_scoring/merge_del_dup_scores.R \
      ${freq_code}.DEL.gene_scores.$best_model.tsv \
      ${freq_code}.DUP.gene_scores.$best_model.tsv \
      ${freq_code}.gene_scores.tsv
    gzip -f ${freq_code}.gene_scores.tsv

    # Plot correlations of raw features vs scores
    mkdir gene_score_corplots
    /opt/rCNV2/analysis/gene_scoring/plot_score_feature_cors.R \
      ${freq_code}.gene_scores.tsv.gz \
      ${raw_gene_features} \
      gene_score_corplots/${freq_code}.gene_scores.raw_feature_cors

    # Make CNV type-dependent truth sets
    for CNV in DEL DUP; do
      case $CNV in
        "DEL")
          cat gene_lists/ClinGen.hc_haploinsufficient.genes.list \
              gene_lists/DDG2P.hc_lof.genes.list \
              gene_lists/cell_essential.genes.list \
              gene_lists/mouse_het_lethal.genes.list \
          | sort -Vk1,1 | uniq \
          > master_union_truth_set.lof.tsv

          # Write truth set input tsv
          for wrapper in 1; do
            echo -e "Union truth set\tmaster_union_truth_set.lof.tsv\tgrey25"
            echo -e "ClinGen dom. HI\tgene_lists/ClinGen.hc_haploinsufficient.genes.list\t#9F2B1C"
            echo -e "DECIPHER dom. LoF\tgene_lists/DDG2P.hc_lof.genes.list\t#D43925"
            echo -e "Cell essential\tgene_lists/cell_essential.genes.list\t#DD6151"
            echo -e "Mouse het. lethal\tgene_lists/mouse_het_lethal.genes.list\t#E5887C"
          done > DEL.roc_truth_sets.tsv
          ;;

        "DUP")
          cat gene_lists/ClinGen.all_triplosensitive.genes.list \
              gene_lists/DDG2P.hc_gof.genes.list \
              gene_lists/COSMIC.hc_oncogenes.genes.list \
          | sort -Vk1,1 | uniq \
          > master_union_truth_set.gof.tsv

          # Write truth set input tsv
          for wrapper in 1; do
            echo -e "Union truth set\tmaster_union_truth_set.gof.tsv\tgrey25"
            echo -e "ClinGen dom. TS\tgene_lists/ClinGen.all_triplosensitive.genes.list\t#1A5985"
            echo -e "DECIPHER dom. GoF\tgene_lists/DDG2P.hc_gof.genes.list\t#2376B2"
            echo -e "COSMIC dom. oncogenes\tgene_lists/COSMIC.hc_oncogenes.genes.list\t#4F91C1"
          done > DUP.roc_truth_sets.tsv
          ;;
      esac    
    done

    # Plot ROC, PRC, and enrichments
    mkdir ${freq_code}_gene_scoring_QC_plots/
    /opt/rCNV2/analysis/gene_scoring/plot_gene_score_qc.R \
      ${freq_code}.gene_scores.tsv.gz \
      DEL.roc_truth_sets.tsv \
      DUP.roc_truth_sets.tsv \
      gene_lists/gold_standard.haplosufficient.genes.list \
      gene_lists/gold_standard.triploinsensitive.genes.list \
      ${freq_code}_gene_scoring_QC_plots/${freq_code}_gene_score_qc

    # Copy all results to Google bucket
    gsutil -m cp \
      ${freq_code}_gene_scoring_model_comparisons/* \
      ${rCNV_bucket}/analysis/gene_scoring/model_comparisons/
    gsutil -m cp \
      gene_score_corplots/* \
      ${rCNV_bucket}/analysis/gene_scoring/plots/
    gsutil -m cp \
      ${freq_code}_gene_scoring_QC_plots/* \
      ${rCNV_bucket}/analysis/gene_scoring/plots/
    gsutil -m cp \
      ${freq_code}.gene_scores.tsv.gz \
      ${rCNV_bucket}/results/gene_scoring/
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    String best_average_auc = read_string("best_average_auc.tsv")
    File final_scores = "${freq_code}.gene_scores.tsv.gz"
  }
}
