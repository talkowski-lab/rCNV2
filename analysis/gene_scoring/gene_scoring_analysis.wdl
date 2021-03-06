######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Gene scoring for haploinsufficiency and triplosensitivity


workflow gene_burden_analysis {
  File training_hpo_list
  File effective_case_sample_sizes
  String prefix
  File metacohort_list
  File metacohort_sample_table
  File gtf
  File gtf_idx
  Int pad_controls
  Int max_cnv_size
  String weight_mode
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  String meta_model_prefix
  Int min_cnvs_per_gene_training
  Float prior_frac
  Float prior_lnor_thresholding_pct
  File training_blacklist
  File gene_features
  File raw_gene_features
  Float max_true_bfdp
  Float min_false_bfdp
  Float elnet_alpha
  Float elnet_l1_l2_mix
  String rCNV_bucket
  File contiglist

  Array[Array[String]] contigs = read_tsv(contiglist)

  Array[String] cnv_types = ["DEL", "DUP"]

  Array[String] models = ["logit", "svm", "randomforest", "lda", "naivebayes", "sgd", "neuralnet"]

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
        weight_mode=weight_mode,
        min_cds_ovr_del=min_cds_ovr_del,
        min_cds_ovr_dup=min_cds_ovr_dup,
        rCNV_bucket=rCNV_bucket
    }
    # Run assocation tests per cohort
    call burden_test as rCNV_burden_test {
      input:
        training_hpo_list=training_hpo_list,
        effective_case_sample_sizes=effective_case_sample_sizes,
        prefix=prefix,
        contig=contig[0],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        gtf=gtf,
        gtf_idx=gtf_idx,
        pad_controls=pad_controls,
        max_cnv_size=max_cnv_size,
        weight_mode=weight_mode,
        min_cds_ovr_del=min_cds_ovr_del,
        min_cds_ovr_dup=min_cds_ovr_dup,
        max_genes_per_cnv=max_genes_per_cnv,
        rCNV_bucket=rCNV_bucket
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
        rCNV_bucket=rCNV_bucket
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
        freq_code="rCNV",
        rCNV_bucket=rCNV_bucket
    }
    call get_underpowered_genes {
      input:
        stats_beds=rCNV_burden_test.stats_beds,
        stats_beds_idxs=rCNV_burden_test.stats_bed_idxs,
        prefix=prefix,
        CNV=CNV,
        metacohort_list=metacohort_list,
        min_cnvs_per_gene_training=min_cnvs_per_gene_training,
        freq_code="rCNV",
        rCNV_bucket=rCNV_bucket
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
      prior_frac=prior_frac,
      prior_lnor_thresholding_pct=prior_lnor_thresholding_pct
  }

  # Score genes for each model
  scatter ( model in models ) {
    call score_genes as score_genes_DEL {
      input:
        CNV="DEL",
        BFDP_stats=calc_priors_bfdp.del_bfdp,
        blacklist=training_blacklist,
        underpowered_genes=get_underpowered_genes.underpowered_genes[0],
        gene_features=gene_features,
        model=model,
        max_true_bfdp=max_true_bfdp,
        min_false_bfdp=min_false_bfdp,
        elnet_alpha=elnet_alpha,
        elnet_l1_l2_mix=elnet_l1_l2_mix,
        freq_code="rCNV",
        rCNV_bucket=rCNV_bucket
    }
    call score_genes as score_genes_DUP {
      input:
        CNV="DUP",
        BFDP_stats=calc_priors_bfdp.dup_bfdp,
        blacklist=training_blacklist,
        underpowered_genes=get_underpowered_genes.underpowered_genes[1],
        gene_features=gene_features,
        model=model,
        max_true_bfdp=max_true_bfdp,
        min_false_bfdp=min_false_bfdp,
        elnet_alpha=elnet_alpha,
        elnet_l1_l2_mix=elnet_l1_l2_mix,
        freq_code="rCNV",
        rCNV_bucket=rCNV_bucket
    }
  }

  # Score with ensemble classifier, and return updated array of all scores + ensemble
  call score_ensemble as score_ensemble_DEL {
    input:
      CNV="DEL",
      scores=score_genes_DEL.scores_tsv,
      freq_code="rCNV",
      rCNV_bucket=rCNV_bucket
  }
  call score_ensemble as score_ensemble_DUP {
    input:
      CNV="DUP",
      scores=score_genes_DUP.scores_tsv,
      freq_code="rCNV",
      rCNV_bucket=rCNV_bucket
  }

  # Determine best model & QC final scores
  call qc_scores {
    input:
      del_scores=score_ensemble_DEL.all_scores,
      dup_scores=score_ensemble_DUP.all_scores,
      models=models,
      raw_gene_features=raw_gene_features,
      freq_code="rCNV",
      rCNV_bucket=rCNV_bucket
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
  String weight_mode
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  String rCNV_bucket

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
          --weight-mode ${weight_mode} \
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
    docker: "talkowski/rcnv@sha256:b88c6669695764493aa778498afc2df9407fe35d3e45aa127922c9701fe8f64c"
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
  File training_hpo_list
  File effective_case_sample_sizes
  String prefix
  String contig
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  File gtf
  File gtf_idx
  Int pad_controls
  Int max_cnv_size
  String weight_mode
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  String rCNV_bucket

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
    tabix ${gtf} ${contig} | bgzip -c > ${contig}.gtf.gz

    # Create string of HPOs to keep
    keep_hpos=$( cat ${training_hpo_list} | paste -s -d\; )

    # Iterate over metacohorts to compute single-cohort stats
    while read meta cohorts; do
      echo $meta

      # Set metacohort-specific parameters
      cnv_bed="cleaned_cnv/$meta.${freq_code}.bed.gz"
      effective_case_n=$( fgrep -w $meta ${effective_case_sample_sizes} | cut -f2 )

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
          --weight-mode ${weight_mode} \
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
        /opt/rCNV2/analysis/genes/gene_burden_test.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --cnv $CNV \
          --effective-case-n $effective_case_n \
          --bgzip \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.${contig}.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.${contig}.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.${contig}.bed.gz"
      done
    done < <( fgrep -v "mega" ${metacohort_list} )
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:b88c6669695764493aa778498afc2df9407fe35d3e45aa127922c9701fe8f64c"
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
    docker: "talkowski/rcnv@sha256:b88c6669695764493aa778498afc2df9407fe35d3e45aa127922c9701fe8f64c"
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
  String freq_code
  String rCNV_bucket

  command <<<
    set -e

    # Make list of stats files to be considered
    find / -name "*.${prefix}.${freq_code}.${CNV}.gene_burden.stats.*.bed.gz" \
    > stats.paths.list
    # Debug: print paths to stdout
    cat stats.paths.list

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
    done < <( fgrep -v "mega" ${metacohort_list} )

    # Make input for meta-analysis
    while read meta cohorts; do
      echo -e "$meta\t$meta.${prefix}.${freq_code}.${CNV}.gene_burden.stats.bed.gz"
    done < <( fgrep -v "mega" ${metacohort_list} ) \
    > ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt
    
    # Run meta-analysis
    /opt/rCNV2/analysis/genes/gene_meta_analysis.R \
      --model ${meta_model_prefix} \
      --p-is-phred \
      --spa \
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
    docker: "talkowski/rcnv@sha256:b88c6669695764493aa778498afc2df9407fe35d3e45aa127922c9701fe8f64c"
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
  Array[Array[File]] stats_beds
  Array[Array[File]] stats_beds_idxs
  String prefix
  String CNV
  File metacohort_list
  Int min_cnvs_per_gene_training
  String freq_code
  String rCNV_bucket

  command <<<
    set -e

    # Make list of stats files to be considered
    find / -name "*.${prefix}.${freq_code}.${CNV}.gene_burden.stats.*.bed.gz" \
    > stats.paths.list
    # Debug: print paths to stdout
    cat stats.paths.list

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
    done < <( fgrep -v "mega" ${metacohort_list} )

    # Make input for meta-analysis
    while read meta cohorts; do
      echo -e "$meta\t$meta.${prefix}.${freq_code}.${CNV}.gene_burden.stats.bed.gz"
    done < <( fgrep -v "mega" ${metacohort_list} ) \
    > ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt

    # Extract list of genes with < min_cnvs_per_gene_training (to be used as blacklist later for training)
    # Note: this cutoff must be pre-determined with eval_or_vs_cnv_counts.R, but is not automated here
    /opt/rCNV2/analysis/gene_scoring/get_underpowered_genes.R \
      --min-cnvs ${min_cnvs_per_gene_training} \
      --gene-counts-out ${prefix}.${freq_code}.${CNV}.counts_per_gene.tsv \
      ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt \
      ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed
    awk -v OFS="\t" '{ print $1, $2, $3, $4 }' \
      ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed \
    | bgzip -c \
    > ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz
    tabix -p bed -f ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz

    # Copy meta-analysis results to Google bucket
    gsutil -m cp \
      ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz \
      ${prefix}.${freq_code}.${CNV}.counts_per_gene.tsv \
      ${rCNV_bucket}/analysis/gene_scoring/data/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:b27de66b70fee3590dbfe965e22456082b9ec735ea404720eeb25958eee9155e"
    preemptible: 1
    memory: "8 GB"
    bootDiskSizeGb: "30"
    disks: "local-disk 50 HDD"
  }

  output {
    File underpowered_genes = "${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz"
    File counts_per_gene = "${prefix}.${freq_code}.${CNV}.counts_per_gene.tsv"
  }
}


# Estimate prior effect sizes and compute BFDPs per gene
task calc_priors_bfdp {
  File del_meta_stats
  File dup_meta_stats
  Array[File] underpowered_genes
  String freq_code
  String rCNV_bucket
  Float prior_frac
  Float prior_lnor_thresholding_pct

  command <<<
    set -e

    # Localize necessary references
    gsutil -m cp \
      ${rCNV_bucket}/analysis/gene_scoring/refs/${freq_code}.gene_scoring.training_gene_blacklist.bed.gz \
      ${rCNV_bucket}/analysis/gene_scoring/gene_lists/*.genes.list \
      ./
    mkdir gene_lists/
    gsutil -m cp \
      gs://rcnv_project/cleaned_data/genes/gene_lists/*genes.list \
      gene_lists/

    # Locate underpowered genes
    find / -name "*.gene_burden.underpowered_genes.bed.gz" | xargs -I {} mv {} ./

    # Create CNV-type-specific gene blacklists
    for CNV in DEL DUP; do
      zcat *.$CNV.gene_burden.underpowered_genes.bed.gz \
        ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz \
      | fgrep -v "#" | cut -f4 | sort -V | uniq \
      > ${freq_code}.$CNV.training_blacklist.genes.list
    done

    # Compute prior effect sizes
    /opt/rCNV2/analysis/gene_scoring/estimate_prior_effect_sizes.R \
      --pct ${prior_lnor_thresholding_pct} \
      ${del_meta_stats} \
      ${dup_meta_stats} \
      ${freq_code}.DEL.training_blacklist.genes.list \
      ${freq_code}.DUP.training_blacklist.genes.list \
      gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
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
    awk -v FS="\t" '{ if ($1=="var0" && $2=="DEL") print $3 }' \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
    > var0_del.tsv
    awk -v FS="\t" '{ if ($1=="var0" && $2=="DUP") print $3 }' \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
    > var0_dup.tsv
    awk -v FS="\t" '{ if ($1=="var1" && $2=="DEL") print $3 }' \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
    > var1_del.tsv
    awk -v FS="\t" '{ if ($1=="var1" && $2=="DUP") print $3 }' \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
    > var1_dup.tsv

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
        --blacklist ${freq_code}.$CNV.training_blacklist.genes.list \
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
    docker: "talkowski/rcnv@sha256:b27de66b70fee3590dbfe965e22456082b9ec735ea404720eeb25958eee9155e"
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
  File blacklist
  File underpowered_genes
  File gene_features
  Float max_true_bfdp
  Float min_false_bfdp
  String model
  Float elnet_alpha
  Float elnet_l1_l2_mix
  String freq_code
  String rCNV_bucket

  command <<<
    set -e

    # Copy necessary references
    gsutil -m cp \
      ${rCNV_bucket}/refs/GRCh37.centromeres_telomeres.bed.gz \
      ${rCNV_bucket}/analysis/gene_scoring/gene_lists/*genes.list \
      ./

    # Merge blacklist and list of underpowered genes
    zcat ${blacklist} ${underpowered_genes} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
    | uniq \
    | cat <( zcat ${blacklist} | grep -e '^#' ) - \
    | bgzip -c  \
    > blacklist_plus_underpowered.bed.gz

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
      --blacklist blacklist_plus_underpowered.bed.gz \
      --true-positives $true_pos \
      --true-negatives $true_neg \
      --model ${model} \
      --max-true-bfdp ${max_true_bfdp} \
      --min-false-bfdp ${min_false_bfdp} \
      --regularization-alpha ${elnet_alpha} \
      --regularization-l1-l2-mix ${elnet_l1_l2_mix} \
      --outfile ${freq_code}.${CNV}.gene_scores.${model}.tsv \
      ${BFDP_stats} \
      ${gene_features}

    # Copy scores to Google bucket
    gsutil -m cp \
      ${freq_code}.${CNV}.gene_scores.${model}.tsv \
      ${rCNV_bucket}/analysis/gene_scoring/all_models/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:4b494ff5dde3eb5a453863f425b4a6d5994bc56b68448332979d26ce92c8ee50"
    preemptible: 1
    memory: "8 GB"
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
    docker: "talkowski/rcnv@sha256:4b494ff5dde3eb5a453863f425b4a6d5994bc56b68448332979d26ce92c8ee50"
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

    # Merge scores from best model (highest harmonic mean AUC)
    sed -n '2p' $compdir/${freq_code}_gene_scoring_model_comparisons.summary_table.tsv \
    | cut -f1 > best_model.tsv
    best_model=$( cat best_model.tsv )
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
    docker: "talkowski/rcnv@sha256:2fdd11e54719ab4c2ec1324e3785e7806c66b6e1cbf894588c286e72c45df75c"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    String best_model = read_string("best_model.tsv")
    File final_scores = "${freq_code}.gene_scores.tsv.gz"
  }
}