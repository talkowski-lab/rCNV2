######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Gene scoring for haploinsufficiency and triplosensitivity


workflow gene_burden_analysis {
  String hpo
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
  Int cen_tel_dist
  Float prior_frac
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
    # Run assocation tests per cohort
    call burden_test as rCNV_burden_test {
      input:
        hpo=hpo,
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

  # Merge association statistics and run meta-analysis
  scatter ( CNV in cnv_types ) {
    call merge_and_meta_analysis as rCNV_meta_analysis {
      input:
        stats_beds=rCNV_burden_test.stats_beds,
        stats_beds_idxs=rCNV_burden_test.stats_bed_idxs,
        hpo=hpo,
        prefix=prefix,
        CNV=CNV,
        metacohort_list=metacohort_list,
        meta_model_prefix=meta_model_prefix,
        min_cnvs_per_gene_training=min_cnvs_per_gene_training,
        freq_code="rCNV",
        rCNV_bucket=rCNV_bucket
    }
  }

  # Make training blacklist, estimate prior effect sizes, and compute BFDPs
  call blacklist_priors_bfdp {
    input:
      del_meta_stats=rCNV_meta_analysis.meta_stats_bed[0],
      dup_meta_stats=rCNV_meta_analysis.meta_stats_bed[1],
      freq_code="rCNV",
      rCNV_bucket=rCNV_bucket,
      cen_tel_dist=cen_tel_dist,
      prior_frac=prior_frac
  }

  # Score genes for each model
  scatter ( model in models ) {
    call score_genes as score_genes_DEL {
      input:
        CNV="DEL",
        BFDP_stats=blacklist_priors_bfdp.del_bfdp,
        blacklist=blacklist_priors_bfdp.training_blacklist,
        underpowered_genes=rCNV_meta_analysis.underpowered_genes[0],
        gene_features=gene_features,
        model=model,
        max_true_bfdpmax_true_bfdp,
        min_false_bfdp=min_false_bfdp,
        elnet_alpha=elnet_alpha,
        elnet_l1_l2_mix=elnet_l1_l2_mix,
        freq_code="rCNV",
        rCNV_bucket=rCNV_bucket
    }
    call score_genes as score_genes_DUP {
      input:
        CNV="DUP",
        BFDP_stats=blacklist_priors_bfdp.dup_bfdp,
        blacklist=blacklist_priors_bfdp.training_blacklist,
        underpowered_genes=rCNV_meta_analysis.underpowered_genes[1],
        gene_features=gene_features,
        model=model,
        max_true_bfdpmax_true_bfdp,
        min_false_bfdp=min_false_bfdp,
        elnet_alpha=elnet_alpha,
        elnet_l1_l2_mix=elnet_l1_l2_mix,
        freq_code="rCNV",
        rCNV_bucket=rCNV_bucket
    }
  }

  # Determine best model & QC final scores
  call qc_scores {
    input:
      del_scores=score_genes_DEL.scores_tsv,
      dup_scores=score_genes_DUP.scores_tsv,
      raw_gene_features=raw_gene_features,
      freq_code="rCNV",
      rCNV_bucket=rCNV_bucket
  }

  output {
    File gene_scores = qc_scores.final_scores
  }
}


# Run burden test for a single chromosome for all cohorts
task burden_test {
  String hpo
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
    tabix ${gtf} ${contig} | bgzip -c > subset.gtf.gz

    # Iterate over metacohorts to compute single-cohort stats
    while read meta cohorts; do
      echo $meta

      # Set metacohort-specific parameters
      cnv_bed="cleaned_cnv/$meta.${freq_code}.bed.gz"

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
          --hpo ${hpo} \
          --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
          --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
          --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
          -z \
          --verbose \
          -o "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.${contig}.bed.gz" \
          "$cnv_bed" \
          subset.gtf.gz
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.${contig}.bed.gz"

        # Perform burden test
        /opt/rCNV2/analysis/genes/gene_burden_test.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --cnv $CNV \
          --case-hpo ${hpo} \
          --bgzip \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.${contig}.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.${contig}.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.${contig}.bed.gz"
      done
    done < <( fgrep -v "mega" ${metacohort_list} )
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:da2df0f4afcfa93d17e27ae5752f1665f0c53c8feed06d6d98d1da53144d8e1f"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] stats_beds = glob("*.gene_burden.stats.${contig}.bed.gz")
    Array[File] stats_bed_idxs = glob("*.gene_burden.stats.${contig}.bed.gz.tbi")
  }
}


# Merge burden tests and perform meta-analysis per cohort per CNV type
task merge_and_meta_analysis {
  Array[Array[File]] stats_beds
  Array[Array[File]] stats_beds_idxs
  String hpo
  String prefix
  String CNV
  File metacohort_list
  String meta_model_prefix
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
    done < ${metacohort_list}

    # Make input for meta-analysis
    while read meta cohorts; do
      echo -e "$meta\t$meta.${prefix}.${freq_code}.${CNV}.gene_burden.stats.bed.gz"
    done < <( fgrep -v mega ${metacohort_list} ) \
    > ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt
    
    # Run meta-analysis
    /opt/rCNV2/analysis/genes/gene_meta_analysis.R \
      --model ${meta_model_prefix} \
      --p-is-phred \
      --spa \
      ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt \
      ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed
    bgzip -f ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed
    tabix -f ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz

    # Extract list of genes with < min_cnvs_per_gene_training (to be used as blacklist later for training)
    /opt/rCNV2/analysis/gene_scoring/get_underpowered_genes.R \
      --min-cnvs ${min_cnvs_per_gene_training} \
      ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt \
      ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed
    bgzip -f ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed
    tabix -f ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz

    # Copy meta-analysis results to Google bucket
    gsutil -m cp \
      ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz* \
      ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.list \
      ${rCNV_bucket}/analysis/gene_scoring/data/
  >>>

  runtime {
    # TODO: update docker
    # docker: "talkowski/rcnv@sha256:da2df0f4afcfa93d17e27ae5752f1665f0c53c8feed06d6d98d1da53144d8e1f"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    File meta_stats_bed = "${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz"
    File meta_stats_bed_idx = "${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz.tbi"
    File underpowered_genes = "${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz"
  }
}


# Determine genes to be blacklisted during model training, estimate prior effect 
# sizes, and compute BFDPs per gene
task blacklist_priors_bfdp {
  File del_meta_stats
  File dup_meta_stats
  String freq_code
  String rCNV_bucket
  Int cen_tel_dist
  Float prior_frac

  command <<<
    set -e

    # Create blacklist: Remove all genes within ±1Mb of a telomere/centromere, 
    # or those within known genomic disorder regions
    gsutil -m cp ${rCNV_bucket}/refs/GRCh37.centromeres_telomeres.bed.gz ./
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs.*.bed.gz refs/
    zcat GRCh37.centromeres_telomeres.bed.gz \
    | awk -v FS="\t" -v OFS="\t" -v d=${cen_tel_dist} \
      '{ if ($2-d<0) print $1, "0", $3+d; else print $1, $2-d, $3+d }' \
    | cat - <( zcat refs/lit_GDs.*.bed.gz | cut -f1-3 | fgrep -v "#" ) \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    > ${freq_code}.gene_scoring.training_mask.bed.gz
    bedtools intersect -u \
      -a ${del_meta_stats} \
      -b ${freq_code}.gene_scoring.training_mask.bed.gz \
    | cut -f1-4 \
    | cat <( echo -e "#chr\tstart\tend\tgene" ) - \
    | bgzip -c \
    > ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz

    # Copy gene lists
    gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/gene_lists ./

    # Define list of high-confidence dosage-sensitive genes
    fgrep -wf gene_lists/DDG2P.hc_lof.genes.list \
      gene_lists/ClinGen.hc_haploinsufficient.genes.list \
    > gold_standard.ad_disease.genes.list
    fgrep -wf gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
      gold_standard.ad_disease.genes.list \
    > gold_standard.haploinsufficient.genes.list

    # Define list of high-confidence dosage-insensitive genes
    cat gene_lists/HP0000118.HPOdb.genes.list \
      gene_lists/DDG2P*.genes.list \
      gene_lists/ClinGen*.genes.list \
    | fgrep -wvf - gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
    > gold_standard.no_disease_assoc.genes.list
    fgrep -wf gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list \
      gold_standard.no_disease_assoc.genes.list \
    > gold_standard.haplosufficient.genes.list

    # Compute prior effect sizes
    /opt/rCNV2/analysis/gene_scoring/estimate_prior_effect_sizes.R \
      ${del_meta_stats} \
      ${dup_meta_stats} \
      ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz \
      gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
      gold_standard.haploinsufficient.genes.list \
      gold_standard.haplosufficient.genes.list \
      ${freq_code}.prior_estimation
    awk -v FS="\t" '{ if ($1=="theta0" && $2=="DEL") print $3 }' \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
    > theta0_del.tsv
    awk -v FS="\t" '{ if ($1=="theta0" && $2=="DUP") print $3 }' \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
    > theta0_dup.tsv
    awk -v FS="\t" '{ if ($1=="theta1" && $2=="DEL") print $3 }' \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
    > theta1.tsv
    awk -v FS="\t" '{ if ($1=="var1" && $2=="DEL") print $3 }' \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
    > var1.tsv

    # Compute BFDP per gene
    for CNV in DEL DUP; do
      # Set CNV-specific variables
      case $CNV in
        "DEL")
          theta0=$( cat theta0_del.tsv )
          statsbed=${del_meta_stats}
          ;;
        "DUP")
          theta0=$( cat theta0_dup.tsv )
          statsbed=${dup_meta_stats}
          ;;
      esac

      # Compute BF & BFDR for all genes
      /opt/rCNV2/analysis/gene_scoring/calc_gene_bfs.py \
        --theta0 $theta0 \
        --theta1 $( cat theta1.tsv ) \
        --var0 $( cat var1.tsv ) \
        --prior ${prior_frac} \
        --blacklist ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz \
        --outfile ${freq_code}.$CNV.gene_abfs.tsv \
        "$statsbed"
    done

    # Copy results to Google bucket
    gsutil -m cp \
      ${freq_code}.*.gene_abfs.tsv \
      ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz \
      ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
      ${rCNV_bucket}/analysis/gene_scoring/data/
    gsutil -m cp \
      gold_standard.*.genes.list \
      ${rCNV_bucket}/analysis/gene_scoring/gene_lists/
    gsutil -m cp \
      *.pdf \
      ${rCNV_bucket}/analysis/gene_scoring/plots/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:0489fe64e3f34b3751a15930c2d85b5b7291fe5cdc0db98f3ba001d06f8f6617"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    File training_blacklist = "${freq_code}.gene_scoring.training_gene_blacklist.bed.gz"
    File priors_tsv = "${freq_code}.prior_estimation.empirical_prior_estimates.tsv"
    Float theta0_del = read_float("theta0_del.tsv")
    Float theta0_dup = read_float("theta0_dup.tsv")
    Float theta1 = read_float("theta1.tsv")
    Float var0 = read_float("var1.tsv")
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

    # Copy centromeres bed
    gsutil -m cp ${rCNV_bucket}/refs/GRCh37.centromeres_telomeres.bed.gz ./

    # Merge blacklist and list of underpowered genes
    zcat ${blacklist} ${underpowered_genes} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
    | uniq \
    | cat <( zcat ${blacklist} | grep -e '^#' ) - \
    | bgzip -c  \
    > blacklist_plus_underpowered.bed.gz

    # Score genes
    /opt/rCNV2/analysis/gene_scoring/score_genes.py \
      --centromeres GRCh37.centromeres_telomeres.bed.gz \
      --blacklist blacklist_plus_underpowered.bed.gz \
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
    docker: "talkowski/rcnv@sha256:ab9eeb96ddc5a72af3c3c67d2ea82bd3410dfe6e4ece81b7660868b1c546846d"
    preemptible: 1
    memory: "8 GB"
    bootDiskSizeGb: "20"
  }

  output {
    File scores_tsv = "${freq_code}.${CNV}.gene_scores.${model}.tsv"
  }
}


# Compare scores between models, determine best model, and QC final set of scores
task qc_scores {
  Array[File] del_scores
  Array[File] dup_scores
  File raw_gene_features
  String freq_code
  String rCNV_bucket

  command <<<
    set -e

    # Localize scores to working directory
    find / -name "*.gene_scores.*.tsv" | xargs -I {} mv {} ./

    # Copy gene lists
    gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/gene_lists ./
    gsutil -m cp ${rCNV_bucket}/analysis/gene_scoring/gene_lists/* ./gene_lists/

    # Evaluate every model to determine best overall predictor
    compdir=${freq_code}_gene_scoring_model_comparisons
    if ! [ -e $compdir ]; then
      mkdir $compdir
    fi
    for CNV in DEL DUP; do
      for wrapper in 1; do
        for model in logit svm randomforest lda naivebayes sgd neuralnet; do
          echo $model
          echo ${freq_code}.$CNV.gene_scores.$model.tsv
        done | paste - -
      done > ${freq_code}.$CNV.model_evaluation.input.tsv
    done
    /opt/rCNV2/analysis/gene_scoring/compare_models.R \
      ${freq_code}.DEL.model_evaluation.input.tsv \
      ${freq_code}.DUP.model_evaluation.input.tsv \
      gene_lists/gold_standard.haploinsufficient.genes.list \
      gene_lists/gold_standard.haplosufficient.genes.list \
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
          # Union (ClinGen HI + DDG2P dominant LoF)
          cat gene_lists/ClinGen.hmc_haploinsufficient.genes.list \
              gene_lists/DDG2P.hmc_lof.genes.list \
          | sort -Vk1,1 | uniq \
          > union_truth_set.lof.tsv

          # Intersection (ClinGen HI + DDG2P dominant LoF)
          fgrep -wf \
            gene_lists/ClinGen.hmc_haploinsufficient.genes.list \
            gene_lists/DDG2P.hmc_lof.genes.list \
          | sort -Vk1,1 | uniq \
          > intersection_truth_set.lof.tsv

          # ClinGen HI alone
          cat gene_lists/ClinGen.all_haploinsufficient.genes.list \
          | sort -Vk1,1 | uniq \
          > clingen_truth_set.lof.tsv

          # DDG2P lof + other alone
          cat gene_lists/DDG2P.all_lof.genes.list \
              gene_lists/DDG2P.all_other.genes.list \
          | sort -Vk1,1 | uniq \
          > ddg2p_truth_set.lof.tsv

          # Combine all truth sets into master union
          cat union_truth_set.lof.tsv \
              clingen_truth_set.lof.tsv \
              ddg2p_truth_set.lof.tsv \
          | sort -Vk1,1 | uniq \
          > master_union_truth_set.lof.tsv

          # Write truth set input tsv
          for wrapper in 1; do
            echo -e "Union truth set\tmaster_union_truth_set.lof.tsv\tgrey25"
            echo -e "ClinGen & DECIPHER (union)\tunion_truth_set.lof.tsv\t#9F2B1C"
            echo -e "ClinGen & DECIPHER (int.)\tintersection_truth_set.lof.tsv\t#D43925"
            echo -e "ClinGen dom. HI\tclingen_truth_set.lof.tsv\t#DD6151"
            echo -e "DECIPHER dom. LoF/unk.\tddg2p_truth_set.lof.tsv\t#E5887C"
          done > DEL.roc_truth_sets.tsv
          ;;

        "DUP")
          # Union (ClinGen HI + DDG2P dominant CG)
          cat gene_lists/ClinGen.hmc_triplosensitive.genes.list \
              gene_lists/DDG2P.hmc_gof.genes.list \
          | sort -Vk1,1 | uniq \
          > union_truth_set.gof.tsv

          # ClinGen triplo alone
          cat gene_lists/ClinGen.all_triplosensitive.genes.list \
          | sort -Vk1,1 | uniq \
          > clingen_truth_set.triplo.tsv

          # DDG2P gof + other alone
          cat gene_lists/DDG2P.all_gof.genes.list \
              gene_lists/DDG2P.all_other.genes.list \
          | sort -Vk1,1 | uniq \
          > ddg2p_truth_set.gof.tsv

          # Combine all truth sets into master union
          cat union_truth_set.gof.tsv \
              clingen_truth_set.triplo.tsv \
              ddg2p_truth_set.gof.tsv \
          | sort -Vk1,1 | uniq \
          > master_union_truth_set.gof.tsv

          # Write truth set input tsv
          for wrapper in 1; do
            echo -e "Union truth set\tmaster_union_truth_set.gof.tsv\tgrey25"
            echo -e "ClinGen & DECIPHER GoF (union)\tunion_truth_set.gof.tsv\t#1A5985"
            echo -e "ClinGen dom. TS\tclingen_truth_set.triplo.tsv\t#2376B2"
            echo -e "DECIPHER dom. GoF/unk.\tddg2p_truth_set.gof.tsv\t#4F91C1"
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
    docker: "talkowski/rcnv@sha256:bde36542e44b22b1cc17940a4faa3ec79a386e8e3d7f7214e0b8d5c6150fed58"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    String best_model = read_string("best_model.tsv")
    File final_scores = "${freq_code}.gene_scores.tsv.gz"
  }
}