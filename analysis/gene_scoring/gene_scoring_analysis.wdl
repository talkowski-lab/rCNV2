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
  Float prior_frac
  String rCNV_bucket
  File contiglist

  Array[Array[String]] contigs = read_tsv(contiglist)

  Array[String] cnv_types = ["DEL", "DUP"]

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
        freq_code="rCNV"
    }
  }

  # Make training blacklist, estimate prior effect sizes, and compute BFDPs
  call blacklist_priors_bfdp {
    input:
      del_meta_stats=rCNV_meta_analysis.meta_stats_bed[0],
      dup_meta_stats=rCNV_meta_analysis.meta_stats_bed[1],


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
  String freq_code

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
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:da2df0f4afcfa93d17e27ae5752f1665f0c53c8feed06d6d98d1da53144d8e1f"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    File meta_stats_bed = "${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz"
    File meta_stats_bed_idx = "${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz.tbi"
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
    zcat GRCh37.centromeres_telomeres.bed.gz \
    | awk -v FS="\t" -v OFS="\t" -v d=${cen_tel_dist} \
      '{ if ($2-d<0) print $1, "0", $3+d; else print $1, $2-d, $3+d }' \
    | cat - <( zcat /opt/rCNV2/refs/UKBB_GD.Owen_2018.bed.gz | cut -f1-3 ) \
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
    theta0_del=$( awk -v FS="\t" '{ if ($1=="theta0" && $2=="DEL") print $3 }' ${freq_code}.prior_estimation.empirical_prior_estimates.tsv )
    theta0_dup=$( awk -v FS="\t" '{ if ($1=="theta0" && $2=="DUP") print $3 }' ${freq_code}.prior_estimation.empirical_prior_estimates.tsv )
    theta1_del=$( awk -v FS="\t" '{ if ($1=="theta1" && $2=="DEL") print $3 }' ${freq_code}.prior_estimation.empirical_prior_estimates.tsv )
    var1=$( awk -v FS="\t" '{ if ($1=="var1" && $2=="DEL") print $3 }' ${freq_code}.prior_estimation.empirical_prior_estimates.tsv )

    # Compute BFDP per gene
    for CNV in DEL DUP; do
      # Set CNV-specific variables
      case $CNV in
        "DEL")
          theta0=${theta0_del}
          statsbed=${del_meta_stats}
          ;;
        "DUP")
          theta0=${theta0_dup}
          statsbed=${dup_meta_stats}
          ;;
      esac

      # Compute BF & BFDR for all genes
      /opt/rCNV2/analysis/gene_scoring/calc_gene_bfs.py \
        --theta0 $theta0 \
        --theta1 ${theta1} \
        --var0 ${var1} \
        --prior ${prior_frac} \
        --blacklist ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz \
        --outfile ${freq_code}.$CNV.gene_abfs.tsv \
        "$statsbed"
    done
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:da2df0f4afcfa93d17e27ae5752f1665f0c53c8feed06d6d98d1da53144d8e1f"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    File training_blacklist = "${freq_code}.gene_scoring.training_gene_blacklist.bed.gz"
    File priors_tsv = "${freq_code}.prior_estimation.empirical_prior_estimates.tsv"
    Float theta0_del = "$theta0_del"
    Float theta0_dup = "$theta0_dup"
    Float theta1 = "$theta1"
    Float var0 = "$var1"
    File del_bfdp = "${freq_code}.DEL.gene_abfs.tsv"
    File dup_bfdp = "${freq_code}.DUP.gene_abfs.tsv"
    Array[File] prior_plots = glob("${freq_code}.prior_estimation*pdf")
  }

  
}