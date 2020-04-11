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


  >>>


  runtime {
    docker: "talkowski/rcnv@sha256:da2df0f4afcfa93d17e27ae5752f1665f0c53c8feed06d6d98d1da53144d8e1f"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    File training_blacklist = "${freq_code}.gene_scoring.training_gene_blacklist.bed.gz"
  }
}