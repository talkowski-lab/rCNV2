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
  Int pad_controls
  Int max_cnv_size
  String weight_mode
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  String meta_model_prefix
  String rCNV_bucket
  File contiglist

  Array[String] contigs = read_tsv(contiglist)

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
        pad_controls=pad_controls,
        max_cnv_size=max_cnv_size,
        weight_mode=weight_mode,
        min_cds_ovr_del=min_cds_ovr_del,
        min_cds_ovr_dup=min_cds_ovr_dup,
        max_genes_per_cnv=max_genes_per_cnv,
        rCNV_bucket=rCNV_bucket
    }
  }

  # # Merge association statistics and run meta-analysis
  # scatter ( cnv in cnv_types ) {}
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
          -o "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz" \
          "$cnv_bed" \
          subset.gtf.gz
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
      done
    done < ${metacohort_list}
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
  }
}