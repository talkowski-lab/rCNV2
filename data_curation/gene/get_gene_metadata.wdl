######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Collect gene-level metadata

workflow get_gene_metadata {
  File gtf
  File ref_fasta
  File ref_fasta_idx
  Float eigenfeatures_min_var_exp
  String gtf_prefix
  File contiglist
  String rCNV_bucket

  Array[Array[String]] contigs = read_tsv(contiglist)

  # Scatter over contigs and gather metadata for all genes per contig
  scatter ( contig in contigs ) {
    call get_genomic_data {
      input:
        gtf=gtf,
        ref_fasta=ref_fasta,
        ref_fasta_idx=ref_fasta_idx,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket
    }
  }

  # Merge across contigs & upload to rCNV bucket
  call merge_metadata as merge_genomic_data {
    input:
      data_shards=get_genomic_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.genomic_features",
      eigen_prefix="genomic_eigenfeature",
      rCNV_bucket=rCNV_bucket
  }

  output {
    File genes_genomic_metadata = merge_genomic_data.merged_data
    File genes_genomic_eigen_metadata = merge_genomic_data.merged_data_eigen
  }
}


# Collect genomic metadata for all genes from a single contig
task get_genomic_data {
  File gtf
  File ref_fasta
  File ref_fasta_idx
  String prefix
  String contig
  String rCNV_bucket

  command <<<
    set -e

    # Download necessary data & references
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/** refs/
    gsutil -m cp ${rCNV_bucket}/refs/** refs/

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz
    samtools faidx ${ref_fasta} ${contig} \
    | bgzip -c \
    > ref.fa.gz
    samtools faidx ref.fa.gz

    # Prep athena input file
    for wrapper in 1; do 
      echo -e "refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz\tcoverage\trepeat_cov"
      echo -e "refs/GRCh37.somatic_hypermutable_sites.bed.gz\tcoverage\thypermutable_cov"
      echo -e "https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign100mer.bw\tmap-mean\talign_100mer"
    done > gene_features.athena_tracklist.tsv

    # Collect genomic metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-genomic \
      --centro-telo-bed refs/GRCh37.centromeres_telomeres.bed.gz \
      --ref-fasta ref.fa.gz \
      --athena-tracks gene_features.athena_tracklist.tsv \
      --outbed ${prefix}.genomic_features.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:3f33790812f5d5a9c27104f437a1e1e01513b0bd89c199c08d94ec6b144533ae"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.genomic_features.bed.gz"
  }
}


# Generic function to merge metadata tables across chromosomes & upload to rCNV_bucket
task merge_metadata {
  Array[File] data_shards
  Float eigenfeatures_min_var_exp
  String prefix
  String eigen_prefix
  String rCNV_bucket

  command <<<
    set -e

    zcat ${data_shards[0]} \
    | head -n1 \
    > header.txt

    zcat ${sep=" " data_shards} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | cat header.txt - \
    | bgzip -c \
    > ${prefix}.bed.gz

    athena eigen-bins \
      --min-variance ${eigenfeatures_min_var_exp} \
      --fill-missing mean \
      --skip-columns 4 \
      --prefix ${eigen_prefix} \
      --bgzip \
      ${prefix}.bed.gz \
      ${prefix}.eigenfeatures.bed.gz

    gsutil -m cp \
      ${prefix}*bed.gz \
      ${rCNV_bucket}/cleaned_data/genes/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:3f33790812f5d5a9c27104f437a1e1e01513b0bd89c199c08d94ec6b144533ae"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File merged_data = "${prefix}.bed.gz"
   File merged_data_eigen = "${prefix}.eigenfeatures.bed.gz"
  }
}

