######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
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
    call get_expression_data {
      input:
        gtf=gtf,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket
    }
    call get_chromatin_data {
      input:
        gtf=gtf,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket
    }
    call get_constraint_data {
      input:
        gtf=gtf,
        ref_fasta=ref_fasta,
        ref_fasta_idx=ref_fasta_idx,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket
    }
    call join_data {
      input:
        genomic_metadata=get_genomic_data.metadata_table,
        expression_metadata=get_expression_data.metadata_table,
        chromatin_metadata=get_chromatin_data.metadata_table,
        constraint_metadata=get_constraint_data.metadata_table,
        prefix="${gtf_prefix}.merged_features.${contig[0]}"
    }
  }

  # Concatenate across contigs & upload to rCNV bucket
  call cat_metadata as merge_genomic_data {
    input:
      data_shards=get_genomic_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.genomic_features",
      eigen_prefix="genomic_eigenfeature",
      rCNV_bucket=rCNV_bucket
  }
  call cat_metadata as merge_expression_data {
    input:
      data_shards=get_expression_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.expression_features",
      eigen_prefix="expression_eigenfeature",
      rCNV_bucket=rCNV_bucket
  }
  call cat_metadata as merge_chromatin_data {
    input:
      data_shards=get_chromatin_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.chromatin_features",
      eigen_prefix="chromatin_eigenfeature",
      rCNV_bucket=rCNV_bucket
  }
  call cat_metadata as merge_constraint_data {
    input:
      data_shards=get_constraint_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.constraint_features",
      eigen_prefix="constraint_eigenfeature",
      rCNV_bucket=rCNV_bucket
  }
  call cat_metadata as merge_joined_data {
    input:
      data_shards=join_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.all_features",
      eigen_prefix="joined_eigenfeature",
      rCNV_bucket=rCNV_bucket
  }

  # Outputs
  output {
    File genes_genomic_metadata = merge_genomic_data.merged_data
    File genes_genomic_eigen_metadata = merge_genomic_data.merged_data_eigen
    File genes_expression_metadata = merge_expression_data.merged_data
    File genes_expression_eigen_metadata = merge_expression_data.merged_data_eigen
    File genes_chromatin_metadata = merge_chromatin_data.merged_data
    File genes_chromatin_eigen_metadata = merge_chromatin_data.merged_data_eigen
    File genes_constraint_metadata = merge_constraint_data.merged_data
    File genes_constraint_eigen_metadata = merge_constraint_data.merged_data_eigen
    File genes_joined_metadata = merge_joined_data.merged_data
    File genes_joined_eigen_metadata = merge_joined_data.merged_data_eigen
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
      --outbed ${prefix}.genomic_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:5239898782e9936bf377373935fce5829907f68276b35e196204ba3f4615496e"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.genomic_features.${contig}.bed.gz"
  }
}


# Collect expression metadata for all genes from a single contig
task get_expression_data {
  File gtf
  String prefix
  String contig
  String rCNV_bucket

  command <<<
    set -e

    # Download precomputed GTEx expression matrices (generated by preprocess_GTEx.py)
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/genes/annotations/gtex_stats \
      ./

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz

    # Collect genomic metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-expression \
      --gtex-medians gtex_stats/*.GTEx_v7_expression_stats.median.tsv.gz \
      --gtex-mads gtex_stats/*.GTEx_v7_expression_stats.mad.tsv.gz \
      --gtex-pca gtex_stats/*.GTEx_v7_expression_stats.pca.tsv.gz \
      --outbed ${prefix}.expression_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    # TODO: UPDATE DOCKER
    # docker: "talkowski/rcnv@sha256:5239898782e9936bf377373935fce5829907f68276b35e196204ba3f4615496e"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.expression_features.${contig}.bed.gz"
  }
}


# Collect chromatin metadata for all genes from a single contig
task get_chromatin_data {
  File gtf
  String prefix
  String contig
  String rCNV_bucket

  command <<<
    set -e

    # Download precomputed Roadmap chromatin data (generated by preprocess_REP.py)
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/genes/annotations/roadmap_stats \
      ./

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz

    # Collect genomic metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-chromatin \
      --roadmap-medians roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats.mean.tsv.gz \
      --roadmap-mads roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats.sd.tsv.gz \
      --roadmap-pca roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats.pca.tsv.gz \
      --outbed ${prefix}.chromatin_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    # TODO: UPDATE DOCKER
    # docker: "talkowski/rcnv@sha256:5239898782e9936bf377373935fce5829907f68276b35e196204ba3f4615496e"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.chromatin_features.${contig}.bed.gz"
  }
}


# Collect constraint metadata for all genes from a single contig
task get_constraint_data {
  File gtf
  File ref_fasta
  File ref_fasta_idx
  String prefix
  String contig
  String rCNV_bucket

  command <<<
    set -e

    # Download necessary data & references
    wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
    wget http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt
    gsutil -m cp ${rCNV_bucket}/cleaned_data/genes/annotations/EDS.Wang_2018.tsv.gz ./
    wget https://storage.googleapis.com/gnomad-public/legacy/exac_browser/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz
    wget https://doi.org/10.1371/journal.pgen.1001154.s002

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz
    samtools faidx ${ref_fasta} ${contig} \
    | bgzip -c \
    > ref.fa.gz
    samtools faidx ref.fa.gz

    # Collect genomic metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-constraint \
      --ref-fasta ref.fa.gz \
      --gnomad-constraint gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
      --exac-cnv forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz \
      --rvis-tsv RVIS_Unpublished_ExACv2_March2017.txt \
      --eds-tsv EDS.Wang_2018.tsv.gz \
      --hi-tsv journal.pgen.1001154.s002 \
      --outbed ${prefix}.constraint_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:9da6bae9884d16b82a7bc702c52a6646529d014b529062b4df19d6f3ee1acc7d"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.constraint_features.${contig}.bed.gz"
  }
}

# Join metadata for multiple conditions
task join_data {
  File genomic_metadata
  File expression_metadata
  File chromatin_metadata
  File constraint_metadata
  String prefix

  command <<<
    set -e

    /opt/rCNV2/data_curation/gene/join_gene_metadata.R \
      ${genomic_metadata} \
      ${expression_metadata} \
      ${chromatin_metadata} \
      ${constraint_metadata} \
    | bgzip -c \
    > ${prefix}.bed.gz
  >>>

  runtime {
    # TODO: UPDATE DOCKER
    # docker: "talkowski/rcnv@sha256:9da6bae9884d16b82a7bc702c52a6646529d014b529062b4df19d6f3ee1acc7d"
    preemptible: 1
    disks: "local-disk 50 SSD"
    bootDiskSizeGb: "20"
  }

  output {
    File metadata_table = "${prefix}.bed.gz"
  }
}


# Generic function to merge metadata tables across chromosomes & upload to rCNV_bucket
task cat_metadata {
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
      ${rCNV_bucket}/cleaned_data/genes/metadata/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:071cc1595d67604f517d31a5743961834d3ba6f9f95151d2a4a1e64708e307a9"
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

