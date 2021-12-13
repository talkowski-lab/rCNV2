######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Generate final plots for large segment permutation analyses

# Necessary to run in GCP because it requires too much memory to plot locally


workflow PlotLargeSegPermAnalyses {
  String prefix
  Int n_seg_perms
  String rCNV_bucket
  String rCNV_docker

  call PlotPermBySize {
    input:
      prefix=prefix,
      n_seg_perms=n_seg_perms,
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker
  }

  call PlotPermByGene {
    input:
      prefix=prefix,
      n_seg_perms=n_seg_perms,
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker
  }

  output {
    File perm_by_size_plots = PlotPermBySize.tarball
    File perm_by_gene_plots = PlotPermByGene.tarball
  }
}


task PlotPermBySize {
  String prefix
  Int n_seg_perms
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -euo pipefail

    # Localize all necessary data
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/genes/gene_lists \
      ${rCNV_bucket}/results/segment_association/* \
      ${rCNV_bucket}/analysis/paper/data/large_segments/${prefix}.master_segments.bed.gz \
      ${rCNV_bucket}/analysis/paper/data/large_segments/permutations/${prefix}*${n_seg_perms}_permuted_segments.bed.gz \
      ./

    # Plot
    mkdir perm_test_plots
    /opt/rCNV2/analysis/paper/plot/large_segments/plot_segment_permutations.R \
      --constrained-genes gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
      rCNV.final_segments.loci.bed.gz \
      ${prefix}.master_segments.bed.gz \
      ${prefix}.${n_seg_perms}_permuted_segments.bed.gz \
      ${prefix}.lit_GDs.${n_seg_perms}_permuted_segments.bed.gz \
      perm_test_plots/ \
      "${prefix}.segment_perms"
    gzip perm_test_plots/${prefix}.segment_perms.permBySize.stats.tsv

    # Copy results to rCNV bucket
    gsutil -m cp -r perm_test_plots/* \
      ${rCNV_bucket}/analysis/paper/plots/large_segments/perm_test_plots/

    # Rename & compress for local review
    mv perm_test_plots PermBySizePlots
    tar -czvf PermBySizePlots.tar.gz PermBySizePlots
  >>>

  output {
    File tarball = "PermBySizePlots.tar.gz"
  }

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "40 GB"
    cpu: 8
  }
}


task PlotPermByGene {
  String prefix
  Int n_seg_perms
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -euo pipefail

    # Localize all necessary data
    gsutil -m cp \
      ${rCNV_bucket}/results/segment_association/* \
      ${rCNV_bucket}/analysis/paper/data/large_segments/${prefix}.master_segments.bed.gz \
      ${rCNV_bucket}/analysis/paper/data/large_segments/permutations/${prefix}*${n_seg_perms}_permuted_segments_bygene.tsv.gz \
      ./

    # Plot
    mkdir perm_test_plots
    /opt/rCNV2/analysis/paper/plot/large_segments/plot_segment_bygene_perms.R \
      rCNV.final_segments.loci.bed.gz \
      ${prefix}.master_segments.bed.gz \
      ${prefix}.${n_seg_perms}_permuted_segments_bygene.tsv.gz \
      ${prefix}.lit_GDs.${n_seg_perms}_permuted_segments_bygene.tsv.gz \
      perm_test_plots/ \
      "${prefix}.segment_perms_bygene"
    gzip perm_test_plots/${prefix}.segment_perms_bygene.permByGene.stats.tsv

    # Copy results to rCNV bucket
    gsutil -m cp -r perm_test_plots/* \
      ${rCNV_bucket}/analysis/paper/plots/large_segments/perm_test_plots/

    # Rename & compress for local review
    mv perm_test_plots PermByGenePlots
    tar -czvf PermByGenePlots.tar.gz PermByGenePlots
  >>>

  output {
    File tarball = "PermByGenePlots.tar.gz"
  }

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "40 GB"
    cpu: 8
  }
}