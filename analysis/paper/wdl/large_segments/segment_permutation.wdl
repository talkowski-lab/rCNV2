######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Parallelize segment permutations for rCNV2 paper analyses


workflow segment_permutation {
  File binned_genome
  Int total_n_perms
  Int perms_per_shard
  String perm_prefix
  String rCNV_bucket

  call perm_prep {
    input:
      binned_genome=binned_genome,
      total_n_perms=total_n_perms,
      perms_per_shard=perms_per_shard
  }

  scatter ( seed in perm_prep.seeds ) {
    call perm_shard {
      input:
        seed=seed,
        n_perms=perms_per_shard,
        whitelist=perm_prep.whitelist,
        perm_prefix="${perm_prefix}.starting_seed_${seed}",
        rCNV_bucket=rCNV_bucket
    }
  }

  call merge_perms {
    input:
      perm_tables=perm_shard.perm_table,
      perm_prefix="${perm_prefix}.${total_n_perms}_permuted_segments",
      rCNV_bucket=rCNV_bucket
  }
}


# Prepare permutation seeds and whitelist
task perm_prep {
  File binned_genome
  Int total_n_perms
  Int perms_per_shard

  command <<<
    set -e

    # Make whitelist
    bedtools merge -i ${binned_genome} | bgzip -c > whitelist.bed.gz

    # Make seeds
    seq 1 ${perms_per_shard} ${total_n_perms} > seeds.txt
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:af27d859af7d4675ae014cb66d1f48bc8ceb7033d4f25677cfc8ec6baac47751"
    preemptible: 1
  }

  output {
    File whitelist = "whitelist.bed.gz"
    Array[String] seeds = read_lines("seeds.txt")
  }
}


# Single shard of segment permutation
task perm_shard {
  String seed
  Int n_perms
  File whitelist
  String perm_prefix
  String rCNV_bucket

  command <<<
    set -e

    # Localize necessary reference files
    mkdir refs
    gsutil -m cp \
      ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
      ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
      ${rCNV_bucket}/analysis/paper/data/large_segments/clustered_nahr_regions.bed.gz \
      ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs.*.bed.gz \
      refs/
    gsutil -m cp \
      ${rCNV_bucket}/results/segment_association/rCNV.final_segments.loci.bed.gz \
      ./

    # Reformat NAHR segments
    cat <( echo -e "#chr\tstart\tend\tnahr_id\tcnv\tn_genes\tgenes" ) \
    <( zcat refs/clustered_nahr_regions.bed.gz | grep -ve '^#' \
       | awk -v FS="\t" -v OFS="\t" \
         '{ print $1, $2, $3, $4"_DEL", "DEL", $5, $6"\n"$1, $2, $3, $4"_DUP", "DUP", $5, $6 }' ) \
    | bgzip -c \
    > clustered_nahr_regions.reformatted.bed.gz

    # Permute segments
    /opt/rCNV2/analysis/paper/scripts/large_segments/shuffle_segs.py \
      --genome refs/GRCh37.autosomes.genome \
      --whitelist ${whitelist} \
      --n-perms ${n_perms} \
      --first-seed ${seed} \
      --outfile ${perm_prefix}.bed.gz \
      --bgzip \
      <( zcat rCNV.final_segments.loci.bed.gz | cut -f1-5,19 )

    # Annotate with genes
    /opt/rCNV2/analysis/sliding_windows/get_genes_per_region.py \
      --bgzip \
      --outbed ${perm_prefix}.w_genes.bed.gz \
      ${perm_prefix}.bed.gz \
      refs/gencode.v19.canonical.pext_filtered.gtf.gz

    # Reformat to match original locus association stats file (add effect sizes, etc)
    /opt/rCNV2/analysis/paper/scripts/large_segments/reformat_shuffled_sig_segs.R \
      --bgzip \
      ${perm_prefix}.w_genes.bed.gz \
      rCNV.final_segments.loci.bed.gz \
      ${perm_prefix}.reformatted.permuted_loci.bed.gz

    # Compile new segment table with existing GDs and NAHR regions
    /opt/rCNV2/analysis/paper/scripts/large_segments/compile_segment_table.py \
      --final-loci ${perm_prefix}.reformatted.permuted_loci.bed.gz \
      --hc-gds refs/lit_GDs.hc.bed.gz \
      --mc-gds refs/lit_GDs.mc.bed.gz \
      --lc-gds refs/lit_GDs.lc.bed.gz \
      --nahr-cnvs clustered_nahr_regions.reformatted.bed.gz \
      --outfile ${perm_prefix}.master_segments.bed.gz \
      --gd-recip "10e-10" \
      --nahr-recip 0.25 \
      --bgzip
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:77adc010fe1501bf69f75988e7251f2a7b1c80fa266a9caea3a4b6f87ea3ff0f"
    preemptible: 1
  }

  output {
    File perm_table = "${perm_prefix}.master_segments.bed.gz"
  }  
}


# Merge permuted results and add permutation index
task merge_perms {
  Array[File] perm_tables
  String perm_prefix
  String rCNV_bucket

  command <<<
    set -e
    
    # Write header and add perm idx as last column
    zcat ${perm_tables[0]} \
    | sed -n '1p' \
    | awk -v OFS="\t" '{ print $0, "perm_idx" }' \
    > header.tsv

    # Merge and sort permutation results
    zcat ${sep=" " perm_tables} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    > perm_res.body.bed

    # Make list of perm idxs
    cut -f4 perm_res.body.bed \
    | cut -f1 -d\_ \
    | sed 's/^perm//g' \
    > perm_idxs.tsv

    # Merge results and remove perm prefix from region IDs
    sed 's/perm[0-9]\+_//g' perm_res.body.bed \
    | paste - perm_idxs.tsv \
    | cat header.tsv - \
    | bgzip -c \
    > ${perm_prefix}.bed.gz

    # Copy to rCNV bucket
    gsutil -m cp \
      ${perm_prefix}.bed.gz \
      ${rCNV_bucket}/analysis/paper/data/large_segments/permutations/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:77adc010fe1501bf69f75988e7251f2a7b1c80fa266a9caea3a4b6f87ea3ff0f"
    preemptible: 1
    disks: "local-disk 200 SSD"
  }

  output {
    File merged_results = "${perm_prefix}.bed.gz"
  }
}
