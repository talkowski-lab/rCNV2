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
    call perm_shard_litGDs {
      input:
        seed=seed,
        n_perms=perms_per_shard,
        whitelist=perm_prep.whitelist,
        perm_prefix="${perm_prefix}.lit_GDs.starting_seed_${seed}",
        rCNV_bucket=rCNV_bucket
    }
  }

  call merge_perms as merge_loci {
    input:
      perm_tables=perm_shard.perm_table,
      perm_prefix="${perm_prefix}.${total_n_perms}_permuted_segments",
      rCNV_bucket=rCNV_bucket
  }
  call merge_perms as merge_lit_gds {
    input:
      perm_tables=perm_shard_litGDs.perm_table,
      perm_prefix="${perm_prefix}.lit_GDs.${total_n_perms}_permuted_segments",
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
    docker: "talkowski/rcnv@sha256:6eb55e14c31a06071dc2c220566cc3665cbfdc6a662efd7c0776f1a1d360f829"
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
    docker: "talkowski/rcnv@sha256:2176cfce20a878f1dd3288e8e534c57689c58663d0354ef68b5c7c18bf34d4b7"
    preemptible: 1
  }

  output {
    File perm_table = "${perm_prefix}.master_segments.bed.gz"
  }  
}


# Single shard of segment permutation for literature-reported GDs
task perm_shard_litGDs {
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

    # Pool & reformat GDs for permutation
    cat <( zcat refs/lit_GDs.hc.bed.gz | sed -n '1p' | cut -f1-5 \
           | paste - <( echo "intervals" ) ) \
        <( zcat refs/lit_GDs.*.bed.gz \
           | fgrep -v "#" | cut -f1-5 \
           | sort -Vk1,1 -k2,2n -k3,3n -k5,5V \
           | awk -v OFS="\t" '{ print $0, $1":"$2"-"$3 }' ) \
    | bgzip -c \
    > lit_GDs.pooled.bed.gz

    # Permute literature GDs
    /opt/rCNV2/analysis/paper/scripts/large_segments/shuffle_segs.py \
      --genome refs/GRCh37.autosomes.genome \
      --whitelist ${whitelist} \
      --n-perms ${n_perms} \
      --first-seed ${seed} \
      --outfile ${perm_prefix}.all.bed.gz \
      --bgzip \
      lit_GDs.pooled.bed.gz

    # Annotate with genes
    /opt/rCNV2/analysis/sliding_windows/get_genes_per_region.py \
      --bgzip \
      --outbed ${perm_prefix}.all.w_genes.bed.gz \
      ${perm_prefix}.all.bed.gz \
      refs/gencode.v19.canonical.pext_filtered.gtf.gz

    # Split out by confidence
    for conf in HC MC LC; do 
      cat <( zcat ${perm_prefix}.all.w_genes.bed.gz | sed -n '1p' \
             | awk '{ print "#"$0 }' ) \
          <( zcat ${perm_prefix}.all.w_genes.bed.gz \
             | awk -v conf=$conf '{ if ($4 ~ "_"conf"_GD_") print $0 }' ) \
      | bgzip -c \
      > ${perm_prefix}.$conf.bed.gz
    done

    # Create dummy of genome-wide significant segments for final segment compilation
    # Uses first row of final segments file with chromosome set to "Z" to ensure
    # no possible overlaps
    cat <( zcat rCNV.final_segments.loci.bed.gz | sed -n '1p' ) \
        <( zcat rCNV.final_segments.loci.bed.gz | sed -n '2p' \
           | awk -v OFS="\t" '{ $1="Z"; print $0 }' ) \
    > dummy.loci.bed

    # Compile new segment table with existing GDs and NAHR regions
    /opt/rCNV2/analysis/paper/scripts/large_segments/compile_segment_table.py \
      --final-loci dummy.loci.bed \
      --hc-gds ${perm_prefix}.HC.bed.gz \
      --mc-gds ${perm_prefix}.MC.bed.gz \
      --lc-gds ${perm_prefix}.LC.bed.gz \
      --nahr-cnvs clustered_nahr_regions.reformatted.bed.gz \
      --outfile ${perm_prefix}.master_segments.w_dummy.bed.gz \
      --gd-recip "10e-10" \
      --nahr-recip 0.25 \
      --bgzip
    zcat ${perm_prefix}.master_segments.w_dummy.bed.gz \
    | grep -ve '^Z' \
    | bgzip -c \
    > ${perm_prefix}.master_segments.bed.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:2176cfce20a878f1dd3288e8e534c57689c58663d0354ef68b5c7c18bf34d4b7"
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
    docker: "talkowski/rcnv@sha256:2176cfce20a878f1dd3288e8e534c57689c58663d0354ef68b5c7c18bf34d4b7"
    preemptible: 1
    disks: "local-disk 200 SSD"
  }

  output {
    File merged_results = "${perm_prefix}.bed.gz"
  }
}
