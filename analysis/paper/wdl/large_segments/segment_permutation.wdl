######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Parallelize segment permutations for rCNV2 paper analyses


workflow segment_permutation {
  File binned_genome
  File phenotype_list
  Int total_n_perms
  Int perms_per_shard
  String perm_prefix
  String rCNV_bucket
  String rCNV_docker

  call perm_prep {
    input:
      binned_genome=binned_genome,
      phenotype_list=phenotype_list,
      total_n_perms=total_n_perms,
      perms_per_shard=perms_per_shard,
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker,
      prefix=perm_prefix
  }

  scatter ( seed in perm_prep.seeds ) {
    call perm_shard {
      input:
        seed=seed,
        n_perms=perms_per_shard,
        whitelist=perm_prep.whitelist,
        blacklist=perm_prep.blacklist,
        del_max_p_bed=perm_prep.del_max_p_bed,
        dup_max_p_bed=perm_prep.dup_max_p_bed,
        perm_prefix="${perm_prefix}.starting_seed_${seed}",
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker
    }
    call perm_shard_litGDs {
      input:
        seed=seed,
        n_perms=perms_per_shard,
        whitelist=perm_prep.whitelist,
        blacklist=perm_prep.blacklist,
        del_max_p_bed=perm_prep.del_max_p_bed,
        dup_max_p_bed=perm_prep.dup_max_p_bed,
        perm_prefix="${perm_prefix}.lit_GDs.starting_seed_${seed}",
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker
    }
  }

  call merge_perms as merge_loci {
    input:
      perm_tables=perm_shard.perm_table,
      perm_prefix="${perm_prefix}.${total_n_perms}_permuted_segments",
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker
  }
  call merge_perms as merge_lit_gds {
    input:
      perm_tables=perm_shard_litGDs.perm_table,
      perm_prefix="${perm_prefix}.lit_GDs.${total_n_perms}_permuted_segments",
      rCNV_bucket=rCNV_bucket,
      rCNV_docker=rCNV_docker
  }
}


# Prepare permutation seeds, whitelist, and BEDs of best P per window
task perm_prep {
  File binned_genome
  File phenotype_list
  Int total_n_perms
  Int perms_per_shard
  String rCNV_bucket
  String rCNV_docker
  String prefix

  command <<<
    set -e

    # Make whitelist after removing untestable bins (those with <2 cohorts for meta-analysis)
    gsutil -m cp \
      ${rCNV_bucket}/analysis/analysis_refs/*.cohort_exclusion.bed.gz \
      ./
    zcat *.cohort_exclusion.bed.gz | sed 's/;/\t/g' \
    | awk -v FS="\t" -v OFS="\t" '{ if (NF<=7) print $1, $2, $3 }' \
    | fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - | bgzip -c > whitelist.bed.gz

    # Invert whitelist as explicit blacklist (pybedtools has some unusual shuffle behavior)
    gsutil -m cat ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
    | awk -v OFS="\t" '{ print $1, 0, 500000000 }' \
    | bedtools subtract -a - -b whitelist.bed.gz \
    | bgzip -c \
    > blacklist.bed.gz

    # Make seeds
    seq 1 ${perms_per_shard} ${total_n_perms} > seeds.txt

    # Localize meta-analysis stats
    mkdir meta_stats/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/sliding_windows/**.rCNV.**.sliding_window.meta_analysis.stats.bed.gz \
      meta_stats/

    # Strip out P-value column per HPO & CNV pair
    mkdir meta_stats/matrices/
    while read nocolon hpo; do
      echo $nocolon
      for cnv in DEL DUP; do
        echo $cnv
        statsfile=meta_stats/$nocolon.rCNV.$cnv.sliding_window.meta_analysis.stats.bed.gz
        idx=$( zcat $statsfile | head -n1 | sed 's/\t/\n/g' \
               | awk -v OFS="\t" '{ print $1, NR }' \
               | fgrep -w meta_neg_log10_p | cut -f2 )
        if [ -z $idx ]; then
          echo "UNABLE TO FIND COLUMN NAMED meta_neg_log10_p IN $statsfile"
          exit
        fi
        zcat $statsfile | sed '1d' | cut -f $idx \
        | cat <( echo -e $nocolon"_"$cnv ) - \
        > meta_stats/matrices/$nocolon.$cnv.meta_neg_log10_p.tsv
      done
    done < ${phenotype_list}

    # Collect bin coordinates
    zcat \
      meta_stats/$( head -n1 ${phenotype_list} | cut -f1 ).rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
    | cut -f1-3 \
    > window_coordinates.bed

    # Make matrices for primary P-values across phenotypes per CNV type 
    for cnv in DEL DUP; do
      paste \
        window_coordinates.bed \
        meta_stats/matrices/*.$cnv.meta_neg_log10_p.tsv \
      | bgzip -c \
      > meta_stats/matrices/${prefix}.$cnv.meta_neg_log10_p.all_hpos.bed.gz
    done

    # Compress to single BED of max P per window
    for cnv in DEL DUP; do
      /opt/rCNV2/analysis/paper/scripts/large_segments/get_best_p_per_window.R \
        meta_stats/matrices/${prefix}.$cnv.meta_neg_log10_p.all_hpos.bed.gz \
        ${prefix}.best_meta_neg_log10_p_per_window.$cnv.bed
      bgzip -f ${prefix}.best_meta_neg_log10_p_per_window.$cnv.bed
    done
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
  }

  output {
    File whitelist = "whitelist.bed.gz"
    File blacklist = "blacklist.bed.gz"
    Array[String] seeds = read_lines("seeds.txt")
    File del_max_p_bed = "${prefix}.best_meta_neg_log10_p_per_window.DEL.bed.gz"
    File dup_max_p_bed = "${prefix}.best_meta_neg_log10_p_per_window.DUP.bed.gz"
  }
}


# Single shard of segment permutation
task perm_shard {
  String seed
  Int n_perms
  File whitelist
  File blacklist
  File del_max_p_bed
  File dup_max_p_bed
  String perm_prefix
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -e

    # Localize necessary reference files
    mkdir refs
    gsutil -m cp \
      ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
      ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
      ${rCNV_bucket}/analysis/paper/data/large_segments/clustered_nahr_regions.bed.gz \
      ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs.*.bed.gz \
      ${rCNV_bucket}/analysis/paper/data/misc/redin_bca_breakpoints.bed.gz \
      ${rCNV_bucket}/analysis/analysis_refs/test_phenotypes.list \
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
      --blacklist ${blacklist} \
      --n-perms ${n_perms} \
      --first-seed ${seed} \
      --outfile ${perm_prefix}.bed.gz \
      --bgzip \
      <( zcat rCNV.final_segments.loci.bed.gz | cut -f1-5,20 )

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

    # Recompute max P-value for all permuted loci 
    # (dummy for regular calc_all_seg_stats.py, for speed)
    echo -e "region_id\thpo\tcnv\tlnor\tlnor_lower\tlnor_upper\tpvalue\tpvalue_secondary" \
    > ${perm_prefix}.final_segments.loci.all_sumstats.tsv
    for cnv in DEL DUP; do
      if [ $cnv == "DEL" ]; then
        best_p_bed=${del_max_p_bed}
      else
        best_p_bed=${dup_max_p_bed}
      fi
      zcat ${perm_prefix}.reformatted.permuted_loci.bed.gz \
      | fgrep -w $cnv | cut -f1-4 | sort -Vk1,1 -k2,2n -k3,3n \
      | bedtools map -g refs/GRCh37.autosomes.genome \
          -c 4 -o max -a - -b $best_p_bed \
      | awk -v FS="\t" -v OFS="\t" -v cnv=$cnv \
        '{ if ($NF==".") best=0; else best=$NF; print $4, "best", cnv, "NA", "NA", "NA", best, "NA" }' \
      >> ${perm_prefix}.final_segments.loci.all_sumstats.tsv
    done

    # Compile new segment table with existing GDs and NAHR regions
    /opt/rCNV2/analysis/paper/scripts/large_segments/compile_segment_table.py \
      --final-loci ${perm_prefix}.reformatted.permuted_loci.bed.gz \
      --hc-gds refs/lit_GDs.hc.bed.gz \
      --mc-gds refs/lit_GDs.mc.bed.gz \
      --lc-gds refs/lit_GDs.lc.bed.gz \
      --nahr-cnvs clustered_nahr_regions.reformatted.bed.gz \
      --bca-tsv refs/redin_bca_breakpoints.bed.gz \
      --meta-sumstats ${perm_prefix}.final_segments.loci.all_sumstats.tsv \
      --gd-recip "10e-10" \
      --nahr-recip 0.5 \
    | awk -v FS="\t" '{ if ($4 !~ /_GD_|nahr_/) print $0 }' \
    | bgzip -c \
    > ${perm_prefix}.master_segments.bed.gz
  >>>

  runtime {
    docker: "${rCNV_docker}"
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
  File blacklist
  File del_max_p_bed
  File dup_max_p_bed
  String perm_prefix
  String rCNV_bucket
  String rCNV_docker

  command <<<
    set -e

    # Localize necessary reference files
    mkdir refs
    gsutil -m cp \
      ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
      ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
      ${rCNV_bucket}/analysis/paper/data/large_segments/clustered_nahr_regions.bed.gz \
      ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs.*.bed.gz \
      ${rCNV_bucket}/analysis/paper/data/misc/redin_bca_breakpoints.bed.gz \
      ${rCNV_bucket}/analysis/analysis_refs/test_phenotypes.list \
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
      --blacklist ${blacklist} \
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

    # Recompute max P-value for all permuted loci 
    # (dummy for regular calc_all_seg_stats.py, for speed)
    echo -e "region_id\thpo\tcnv\tlnor\tlnor_lower\tlnor_upper\tpvalue\tpvalue_secondary" \
    > ${perm_prefix}.permuted_gds.all_sumstats.tsv
    for cnv in DEL DUP; do
      if [ $cnv == "DEL" ]; then
        best_p_bed=${del_max_p_bed}
      else
        best_p_bed=${dup_max_p_bed}
      fi
      zcat ${perm_prefix}.all.w_genes.bed.gz \
      | fgrep -w $cnv | cut -f1-4 | sort -Vk1,1 -k2,2n -k3,3n \
      | bedtools map -g refs/GRCh37.autosomes.genome \
          -c 4 -o max -a - -b $best_p_bed \
      | awk -v FS="\t" -v OFS="\t" -v cnv=$cnv \
        '{ if ($NF==".") best=0; else best=$NF; print $4, "best", cnv, "NA", "NA", "NA", best, "NA" }' \
      >> ${perm_prefix}.permuted_gds.all_sumstats.tsv
    done

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
      --keep-all-gds \
      --nahr-cnvs clustered_nahr_regions.reformatted.bed.gz \
      --bca-tsv refs/redin_bca_breakpoints.bed.gz \
      --meta-sumstats ${perm_prefix}.permuted_gds.all_sumstats.tsv \
      --gd-recip "10e-10" \
      --nahr-recip 0.5 \
    | awk -v FS="\t" '{ if ($4 !~ /^nahr_/) print $0 }' \
    | grep -ve '^Z' \
    | bgzip -c \
    > ${perm_prefix}.master_segments.bed.gz
  >>>

  runtime {
    docker: "${rCNV_docker}"
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
  String rCNV_docker

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
    docker: "${rCNV_docker}"
    preemptible: 1
    disks: "local-disk 200 SSD"
  }

  output {
    File merged_results = "${perm_prefix}.bed.gz"
  }
}
