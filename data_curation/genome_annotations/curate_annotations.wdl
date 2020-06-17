######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Curate genome annotation tracks, conduct global burden testing, and create cis-regulatory blocks

workflow curate_annotations {
	File chromhmm_tracklist
  Int tracks_per_shard
  Int min_element_size
  Int max_element_size
  Float fdr_cutoff
  String rCNV_bucket
  String prefix

  # Shuffle & shard tracklists
  call shard_tracklist as shard_tracklist_chromhmm {
    input:
      tracklist=chromhmm_tracklist,
      tracks_per_shard=tracks_per_shard,
      prefix="${prefix}.chromhmm_shards"
  }

  # Scatter over sharded tracklists and curate tracks
  scatter ( shard in shard_tracklist_chromhmm.shards ) {
    call curate_and_burden_shard as curate_and_burden_shard_chromhmm {
      input:
        tracklist=shard,
        min_element_size=min_element_size,
        max_element_size=max_element_size,
        rCNV_bucket=rCNV_bucket,
        prefix="${prefix}.chromhmm_shard",
        keep_tracks="FALSE"
    }
  }

  # Collect shards for each tracklist
  call merge_shards as merge_chromhmm {
    input:
      stat_shards=curate_and_burden_shard_chromhmm.stats,
      prefix="${prefix}.chromhmm"
  }

  # Merge stats across all tracklists and compute meta-analysis burden stats
  call merge_and_meta {
    input:
      stats=[merge_chromhmm.merged_stats],
      fdr_cutoff=fdr_cutoff,
      rCNV_bucket=rCNV_bucket,
      prefix=prefix
  }

  # # Re-shard all significant tracks for final curation
  # call shard_tracklist as shard_tracklist_signif {
  #   input:
  #     tracklist=merge_shards_and_meta.signif_tracklist,
  #     tracks_per_shard=tracks_per_shard,
  #     prefix="${prefix}.signif_tracks"
  # }

  # # Scatter over sharded significant tracklist and curate tracks
  # scatter ( shard in shard_tracklist_signif.shards ) {
  #   call curate_and_burden_shard as curate_and_burden_shard_signif {
  #     input:
  #       tracklist=shard,
  #       min_element_size=min_element_size,
  #       max_element_size=max_element_size,
  #       rCNV_bucket=rCNV_bucket,
  #       prefix="${prefix}.signif_shard",
  #       keep_tracks="TRUE"
  #   }
  # }

  # Cluster final cis-regulatory blocks from significant tracks
  # TBD

  # Final outputs
  output {
    File stats_all_tracks = merge_shards_and_meta.meta_stats
    File signif_tracklist = merge_shards_and_meta.signif_tracklist
  }
}


# Shuffle & shard a list of track paths
task shard_tracklist {
  File tracklist
  Int tracks_per_shard
  String prefix

  command <<<
    set -e 

    /opt/rCNV2/utils/evenSplitter.R \
      -L ${tracks_per_shard} \
      --shuffle \
      ${tracklist} \
      ${prefix}
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:2942c7386b43479d02b29506ad1f28fdcff17bdf8b279f2e233be0c4d2cd50fa"
    preemptible: 1
  }

  output {
    Array[File] shards = glob("${prefix}*")
  }
}


# Download & curate tracks, and count CNVs for each track (optional)
task curate_and_burden_shard {
  File tracklist
  Int min_element_size
  Int max_element_size
  String rCNV_bucket
  String prefix
  String keep_tracks

  command <<<
    set -e

    # Download necessary reference files
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/refs/** refs/
    gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/GRCh37.genome refs/
    gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/rCNV_metacohort* refs/

    # Curate all tracks in tracklist
    while read path; do
      echo -e "Curating $path"
      /opt/rCNV2/data_curation/genome_annotations/curate_track.py \
        --genome refs/GRCh37.genome \
        --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
        --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
        --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
        --min-size ${min_element_size} \
        --max-size ${max_element_size} \
        --stats \
        $path
    done < ${tracklist}

    # Collect table of curated stats
    cat *.curated.stats.tsv | fgrep -v "#" | sort -Vk1,1 | uniq \
    | awk -v FS="\t" -v OFS="\t" '{ print $0, $1".curated.bed.gz" }' \
    | cat <( cat *.curated.stats.tsv | grep -e '^#' | sed -n '1p' \
             | awk -v OFS="\t" '{ print $0, "local_path" }' ) \
          - \
    > ${prefix}.stats.tsv

    # Final behavior depends on value of keep_tracks
    if [ ${keep_tracks} == "TRUE" ]; then
      # Copy all curated tracks to final gs:// bucket
      gsutil -m cp \
        *.curated.bed.gz \
        ${rCNV_bucket}/cleaned_data/genome_annotations/significant_tracks/

      # Dummy output
      touch ${prefix}.stats.with_counts.tsv
      gzip -f ${prefix}.stats.with_counts.tsv
    else
      # Localize noncoding CNV data
      mkdir cnvs/
      gsutil -m cp ${rCNV_bucket}/cleaned_data/cnv/noncoding/** cnvs/

      # Compute CNV counts per annotation track per metacohort per CNV type
      fgrep -v mega refs/rCNV_metacohort_list.txt \
      | awk -v OFS="\t" '{ print $1, "cnvs/"$1".rCNV.strict_noncoding.bed.gz" }' \
      > metacohorts.cnv_paths.tsv
      /opt/rCNV2/data_curation/genome_annotations/count_cnvs_per_track.py \
        --cohorts metacohorts.cnv_paths.tsv \
        --track-stats ${prefix}.stats.tsv \
        --outfile ${prefix}.stats.with_counts.tsv.gz \
        --gzip 
    fi
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:2942c7386b43479d02b29506ad1f28fdcff17bdf8b279f2e233be0c4d2cd50fa"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }

  output {
    File stats = "${prefix}.stats.with_counts.tsv.gz"
  }
}


# Merge stats across shards
task merge_shards {
  Array[File] stat_shards
  String prefix


  command <<<
    set -e

    # Merge all stats
    zcat ${stat_shards[0]} | sed -n '1p' > header.tsv
    zcat ${sep=" " stat_shards} | grep -ve '^#' | cat header.tsv - | gzip -c \
    > ${prefix}.merged_stats.with_counts.tsv.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:2942c7386b43479d02b29506ad1f28fdcff17bdf8b279f2e233be0c4d2cd50fa"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }

  output {
    File merged_stats = "${prefix}.merged_stats.with_counts.tsv.gz"
  }
}


# Combine track stat across tracklists and conduct unified meta-analysis burden test
task merge_shards_and_meta {
  Array[File] stats
  Float fdr_cutoff
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Merge all stats
    zcat ${stats[0]} | sed -n '1p' > header.tsv
    zcat ${sep=" " stats} | grep -ve '^#' | cat header.tsv - | gzip -c \
    > ${prefix}.merged_stats.with_counts.tsv.gz

    /opt/rCNV2/data_curation/genome_annotations/trackwise_cnv_burden_meta_analysis.R \
      --model "fe" \
      --spa \
      ${prefix}.merged_stats.with_counts.tsv.gz \
      ${prefix}.burden_stats.tsv
    gzip -f ${prefix}.burden_stats.tsv
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:2942c7386b43479d02b29506ad1f28fdcff17bdf8b279f2e233be0c4d2cd50fa"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }

  output {
    File meta_stats = "${prefix}.burden_stats.tsv.gz"
  }

}
