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
  String case_hpo
  Float min_element_overlap
  Float p_cutoff
  Float min_prop_tracks_per_crb
  Int clustering_neighborhood_dist
  Int min_crb_separation
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
    call curate_and_burden as curate_and_burden_chromhmm {
      input:
        tracklist=shard,
        min_element_size=min_element_size,
        max_element_size=max_element_size,
        case_hpo=case_hpo,
        min_element_overlap=min_element_overlap,
        rCNV_bucket=rCNV_bucket,
        prefix="${prefix}.chromhmm_shard"
    }
  }

  # Collect shards for each tracklist
  call merge_shards as merge_chromhmm {
    input:
      stat_shards=curate_and_burden_chromhmm.stats,
      prefix="${prefix}.chromhmm"
  }

  # Merge stats across all tracklists
  call merge_shards as merge_all {
    input:
      stat_shards=[merge_chromhmm.merged_stats],
      prefix="${prefix}.all"
  }

  # Merge all tracklists
  call merge_tracklists {
    input:
      tracklists=[chromhmm_tracklist],
      prefix="${prefix}"
  }

  # Merge stats across all tracklists and compute meta-analysis burden stats
  call meta_burden_test {
    input:
      stats=merge_all.merged_stats,
      p_cutoff=p_cutoff,
      merged_tracklist=merge_tracklists.merged_tracklist,
      rCNV_bucket=rCNV_bucket,
      prefix=prefix,
      clear_sig="TRUE"
  }

  # Re-shard all significant tracks for final curation
  call shard_tracklist as shard_tracklist_signif {
    input:
      tracklist=meta_burden_test.signif_tracklist,
      tracks_per_shard=tracks_per_shard,
      prefix="${prefix}.signif_tracks"
  }

  # Scatter over sharded significant tracklist and curate tracks
  scatter ( shard in shard_tracklist_signif.shards ) {
    call curate_only as curate_only_signif {
      input:
        tracklist=shard,
        min_element_size=min_element_size,
        max_element_size=max_element_size,
        rCNV_bucket=rCNV_bucket,
        prefix="${prefix}.signif_shard"
    }
  }

  call cluster_elements {
    input:
      completion_tokens=curate_only_signif.completion_token,
      min_prop_tracks_per_crb=min_prop_tracks_per_crb,
      clustering_neighborhood_dist=clustering_neighborhood_dist,
      min_crb_separation=min_crb_separation,
      rCNV_bucket=rCNV_bucket,
      prefix=prefix
  }
  
  # Final outputs
  output {
    File final_crbs = cluster_elements.final_crbs
    File final_crb_elements = cluster_elements.final_crb_elements
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
    docker: "talkowski/rcnv@sha256:7bc27c40820b1f3fe66beb438ef543b8e247bc8d040629b36701c75cd1ada90b"
    preemptible: 1
  }

  output {
    Array[File] shards = glob("${prefix}*")
  }
}


# Download & curate tracks (no burden test)
task curate_only {
  File tracklist
  Int min_element_size
  Int max_element_size
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Download necessary reference files
    mkdir refs/
    gsutil -m cp \
      ${rCNV_bucket}/refs/** \
      ${rCNV_bucket}/analysis/analysis_refs/GRCh37.genome \
      ${rCNV_bucket}/analysis/analysis_refs/rCNV_metacohort* \
      ${rCNV_bucket}/analysis/analysis_refs/HPOs_by_metacohort.table.tsv \
      refs/

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

    # Copy all curated tracks to final gs:// bucket
    gsutil -m cp \
      *.curated.bed.gz \
      ${rCNV_bucket}/cleaned_data/genome_annotations/significant_tracks/

    # Dummy output
    touch completion.txt
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:e2297640f2734f6374832df1ae5388df5bc0ccb692819be3442f4d9a319ab299"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }

  output {
    File completion_token = "completion.txt"
  }
}


# Download & curate tracks, and count CNVs for each track
task curate_and_burden {
  File tracklist
  Int min_element_size
  Int max_element_size
  Float min_element_overlap
  String case_hpo
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Download necessary reference files
    mkdir refs/
    gsutil -m cp \
      ${rCNV_bucket}/refs/** \
      ${rCNV_bucket}/analysis/analysis_refs/GRCh37.genome \
      ${rCNV_bucket}/analysis/analysis_refs/rCNV_metacohort* \
      ${rCNV_bucket}/analysis/analysis_refs/HPOs_by_metacohort.table.tsv \
      refs/

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

    # Localize noncoding CNV data
    mkdir cnvs/
    gsutil -m cp ${rCNV_bucket}/cleaned_data/cnv/noncoding/** cnvs/

    # Compute CNV counts per annotation track per metacohort per CNV type
    while read cohort; do
      for dummy in 1; do
        echo $cohort
        cidx=$( sed -n '1p' refs/HPOs_by_metacohort.table.tsv \
                | sed 's/\t/\n/g' \
                | awk -v cohort=$cohort '{ if ($1==cohort) print NR }' )
        fgrep -w ${case_hpo} refs/HPOs_by_metacohort.table.tsv | cut -f$cidx
        fgrep -w "HEALTHY_CONTROL" refs/HPOs_by_metacohort.table.tsv | cut -f$cidx
        echo -e "cnvs/$cohort.rCNV.strict_noncoding.bed.gz"
      done | paste -s
    done < <( fgrep -v mega refs/rCNV_metacohort_list.txt | cut -f1 ) \
    > metacohorts.input.tsv
    /opt/rCNV2/data_curation/genome_annotations/count_cnvs_per_track.py \
      --cohorts metacohorts.input.tsv \
      --track-stats ${prefix}.stats.tsv \
      --frac-overlap ${min_element_overlap} \
      --case-hpo ${case_hpo} \
      --norm-by-samplesize \
      --outfile ${prefix}.stats.with_counts.tsv.gz \
      --gzip
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:e2297640f2734f6374832df1ae5388df5bc0ccb692819be3442f4d9a319ab299"
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
    docker: "talkowski/rcnv@sha256:e2297640f2734f6374832df1ae5388df5bc0ccb692819be3442f4d9a319ab299"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }

  output {
    File merged_stats = "${prefix}.merged_stats.with_counts.tsv.gz"
  }
}


# Merge all tracklists
task merge_tracklists {
  Array[File] tracklists
  String prefix

  command <<<
    set -e 

    cat ${sep=" " tracklists} > ${prefix}.all_tracks.list
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:e2297640f2734f6374832df1ae5388df5bc0ccb692819be3442f4d9a319ab299"
    preemptible: 1
  }

  output {
    File merged_tracklist = "${prefix}.all_tracks.list"
  }
}

# Conduct meta-analysis burden test for all tracks
task meta_burden_test {
  File stats
  Float p_cutoff
  File merged_tracklist
  String rCNV_bucket
  String prefix
  String clear_sig

  command <<<
    set -e

    # Clear all old significant tracks
    if [ ${clear_sig} == "TRUE" ]; then
      gsutil -m rm ${rCNV_bucket}/cleaned_data/genome_annotations/significant_tracks/*.curated.bed.gz || true
    fi

    # Perform burden analysis
    /opt/rCNV2/data_curation/genome_annotations/trackwise_cnv_burden_meta_analysis.R \
      --p-cutoff ${p_cutoff} \
      --signif-tracks ${prefix}.signif_paths_and_tracks.list \
      ${stats} \
      ${merged_tracklist} \
      ${prefix}.burden_stats.tsv

    # Extract significant tracks
    gzip -f ${prefix}.burden_stats.tsv
    cut -f1 ${prefix}.signif_paths_and_tracks.list \
    > ${prefix}.signif_tracks.list
    cut -f2 ${prefix}.signif_paths_and_tracks.list \
    > ${prefix}.signif_tracknames.list

    # Copy final stats to gs:// bucket
    gsutil -m cp \
      ${prefix}.burden_stats.tsv.gz \
      ${rCNV_bucket}/cleaned_data/genome_annotations/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:e2297640f2734f6374832df1ae5388df5bc0ccb692819be3442f4d9a319ab299"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }

  output {
    File meta_stats = "${prefix}.burden_stats.tsv.gz"
    File signif_tracklist = "${prefix}.signif_tracks.list"
    File signif_tracknames = "${prefix}.signif_tracknames.list"
  }
}


# Cluster elements from significant tracks into CRBs
task cluster_elements {
  Array[File] completion_tokens
  Float min_prop_tracks_per_crb
  Int clustering_neighborhood_dist
  Int min_crb_separation
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Copy necessary reference files
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/GRCh37.genome refs/

    # Copy & index all final curated significant tracks locally
    mkdir sig_tracks/
    gsutil -m cp \
      ${rCNV_bucket}/cleaned_data/genome_annotations/significant_tracks/*.curated.bed.gz \
      sig_tracks/
    find sig_tracks/ -name "*.curated.bed.gz" \
    | xargs -I {} tabix -f {}

    # Subset genome file to autosomes
    grep -e '^[1-9]' refs/GRCh37.genome \
    > autosomes.genome

    # Cluster significant tracks into CRBs
    /opt/rCNV2/data_curation/genome_annotations/build_crbs.py \
      --genome autosomes.genome \
      --prop-min-elements ${min_prop_tracks_per_crb} \
      --neighborhood-dist ${clustering_neighborhood_dist} \
      --min-crb-separation ${min_crb_separation} \
      --crb-prefix "${prefix}_CRB" \
      --crb-outbed ${prefix}.crbs.bed.gz \
      --element-outbed ${prefix}.crb_elements.bed.gz \
      --bgzip \
      sig_tracks/*.curated.bed.gz

    # Index CRBs and elements
    tabix -f ${prefix}.crbs.bed.gz
    tabix -f ${prefix}.crb_elements.bed.gz

    # Copy final CRBs and elements to gs:// bucket
    gsutil -m cp \
      ${prefix}.crbs.bed.gz* \
      ${prefix}.crb_elements.bed.gz* \
      ${rCNV_bucket}/cleaned_data/genome_annotations/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:b0345a2b594dc51b0be4f75466be5f253bf6b024d707272ec9d44b6f5345912b"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }

  output {
    File final_crbs = "${prefix}.crbs.bed.gz"
    File final_crb_elements = "${prefix}.crb_elements.bed.gz"
  }
}

