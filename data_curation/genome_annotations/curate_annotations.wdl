######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Curate genome annotation tracks, conduct global burden testing, and create cis-regulatory blocks

workflow curate_annotations {
	File chromhmm_tracklist
  File encode_dnaaccessibility_tracklist
  File encode_histone_mods_tracklist
  File encode_tfbs_tracklist
  File encode_transcription_tracklist
  File enhancer_databases_tracklist
  File misc_tracklist
  Int tracks_per_shard
  Int min_element_size
  Int max_element_size
  String case_hpo
  Float min_element_overlap
  Float p_cutoff
  File contiglist
  Float min_prop_tracks_per_crb
  Float min_prop_track_representation
  Int clustering_neighborhood_dist
  Int min_crb_separation
  Int max_crb_size
  Int crb_blacklist_buffer
  Int crb_blacklist_buffer_min_size
  String rCNV_bucket
  String prefix

  Array[Array[String]] contigs = read_tsv(contiglist)

  # Process Roadmap ChromHMM tracks
  call shard_tracklist as shard_tracklist_chromhmm {
    input:
      tracklist=chromhmm_tracklist,
      tracks_per_shard=tracks_per_shard,
      prefix="${prefix}.chromhmm_shards"
  }
  scatter ( shard in shard_tracklist_chromhmm.shards ) {
    call curate_and_burden as curate_and_burden_chromhmm {
      input:
        tracklist=shard,
        min_element_size=min_element_size,
        max_element_size=max_element_size,
        case_hpo=case_hpo,
        min_element_overlap=min_element_overlap,
        rCNV_bucket=rCNV_bucket,
        prefix="${prefix}.chromhmm_shard",
        track_prefix=""
    }
  }
  call merge_shards as merge_chromhmm {
    input:
      stat_shards=curate_and_burden_chromhmm.stats,
      prefix="${prefix}.chromhmm"
  }

  # Process ENCODE DNA accessibility tracks
  call shard_tracklist as shard_tracklist_encode_dnaaccessibility {
    input:
      tracklist=encode_dnaaccessibility_tracklist,
      tracks_per_shard=tracks_per_shard,
      prefix="${prefix}.encode_dnaaccessibility_shards"
  }
  scatter ( shard in shard_tracklist_encode_dnaaccessibility.shards ) {
    call curate_and_burden as curate_and_burden_encode_dnaaccessibility {
      input:
        tracklist=shard,
        min_element_size=min_element_size,
        max_element_size=max_element_size,
        case_hpo=case_hpo,
        min_element_overlap=min_element_overlap,
        rCNV_bucket=rCNV_bucket,
        prefix="${prefix}.encode_dnaaccessibility_shard",
        track_prefix="encode_dnaaccessibility"
    }
  }
  call merge_shards as merge_encode_dnaaccessibility {
    input:
      stat_shards=curate_and_burden_encode_dnaaccessibility.stats,
      prefix="${prefix}.encode_dnaaccessibility"
  }

  # Process ENCODE histone modification tracks
  call shard_tracklist as shard_tracklist_encode_histone_mods {
    input:
      tracklist=encode_histone_mods_tracklist,
      tracks_per_shard=tracks_per_shard,
      prefix="${prefix}.encode_histone_mods_shards"
  }
  scatter ( shard in shard_tracklist_encode_histone_mods.shards ) {
    call curate_and_burden as curate_and_burden_encode_histone_mods {
      input:
        tracklist=shard,
        min_element_size=min_element_size,
        max_element_size=max_element_size,
        case_hpo=case_hpo,
        min_element_overlap=min_element_overlap,
        rCNV_bucket=rCNV_bucket,
        prefix="${prefix}.encode_histone_mods_shard",
        track_prefix="encode_histone_mods"
    }
  }
  call merge_shards as merge_encode_histone_mods {
    input:
      stat_shards=curate_and_burden_encode_histone_mods.stats,
      prefix="${prefix}.encode_histone_mods"
  }

  # Process ENCODE TFBS tracks
  call shard_tracklist as shard_tracklist_encode_tfbs {
    input:
      tracklist=encode_tfbs_tracklist,
      tracks_per_shard=tracks_per_shard,
      prefix="${prefix}.encode_tfbs_shards"
  }
  scatter ( shard in shard_tracklist_encode_tfbs.shards ) {
    call curate_and_burden as curate_and_burden_encode_tfbs {
      input:
        tracklist=shard,
        min_element_size=min_element_size,
        max_element_size=max_element_size,
        case_hpo=case_hpo,
        min_element_overlap=min_element_overlap,
        rCNV_bucket=rCNV_bucket,
        prefix="${prefix}.encode_tfbs_shard",
        track_prefix="encode_tfbs"
    }
  }
  call merge_shards as merge_encode_tfbs {
    input:
      stat_shards=curate_and_burden_encode_tfbs.stats,
      prefix="${prefix}.encode_tfbs"
  }

  # Process ENCODE transcription tracks
  call shard_tracklist as shard_tracklist_encode_transcription {
    input:
      tracklist=encode_transcription_tracklist,
      tracks_per_shard=tracks_per_shard,
      prefix="${prefix}.encode_transcription_shards"
  }
  scatter ( shard in shard_tracklist_encode_transcription.shards ) {
    call curate_and_burden as curate_and_burden_encode_transcription {
      input:
        tracklist=shard,
        min_element_size=min_element_size,
        max_element_size=max_element_size,
        case_hpo=case_hpo,
        min_element_overlap=min_element_overlap,
        rCNV_bucket=rCNV_bucket,
        prefix="${prefix}.encode_transcription_shard",
        track_prefix="encode_transcription"
    }
  }
  call merge_shards as merge_encode_transcription {
    input:
      stat_shards=curate_and_burden_encode_transcription.stats,
      prefix="${prefix}.encode_transcription"
  }

  # Process various enhancer database tracks
  call shard_tracklist as shard_tracklist_enhancer_databases {
    input:
      tracklist=enhancer_databases_tracklist,
      tracks_per_shard=tracks_per_shard,
      prefix="${prefix}.enhancer_databases_shards"
  }
  scatter ( shard in shard_tracklist_enhancer_databases.shards ) {
    call curate_and_burden as curate_and_burden_enhancer_databases {
      input:
        tracklist=shard,
        min_element_size=min_element_size,
        max_element_size=max_element_size,
        case_hpo=case_hpo,
        min_element_overlap=min_element_overlap,
        rCNV_bucket=rCNV_bucket,
        prefix="${prefix}.enhancer_databases_shard",
        track_prefix="enhancer_databases"
    }
  }
  call merge_shards as merge_enhancer_databases {
    input:
      stat_shards=curate_and_burden_enhancer_databases.stats,
      prefix="${prefix}.enhancer_databases"
  }

  # Process miscellaneous tracks
  call shard_tracklist as shard_tracklist_misc {
    input:
      tracklist=misc_tracklist,
      tracks_per_shard=tracks_per_shard,
      prefix="${prefix}.misc_shards"
  }
  scatter ( shard in shard_tracklist_misc.shards ) {
    call curate_and_burden as curate_and_burden_misc {
      input:
        tracklist=shard,
        min_element_size=min_element_size,
        max_element_size=max_element_size,
        case_hpo=case_hpo,
        min_element_overlap=min_element_overlap,
        rCNV_bucket=rCNV_bucket,
        prefix="${prefix}.misc_shard",
        track_prefix="misc"
    }
  }
  call merge_shards as merge_misc {
    input:
      stat_shards=curate_and_burden_misc.stats,
      prefix="${prefix}.misc"
  }

  # Merge stats across all tracklists
  call merge_shards as merge_all {
    input:
      stat_shards=[merge_chromhmm.merged_stats, merge_encode_dnaaccessibility.merged_stats, merge_encode_histone_mods.merged_stats, merge_encode_tfbs.merged_stats, merge_encode_transcription.merged_stats, merge_enhancer_databases.merged_stats, merge_misc.merged_stats],
      prefix="${prefix}.all"
  }

  # Merge stats across all tracklists and compute meta-analysis burden stats
  call meta_burden_test {
    input:
      stats=merge_all.merged_stats,
      p_cutoff=p_cutoff,
      rCNV_bucket=rCNV_bucket,
      prefix=prefix,
      clear_sig="TRUE"
  }

  # Re-shard all significant tracks for final curation
  call shard_tracklist as shard_tracklist_signif {
    input:
      tracklist=meta_burden_test.signif_tracks,
      tracks_per_shard=tracks_per_shard,
      prefix="${prefix}.signif_tracks"
  }

  # Scatter over sharded significant tracklist and curate tracks
  scatter ( shard in shard_tracklist_signif.shards ) {
    call curate_only as curate_only_signif {
      input:
        tracknames_and_paths=shard,
        min_element_size=min_element_size,
        max_element_size=max_element_size,
        rCNV_bucket=rCNV_bucket,
        prefix="${prefix}.signif_shard"
    }
  }

  # Cluster CRBs, parallelized per chromosome
  scatter ( contig in contigs ) {
    call cluster_elements {
    input:
      completion_tokens=curate_only_signif.completion_token,
      min_prop_tracks_per_crb=min_prop_tracks_per_crb,
      min_prop_track_representation=min_prop_track_representation,
      clustering_neighborhood_dist=clustering_neighborhood_dist,
      min_crb_separation=min_crb_separation,
      max_crb_size=max_crb_size,
      blacklist_buffer=crb_blacklist_buffer,
      blacklist_buffer_min_size=crb_blacklist_buffer_min_size,
      contig=contig[0],
      rCNV_bucket=rCNV_bucket,
      prefix=prefix
    }
  }
  call merge_final_beds as merge_crbs {
    input:
      beds=cluster_elements.final_crbs,
      outfile_prefix="${prefix}.crbs",
      out_bucket="${rCNV_bucket}/cleaned_data/genome_annotations/"
  }
  call merge_final_beds as merge_crb_elements {
    input:
      beds=cluster_elements.final_crb_elements,
      outfile_prefix="${prefix}.crb_elements",
      out_bucket="${rCNV_bucket}/cleaned_data/genome_annotations/"
  }
  call merge_final_beds as merge_crb_whitelists {
    input:
      beds=cluster_elements.crb_whitelist,
      outfile_prefix="${prefix}.crb_whitelist",
      out_bucket="${rCNV_bucket}/cleaned_data/genome_annotations/"
  }
  
  
  # Final outputs
  output {
    File final_crbs = merge_crbs.merged_bed
    File final_crb_elements = merge_crb_elements.merged_bed
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
    docker: "talkowski/rcnv@sha256:90f47bc9fe55578e830ec617328188267e5449104a3b0b8d36fdf9d9a9fd37a0"
    preemptible: 1
  }

  output {
    Array[File] shards = glob("${prefix}*")
  }
}


# Download & curate tracks (no burden test)
task curate_only {
  File tracknames_and_paths
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
    while read trackname path; do
      echo -e "Curating $path"
      /opt/rCNV2/data_curation/genome_annotations/curate_track.py \
        --genome refs/GRCh37.genome \
        --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
        --blacklist refs/GRCh37.somatic_hypermutable_sites.200kb_clustered.bed.gz \
        --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
        --min-size ${min_element_size} \
        --max-size ${max_element_size} \
        --stats \
        --prefix "$tprefix" \
        "$path"
    done < ${tracknames_and_paths}

    # Copy all curated tracks to final gs:// bucket
    gsutil -m cp \
      *.curated.bed.gz \
      ${rCNV_bucket}/cleaned_data/genome_annotations/significant_tracks/

    # Dummy output
    touch completion.txt
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:90f47bc9fe55578e830ec617328188267e5449104a3b0b8d36fdf9d9a9fd37a0"
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
  String track_prefix

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

    # Curate all tracks in tracklist, using second column as track-specific prefix (if provided)
    while IFS=$'\t' read path tprefix; do
      if [ -z $tprefix ]; then
        tprefix=${track_prefix}
      fi
      echo -e "Curating $path"
      /opt/rCNV2/data_curation/genome_annotations/curate_track.py \
        --genome refs/GRCh37.genome \
        --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
        --blacklist refs/GRCh37.somatic_hypermutable_sites.200kb_clustered.bed.gz \
        --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
        --min-size ${min_element_size} \
        --max-size ${max_element_size} \
        --stats \
        --prefix "$tprefix" \
        "$path"
    done < ${tracklist}
    export IFS=$' \t\n'

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
        echo -e "cnvs/$cohort.rCNV.loose_noncoding.bed.gz"
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
    docker: "talkowski/rcnv@sha256:90f47bc9fe55578e830ec617328188267e5449104a3b0b8d36fdf9d9a9fd37a0"
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
    docker: "talkowski/rcnv@sha256:90f47bc9fe55578e830ec617328188267e5449104a3b0b8d36fdf9d9a9fd37a0"
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
    docker: "talkowski/rcnv@sha256:90f47bc9fe55578e830ec617328188267e5449104a3b0b8d36fdf9d9a9fd37a0"
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
      --cutoff ${p_cutoff} \
      --signif-tracks ${prefix}.signif_paths_and_tracks.list \
      --volcano ${prefix}.track_volcanos.png \
      ${stats} \
      ${prefix}.burden_stats.tsv
    gzip -f ${prefix}.burden_stats.tsv

    # Copy final stats to gs:// bucket
    gsutil -m cp \
      ${prefix}.burden_stats.tsv.gz \
      ${rCNV_bucket}/cleaned_data/genome_annotations/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:90f47bc9fe55578e830ec617328188267e5449104a3b0b8d36fdf9d9a9fd37a0"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }

  output {
    File meta_stats = "${prefix}.burden_stats.tsv.gz"
    File signif_tracks = "${prefix}.signif_paths_and_tracks.list"
    File volcano_plots = "${prefix}.track_volcanos.png"
  }
}


# Cluster elements from significant tracks into CRBs (per chromosome)
task cluster_elements {
  Array[File] completion_tokens
  Float min_prop_tracks_per_crb
  Float min_prop_track_representation
  Int clustering_neighborhood_dist
  Int min_crb_separation
  Int max_crb_size
  Int blacklist_buffer
  Int blacklist_buffer_min_size
  String contig
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Copy necessary reference files
    mkdir refs/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/analysis_refs/GRCh37.genome \
      ${rCNV_bucket}/refs/*.bed.gz \
      refs/

    # Copy & index all final curated significant tracks locally
    mkdir sig_tracks/
    gsutil -m cp \
      ${rCNV_bucket}/cleaned_data/genome_annotations/significant_tracks/*.curated.bed.gz \
      sig_tracks/
    find sig_tracks/ -name "*.curated.bed.gz" \
    | xargs -I {} tabix -f {}

    # Subset genome file to contig of interest
    awk -v FS="\t" -v OFS="\t" -v contig=${contig} \
      '{ if ($1==contig) print $0 }' refs/GRCh37.genome \
    > ${contig}.genome

    # Create whitelist for CRB clustering
    gsutil -m cp \
      ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
      ${rCNV_bucket}/cleaned_data/genes/gene_lists/gnomad.v2.1.1.likely_unconstrained.genes.list \
      ./
    tabix gencode.v19.canonical.pext_filtered.gtf.gz ${contig} \
    | gzip -c \
    > ${contig}.gtf.gz
    /opt/rCNV2/data_curation/genome_annotations/make_crb_whitelist.py \
      --bgzip \
      ${contig}.gtf.gz \
      gnomad.v2.1.1.likely_unconstrained.genes.list \
      ${contig}.genome \
      ${contig}.crb_whitelist.bed.gz

    # Cluster significant tracks into CRBs
    /opt/rCNV2/data_curation/genome_annotations/build_crbs.py \
      --genome ${contig}.genome \
      --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
      --blacklist refs/GRCh37.somatic_hypermutable_sites.200kb_clustered.bed.gz \
      --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
      --blacklist-buffer ${blacklist_buffer} \
      --blacklist-buffer-min-size ${blacklist_buffer_min_size} \
      --whitelist ${contig}.crb_whitelist.bed.gz \
      --prop-min-elements ${min_prop_tracks_per_crb} \
      --prop-min-tracks ${min_prop_track_representation} \
      --neighborhood-dist ${clustering_neighborhood_dist} \
      --min-crb-separation ${min_crb_separation} \
      --max-crb-size ${max_crb_size} \
      --crb-prefix "${prefix}_CRB" \
      --crb-outbed ${prefix}.${contig}.crbs.bed.gz \
      --element-outbed ${prefix}.${contig}.crb_elements.bed.gz \
      --bgzip \
      sig_tracks/*.curated.bed.gz

    # Index CRBs and elements
    tabix -f ${prefix}.${contig}.crbs.bed.gz
    tabix -f ${prefix}.${contig}.crb_elements.bed.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:90f47bc9fe55578e830ec617328188267e5449104a3b0b8d36fdf9d9a9fd37a0"
    preemptible: 1
    memory: "15 GB"
    bootDiskSizeGb: "30"
    disks: "local-disk 100 HDD"
  }

  output {
    File final_crbs = "${prefix}.${contig}.crbs.bed.gz"
    File final_crb_elements = "${prefix}.${contig}.crb_elements.bed.gz"
    File crb_whitelist = "${contig}.crb_whitelist.bed.gz"
  }
}


# Merge array of BEDs (with header) and copy to gs:// bucket
task merge_final_beds {
  Array[File] beds
  String outfile_prefix
  String out_bucket

  command <<<
    set -e

    # Merge all stats
    zcat ${beds[0]} | sed -n '1p' | grep -e '^#' > ${outfile_prefix}.bed || true
    zcat ${sep=" " beds} | grep -ve '^#' \
    >> ${outfile_prefix}.bed
    bgzip -f ${outfile_prefix}.bed
    tabix -f ${outfile_prefix}.bed.gz

    # Copy merged BED gs:// bucket
    gsutil -m cp \
      ${outfile_prefix}.bed.gz* \
      ${out_bucket}/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:90f47bc9fe55578e830ec617328188267e5449104a3b0b8d36fdf9d9a9fd37a0"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 50 HDD"
  }

  output {
    File merged_bed = "${outfile_prefix}.bed.gz"
  }
}

