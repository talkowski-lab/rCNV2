#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Curate all genome annotation tracks used for noncoding association testing


# Launch docker image
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global paramters
export rCNV_bucket="gs://rcnv_project"


# Download necessary reference files
mkdir refs/
gsutil -m cp ${rCNV_bucket}/refs/** refs/
gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/GRCh37.genome refs/
gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/rCNV_metacohort* refs/
mkdir cnvs/
gsutil -m cp ${rCNV_bucket}/cleaned_data/cnv/noncoding/** cnvs/


# Preprocess Roadmap Epigenomics ChromHMM states
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/all.mnemonics.bedFiles.tgz
mkdir roadmap_raw/
tar -xzvf all.mnemonics.bedFiles.tgz -C roadmap_raw/
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/labelmap_18_core_K27ac.tab
gsutil -m cp ${rCNV_bucket}/refs/REP_sample_manifest.tsv ./
mkdir chromhmm_beds/
/opt/rCNV2/data_curation/genome_annotations/preprocess_chromhmm.py \
  --state-manifest labelmap_18_core_K27ac.tab \
  --sample-manifest REP_sample_manifest.tsv \
  --prefix chromhmm_beds/roadmap_chromhmm \
  roadmap_raw/
# Copy ChromHMM beds (and tracklist) to gs:// bucket (note: requires permissions)
gsutil -m cp -r \
  chromhmm_beds \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find chromhmm_beds/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
> chromhmm_tracks.gs_paths.list
gsutil -m cp \
  chromhmm_tracks.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate all annotations in an arbitrary input list of paths
tracklist="test.annotations.list"
min_element_size=10
max_element_size=100000
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


# Compute CNV counts per annotation track per metacohort per CNV type
cat *.curated.stats.tsv | fgrep -v "#" | sort -Vk1,1 | uniq \
| awk -v FS="\t" -v OFS="\t" '{ print $0, $1".curated.bed.gz " }' \
| cat <( cat *.curated.stats.tsv | grep -e '^#' | sed -n '1p' \
         | awk -v OFS="\t" '{ print $0, "local_path" }' ) \
      - \
> all_tracks.stats.tsv
fgrep -v mega refs/rCNV_metacohort_list.txt \
| awk -v OFS="\t" '{ print $1, "cnvs/"$1".rCNV.strict_noncoding.bed.gz" }' \
> metacohorts.cnv_paths.tsv
/opt/rCNV2/data_curation/genome_annotations/count_cnvs_per_track.py \
  --cohorts metacohorts.cnv_paths.tsv \
  --track-stats all_tracks.stats.tsv \
  --outfile all_tracks.stats.with_counts.tsv.gz \
  --gzip 


# Burden meta-analysis of cohort CNV counts for each track
/opt/rCNV2/data_curation/genome_annotations/trackwise_cnv_burden_meta_analysis.R \
  --model "fe" \
  --spa \
  all_tracks.stats.with_counts.tsv.gz \
  all_tracks.burden_stats.tsv
gzip -f all_tracks.burden_stats.tsv

