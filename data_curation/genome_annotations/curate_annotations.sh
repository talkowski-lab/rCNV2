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
gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/HPOs_by_metacohort.table.tsv refs/
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


# Development parameters for curate_annotations.wdl
prefix="all_tracks"
tracklist="test.annotations.list"
min_element_size=5
max_element_size=200000
case_hpo="HP:0000118"
min_element_overlap=1.0
p_cutoff=0.05


# Curate all annotations in an arbitrary input list of paths
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
> ${prefix}.stats.tsv
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
  --outfile ${prefix}.stats.with_counts.tsv.gz \
  --gzip 


# # Dev code:
# gsutil -m cp \
#   gs://fc-cc4e446a-02a8-44af-bb70-2ca6013099b4/83202b3d-9437-41aa-a82b-2d0c933d02d7/curate_annotations/84695746-fc0d-4781-9a64-b2b71323ea07/call-merge_chromhmm/rCNV.chromhmm.merged_stats.with_counts.tsv.gz \
#   gs://fc-cc4e446a-02a8-44af-bb70-2ca6013099b4/83202b3d-9437-41aa-a82b-2d0c933d02d7/curate_annotations/84695746-fc0d-4781-9a64-b2b71323ea07/call-merge_tracklists/cacheCopy/rCNV.all_tracks.list \
#   ./
# stats=rCNV.chromhmm.merged_stats.with_counts.tsv.gz
# merged_tracklist=rCNV.all_tracks.list


# Burden meta-analysis of cohort CNV counts for each track
/opt/rCNV2/data_curation/genome_annotations/trackwise_cnv_burden_meta_analysis.R \
  --p-cutoff ${p_cutoff} \
  --signif-tracks ${prefix}.signif_paths_and_tracks.list \
  ${stats} \
  ${merged_tracklist} \
  ${prefix}.burden_stats.tsv
gzip -f ${prefix}.burden_stats.tsv
cut -f1 ${prefix}.signif_paths_and_tracks.list \
> ${prefix}.signif_tracks.list
cut -f2 ${prefix}.signif_paths_and_tracks.list \
> ${prefix}.signif_tracknames.list


# Dev code:
min_prop_tracks_per_crb=0.1
clustering_neighborhood_dist=10000
clustering_neighborhood_dist=10000


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

