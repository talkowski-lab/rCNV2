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


# Preprocess ENCODE TAD boundaries
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genome_annotations/manifests/encode.hic_tads.manifest.tsv.gz \
  ./
mkdir encode_tad_boundaries/
while read trackname path; do
  wget -O $trackname.raw.bed.gz $path
  bedtools flank \
    -i $trackname.raw.bed.gz \
    -g <( awk '{ print "chr"$0 }' refs/GRCh37.genome ) \
    -b 5000 \
  | cut -f1-3 \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bgzip -c \
  > encode_tad_boundaries/$trackname.tad_boundaries.bed.gz
done < <( zcat encode.hic_tads.manifest.tsv.gz \
          | awk -v FS="\t" -v OFS="\t" '{ print $1, $NF }' \
          | sed '1d' )
# Copy ENCODE beds (and tracklist) to gs:// bucket (note: requires permissions)
gsutil -m cp -r \
  encode_tad_boundaries \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find encode_tad_boundaries/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} -v OFS="\t" \
  '{ print gs"/cleaned_data/genome_annotations/"$1, "encode_tad_boundaries" }' \
> encode.tad_boundaries.gs_paths.list
gsutil -m cp \
  encode.tad_boundaries.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate HACER enhancers
mkdir hacer_enhancers
gsutil -m cp \
  ${rCNV_bucket}/raw_data/genome_annotations/HACER.GROseq.manifest.tsv \
  ${rCNV_bucket}/raw_data/genome_annotations/HACER.PROseq.manifest.tsv \
  ${rCNV_bucket}/raw_data/genome_annotations/HACER.CAGE.manifest.tsv \
  ./
while read sample url; do
  echo $sample
  if [ $( echo $url | grep -e '\.zip$' | wc -l ) -gt 0 ]; then
    wget $url
    if [ -e $sample/ ]; then
      rm -r $sample/
    fi
    mkdir $sample/
    unzip $( basename $url )  -d $sample
    for file in $sample/*; do
      tname=$( basename $file | sed 's/.txt//g' | cut -f2- -d\- | paste -s -d \- )
      sed '1d' $file | cut -f2-4 | bgzip -c > hacer_enhancers/$tname.bed.gz
    done
  else
    wget $url -q -O - | sed '1d' | cut -f2-4 | bgzip -c > hacer_enhancers/$sample.bed.gz
  fi
done < <( cat HACER.GROseq.manifest.tsv HACER.PROseq.manifest.tsv )
gsutil -m cp -r \
  hacer_enhancers \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find hacer_enhancers/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
| cat - <( awk -v FS="\t" '{ print $2 }' HACER.CAGE.manifest.tsv ) \
| awk -v OFS="\t" '{ print $1, "hacer_enh" }' \
> HACER.track_urls.list
gsutil -m cp \
  HACER.track_urls.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate super enhancers and typical enhancers per sample from SEdb
gsutil cp ${rCNV_bucket}/raw_data/genome_annotations/SEdb.sample_manifest.tsv ./
mkdir SEdb_super_enhancers/
/opt/rCNV2/data_curation/genome_annotations/preprocess_SEdb.py \
  --format 'SE' \
  --outdir SEdb_super_enhancers/ \
  http://www.licpathway.net/sedb/download/package/SE_package.bed \
  SEdb.sample_manifest.tsv
mkdir SEdb_enhancers
/opt/rCNV2/data_curation/genome_annotations/preprocess_SEdb.py \
  --format 'TE' \
  --outdir SEdb_enhancers/ \
  http://www.licpathway.net/sedb/download/package/TE_package.bed \
  SEdb.sample_manifest.tsv
gsutil -m cp -r \
  SEdb_super_enhancers \
  SEdb_enhancers \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find SEdb_super_enhancers/ -name "*bed.gz" \
|  awk -v gs=${rCNV_bucket} '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
> SEdb.track_gspaths.list
find SEdb_enhancers/ -name "*bed.gz" \
|  awk -v gs=${rCNV_bucket} '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
>> SEdb.track_gspaths.list
gsutil -m cp \
  SEdb.track_gspaths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate super enhancers from dbSUPER
wget http://asntech.org/dbsuper/data/bed/hg19/all_hg19_bed.zip
mkdir dbSUPER_super_enhancers/
unzip all_hg19_bed.zip
for file in all_hg19_bed/*.bed; do
  newname="$( basename "$file" | sed -e 's/\.bed//g' -e 's/\ /_/g' ).bed"
  mv "$file" dbSUPER_super_enhancers/$newname
  bgzip -f dbSUPER_super_enhancers/$newname
done
gsutil -m cp -r \
  dbSUPER_super_enhancers \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find dbSUPER_super_enhancers/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} -v OFS="\t" \
  '{ print gs"/cleaned_data/genome_annotations/"$1, "dbSUPER" }' \
> dbSUPER.gs_paths.list
gsutil -m cp \
  dbSUPER.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate VISTA validated enhancers
gsutil -m cp ${rCNV_bucket}/raw_data/genome_annotations/vista_raw.fa.gz ./
zcat vista_raw.fa.gz \
| fgrep positive | fgrep Human \
| awk -v FS="|" '{ print $2 }' \
| sed -e 's/:\|-/\t/g' \
| sort -Vk1,1 -k2,2n -k3,3n \
| bgzip -c \
> vista_enhancers.bed.gz
gsutil -m cp vista_enhancers.bed.gz \
  ${rCNV_bucket}/cleaned_data/genome_annotations/misc_tracks/


# Curate DENdb enhancers split by cell type
wget https://www.cbrc.kaust.edu.sa/dendb/src/enhancers.csv.zip
unzip enhancers.csv.zip
mkdir DENdb_enhancers
/opt/rCNV2/data_curation/genome_annotations/preprocess_DENdb.py \
  enhancers.csv \
  DENdb_enhancers/
gsutil -m cp -r \
  DENdb_enhancers \
  ${rCNV_bucket}/cleaned_data/genome_annotations/  
find DENdb_enhancers/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} -v OFS="\t" \
  '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
> DENdb.gs_paths.list
gsutil -m cp \
  DENdb.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate SEA super enhancers split by cell type
wget http://113.6.251.5:8080/SEA3/download/SEA00101.bed
awk -v FS="\t" -v OFS="\t" '{ print $2, $3, $4, $7 }' SEA00101.bed \
| sed 's/\ /_/g' | bgzip -c \
> SEA.hg38.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
liftOver -minMatch=0.75 -bedPlus=4 \
  SEA.hg38.gz \
  hg38ToHg19.over.chain.gz \
  SEA.hg19.bed \
  SEA.hg38_liftFailed.txt
mkdir SEA_super_enhancers/
/opt/rCNV2/data_curation/genome_annotations/preprocess_SEA.py \
  SEA.hg19.bed \
  SEA_super_enhancers/
gsutil -m cp -r \
  SEA_super_enhancers \
  ${rCNV_bucket}/cleaned_data/genome_annotations/  
find SEA_super_enhancers/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} -v OFS="\t" \
  '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
> SEA.gs_paths.list
gsutil -m cp \
  SEA.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate various PsychENCODE datasets
wget http://resource.psychencode.org/Datasets/Derived/DER-03a_hg19_PEC_enhancers.bed
sed '1d' DER-03a_hg19_PEC_enhancers.bed \
| awk -v OFS="\t" '{ print $1, $2, $3 }' \
| sort -Vk1,1 -k2,2n -k3,3n \
| bgzip -c \
> psychencode.enhancers.bed.gz
wget http://resource.psychencode.org/Datasets/Derived/DER-18_TAD_adultbrain.bed
gsutil -m cp psychencode.enhancers.bed.gz \
  ${rCNV_bucket}/cleaned_data/genome_annotations/misc_tracks/



# Combine all enhancer tracks into single tracklist for sharding on FireCloud
gsutil cat \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/EnhancerAtlas2.track_urls.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/HACER.track_urls.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/SEdb.track_gspaths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/dbSUPER.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/DENdb.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/SEA.gs_paths.list \
> enhancer_databases.track_urls.list
echo -e \
"${rCNV_bucket}/cleaned_data/genome_annotations/misc_tracks/vista_enhancers.bed.gz
${rCNV_bucket}/cleaned_data/genome_annotations/misc_tracks/psychencode.enhancers.bed.gz
http://resource.psychencode.org/Datasets/Derived/DER-03b_hg19_high_confidence_PEC_enhancers.bed\tpsychencode" \
>> enhancer_databases.track_urls.list
gsutil -m cp \
  enhancer_databases.track_urls.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Combine all miscellaneous tracks into single tracklist for sharding on FireCloud
gsutil cat \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/encode.tad_boundaries.gs_paths.list \
> misc_genome_annotations.track_urls.list
echo -e \
"http://resource.psychencode.org/Datasets/Derived/DER-05_PFC_H3K27ac_peaks.bed\tpsychencode
http://resource.psychencode.org/Datasets/Derived/DER-06_TC_H3K27ac_peaks.bed\tpsychencode
http://resource.psychencode.org/Datasets/Derived/DER-07_CBC_H3K27ac_peaks.bed\tpsychencode
http://resource.psychencode.org/Datasets/Pipeline/TARs/PIP-04_all_TARs.70pc.active.bed\tpsychencode"
gsutil -m cp \
  misc_genome_annotations.track_urls.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Development parameters for curate_annotations.wdl
prefix="all_tracks"
tracklist="test.annotations.list"
min_element_size=5
max_element_size=200000
case_hpo="HP:0000118"
min_element_overlap=1.0
p_cutoff=0.05
track_prefix="encode_tfbs"

# Curate all annotations in an arbitrary input list of paths
while IFS=$'\t' read path tprefix; do
  if [ -z $tprefix ]; then
    tprefix=${track_prefix}
  fi
  echo -e "Curating $path"
  /opt/rCNV2/data_curation/genome_annotations/curate_track.py \
    --genome refs/GRCh37.genome \
    --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
    --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
    --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
    --min-size ${min_element_size} \
    --max-size ${max_element_size} \
    --stats \
    --prefix "$tprefix" \
    "$path"
done < ${tracklist}
export IFS=$' \t\n'


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
  ${prefix}.burden_stats.tsv
gzip -f ${prefix}.burden_stats.tsv
cut -f1 ${prefix}.signif_paths_and_tracks.list \
> ${prefix}.signif_tracks.list
cut -f2 ${prefix}.signif_paths_and_tracks.list \
> ${prefix}.signif_tracknames.list


# # Dev code:
# min_prop_tracks_per_crb=0.1
# clustering_neighborhood_dist=10000
# min_crb_separation=10000


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
  --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
  --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
  --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
  --prop-min-elements ${min_prop_tracks_per_crb} \
  --neighborhood-dist ${clustering_neighborhood_dist} \
  --min-crb-separation ${min_crb_separation} \
  --crb-prefix "${prefix}_CRB" \
  --crb-outbed ${prefix}.crbs.bed.gz \
  --element-outbed ${prefix}.crb_elements.bed.gz \
  --bgzip \
  sig_tracks/*.curated.bed.gz

