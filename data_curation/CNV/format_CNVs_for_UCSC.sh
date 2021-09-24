#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Reformat all curated CNV callsets as tracks appropriate for UCSC Genome Browser


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv


# Copy all raw and filtered CNV data from Google Bucket (note: requires permissions)
gcloud auth login
mkdir ./cleaned_cnv
gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* ./cleaned_cnv/
mkdir refs/
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/


# Function to write track header
function track_header {
  cohort=$1
  pheno=$2
  CNV=$3
  color=$4

  echo -e "track name=\"$cohort $pheno $CNV\" visibility=4 color=$color itemRgb=On"
  echo -e "#chrom chromStart chromEnd name"
}


# Format each metacohort as a separate track, but store them in the same file
for freq in rCNV; do
  for CNV in DEL DUP; do
    case $CNV in
      DEL)
        case_color="200,0,0,"
        control_color="150,85,85,"
        ;;
      DUP)
        case_color="0,0,200,"
        control_color="85,85,150,"
        ;;
    esac
    for wrapper in 1; do
      # Case CNVs
      while read meta cohorts; do
        track_header $meta case $CNV $case_color
        zcat cleaned_cnv/$meta.$freq.bed.gz \
        | sed '1d' \
        | awk -v OFS="\t" -v CNV=$CNV '{ if ($NF !~ "HEALTHY_CONTROL" && $5==CNV) \
            print "chr"$1, $2, $3, $4"__"$NF }' \
        | bedtools subtract \
          -a - \
          -b <( awk -v OFS="\t" '{ print "chr"$1, $2, 1000000000 }' refs/GRCh37.genome )
      done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt )
      # Control CNVs
      while read meta cohorts; do
        track_header $meta control $CNV $control_color
        zcat cleaned_cnv/$meta.$freq.bed.gz \
        | sed '1d' \
        | awk -v OFS="\t" -v CNV=$CNV '{ if ($NF=="HEALTHY_CONTROL" && $5==CNV) \
            print "chr"$1, $2, $3, $4"__"$NF }' \
        | bedtools subtract \
          -a - \
          -b <( awk -v OFS="\t" '{ print "chr"$1, $2, 1000000000 }' refs/GRCh37.genome )
      done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt )
    done \
    | gzip -c \
    > ${freq}_UCSC.all_metacohorts.$CNV.tsv.gz
  done
done


# Copy all UCSC tracks to Google Bucket (note: requires permissions)
gsutil -m cp ./*_UCSC.*.tsv.gz \
  gs://rcnv_project/cleaned_data/ucsc_cnv/

