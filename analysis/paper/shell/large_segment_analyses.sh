#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for secondary analyses of large rCNV-associated segments


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"


# Localize all analysis refs, sliding window meta-analysis stats, and large segment results
mkdir refs/ 
gsutil -m cp \
   ${rCNV_bucket}/analysis/analysis_refs/** \
   ${rCNV_bucket}/analysis/paper/data/hpo/${prefix}.reordered_hpos.txt \
   refs/
mkdir meta_stats/
gsutil -m cp \
  ${rCNV_bucket}/analysis/sliding_windows/**.rCNV.**.sliding_window.meta_analysis.stats.bed.gz \
  meta_stats/
gsutil -m cp \
  ${rCNV_bucket}/results/segment_association/* \
  ./
gsutil -m cp \
  ${rCNV_bucket}/analysis/paper/data/large_segments/clustered_nahr_regions.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs.*.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/wgs_common_cnvs.*.bed.gz \
  refs/


# Tabix all meta-analysis stats
find meta_stats/ -name "*meta_analysis.stats.bed.gz" | xargs -I {} tabix -f {}


# Compute effect size and max P-value per phenotype per final segment
/opt/rCNV2/analysis/paper/scripts/large_segments/calc_all_seg_stats.py \
  -o ${prefix}.final_segments.loci.all_sumstats.tsv \
  rCNV.final_segments.loci.bed.gz \
  refs/test_phenotypes.list \
  meta_stats
gzip -f ${prefix}.final_segments.loci.all_sumstats.tsv
gsutil -m cp \
  ${prefix}.final_segments.loci.all_sumstats.tsv.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/


# Build master BED of all regions, including final segments, known GDs, and 
# other predicted NAHR-mediated CNVs
cat <( echo -e "#chr\tstart\tend\tnahr_id\tcnv\tn_genes\tgenes" ) \
    <( zcat refs/clustered_nahr_regions.bed.gz | grep -ve '^#' \
       | awk -v FS="\t" -v OFS="\t" \
         '{ print $1, $2, $3, $4"_DEL", "DEL", $5, $6"\n"$1, $2, $3, $4"_DUP", "DUP", $5, $6 }' ) \
| bgzip -c \
> clustered_nahr_regions.reformatted.bed.gz
/opt/rCNV2/analysis/paper/scripts/large_segments/compile_segment_table.py \
  --final-loci rCNV.final_segments.loci.bed.gz \
  --hc-gds refs/lit_GDs.hc.bed.gz \
  --lc-gds refs/lit_GDs.lc.bed.gz \
  --nahr-cnvs clustered_nahr_regions.reformatted.bed.gz \
  --outfile ${prefix}.master_segments.bed.gz \
  --common-dels refs/wgs_common_cnvs.DEL.bed.gz \
  --common-dups refs/wgs_common_cnvs.DUP.bed.gz \
  --common-cnv-cov 0.5 \
  --gd-recip "10e-10" \
  --nahr-recip 0.2 \
  --bgzip


# Collapse overlapping DEL/DUP segments for sake of plotting
while read intervals rid; do
  echo "$intervals" | sed -e 's/\;/\n/g' -e 's/\:\|\-/\t/g' \
  | awk -v OFS="\t" -v rid=$rid '{ print $0, rid }'
done < <( zcat rCNV.final_segments.loci.bed.gz | grep -ve '^#' \
          | awk -v FS="\t" -v OFS="\t" '{ print $(NF-3), $4 }' ) \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| bedtools merge -i - -d 1000000 -c 4 -o distinct \
| cut -f4 \
> locus_clusters.txt


# Plot master grid summarizing segment association across all phenotypes
/opt/rCNV2/analysis/paper/plot/large_segments/plot_association_grid.R \
  --clusters locus_clusters.txt \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.final_segments.loci.all_sumstats.tsv.gz \
  refs/${prefix}.reordered_hpos.txt \
  ${prefix}.large_segments.association_grid.pdf
gsutil -m cp \
  ${prefix}.large_segments.association_grid.pdf \
  ${rCNV_bucket}/analysis/paper/plots/large_segments/


