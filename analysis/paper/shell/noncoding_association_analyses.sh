#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for secondary analyses of noncoding association for rCNV2 manuscript


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d2"


# Download necessary data (note: requires permissions)
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/cnv ./
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genome_annotations/rCNV.burden_stats.tsv.gz \
  ${rCNV_bucket}/cleaned_data/genome_annotations/rCNV.crbs.bed.gz \
  ${rCNV_bucket}/analysis/crb_burden/other_data/unconstrained_cnv_counts.*.tsv.gz \
  ${rCNV_bucket}/results/segment_association/* \
  ./
mkdir refs/
gsutil -m cp -r \
  ${rCNV_bucket}/analysis/analysis_refs/* \
  ${rCNV_bucket}/refs/REP_state_manifest.tsv \
  refs/
mkdir meta_stats/
gsutil -m cp \
  gs://rcnv_project/analysis/crb_burden/**rCNV.loose_noncoding.*.crb_burden.meta_analysis.stats.bed.gz \
  meta_stats/


# Plot observed effect size distributions split by gene set membership
# Used to justify inclusion of some coding effects in noncoding association test
total_dev=$( fgrep -v "#" refs/rCNV2.hpos_by_severity.developmental.counts.tsv \
             | awk '{ sum+=$2 }END{ print sum }' )
echo -e "DEVELOPMENTAL\tStrong-effect HPOs\t$total_dev" \
| paste - <( fgrep -v "#" refs/rCNV2.hpos_by_severity.developmental.counts.tsv \
             | cut -f2 | paste -s ) \
| paste - <( echo $total_dev ) | cat refs/HPOs_by_metacohort.table.tsv - \
> HPOs_by_metacohort.w_DEV.table.tsv
/opt/rCNV2/analysis/paper/plot/noncoding_association/plot_unconstrained_effect_sizes.R \
  unconstrained_cnv_counts.DEL.tsv.gz \
  unconstrained_cnv_counts.DUP.tsv.gz \
  HPOs_by_metacohort.w_DEV.table.tsv \
  ${prefix}


# Make supplementary table of noncoding track stats
/opt/rCNV2/analysis/paper/scripts/noncoding_association/format_track_stats_table.R \
  --saddlepoint-adj \
  rCNV.burden_stats.tsv.gz \
  ${prefix}


# Plot annotation track distributions & stats
/opt/rCNV2/analysis/paper/plot/noncoding_association/plot_track_stats.R \
  rCNV.burden_stats.tsv.gz \
  ${prefix}


# Volcano plots of DEL & DUP track-level burden tests
/opt/rCNV2/analysis/paper/plot/noncoding_association/plot_volcanos.R \
  --saddlepoint-adj \
  rCNV.burden_stats.tsv.gz \
  ${prefix}


# Plot ChromHMM enrichments as positive controls
/opt/rCNV2/analysis/paper/plot/noncoding_association/chromhmm_enrichments.R \
  --saddlepoint-adj \
  rCNV.burden_stats.tsv.gz \
  refs/REP_state_manifest.tsv \
  ${prefix}


# Plot CRB distributions & stats
/opt/rCNV2/analysis/paper/plot/noncoding_association/plot_crb_stats.R \
  rCNV.crbs.bed.gz \
  ${prefix}


# Get all significant CRBs
for cnv in DEL DUP; do
  # Make input for get_significant_genes_v2.py
  while read nocolon hpo; do
    echo $hpo
    echo -e "meta_stats/$nocolon.rCNV.loose_noncoding.$cnv.crb_burden.meta_analysis.stats.bed.gz"
    head -n1 refs/crb_burden.rCNV.loose_noncoding.DEL.bonferroni_pval.hpo_cutoffs.tsv | cut -f2
  done < refs/test_phenotypes.list | paste - - - > $cnv.stats.list

  # Extract all significant CRBs
  # Note: can use same script as for extracting significant genes because format of sumstats is identical
  /opt/rCNV2/analysis/genes/get_significant_genes_v2.py \
    $cnv.stats.list \
    --secondary-p-cutoff 0.05 \
    --min-nominal 2 \
    --secondary-or-nominal \
    --fdr-q-cutoff 0.01 \
    --secondary-for-fdr \
    --outfile $cnv.all_sig_crbs.tsv
done


# Check if any CRBs don't overlap with significant large segments
for CNV in DEL DUP; do
  zcat rCNV.crbs.bed.gz | fgrep -v "#" \
  | awk -v OFS="\t" '{ print $4, $1, $2, $3 }' | sort -k1,1 \
  | join -j 1 -t $'\t' <( cat $CNV.all_sig_crbs.tsv | fgrep -v "#" | sort -k1,1 ) - \
  | awk -v OFS="\t" -v CNV=$CNV '{ print $3, $4, $5, $1, $2, CNV }' \
  # | bedtools intersect -v -a - \
  #   -b <( zcat rCNV.final_segments.loci.bed.gz | fgrep -w $CNV )
done

# Copy all plots to final gs:// directory
gsutil -m cp -r \
  ${prefix}.cnv_lnORs_by_genic_context.pdf \
  ${prefix}.annotation_burden_stats.tsv.gz \
  ${prefix}.track_stats.*.pdf \
  ${prefix}.*volcano*pdf \
  ${prefix}.chromhmm_enrichment.*.pdf \
  ${prefix}.crb_stats.*.pdf \
  ${rCNV_bucket}/analysis/paper/plots/noncoding_association/
