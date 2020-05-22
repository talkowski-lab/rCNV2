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
  ${rCNV_bucket}/cleaned_data/binned_genome/GRCh37.200kb_bins_10kb_steps.raw.bed.gz \
  ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
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
  ${rCNV_bucket}/analysis/paper/data/large_segments/rCNV2_common_cnvs.*.bed.gz \
  refs/
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  ./
gsutil -m cp \
  ${rCNV_bucket}/analysis/paper/data/hpo/rCNV2_analysis_d1.hpo_jaccard_matrix.tsv \
  ./


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
af_suffix="01pct"
for CNV in DEL DUP; do
  zcat \
    refs/wgs_common_cnvs.$CNV.$af_suffix.bed.gz \
    refs/rCNV2_common_cnvs.$CNV.$af_suffix.bed.gz \
  | grep -ve '^#' \
  | cut -f1-3 \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i -\
  | bgzip -c \
  > combined_common_cnvs.$CNV.$af_suffix.bed.gz
done
TAB=$( printf '\t' )
cat << EOF > genelists_to_annotate.tsv
gnomAD_constrained${TAB}gene_lists/gnomad.v2.1.1.lof_constrained.genes.list
gnomAD_tolerant${TAB}gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list
CLinGen_HI${TAB}gene_lists/ClinGen.hmc_haploinsufficient.genes.list
CLinGen_TS${TAB}gene_lists/ClinGen.hmc_triplosensitive.genes.list
DECIPHER_LoF${TAB}gene_lists/DDG2P.hmc_lof.genes.list
DECIPHER_GoF${TAB}gene_lists/DDG2P.hmc_gof.genes.list
OMIM${TAB}gene_lists/HP0000118.HPOdb.genes.list
EOF
while read nocolon hpo; do
  echo -e "${hpo}\tgene_lists/${nocolon}.HPOdb.genes.list"
done < refs/test_phenotypes.list \
> hpo_genelists.tsv
/opt/rCNV2/analysis/paper/scripts/large_segments/compile_segment_table.py \
  --final-loci rCNV.final_segments.loci.bed.gz \
  --hc-gds refs/lit_GDs.hc.bed.gz \
  --mc-gds refs/lit_GDs.mc.bed.gz \
  --lc-gds refs/lit_GDs.lc.bed.gz \
  --nahr-cnvs clustered_nahr_regions.reformatted.bed.gz \
  --outfile ${prefix}.master_segments.bed.gz \
  --common-dels combined_common_cnvs.DEL.$af_suffix.bed.gz \
  --common-dups combined_common_cnvs.DUP.$af_suffix.bed.gz \
  --common-cnv-cov 0.5 \
  --hpo-jaccard-matrix rCNV2_analysis_d1.hpo_jaccard_matrix.tsv \
  --min-jaccard-sum 1.0 \
  --genelists genelists_to_annotate.tsv \
  --hpo-genelists hpo_genelists.tsv \
  --gd-recip "10e-10" \
  --nahr-recip 0.25 \
  --bgzip
gsutil -m cp \
  ${prefix}.master_segments.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/


# Run segment permutation tests
# Note: in practice, this is parallelized in the cloud using segment_permutation.wdl
# The code to execute these permutation tests is contained elsewhere
# Copy results of segment permutation tests (note: requires permissions)
n_seg_perms=10000
gsutil -m cp \
  ${rCNV_bucket}/analysis/paper/data/large_segments/permutations/${prefix}.${n_seg_perms}_permuted_segments.bed.gz \
  ./


# Plot segment permutation results
if [ -e perm_test_plots ]; then
  rm -rf perm_test_plots
fi
mkdir perm_test_plots
/opt/rCNV2/analysis/paper/plot/large_segments/plot_segment_permutations.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.master_segments.bed.gz \
  ${prefix}.${n_seg_perms}_permuted_segments.bed.gz \
  perm_test_plots/${prefix}


# # DEV: Run gene set permutation tests
# n_perms=50
# seed=51
# perm_prefix="gene_perm_test"
# # Download & format gene coordinates (note: requires permissions)
# gsutil -m cp \
#   ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
#   refs/
# zcat refs/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
# | cut -f1-4 \
# | bgzip -c \
# > refs/gencode.v19.canonical.pext_filtered.bed.gz
# # Permute gene blocks
# /opt/rCNV2/analysis/paper/scripts/gene_association/shuffle_gene_blocks.py \
#   --genome refs/GRCh37.autosomes.genome \
#   --n-perms ${n_perms} \
#   --first-seed ${seed} \
#   --outfile ${perm_prefix}.tsv.gz \
#   --gzip \
#   <( zcat rCNV.final_segments.loci.bed.gz \
#      | awk -v FS="\t" -v OFS="\t" '{ print $4, $5, $NF }' ) \
#   refs/gencode.v19.canonical.pext_filtered.bed.gz
# # Build necessary inputs for annotation
# echo -e "gnomAD_constrained\tgene_lists/gnomad.v2.1.1.lof_constrained.genes.list" > genelist_to_annotate.tsv
# echo -e "gnomAD_tolerant\tgene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list" >> genelist_to_annotate.tsv
# echo -e "CLinGen_HI\tgene_lists/ClinGen.hmc_haploinsufficient.genes.list" >> genelist_to_annotate.tsv
# echo -e "CLinGen_TS\tgene_lists/ClinGen.hmc_triplosensitive.genes.list" >> genelist_to_annotate.tsv
# echo -e "DECIPHER_LoF\tgene_lists/DDG2P.hmc_lof.genes.list" >> genelist_to_annotate.tsv
# echo -e "DECIPHER_GoF\tgene_lists/DDG2P.hmc_gof.genes.list" >> genelist_to_annotate.tsv
# echo -e "OMIM\tgene_lists/HP0000118.HPOdb.genes.list" >> genelist_to_annotate.tsv
# while read nocolon hpo; do
#   echo -e "$hpo\tgene_lists/$nocolon.HPOdb.genes.list"
# done < refs/test_phenotypes.list \
# > hpo_genelists.tsv
# zcat rCNV.final_segments.loci.bed.gz \
# | grep -ve '^#' \
# | awk -v FS="\t" -v OFS="\t" '{ print $4, $15 }' \
# > segment_hpos.tsv
# # Annotate
# /opt/rCNV2/analysis/paper/scripts/large_segments/annotate_shuffled_seg_gene_blocks.py \
#   --gene-sets genelists_to_annotate.forblocks.tsv \
#   --hpo-genelists hpo_genelists.tsv \
#   --segment-hpos segment_hpos.tsv \
#   --outfile ${perm_prefix}.annotated.tsv.gz \
#   --gzip \
#   ${perm_prefix}.tsv.gz
# Run segment permutation tests while matching on number of genes per segment
# Note: in practice, this is parallelized in the cloud using segment_permutation_bygene.wdl
# The code to execute these permutation tests is contained elsewhere
# Copy results of segment permutation tests (note: requires permissions)
n_seg_perms=10000
gsutil -m cp \
  ${rCNV_bucket}/analysis/paper/data/large_segments/permutations/${prefix}.${n_seg_perms}_permuted_segments_bygene.tsv.gz \
  ./


# Plot segment permutation results
if ! [ -e perm_test_plots ]; then
  mkdir perm_test_plots
fi
/opt/rCNV2/analysis/paper/plot/large_segments/plot_segment_bygene_perms.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.master_segments.bed.gz \
  ${prefix}.${n_seg_perms}_permuted_segments_bygene.tsv.gz \
  perm_test_plots/${prefix}


# Plot effect size covariates
if [ -e effect_size_plots ]; then
  rm -rf effect_size_plots
fi
mkdir effect_size_plots
/opt/rCNV2/analysis/paper/plot/large_segments/plot_effect_sizes.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.master_segments.bed.gz \
  effect_size_plots/${prefix}


# Plot pleiotropy covariates
if [ -e pleiotropy_plots ]; then
  rm -rf pleiotropy_plots
fi
mkdir pleiotropy_plots
/opt/rCNV2/analysis/paper/plot/large_segments/plot_pleiotropy.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.master_segments.bed.gz \
  pleiotropy_plots/${prefix}


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


