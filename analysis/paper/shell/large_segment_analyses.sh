#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for secondary analyses of large rCNV-associated segments


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d2"


# Localize all analysis refs, sliding window meta-analysis stats, and large segment results
mkdir refs/ 
gsutil -m cp \
  ${rCNV_bucket}/analysis/analysis_refs/** \
  ${rCNV_bucket}/analysis/paper/data/hpo/${prefix}.reordered_hpos.txt \
  ${rCNV_bucket}/cleaned_data/binned_genome/GRCh37.200kb_bins_10kb_steps.raw.bed.gz \
  ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  ${rCNV_bucket}/analysis/paper/data/misc/*_dnm_counts*tsv.gz \
  ${rCNV_bucket}/analysis/paper/data/misc/*.gw_sig.genes.list \
  ${rCNV_bucket}/analysis/paper/data/misc/redin_bca_breakpoints.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/misc/gene_mutation_rates.tsv.gz \
  ${rCNV_bucket}/cleaned_data/genes/annotations/gtex_stats/gencode.v19.canonical.pext_filtered.GTEx_v7_expression_stats.median.tsv.gz \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  ${rCNV_bucket}/cleaned_data/cnv/mega.rCNV.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/clustered_nahr_regions.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/loose_unclustered_nahr_regions.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/${prefix}.correct_nahr_labels.tsv \
  ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs.*.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/wgs_common_cnvs.*.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/rCNV2_common_cnvs.*.bed.gz \
  refs/
wget -P refs/ https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
mkdir meta_stats/
gsutil -m cp \
  ${rCNV_bucket}/analysis/sliding_windows/**.rCNV.**.sliding_window.meta_analysis.stats.bed.gz \
  meta_stats/
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  ${rCNV_bucket}/analysis/paper/data/hpo/rCNV2_analysis_d2.hpo_jaccard_matrix.tsv \
  ${rCNV_bucket}/results/segment_association/* \
  ./


# Tabix all meta-analysis stats
find meta_stats/ -name "*meta_analysis.stats.bed.gz" | xargs -I {} tabix -f {}


# Get genome-wide significance thresholds
example_hpo="HP0012759"
del_cutoff=$( awk -v FS="\t" -v hpo=${example_hpo} '{ if ($1==hpo) print $2 }' \
              refs/sliding_window.rCNV.DEL.bonferroni_pval.hpo_cutoffs.tsv )
dup_cutoff=$( awk -v FS="\t" -v hpo=${example_hpo} '{ if ($1==hpo) print $2 }' \
              refs/sliding_window.rCNV.DUP.bonferroni_pval.hpo_cutoffs.tsv )


# Compute effect size and max P-value per phenotype per gw-sig segment
/opt/rCNV2/analysis/paper/scripts/large_segments/calc_all_seg_stats.py \
  -o ${prefix}.final_segments.loci.all_sumstats.tsv \
  rCNV.final_segments.loci.bed.gz \
  refs/test_phenotypes.list \
  meta_stats
gzip -f ${prefix}.final_segments.loci.all_sumstats.tsv
gsutil -m cp \
  ${prefix}.final_segments.loci.all_sumstats.tsv.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/


# Compute effect size and max P-value per phenotype per literature-based segment
zcat refs/lit_GDs.*.bed.gz | fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n \
| awk -v OFS="\t" '{ print $0, $1":"$2"-"$3 }' \
| cat <( zcat refs/lit_GDs.hc.bed.gz | grep -e '^#' | paste - <( echo "cred_interval_coords") ) - \
| bgzip -c > all_gds.bed.gz
/opt/rCNV2/analysis/paper/scripts/large_segments/calc_all_seg_stats.py \
  -o ${prefix}.lit_gds.all_sumstats.tsv \
  all_gds.bed.gz \
  refs/test_phenotypes.list \
  meta_stats
gzip -f ${prefix}.lit_gds.all_sumstats.tsv
gsutil -m cp \
  ${prefix}.lit_gds.all_sumstats.tsv.gz \
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
cat << EOF > dnm_counts_to_annotate.tsv
ASC${TAB}refs/fu_asc_spark_dnm_counts.tsv.gz${TAB}refs/ASC_2021.gw_sig.genes.list
ASC_unaffected${TAB}refs/fu_asc_spark_dnm_counts.unaffecteds.tsv.gz${TAB}refs/ASC_2021.gw_sig.genes.list
DDD${TAB}refs/ddd_dnm_counts.tsv.gz${TAB}refs/DDD_2020.gw_sig.genes.list
EOF
cat \
  <( zcat ${prefix}.final_segments.loci.all_sumstats.tsv.gz ) \
  <( zcat ${prefix}.lit_gds.all_sumstats.tsv.gz | grep -ve '^region_id' ) \
> pooled_sumstats.tsv
/opt/rCNV2/analysis/paper/scripts/large_segments/compile_segment_table.py \
  --final-loci rCNV.final_segments.loci.bed.gz \
  --final-associations rCNV.final_segments.associations.bed.gz \
  --hc-gds refs/lit_GDs.hc.bed.gz \
  --mc-gds refs/lit_GDs.mc.bed.gz \
  --lc-gds refs/lit_GDs.lc.bed.gz \
  --nahr-cnvs clustered_nahr_regions.reformatted.bed.gz \
  --outfile ${prefix}.master_segments.pre_nahr_polishing.bed.gz \
  --loose-nahr-mask refs/loose_unclustered_nahr_regions.bed.gz \
  --cnv-bed refs/mega.rCNV.bed.gz \
  --genome-file refs/GRCh37.genome \
  --common-dels combined_common_cnvs.DEL.$af_suffix.bed.gz \
  --common-dups combined_common_cnvs.DUP.$af_suffix.bed.gz \
  --common-cnv-cov 0.5 \
  --hpo-jaccard-matrix rCNV2_analysis_d2.hpo_jaccard_matrix.tsv \
  --min-jaccard-sum 1.0 \
  --genelists genelists_to_annotate.tsv \
  --hpo-genelists hpo_genelists.tsv \
  --gnomad-constraint-tsv refs/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz  \
  --dnm-tsvs dnm_counts_to_annotate.tsv \
  --snv-mus refs/gene_mutation_rates.tsv.gz \
  --bca-tsv refs/redin_bca_breakpoints.bed.gz \
  --gtex-matrix refs/gencode.v19.canonical.pext_filtered.GTEx_v7_expression_stats.median.tsv.gz \
  --meta-sumstats pooled_sumstats.tsv \
  --neuro-hpos refs/neuro_hpos.list \
  --dev-hpos refs/rCNV2.hpos_by_severity.developmental.list \
  --gd-recip "10e-10" \
  --nahr-recip 0.5 \
  --bgzip
# Polish NAHR labels following manual review (note: requires refs/${prefix}.correct_NAHR_labels.tsv)
/opt/rCNV2/analysis/paper/scripts/large_segments/polish_nahr_labels.R \
  ${prefix}.master_segments.pre_nahr_polishing.bed.gz \
  refs/${prefix}.correct_nahr_labels.tsv \
  ${prefix}.master_segments.bed
bgzip -f ${prefix}.master_segments.bed
gsutil -m cp \
  ${prefix}.master_segments.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/


# Plot basic segment distributions
if [ -e basic_distribs ]; then
  rm -rf basic_distribs
fi
mkdir basic_distribs
/opt/rCNV2/analysis/paper/plot/large_segments/plot_basic_segment_distribs.R \
  --gw-sig $del_cutoff \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.master_segments.bed.gz \
  basic_distribs/${prefix}
/opt/rCNV2/analysis/paper/plot/large_segments/plot_seg_attributes_vs_ngenes.R \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.master_segments.bed.gz \
  basic_distribs/${prefix}
# Get range of effect sizes across all g-w sig associations
zcat rCNV.final_segments.associations.bed.gz | fgrep -v "#" | cut -f10 | sort -nk1,1 | head -n1
zcat rCNV.final_segments.associations.bed.gz | fgrep -v "#" | cut -f11 | sort -nk1,1 | tail -n1
zcat rCNV.final_segments.associations.bed.gz | fgrep -v "#" | cut -f11 \
| awk '{ sum+=$1 }END{ print sum/NR }'


# Run segment permutation tests
# Note: in practice, this is parallelized in the cloud using segment_permutation.wdl
# The code to execute these permutation tests is contained elsewhere
# Copy results of segment permutation tests (note: requires permissions)
n_seg_perms=100000
gsutil -m cp \
  ${rCNV_bucket}/analysis/paper/data/large_segments/permutations/${prefix}*${n_seg_perms}_permuted_segments.bed.gz \
  ./


# Plot segment permutation results
# Note: depending on the number of permutations, the Docker image RAM may need 
# to be increased (e.g., 2GB for 100k permutations is insufficient, but 8GB is enough)
if [ -e perm_test_plots ]; then
  rm -rf perm_test_plots
fi
mkdir perm_test_plots
/opt/rCNV2/analysis/paper/plot/large_segments/plot_segment_permutations.R \
  --constrained-genes gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.master_segments.bed.gz \
  ${prefix}.${n_seg_perms}_permuted_segments.bed.gz \
  ${prefix}.lit_GDs.${n_seg_perms}_permuted_segments.bed.gz \
  perm_test_plots/ \
  "${prefix}.segment_perms"
gzip perm_test_plots/${prefix}.segment_perms.permBySize.stats.tsv


# Run segment permutation tests while matching on number of genes per segment
# Note: in practice, this is parallelized in the cloud using segment_permutation_bygene.wdl
# The code to execute these permutation tests is contained elsewhere
# Copy results of segment permutation tests (note: requires permissions)
n_seg_perms=100000
gsutil -m cp \
  ${rCNV_bucket}/analysis/paper/data/large_segments/permutations/${prefix}*${n_seg_perms}_permuted_segments_bygene.tsv.gz \
  ./


# Plot segment permutation results while matching on number of genes per segment
# Note: depending on the number of permutations, the Docker image RAM may need 
# to be increased (e.g., 2GB for 100k permutations is insufficient)
# TODO: DEBUG THIS
if ! [ -e perm_test_plots ]; then
  mkdir perm_test_plots
fi
/opt/rCNV2/analysis/paper/plot/large_segments/plot_segment_bygene_perms.R \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.master_segments.bed.gz \
  ${prefix}.${n_seg_perms}_permuted_segments_bygene.tsv.gz \
  ${prefix}.lit_GDs.${n_seg_perms}_permuted_segments_bygene.tsv.gz \
  perm_test_plots/ \
  "${prefix}.segment_perms_bygene"
  

# Plot effect size covariates
if [ -e effect_size_plots ]; then
  rm -rf effect_size_plots
fi
mkdir effect_size_plots
/opt/rCNV2/analysis/paper/plot/large_segments/plot_effect_sizes.R \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.master_segments.bed.gz \
  effect_size_plots/${prefix}


# Plot pleiotropy covariates
if [ -e pleiotropy_plots ]; then
  rm -rf pleiotropy_plots
fi
mkdir pleiotropy_plots
/opt/rCNV2/analysis/paper/plot/large_segments/plot_pleiotropy.R \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.master_segments.bed.gz \
  pleiotropy_plots/${prefix}


# Plot distributions of de novo PTVs & missense across segments
if [ -e dnm_distributions ]; then
  rm -rf dnm_distributions
fi
mkdir dnm_distributions
/opt/rCNV2/analysis/paper/plot/large_segments/plot_dnm_comparisons.R \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.master_segments.bed.gz \
  refs/ddd_dnm_counts.tsv.gz \
  refs/asc_dnm_counts.tsv.gz \
  refs/asc_dnm_counts.unaffecteds.tsv.gz \
  refs/gene_mutation_rates.tsv.gz \
  dnm_distributions/${prefix}


# Generate mini Miami plot of example phenotype for main figure panel
if [ -e assoc_stat_plots ]; then
  rm -rf assoc_stat_plots
fi
mkdir assoc_stat_plots
/opt/rCNV2/analysis/paper/plot/large_segments/plot_example_miami.R \
  --del-cutoff ${del_cutoff} \
  --dup-cutoff ${dup_cutoff} \
  --xaxis-label "200kb sliding windows in 10kb steps" \
  meta_stats/${example_hpo}.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
  meta_stats/${example_hpo}.rCNV.DUP.sliding_window.meta_analysis.stats.bed.gz \
  assoc_stat_plots/${prefix}.example_miami.png


# Plot results from genome-wide permutation testing to establish FDR
if ! [ -e assoc_stat_plots ]; then
  mkdir assoc_stat_plots
fi
if ! [ -e meta_stats/perm_res ]; then
  mkdir meta_stats/perm_res
fi
# Gather all permutation P-values 
# (Note: in practice, this has been parallelized in FireCloud with collect_permuted_fdrs.sliding_windows.wdl)
# (Note 2: this is a _different_ WDL than collect_permuted_meta_p_matrices.sliding_windows.wdl)
gsutil -m cp \
  ${rCNV_bucket}/analysis/sliding_windows/permuted_pvalue_matrices/*.permuted_fdrs.tsv.gz \
  meta_stats/perm_res/
# Make plots of permuted P-values vs empirical FDR target
/opt/rCNV2/analysis/paper/plot/large_segments/plot_permuted_fdr.R \
  --fdr-target "${del_cutoff}" \
  meta_stats/perm_res/${prefix}.rCNV.DEL.sliding_window.meta_analysis.stats.permuted_fdrs.tsv.gz \
  meta_stats/perm_res/${prefix}.rCNV.DUP.sliding_window.meta_analysis.stats.permuted_fdrs.tsv.gz \
  refs/HPOs_by_metacohort.table.tsv \
  assoc_stat_plots/${prefix}


# Plot correlation of primary & secondary P-values for all phenotypes & CNV classes
if ! [ -e assoc_stat_plots ]; then
  mkdir assoc_stat_plots
fi
if ! [ -e meta_stats/matrices ]; then
  mkdir meta_stats/matrices
fi
# Strip out P-value columns per HPO & CNV pair
while read nocolon hpo; do
  echo $nocolon
  for cnv in DEL DUP; do
    echo $cnv
    statsfile=meta_stats/$nocolon.rCNV.$cnv.sliding_window.meta_analysis.stats.bed.gz
    for column in meta_neg_log10_p meta_neg_log10_p_secondary; do
      echo $column
      idx=$( zcat $statsfile | head -n1 | sed 's/\t/\n/g' \
             | awk -v column=$column '{ if ($1==column) print NR }' )
      zcat $statsfile | sed '1d' | cut -f$idx \
      | cat <( echo -e "${nocolon}_${cnv}" ) - \
      > meta_stats/matrices/$nocolon.$cnv.$column.tsv
    done
  done
done < refs/test_phenotypes.list
# Collect bin coordinates
zcat \
  meta_stats/$( head -n1 refs/test_phenotypes.list | cut -f1 ).rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz \
| cut -f1-3 \
> window_coordinates.bed
# Make matrices for primary and secondary P-values across phenotypes per CNV type 
for cnv in DEL DUP; do
  for column in meta_neg_log10_p meta_neg_log10_p_secondary; do
    paste \
      window_coordinates.bed \
      meta_stats/matrices/*.$cnv.$column.tsv \
    | bgzip -c \
    > meta_stats/matrices/${prefix}.$cnv.$column.all_hpos.bed.gz
  done
done
# Generate plots
/opt/rCNV2/analysis/paper/plot/large_segments/plot_sliding_window_pval_distribs.R \
  --del-cutoff ${del_cutoff} \
  --dup-cutoff ${dup_cutoff} \
  --del-nomsig-bed ./nomsig_windows.DEL.bed \
  --dup-nomsig-bed ./nomsig_windows.DUP.bed \
  meta_stats/matrices/${prefix}.DEL.meta_neg_log10_p.all_hpos.bed.gz \
  meta_stats/matrices/${prefix}.DUP.meta_neg_log10_p.all_hpos.bed.gz \
  meta_stats/matrices/${prefix}.DEL.meta_neg_log10_p_secondary.all_hpos.bed.gz \
  meta_stats/matrices/${prefix}.DUP.meta_neg_log10_p_secondary.all_hpos.bed.gz \
  refs/${prefix}.reordered_hpos.txt \
  refs/HPOs_by_metacohort.table.tsv \
  assoc_stat_plots/${prefix}
# Calculate fraction of genome with nominal association with at least one phenotype
searchspace=$( zcat meta_stats/matrices/${prefix}.DEL.meta_neg_log10_p.all_hpos.bed.gz \
               | cut -f1-3 | fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n \
               | bedtools merge -i - | awk '{ sum+=$3-$2 }END{ print sum }' )
zcat ./nomsig_windows.*.bed.gz \
| sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - \
| awk -v denom=$searchspace '{ sum+=$3-$2 }END{ print sum/denom }'


# Plot correlation of sizes for original significant regions and fine-mapped reigons
if ! [ -e assoc_stat_plots ]; then
  mkdir assoc_stat_plots
fi
if ! [ -e meta_stats/sigbins ]; then
  mkdir meta_stats/sigbins
fi
# Extract significant windows per phenotype from summary statistics
while read nocolon hpo; do
  echo $nocolon
  for cnv in DEL DUP; do
    echo $cnv
    case $cnv in
      DEL)
        primary_cutoff=$del_cutoff
        ;;
      DUP)
        primary_cutoff=$dup_cutoff
        ;;
      esac
    /opt/rCNV2/analysis/paper/scripts/large_segments/get_sig_windows.R \
      --primary-p-cutoff $primary_cutoff \
      --secondary-p-cutoff 0.05 \
      --min-nominal 2 \
      --secondary-or-nominal \
      --fdr-q-cutoff 0.05 \
      meta_stats/$nocolon.rCNV.$cnv.sliding_window.meta_analysis.stats.bed.gz \
      meta_stats/sigbins/$nocolon.$cnv.sig_windows.bed
    # fgrep -v "#" meta_stats/sigbins/$nocolon.$cnv.sig_windows.bed \
    # | sort -Vk1,1 -k2,2n -k3,3n \
    # | bedtools merge -d 1000000 -i - \
    # | bgzip -c \
    # > meta_stats/sigbins/$nocolon.$cnv.sig_windows.merged.bed.gz
  done
done < refs/test_phenotypes.list
# Cluster all significant windows to determine maximal segment size
for cnv in DEL DUP; do
  cat meta_stats/sigbins/*.$cnv.sig_windows.bed \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -d 200000 -i - \
  | bgzip -c \
  > meta_stats/sigbins/all_HPOs.$cnv.sig_windows.merged.bed.gz
done
# Match final significant refined segments with their original significant windows
# TODO: UPDATE THIS TO MATCH ONE SEGMENT TO ONE ORIGINAL CANDIDATE REGION
# while read nocolon hpo; do
#   for cnv in DEL DUP; do
#     zcat rCNV.final_segments.associations.bed.gz \
#     | fgrep -w $hpo \
#     | fgrep -w $cnv \
#     | cut -f1-4 \
#     | bedtools intersect -wa -wb -a - \
#       -b meta_stats/sigbins/$nocolon.$cnv.sig_windows.merged.bed.gz \
#     | awk -v FS="\t" -v OFS="\t" -v hpo=$hpo -v cnv=$cnv \
#       '{ print $4, cnv, hpo, $7-$6, $3-$2 }'
#   done
# done < refs/test_phenotypes.list \
# | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V \
# | cat <( echo -e "region_id\tcnv\thpo\toriginal_size\trefined_size" ) - \
# | gzip -c \
# > ${prefix}.associations.old_vs_new_size.tsv.gz
# # Scatterplot of original segment size (total sig bp) & finemapped size
# /opt/rCNV2/analysis/paper/plot/large_segments/plot_orig_vs_refined_assoc_sizes.R \
#   ${prefix}.associations.old_vs_new_size.tsv.gz \
#   assoc_stat_plots/${prefix}


# Plot master grid summarizing segment association across all phenotypes
if [ -e association_grid ]; then
  rm -rf association_grid
fi
mkdir association_grid
# Collapse overlapping DEL/DUP segments for sake of plotting
while read intervals rid; do
  echo "$intervals" | sed -e 's/\;/\n/g' -e 's/\:\|\-/\t/g' \
  | awk -v OFS="\t" -v rid=$rid '{ print $0, rid }' \
  | bedtools merge -i - -c 4 -o distinct -d 10000000
done < <( zcat rCNV.final_segments.loci.bed.gz | grep -ve '^#' \
          | awk -v FS="\t" -v OFS="\t" '{ print $(NF-3), $4 }' ) \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| bedtools merge -i - -d 200000 -c 4 -o distinct \
| cut -f4 \
| uniq \
> locus_clusters.gw_sig.txt
while read intervals rid; do
  echo "$intervals" | sed -e 's/\;/\n/g' -e 's/\:\|\-/\t/g' \
  | awk -v OFS="\t" -v rid=$rid '{ print $0, rid }' \
  | bedtools merge -i - -c 4 -o distinct -d 10000000
done < <( zcat all_gds.bed.gz | grep -ve '^#' \
          | awk -v FS="\t" -v OFS="\t" '{ print $8, $4 }' ) \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| bedtools merge -i - -d 200000 -c 4 -o distinct \
| cut -f4 \
| uniq \
> locus_clusters.all_gds.txt
# Merge all sumstats for gw-sig and lit GDs
cat <( zcat ${prefix}.final_segments.loci.all_sumstats.tsv.gz ) \
    <( zcat ${prefix}.lit_gds.all_sumstats.tsv.gz | grep -ve '^#' ) \
| gzip -c \
> ${prefix}.all_segs.all_sumstats.tsv.gz
# Generate master locus grid plot
/opt/rCNV2/analysis/paper/plot/large_segments/plot_association_grid.R \
  --gw-clusters locus_clusters.gw_sig.txt \
  --lit-clusters locus_clusters.all_gds.txt \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.final_segments.loci.bed.gz \
  ${prefix}.master_segments.bed.gz \
  ${prefix}.all_segs.all_sumstats.tsv.gz \
  refs/${prefix}.reordered_hpos.txt \
  association_grid/${prefix}


# Copy all plots to gs:// bucket (note: requires permissions)
gsutil -m cp -r \
  basic_distribs \
  perm_test_plots \
  effect_size_plots \
  pleiotropy_plots \
  dnm_distributions \
  assoc_stat_plots \
  association_grid \
  ${rCNV_bucket}/analysis/paper/plots/large_segments/


