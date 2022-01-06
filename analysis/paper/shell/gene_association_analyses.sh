#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Main code block for secondary analyses of gene association & fine mapping


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d2"
export finemap_dist=1000000
alias addcom="sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta'"


# Download necessary data (note: requires permissions)
gsutil -m cp \
  ${rCNV_bucket}/results/gene_association/* \
  ${rCNV_bucket}/results/segment_association/* \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/rCNV.*.gene_fine_mapping.gene_stats.naive_priors.tsv \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/rCNV.*.gene_fine_mapping.gene_stats.genetics_only.tsv \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/rCNV.*.gene_fine_mapping.gene_stats.merged_no_variation_features.tsv \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/rCNV.*.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv \
  ./
mkdir meta_stats/
gsutil -m cp -r \
  ${rCNV_bucket}/analysis/gene_burden/**meta_analysis.stats.bed.gz \
  meta_stats/
mkdir refs/
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features*bed.gz \
  ${rCNV_bucket}/analysis/paper/data/misc/gene_feature_metadata.tsv \
  ${rCNV_bucket}/analysis/analysis_refs/** \
  ${rCNV_bucket}/analysis/paper/data/hpo/${prefix}.reordered_hpos.txt \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  refs/


# Tabix all meta-analysis stats
find meta_stats/ -name "*meta_analysis.stats.bed.gz" | xargs -I {} tabix -f {}


# Get genome-wide significance thresholds
example_hpo="HP0012759"
del_cutoff=$( awk -v FS="\t" -v hpo=${example_hpo} '{ if ($1==hpo) print $2 }' \
              refs/gene_burden.rCNV.DEL.bonferroni_pval.hpo_cutoffs.tsv )
dup_cutoff=$( awk -v FS="\t" -v hpo=${example_hpo} '{ if ($1==hpo) print $2 }' \
              refs/gene_burden.rCNV.DUP.bonferroni_pval.hpo_cutoffs.tsv )


# Plot correlation heatmap of gene features
if ! [ -e feature_heatmap ]; then
  mkdir feature_heatmap
fi
/opt/rCNV2/analysis/paper/plot/gene_association/plot_gene_feature_cor_matrix.R \
	refs/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz \
  refs/gencode.v19.canonical.pext_filtered.all_features.no_variation.eigenfeatures.bed.gz \
  refs/gene_feature_metadata.tsv \
  feature_heatmap/${prefix}


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
    statsfile=meta_stats/$nocolon.rCNV.$cnv.gene_burden.meta_analysis.stats.bed.gz
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
  meta_stats/$( head -n1 refs/test_phenotypes.list | cut -f1 ).rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz \
| cut -f1-3 \
> gene_coordinates.bed
# Make matrices for primary and secondary P-values across phenotypes per CNV type 
for cnv in DEL DUP; do
  for column in meta_neg_log10_p meta_neg_log10_p_secondary; do
    paste \
      gene_coordinates.bed \
      meta_stats/matrices/*.$cnv.$column.tsv \
    | bgzip -c \
    > meta_stats/matrices/${prefix}.$cnv.$column.all_hpos.bed.gz
  done
done
# Generate plots
/opt/rCNV2/analysis/paper/plot/large_segments/plot_sliding_window_pval_distribs.R \
  --del-cutoff ${del_cutoff} \
  --dup-cutoff ${dup_cutoff} \
  --del-nomsig-bed ./nomsig_genes.DEL.bed \
  --dup-nomsig-bed ./nomsig_genes.DUP.bed \
  meta_stats/matrices/${prefix}.DEL.meta_neg_log10_p.all_hpos.bed.gz \
  meta_stats/matrices/${prefix}.DUP.meta_neg_log10_p.all_hpos.bed.gz \
  meta_stats/matrices/${prefix}.DEL.meta_neg_log10_p_secondary.all_hpos.bed.gz \
  meta_stats/matrices/${prefix}.DUP.meta_neg_log10_p_secondary.all_hpos.bed.gz \
  refs/${prefix}.reordered_hpos.txt \
  refs/HPOs_by_metacohort.table.tsv \
  assoc_stat_plots/${prefix}


# Generate mini Miami plot of example phenotype for supplementary figure panel
if ! [ -e assoc_stat_plots ]; then
  mkdir assoc_stat_plots
fi
mkdir assoc_stat_plots
ngenes=$( zcat meta_stats/${example_hpo}.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz \
          | cut -f1 | fgrep -v "#" | wc -l | addcom )
/opt/rCNV2/analysis/paper/plot/large_segments/plot_example_miami.R \
  --del-cutoff ${del_cutoff} \
  --dup-cutoff ${dup_cutoff} \
  --cutoff-label "Exome-wide significance" \
  --xaxis-label "$ngenes protein-coding genes" \
  meta_stats/${example_hpo}.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz \
  meta_stats/${example_hpo}.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz \
  assoc_stat_plots/${prefix}.example_miami.png


# # Collect list of genes theoretically eligible to be in each credible set
# while read chrom start end csID; do
#   echo -e "$chrom\t$(( $start - $finemap_dist ))\t$(( $end + $finemap_dist ))" \
#   | awk -v OFS="\t" '{ if ($2<0) $2=0; print $1, $2, $3 }' \
#   | bedtools intersect -u -wa \
#     -a refs/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
#     -b - \
#   | cut -f4 | sort | uniq | paste -s -d\; \
#   | awk -v OFS="\t" -v csID=$csID '{ print csID, $1 }'
# done < <( zcat rCNV.final_genes.credible_sets.bed.gz | fgrep -v "#" | cut -f1-4 ) \
# | cat <( echo -e "#credible_set_id\teligible_genes" ) - \
# > credible_set.eligible_genes.tsv


# Plot fine-mapping descriptive panels (number of genes per block, distribution of PIPs, etc)
head -n1 rCNV.DEL.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv \
| awk -v OFS="\t" '{ print $0, "cnv" }' > pip_header.tsv
for CNV in DEL DUP; do
  fgrep -v "#" rCNV.$CNV.gene_fine_mapping.gene_stats.naive_priors.tsv \
  | awk -v OFS="\t" -v CNV=$CNV '{ print $0, CNV }'
done  \
| cat pip_header.tsv - > all_PIPs.prior.tsv
for CNV in DEL DUP; do
  fgrep -v "#" rCNV.$CNV.gene_fine_mapping.gene_stats.genetics_only.tsv \
  | awk -v OFS="\t" -v CNV=$CNV '{ print $0, CNV }'
done  \
| cat pip_header.tsv - > all_PIPs.posterior.tsv
for CNV in DEL DUP; do 
  fgrep -v "#" rCNV.$CNV.gene_fine_mapping.gene_stats.merged_no_variation_features.tsv \
  | awk -v OFS="\t" -v CNV=$CNV '{ print $0, CNV }'
done \
| cat pip_header.tsv - > all_PIPs.full_model.tsv
if ! [ -e finemapping_distribs ]; then
  mkdir finemapping_distribs
fi
/opt/rCNV2/analysis/paper/plot/gene_association/plot_finemapping_distribs.R \
  rCNV.final_genes.credible_sets.bed.gz \
  rCNV.final_genes.associations.bed.gz \
  all_PIPs.prior.tsv \
  all_PIPs.posterior.tsv \
  all_PIPs.full_model.tsv \
  finemapping_distribs/${prefix}


# Plot annotated 2x2 tables of fine-mapped genes
## TODO: ADD SPLITS BY SIGNIFICANCE & ADULT/DEV
if ! [ -e finemapped_gene_grids ]; then
  mkdir finemapped_gene_grids
fi
while read nocolon hpo; do
  echo -e "$hpo\trefs/gene_lists/$nocolon.HPOdb.genes.list"
done < refs/test_phenotypes.list \
> omim.gene_lists.tsv
/opt/rCNV2/analysis/paper/plot/gene_association/plot_finemapped_gene_table.R \
  rCNV.final_genes.credible_sets.bed.gz \
  rCNV.final_genes.associations.bed.gz \
  refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
  refs/gene_lists/gnomad.v2.1.1.mis_constrained.genes.list \
  omim.gene_lists.tsv \
  finemapped_gene_grids/${prefix}


# Plot gene set enrichments for fine-mapped genes vs. various gene metadata
echo -e "lof_constrained\trefs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list" > enrichment.genelists.tsv
echo -e "mis_constrained\trefs/gene_lists/gnomad.v2.1.1.mis_constrained.genes.list" >> enrichment.genelists.tsv
cat refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
    refs/gene_lists/gnomad.v2.1.1.mis_constrained.genes.list \
| sort | uniq \
> gnomad.v2.1.1.all_constrained.genes.list
echo -e "all_constrained\tgnomad.v2.1.1.all_constrained.genes.list" >> enrichment.genelists.tsv
cat refs/gene_lists/DDG2P.*.genes.list | sort | uniq > ddg2p.all_dom.genes.list
echo -e "ddg2p\tddg2p.all_dom.genes.list" >> enrichment.genelists.tsv
echo -e "omim\trefs/gene_lists/HP0000118.HPOdb.genes.list" >> enrichment.genelists.tsv
if ! [ -e finemapping_distribs ]; then
  mkdir finemapping_distribs
fi
/opt/rCNV2/analysis/paper/plot/gene_association/plot_finemapped_enrichments.R \
  rCNV.final_genes.credible_sets.bed.gz \
  rCNV.final_genes.associations.bed.gz \
  refs/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
  enrichment.genelists.tsv \
  finemapping_distribs/${prefix}


# Calculate fraction of credible sets not overlapping significant large segment
for CNV in DEL DUP; do
  zcat rCNV.final_segments.loci.bed.gz \
  | fgrep -w $CNV | cut -f22 | sed 's/\;/\n/g' \
  | sort | uniq \
  | fgrep -wvf - <( zcat rCNV.final_genes.credible_sets.bed.gz ) \
  | fgrep -w $CNV
done | wc -l


# Count number of significant large segments with at least one confident fine-mapped gene
for CNV in DEL DUP; do
  zcat rCNV.final_segments.loci.bed.gz \
  | fgrep -w $CNV | cut -f22 \
  | fgrep -wf <( zcat rCNV.final_genes.genes.bed.gz | fgrep -w $CNV | cut -f4 )
done | wc -l


# Count number of fine-mapped genes also present in OMIM
zcat rCNV.final_genes.genes.bed.gz \
| fgrep -v "#" | cut -f4 | sort | uniq \
| fgrep -wf - refs/gene_lists/HP0000118.HPOdb.genes.list \
| wc -l


# Count number of fine-mapped triplosensitive genes not present in ClinGen TS or DDG2P GoF/other (for abstract)
zcat rCNV.final_genes.genes.bed.gz \
| fgrep -v "#" | fgrep -w DUP | cut -f4 | sort | uniq \
| fgrep -wvf <( cat refs/gene_lists/ClinGen.all_triplosensitive.genes.list \
                    refs/gene_lists/DDG2P.all_gof.genes.list \
                    refs/gene_lists/DDG2P.all_other.genes.list ) \
| wc -l
# Further count number of fine-mapped novel triplo genes with established LoF mechanism in ClinGen/DDG2P
zcat rCNV.final_genes.genes.bed.gz \
| fgrep -v "#" | fgrep -w DUP | cut -f4 | sort | uniq \
| fgrep -wvf <( cat refs/gene_lists/ClinGen.all_triplosensitive.genes.list \
                    refs/gene_lists/DDG2P.all_gof.genes.list \
                    refs/gene_lists/DDG2P.all_other.genes.list ) \
| fgrep -wf <( cat refs/gene_lists/ClinGen.all_haploinsufficient.genes.list \
                   refs/gene_lists/DDG2P.all_lof.genes.list ) \
| wc -l


# Gather mean odds ratio across all gw sig large segments and credsets
for wrapper in 1; do
  zcat rCNV.final_segments.associations.bed.gz \
  | fgrep -v "#" | cut -f10
  for CNV in DEL DUP; do
    zcat rCNV.final_segments.loci.bed.gz \
    | fgrep -w $CNV | cut -f22 | sed 's/\;/\n/g' \
    | sort | uniq \
    | fgrep -wvf - <( zcat rCNV.final_genes.credible_sets.bed.gz ) \
    | fgrep -w $CNV
  done | cut -f9
done | awk '{ sum+=$1 }END{ print sum/NR }'


# Gather mean & range of odds ratios across all credsets
zcat rCNV.final_genes.credible_sets.bed.gz | fgrep -v "#" | cut -f9 \
| awk '{ sum+=$1 }END{ print sum/NR }'
zcat rCNV.final_genes.credible_sets.bed.gz | fgrep -v "#" | cut -f9 | sort -nk1,1 | head -n1
zcat rCNV.final_genes.credible_sets.bed.gz | fgrep -v "#" | cut -f9 | sort -nk1,1 | tail -n1

# # Count number of haploinsufficient genes that qualify as mechanism expansion per DECIPHER + DDG2P
# zcat rCNV.final_genes.genes.bed.gz \
# | fgrep -v "#" | fgrep -w DEL | cut -f4 | sort | uniq \
# | fgrep -wf <( cat refs/gene_lists/ClinGen.all_triplosensitive.genes.list \
#                    refs/gene_lists/DDG2P.all_gof.genes.list \
#                    refs/gene_lists/DDG2P.all_other.genes.list ) \
# | fgrep -wvf <( cat refs/gene_lists/ClinGen.all_haploinsufficient.genes.list \
#                    refs/gene_lists/DDG2P.all_lof.genes.list ) \
# | wc -l


# Count number of constrained fine-mapped genes that aren't present in OMIM
zcat rCNV.final_genes.genes.bed.gz \
| fgrep -v "#" | cut -f4 | sort | uniq \
| fgrep -wf - <( cat refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
                     refs/gene_lists/gnomad.v2.1.1.mis_constrained.genes.list ) \
| fgrep -wvf refs/gene_lists/HP0000118.HPOdb.genes.list \
| wc -l


# Copy all plots to final gs:// directory
gsutil -m cp -r \
  feature_heatmap \
  assoc_stat_plots \
  finemapping_distribs \
  ${prefix}.*finemapped_genes_grid*pdf \
  ${rCNV_bucket}/analysis/paper/plots/gene_association/

