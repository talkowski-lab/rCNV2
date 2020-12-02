#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for secondary analyses of gene association & fine mapping


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"
export finemap_dist=1000000


# Download necessary data (note: requires permissions)
gsutil -m cp \
  ${rCNV_bucket}/results/gene_association/* \
  ${rCNV_bucket}/results/segment_association/* \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/rCNV.*.gene_fine_mapping.gene_stats.merged_no_variation_features.all_genes_from_blocks.tsv \
  ./
# mkdir meta_stats/
# gsutil -m cp -r \
#   ${rCNV_bucket}/analysis/gene_burden/**meta_analysis.stats.bed.gz \
#   meta_stats/
mkdir refs/
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features*bed.gz \
  ${rCNV_bucket}/analysis/paper/data/misc/gene_feature_metadata.tsv \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  ${rCNV_bucket}/analysis/analysis_refs/test_phenotypes.list \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  refs/


# Plot correlation heatmap of gene features
if ! [ -e feature_heatmap ]; then
  mkdir feature_heatmap
fi
/opt/rCNV2/analysis/paper/plot/gene_association/plot_gene_feature_cor_matrix.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
	refs/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz \
  refs/gencode.v19.canonical.pext_filtered.all_features.no_variation.eigenfeatures.bed.gz \
  refs/gene_feature_metadata.tsv \
  feature_heatmap/${prefix}


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


# Plot fine-mapping descriptive panels (number of genes per block, distribution of PIP, etc)
for CNV in DEL DUP; do 
  fgrep -v "#" rCNV.$CNV.gene_fine_mapping.gene_stats.merged_no_variation_features.all_genes_from_blocks.tsv \
  | awk -v OFS="\t" -v CNV=$CNV '{ print $0, CNV }'
done \
| cat <( head -n1 rCNV.DEL.gene_fine_mapping.gene_stats.merged_no_variation_features.all_genes_from_blocks.tsv \
         | awk -v OFS="\t" '{ print $0, "cnv" }' ) - \
> all_genes_from_fine_mapping.tsv
if ! [ -e finemapping_distribs ]; then
  mkdir finemapping_distribs
fi
/opt/rCNV2/analysis/paper/plot/gene_association/plot_finemapping_distribs.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.final_genes.credible_sets.bed.gz \
  rCNV.final_genes.associations.bed.gz \
  all_genes_from_fine_mapping.tsv \
  finemapping_distribs/${prefix}


# Plot master annotated 2x2 table of fine-mapped genes
while read nocolon hpo; do
  echo -e "$hpo\trefs/gene_lists/$nocolon.HPOdb.genes.list"
done < refs/test_phenotypes.list \
> omim.gene_lists.tsv
/opt/rCNV2/analysis/paper/plot/gene_association/plot_finemapped_gene_table.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.final_genes.credible_sets.bed.gz \
  rCNV.final_genes.associations.bed.gz \
  refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
  refs/gene_lists/gnomad.v2.1.1.mis_constrained.genes.list \
  omim.gene_lists.tsv \
  ${prefix}


# Plot gene set enrichments for fine-mapped genes vs. various gene metadata
echo -e "lof_constrained\trefs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list" > enrichment.genelists.tsv
cat refs/gene_lists/DDG2P.*.genes.list | sort | uniq > ddg2p.all_dom.genes.list
echo -e "ddg2p\tddg2p.all_dom.genes.list" >> enrichment.genelists.tsv
echo -e "omim\trefs/gene_lists/HP0000118.HPOdb.genes.list" >> enrichment.genelists.tsv
if ! [ -e finemapping_distribs ]; then
  mkdir finemapping_distribs
fi
/opt/rCNV2/analysis/paper/plot/gene_association/plot_finemapped_enrichments.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
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


# Count number of fine-mapped duplication genes not present in OMIM, DECIPHER, or DDG2P (for abstract)
zcat rCNV.final_genes.genes.bed.gz \
| fgrep -v "#" | fgrep -w DUP | cut -f4 | sort | uniq \
| fgrep -wvf <( cat refs/gene_lists/HP0000118.HPOdb.genes.list \
                    refs/gene_lists/ClinGen.all_triplosensitive.genes.list \
                    refs/gene_lists/DDG2P.all_gof.genes.list ) \
| wc -l


# Count number of triplosensitive genes that qualify as mechanism expansion per DECIPHER + DDG2P
zcat rCNV.final_genes.genes.bed.gz \
| fgrep -v "#" | fgrep -w DUP | cut -f4 | sort | uniq \
| fgrep -wvf <( cat refs/gene_lists/ClinGen.all_triplosensitive.genes.list \
                    refs/gene_lists/DDG2P.all_gof.genes.list ) \
| fgrep -wf <( cat refs/gene_lists/ClinGen.all_haploinsufficient.genes.list \
                   refs/gene_lists/DDG2P.all_lof.genes.list ) \
| wc -l


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
  finemapping_distribs \
  ${prefix}.*finemapped_genes_grid*pdf \
  ${rCNV_bucket}/analysis/paper/plots/gene_association/

