#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Prepare release files for Zenodo deposit


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="Collins_rCNV_2022"


# Make parent directory for all files to upload
mkdir ${prefix}_release_files


# Download and format sliding window summary stats
mkdir ${prefix}.sliding_window_sumstats
gsutil -m cp -r \
  ${rCNV_bucket}/analysis/sliding_windows/**meta_analysis.stats.bed.gz \
  ${prefix}.sliding_window_sumstats/
find ${prefix}.sliding_window_sumstats/*meta_analysis.stats.bed.gz \
| xargs -I {} tabix {}
cat << EOF > ${prefix}.sliding_window_sumstats/README
+----------------------------------------------------+
| Sliding window rCNV association summary statistics |
+----------------------------------------------------+
Contact: Ryan Collins <rlcollins@g.harvard.edu>

This directory contains rCNV association summary statistics for genome-wide sliding 
window meta-analyses for 54 phenotypes as described in Collins et al., 2022.

For each phenotype, we provide one bgzipped BED file each for deletions (DEL) and
duplications (DUP). Each BED file is also provided with a tabix index (*.tbi).

Each BED file contains summary statistics for 267,237 autosomal 200kb windows in 
10kb steps. Each window has the following 20 columns:

1. chr: window chromosome (in GRCh37 coordinates)
2. start: start position of window (in GRCh37 coordinates)
3. end: end position of window (in GRCh37 coordinates)
4. n_nominal_cohorts: number of cohorts in meta-analysis that achieved nominal
   (P<0.05) significance in single-cohort Fisher's exact tests
5. top_cohort: cohort with most significant P-value from single-cohort Fisher's
   exact tests (as in column 4)
6. cohorts_excluded_from_meta: cohorts excluded from meta-analysis due to either 
   having inadequate microarray probe density in control samples or having 
   too few (<300) case samples for the phenotype in question
7. case_freq: simple pooled mean of case CNV frequency across all cohorts
8. control_freq: simple pooled mean of control CNV frequency across all cohorts
9. meta_lnOR: meta-analysis log-odds ratio
10. meta_lnOR_lower: lower 95% confidence bound on log-odds ratio
11. meta_lnOR_upper: upper 95% confidence bound on log-odds ratio
12. meta_z: meta-analysis Z-score
13. meta_neg_log10_p: meta-analysis P-value scaled as -log10(P)
14. meta_neg_log10_fdr_q: Benjamini-Hochberg corrected false discovery rate Q value
15. meta_lnOR_secondary: same as column 9, but after excluding the cohort listed 
    in column 5 from the meta-analysis
16. meta_lnOR_lower_secondary: same as column 10, but after excluding the cohort 
    listed in column 5 from the meta-analysis
17. meta_lnOR_upper_secondary: same as column 11, but after excluding the cohort 
    listed in column 5 from the meta-analysis
18. meta_z_secondary: same as column 12, but after excluding the cohort listed 
    in column 5 from the meta-analysis
19. meta_neg_log10_p_secondary: same as column 13, but after excluding the cohort 
    listed in column 5 from the meta-analysis
20. meta_neg_log10_fdr_q_secondary: same as column 14, but after excluding the 
    cohort listed in column 5 from the meta-analysis
EOF
tar -czvf \
  ${prefix}_release_files/${prefix}.sliding_window_sumstats.tar.gz \
  ${prefix}.sliding_window_sumstats


# Download and format gene-based summary stats
mkdir ${prefix}.gene_association_sumstats
gsutil -m cp -r \
  ${rCNV_bucket}/analysis/gene_burden/**meta_analysis.stats.bed.gz \
  ${prefix}.gene_association_sumstats/
for ofile in ${prefix}.gene_association_sumstats/*gene_burden.meta_analysis.stats.bed.gz; do
  nfile="${prefix}.gene_association_sumstats/$( basename $ofile | sed 's/burden/association/g' )"
  mv -v $ofile $nfile
done
find ${prefix}.gene_association_sumstats/*meta_analysis.stats.bed.gz \
| xargs -I {} tabix {}
cat << EOF > ${prefix}.gene_association_sumstats/README
+------------------------------------------------+
| Gene-based rCNV association summary statistics |
+------------------------------------------------+
Contact: Ryan Collins <rlcollins@g.harvard.edu>

This directory contains rCNV association summary statistics for exome-wide 
gene-based meta-analyses for 54 phenotypes as described in Collins et al., 2022.

For each phenotype, we provide one bgzipped BED file each for deletions (DEL) and
duplications (DUP). Each BED file is also provided with a tabix index (*.tbi).

Each BED file contains summary statistics for 17,263 autosomal protein-coding
genes from Gencode v19. Each gene has the following 21 columns:

1. chr: chromosome of gene (in GRCh37 coordinates)
2. start: start position of canonical transcript (in GRCh37 coordinates)
3. end: end position of canonical transcript (in GRCh37 coordinates)
4. gene: gene symbol per Gencode v19
5. n_nominal_cohorts: number of cohorts in meta-analysis that achieved nominal
   (P<0.05) significance in single-cohort Fisher's exact tests
6. top_cohort: cohort with most significant P-value from single-cohort Fisher's
   exact tests (as in column 4)
7. cohorts_excluded_from_meta: cohorts excluded from meta-analysis due to either 
   having inadequate microarray probe density in control samples or having 
   too few (<300) case samples for the phenotype in question
8. case_freq: simple pooled mean of case CNV frequency across all cohorts
9. control_freq: simple pooled mean of control CNV frequency across all cohorts
10. meta_lnOR: meta-analysis log-odds ratio
11. meta_lnOR_lower: lower 95% confidence bound on log-odds ratio
12. meta_lnOR_upper: upper 95% confidence bound on log-odds ratio
13. meta_z: meta-analysis Z-score
14. meta_neg_log10_p: meta-analysis P-value scaled as -log10(P)
15. meta_neg_log10_fdr_q: Benjamini-Hochberg corrected false discovery rate Q value
16. meta_lnOR_secondary: same as column 10, but after excluding the cohort listed 
    in column 5 from the meta-analysis
17. meta_lnOR_lower_secondary: same as column 11, but after excluding the cohort 
    listed in column 5 from the meta-analysis
18. meta_lnOR_upper_secondary: same as column 12, but after excluding the cohort 
    listed in column 5 from the meta-analysis
19. meta_z_secondary: same as column 13, but after excluding the cohort listed 
    in column 5 from the meta-analysis
20. meta_neg_log10_p_secondary: same as column 14, but after excluding the cohort 
    listed in column 5 from the meta-analysis
21. meta_neg_log10_fdr_q_secondary: same as column 15, but after excluding the 
    cohort listed in column 5 from the meta-analysis
EOF
tar -czvf \
  ${prefix}_release_files/${prefix}.gene_association_sumstats.tar.gz \
  ${prefix}.gene_association_sumstats


# Download and format gene metadata
mkdir ${prefix}.gene_features_matrix
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.no_variation*bed.gz \
  ${prefix}.gene_features_matrix/
mv \
  ${prefix}.gene_features_matrix/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz \
  ${prefix}.gene_features_matrix/${prefix}.raw_gene_features.bed.gz
mv \
  ${prefix}.gene_features_matrix/gencode.v19.canonical.pext_filtered.all_features.no_variation.eigenfeatures.bed.gz \
  ${prefix}.gene_features_matrix/${prefix}.gene_eigenfeatures.bed.gz
find ${prefix}.gene_features_matrix/*bed.gz | xargs -I {} tabix {}
cat << EOF > ${prefix}.gene_features_matrix/README
+-----------------------------+
| Gene-level feature matrices |
+-----------------------------+
Contact: Ryan Collins <rlcollins@g.harvard.edu>

This directory contains gene-level features in two formats for 18,641 autosomal 
protein-coding genes as described in Collins et al., 2022.

All gene symbols and coordinates correspond to the Ensembl-defined canonical 
transcript present in the Gencode v19 GTF.

These features are provided in two formats:

1. Collins_rCNV_2022.raw_gene_features.bed.gz: this file contains the raw feature
   values for each gene prior to transformation/normalization and decomposition.

2. Collins_rCNV_2022.gene_eigenfeatures.bed.gz: this file contains the features 
   from [1] after transformation/normalization and PCA-based decomposition.

Please refer to Table S7 from Collins et al. for a description of each column
name, source, and criteria.
EOF
tar -czvf \
  ${prefix}_release_files/${prefix}.gene_features_matrix.tar.gz \
  ${prefix}.gene_features_matrix


# Download and rename gene scores
gsutil -m cp \
  ${rCNV_bucket}/results/gene_scoring/rCNV.gene_scores.tsv.gz \
  ${prefix}_release_files/${prefix}.dosage_sensitivity_scores.tsv.gz


# Write main README
cat << EOF > ${prefix}_release_files/README
+--------+
| README |
+--------+
Contact: Ryan Collins <rlcollins@g.harvard.edu>

This repository contains data from Collins et al. (2022), including:

1. Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz:
   This file contains predicted probabilities of haploinsufficiency (pHaplo) and
   triplosensitivity (pTriplo) for 18,641 autosomal protein-coding genes as 
   defined in Gencode v19.

2. Collins_rCNV_2022.sliding_window_sumstats.tar.gz:
   This compressed directory contains rCNV association summary statistics for 54
   phenotypes from genome-wide sliding window meta-analyses. Please refer to the
   README file included in this compressed directory for more details.

3. Collins_rCNV_2022.gene_association_sumstats.tar.gz:
   This compressed directory contains rCNV association summary statistics for 54
   phenotypes from exome-wide gene-based meta-analyses. Please refer to the
   README file included in this compressed directory for more details. 

4. Collins_rCNV_2022.gene_features_matrix.tar.gz:
   This compressed directory contains gene-level feature annotations for 145 
   features and 18,641 autosomal protein-coding genes. Please refer to the
   README file included in this compressed directory for more details. 

Smaller data files have been provided as supplemental tables alongside the
publication online.

Please also refer to the original publication for details on data sources, study
design, methods, and other analyses.
EOF


# Copy to Google Bucket
gsutil -m cp -r ${prefix}_release_files ${rCNV_bucket}/public/

