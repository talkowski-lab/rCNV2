#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Curate & annotate genes used in subsequent analyses


# Launch docker image
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"


# Add helper alias for formatting long integers
alias addcom="sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta'"


# Download analysis references
mkdir refs/
gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/** refs/
gsutil -m cp ${rCNV_bucket}/refs/** refs/
gsutil -m cp ${rCNV_bucket}/analysis/paper/data/misc/** refs/


# Download all gene data from Gencode v19
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_translations.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.metadata.Transcript_source.gz


# Build table of canonical protein-coding transcripts from Gencode
/opt/rCNV2/data_curation/gene/get_canonical_transcripts.py \
  gencode.v19.annotation.gtf.gz \
  gencode.v19.pc_translations.fa.gz \
  gencode.v19.pc_transcripts.fa.gz \
  gencode.v19.metadata.Transcript_source.gz \
  gencode.v19.canonical.tsv
gzip -f gencode.v19.canonical.tsv


# Subset GTF to autosomal canonical protein-coding transcripts
cat \
  <( zcat gencode.v19.canonical.tsv.gz \
     | cut -f2 | sed '1d' \
     | fgrep -wf - <( zcat gencode.v19.annotation.gtf.gz ) ) \
  <( zcat gencode.v19.canonical.tsv.gz \
     | cut -f1 | sed '1d' \
     | fgrep -wf - <( zcat gencode.v19.annotation.gtf.gz ) \
     | awk '($3=="gene")' ) \
| sed -e 's/^chr//' \
| grep -e '^[1-9]' \
| sort -k10,10 \
| /opt/rCNV2/data_curation/gene/filter_UTRs.py stdin stdout \
| sort -k1,1V -k4,4n \
| bgzip -c \
> gencode.v19.canonical.gtf.gz
tabix -f gencode.v19.canonical.gtf.gz


# Apply pext filter to exons from canonical transcripts
# Note: pext pre-formatting takes several hours locally; was parallelized in cloud 
# using process_pext.wdl, and can be downloaded (with permissions) at 
# ${rCNV_bucket}/cleaned_data/genes/annotations/gnomad.v2.1.1.pext.bed.gz
gsutil -m cp \
  gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.022719.tsv.bgz \
  ./
tabix -s 1 -b 2 -e 2 -S 1 \
  all.possible.snvs.tx_annotated.022719.tsv.bgz
/opt/rCNV2/data_curation/gene/process_pext.py \
  --pan-tissue \
  all.possible.snvs.tx_annotated.022719.tsv.bgz \
| sort -Vk1,1 -k2,2n -k3,3n \
| bgzip -c \
> gnomad.v2.1.1.pext.bed.gz
tabix -f gnomad.v2.1.1.pext.bed.gz
/opt/rCNV2/data_curation/gene/apply_pext_filter.py \
  --min-pext 0.2 \
  -o gencode.v19.canonical.pext_filtered.gtf.gz \
  --bgzip \
  --lost-genes genes_lost_during_pext_filtering.genes.list \
  gencode.v19.canonical.gtf.gz \
  gnomad.v2.1.1.pext.bed.gz
tabix -f gencode.v19.canonical.pext_filtered.gtf.gz


# Make gene list of all canonical autosomal genes used in analysis
for gtf in gencode.v19.canonical gencode.v19.canonical.pext_filtered; do
  zcat ${gtf}.gtf.gz \
  | awk -v FS="\t" '{ if ($3=="transcript") print $9 }' \
  | sed 's/;/\n/g' | fgrep "gene_name" \
  | sed 's/gene_name\ //g' | tr -d '"' \
  | awk '{ print $1 }' | sort -Vk1,1 | uniq \
  > ${gtf}.genes.list
done
zcat gencode.v19.canonical.tsv.gz \
| fgrep -wf gencode.v19.canonical.pext_filtered.genes.list \
| gzip -c \
> gencode.v19.canonical.pext_filtered.tsv.gz


# Preprocess GTEx expression matrix to compute summary data
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt
mkdir gtex_stats/
/opt/rCNV2/data_curation/gene/preprocess_GTEx.py \
  --gzip \
  --prefix gtex_stats/gencode.v19.canonical.pext_filtered.GTEx_v7_expression_stats \
  --n-pcs 20 \
  gencode.v19.canonical.pext_filtered.tsv.gz \
  GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz \
  GTEx_v7_Annotations_SampleAttributesDS.txt
/opt/rCNV2/data_curation/gene/get_variable_expressors.py \
  --min-mean-tpm 5 \
  --min-sample-tpm 1 \
  --min-tissues 3 \
  --variable-tpm-cv 1 \
  --invariant-tpm-cv 0.3 \
  --min-variable-prop-low 0.005 \
  --max-invariant-prop-low 0.001 \
  --min-variable-prop-high 0.05 \
  --max-invariant-prop-high 0.01 \
  --prefix gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors \
  gencode.v19.canonical.pext_filtered.tsv.gz \
  GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz \
  GTEx_v7_Annotations_SampleAttributesDS.txt
# Copy precomputed GTEx summary data to rCNV bucket (note: requires permissions)
gsutil -m cp -r gtex_stats \
  ${rCNV_bucket}/cleaned_data/genes/annotations/
gsutil -m cp \
  gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors*.genes.list \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists/


# Preprocess Epigenome Roadmap ChromHMM data to compute summaries
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/all.mnemonics.bedFiles.tgz
tar -xzvf all.mnemonics.bedFiles.tgz
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/labelmap_18_core_K27ac.tab
gsutil -m cp ${rCNV_bucket}/refs/REP_sample_manifest.tsv ./
mkdir roadmap_stats/
/opt/rCNV2/data_curation/gene/preprocess_REP.py \
  --state-manifest labelmap_18_core_K27ac.tab \
  --sample-manifest REP_sample_manifest.tsv \
  --n-pcs 20 \
  --prefix roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats \
  --gzip \
  gencode.v19.canonical.pext_filtered.gtf.gz \
  ./
# Copy precomputed Roadmap ChromHMM summary data to rCNV bucket (note: requires permissions)
gsutil -m cp -r roadmap_stats \
  ${rCNV_bucket}/cleaned_data/genes/annotations/
# Print HTML-formatted rows for README from ChromHMM data (helper function)
gsutil -m cp ${rCNV_bucket}/refs/REP_state_manifest.tsv ./
/opt/rCNV2/data_curation/gene/print_chromhmm_readme_rows.py REP_state_manifest.tsv  


# Gather per-gene metadata (genomic)
# Note: this is parallelized in FireCloud with get_gene_metadata.wdl
wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
bgzip Homo_sapiens.GRCh37.dna.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
ref_fasta=Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
for wrapper in 1; do 
  echo -e "refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz\tcoverage\trepeat_cov"
  echo -e "https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign100mer.bw\tmap-mean\talign_100mer"
done > gene_features.athena_tracklist.tsv
/opt/rCNV2/data_curation/gene/get_gene_features.py \
  --get-genomic \
  --centro-telo-bed refs/GRCh37.centromeres_telomeres.bed.gz \
  --ref-fasta ${ref_fasta} \
  --athena-tracks gene_features.athena_tracklist.tsv \
  --gnomad-constraint gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
  --no-scaling \
  --outbed gencode.v19.canonical.pext_filtered.genomic_features.bed.gz \
  --bgzip \
  gencode.v19.canonical.pext_filtered.gtf.gz


# Gather per-gene metadata (expression)
/opt/rCNV2/data_curation/gene/get_gene_features.py \
  --get-expression \
  --gtex-medians gtex_stats/gencode.v19.canonical.pext_filtered.GTEx_v7_expression_stats.median.tsv.gz \
  --gtex-mads gtex_stats/gencode.v19.canonical.pext_filtered.GTEx_v7_expression_stats.mad.tsv.gz \
  --gtex-pca gtex_stats/gencode.v19.canonical.pext_filtered.GTEx_v7_expression_stats.pca.tsv.gz \
  --outbed gencode.v19.canonical.pext_filtered.expression_features.bed.gz \
  --bgzip \
  gencode.v19.canonical.pext_filtered.gtf.gz


# Gather per-gene metadata (chromatin)
/opt/rCNV2/data_curation/gene/get_gene_features.py \
  --get-chromatin \
  --roadmap-means roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats.mean.tsv.gz \
  --roadmap-sds roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats.sd.tsv.gz \
  --roadmap-pca roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats.pca.tsv.gz \
  --outbed gencode.v19.canonical.pext_filtered.chromatin_features.bed.gz \
  --bgzip \
  gencode.v19.canonical.pext_filtered.gtf.gz


# Gather per-gene constraint metadata
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
wget http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt
gsutil -m cp ${rCNV_bucket}/cleaned_data/genes/annotations/EDS.Wang_2018.tsv.gz ./
wget https://storage.googleapis.com/gnomad-public/legacy/exac_browser/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz
wget https://doi.org/10.1371/journal.pgen.1001154.s002
/opt/rCNV2/data_curation/gene/get_gene_features.py \
  --get-constraint \
  --ref-fasta ${ref_fasta} \
  --gnomad-constraint gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
  --exac-cnv forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz \
  --rvis-tsv RVIS_Unpublished_ExACv2_March2017.txt \
  --eds-tsv EDS.Wang_2018.tsv.gz \
  --hi-tsv journal.pgen.1001154.s002 \
  --outbed gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  --bgzip \
  gencode.v19.canonical.pext_filtered.gtf.gz


# Gather per-gene metadata (variation)
/opt/rCNV2/data_curation/gene/get_gene_features.py \
  --get-variation \
  --gnomad-constraint gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
  --ddd-dnms refs/ddd_dnm_counts.tsv.gz \
  --asc-dnms refs/asc_dnm_counts.tsv.gz \
  --asc-unaffected-dnms refs/asc_dnm_counts.unaffecteds.tsv.gz \
  --gnomad-svs refs/gnomad_sv_nonneuro_counts.tsv.gz \
  --redin-bcas refs/redin_bca_counts.tsv.gz \
  --outbed gencode.v19.canonical.pext_filtered.variation_features.bed.gz \
  --bgzip \
  gencode.v19.canonical.pext_filtered.gtf.gz


# Join all per-gene metadata
/opt/rCNV2/data_curation/gene/join_gene_metadata.R \
  gencode.v19.canonical.pext_filtered.genomic_features.bed.gz \
  gencode.v19.canonical.pext_filtered.expression_features.bed.gz \
  gencode.v19.canonical.pext_filtered.chromatin_features.bed.gz \
  gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
| bgzip -c \
> gencode.v19.canonical.pext_filtered.all_features.bed.gz


# Copy canonical gene metadata to rCNV bucket (note: requires permissions)
gsutil -m cp gencode.v19.canonical*.gz* ${rCNV_bucket}/cleaned_data/genes/
gsutil -m cp gencode.v19.canonical*genes.list \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists/
gsutil -m cp genes_lost_during_pext_filtering.genes.list \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists/


# Generate BED file of pLoF constrained genes and mutationally tolerant genes from gnomAD
if ! [ -e gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz ]; then
  wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
fi
/opt/rCNV2/data_curation/gene/get_gnomad_genelists.R \
  gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
  gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
  "gnomad.v2.1.1"


# Generate lists of genes based on disruptions in gnomAD-SV
if ! [ -e gencode.v19.canonical.pext_filtered.variation_features.bed.gz ]; then
  gsutil -m cp \
    ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.variation_features.bed.gz \
    ./
fi
/opt/rCNV2/data_curation/gene/get_gnomad-sv_genelists.R \
  gencode.v19.canonical.pext_filtered.variation_features.bed.gz \
  gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
  "gnomad_sv.v2.1.nonneuro"


# Standardize gene lists curated in Karczewski et al. (2020)
git clone https://github.com/macarthur-lab/gnomad_lof.git
wget http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt
fgrep -wf gencode.v19.canonical.pext_filtered.genes.list \
  gnomad_lof/R/ko_gene_lists/list_CEGv2.tsv \
> cell_essential.genes.list
fgrep -wf gencode.v19.canonical.pext_filtered.genes.list \
  gnomad_lof/R/ko_gene_lists/list_NEGv1.tsv \
> cell_nonessential.genes.list
fgrep -wf gnomad_lof/R/ko_gene_lists/list_mouse_het_lethal_genes.tsv \
  HMD_HumanPhenotype.rpt \
| cut -f1 \
| fgrep -wf gencode.v19.canonical.pext_filtered.genes.list \
> mouse_het_lethal.genes.list
fgrep -wf gnomad_lof/R/ko_gene_lists/list_mouse_viable_genes.tsv \
  HMD_HumanPhenotype.rpt \
| cut -f1 \
| fgrep -wf gencode.v19.canonical.pext_filtered.genes.list \
> mouse_dispensable.genes.list


# Copy gnomAD & gnomAD-SV gene files to rCNV bucket (note: requires permissions)
gsutil -m cp gencode.v19.canonical.pext_filtered.constrained.bed.gz \
  ${rCNV_bucket}/analysis/analysis_refs/
gsutil -m cp gnomad.v2.1.1.*.genes.list \
  gnomad_sv.v2.1.nonneuro.*.genes.list \
  cell_*.genes.list \
  mouse_*.genes.list \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists/


# Collate list of olfactory receptor genes and copy to rCNV bucket (note: requires permissions)
gsutil -m cat ${rCNV_bucket}/raw_data/other/HUGO.olfactory_receptors.genes.list \
| fgrep -wf gencode.v19.canonical.pext_filtered.genes.list \
> olfactory_receptors.genes.list
gsutil -m cp \
  olfactory_receptors.genes.list \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists/


# Generate gene list of all genes associated with each HPO term
wget http://compbio.charite.de/jenkins/job/hpo.annotations/lastSuccessfulBuild/artifact/util/annotation/phenotype_to_genes.txt
while read pheno hpo; do
  awk -v hpo=${hpo} -v FS="\t" '{ if ($1==hpo) print $4 }' \
    phenotype_to_genes.txt \
  | fgrep -wf gencode.v19.canonical.pext_filtered.genes.list \
  | sort -Vk1,1 | uniq \
  > $pheno.HPOdb.genes.list
done < refs/test_phenotypes.list
# For HP0000118 & UNKNOWN, take union of all genes associated with any HPO term
sed '1d' phenotype_to_genes.txt \
| awk -v FS="\t" '{ print $4 }' \
| fgrep -wf gencode.v19.canonical.pext_filtered.genes.list \
| sort -Vk1,1 | uniq \
> HP0000118.HPOdb.genes.list
cp HP0000118.HPOdb.genes.list UNKNOWN.HPOdb.genes.list


# Generate gene list of all constrained genes associated with each HPO term
while read pheno hpo; do
  fgrep -wf gnomad.v2.1.1.lof_constrained.genes.list \
    $pheno.HPOdb.genes.list \
  > $pheno.HPOdb.constrained.genes.list
done < refs/test_phenotypes.list


# Copy HPO-based gene lists to rCNV bucket (note: requires permissions)
gsutil -m cp *.HPOdb.*genes.list ${rCNV_bucket}/cleaned_data/genes/gene_lists/


# Generate ClinGen dosage sensitive gene lists & copy to rCNV bucket (note: requires permissions)
wget ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv
/opt/rCNV2/data_curation/gene/parse_clingen_genes.R \
  ClinGen_gene_curation_list_GRCh37.tsv \
  gencode.v19.canonical.pext_filtered.genes.list \
  ./
gsutil -m cp ClinGen.*.genes.list ${rCNV_bucket}/cleaned_data/genes/gene_lists/  


# Generate DECIPHER/DDG2P gene lists & copy to rCNV bucket (note: requires permissions)
wget http://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz
/opt/rCNV2/data_curation/gene/parse_ddg2p_genes.R \
  DDG2P.csv.gz \
  gencode.v19.canonical.pext_filtered.genes.list \
  ./
gsutil -m cp DDG2P.*.genes.list \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists/
gsutil -m cp DDG2P.*.genes_with_hpos.tsv \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists/ddg2p_with_hpos/


# Generate COSMIC cancer census gene lists & copy to rCNV bucket (note: requires permissions)
gsutil -m cp ${rCNV_bucket}/raw_data/other/cancer_gene_census.csv ./
/opt/rCNV2/data_curation/gene/parse_cosmic_genes.R \
  cancer_gene_census.csv \
  gencode.v19.canonical.pext_filtered.genes.list \
  ./
gsutil -m cp COSMIC.*.genes.list \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists/


# Generate HTML table of gene lists for README
# Header
echo "| Gene set | Genes | Filename prefix | Source | Description |" \
> genelist_table.html.txt
echo "| :--- | ---: | :--- | :--- | :---- | " \
>> genelist_table.html.txt
# All genes
for wrapper in 1; do
  echo "All genes"
  wc -l gencode.v19.canonical.pext_filtered.genes.list \
  | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
  | sed 's/\.genes\.list//g' | addcom
  echo "Gencode v19 [Harrow _et al._, _Genome Res._, 2012](https://www.ncbi.nlm.nih.gov/pubmed/22955987)"
  echo "Canonical transcripts from autosomal, protein-coding genes"
done | paste -s \
| sed 's/\t/\ \|\ /g' \
| sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
>> genelist_table.html.txt
# Constrained genes
for wrapper in 1; do
  echo "LoF-constrained genes"
  wc -l gnomad.v2.1.1.lof_constrained.genes.list \
  | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
  | sed 's/\.genes\.list//g' | addcom
  echo "gnomAD v2.1.1 [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7)"
  echo "pLI ≥ 0.9 or in the first LOEUF sextile"
done | paste -s \
| sed 's/\t/\ \|\ /g' \
| sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
>> genelist_table.html.txt
for wrapper in 1; do
  echo "LoF-constrained (strict) genes"
  wc -l gnomad.v2.1.1.lof_constrained_strict.genes.list \
  | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
  | sed 's/\.genes\.list//g' | addcom
  echo "gnomAD v2.1.1 [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7)"
  echo "pLI ≥ 0.9 and in the first LOEUF sextile"
done | paste -s \
| sed 's/\t/\ \|\ /g' \
| sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
>> genelist_table.html.txt
# Missense constrained genes
for wrapper in 1; do
  echo "Missense-constrained genes"
  wc -l gnomad.v2.1.1.mis_constrained.genes.list \
  | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
  | sed 's/\.genes\.list//g' | addcom
  echo "gnomAD v2.1.1 [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7)"
  echo "Missense Z ≥ 3 or in the first MOEUF sextile"
done | paste -s \
| sed 's/\t/\ \|\ /g' \
| sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
>> genelist_table.html.txt
# Likely unconstrained genes
for wrapper in 1; do
  echo "Likely unconstrained genes"
  wc -l gnomad.v2.1.1.likely_unconstrained.genes.list \
  | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
  | sed 's/\.genes\.list//g' | addcom
  echo "gnomAD v2.1.1 [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7)"
  echo "LOEUF ≥ 1, mis. OEUF ≥ 1, synonymous Z-score ~ [-3, 3], pLI ≤ 0.1, LoF O/E in upper 50% of all genes, mis. OE in upper 50% of all genes, observed LoF > 0, observed mis. > 0"
done | paste -s \
| sed 's/\t/\ \|\ /g' \
| sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
>> genelist_table.html.txt
# Mutation tolerant genes
for wrapper in 1; do
  echo "Mutation-tolerant genes"
  wc -l gnomad.v2.1.1.mutation_tolerant.genes.list \
  | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
  | sed 's/\.genes\.list//g' | addcom
  echo "gnomAD v2.1.1 [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7)"
  echo "pLI ≥ 0.01, the last third of LOEUF, missense Z-score ≤ 0, missense OEUF ≥ 1, and synonymous Z-score ~ (-3, 3)"
done | paste -s \
| sed 's/\t/\ \|\ /g' \
| sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
>> genelist_table.html.txt
# Cell essential & nonessential genes
for categ in essential nonessential; do
  for wrapper in 1; do
    echo -e "Cell $categ genes"
    wc -l cell_${categ}.genes.list \
    | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
    | sed 's/\.genes\.list//g' | addcom
    echo "[Hart _et al._, _G3_, 2017](https://www.g3journal.org/content/7/8/2719)"
    echo "Genes $categ in human cell lines as determined by CRISPR/Cas9 screens. Curated in [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7)"
  done | paste -s \
  | sed 's/\t/\ \|\ /g' \
  | sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
  >> genelist_table.html.txt
done
# Mouse het LoF lethal genes
for wrapper in 1; do
  echo "Mouse heterozygous LoF lethal"
  wc -l mouse_het_lethal.genes.list \
  | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
  | sed 's/\.genes\.list//g' | addcom
  echo "[Motenko _et al._, _Mamm. Genome_, 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4534495/)"
  echo "Genes lethal in mouse models when heterozygously inactivated in mice. Curated in [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7)"
done | paste -s \
| sed 's/\t/\ \|\ /g' \
| sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
>> genelist_table.html.txt
# Mouse dispensable genes
for wrapper in 1; do
  echo "Mouse dispensable"
  wc -l mouse_dispensable.genes.list \
  | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
  | sed 's/\.genes\.list//g' | addcom
  echo "[Motenko _et al._, _Mamm. Genome_, 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4534495/)"
  echo "Genes with no phenotype reported when heterozygously inactivated in mice. Curated in [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7)"
done | paste -s \
| sed 's/\t/\ \|\ /g' \
| sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
>> genelist_table.html.txt
# Olfactory receptor genes
for wrapper in 1; do
  echo "Olfactory receptors"
  wc -l olfactory_receptors.genes.list \
  | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
  | sed 's/\.genes\.list//g' | addcom
  echo "[Braschi _et al._, _Nucleic Acids Res_, 2019](https://pubmed.ncbi.nlm.nih.gov/30304474/)"
  echo "Genes from any HUGO-recognized family of olfactory receptor genes."
done | paste -s \
| sed 's/\t/\ \|\ /g' \
| sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
>> genelist_table.html.txt
# GTEx expression outlier-defined gene sets
for dir in low high; do
  for tol in variable invariant; do
    for wrapper in 1; do
      echo "$( echo "${dir}" | cut -c1 | tr 'a-z' 'A-Z' )$( echo ${dir} | cut -c2- ) expression-${tol} genes"
      cat gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.${dir}_expression_${tol}.genes.list | wc -l | addcom
      echo "\`gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.${dir}_expression_${tol}\`"
      echo "GTex v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597)"
      cbase="≥ 5 mean TPM and 25<sup>th</sup> percentile ≥ 0 in ≥ 3 tissues"
      if [ $dir == low ]; then
        if [ $tol == variable ]; then
          echo -e "${cbase}; mean ≥ 0.5% of samples are low expression outliers; mean TPM CV ≥ 0.5"
        else
          echo -e "${cbase}; mean ≤ 0.1% of samples are low expression outliers; mean TPM CV ≤ 0.3; no samples with TPM < 1 in any tissue"
        fi
      else
        if [ $tol == variable ]; then
          echo -e "${cbase}; mean ≥ 5% of samples are high expression outliers; mean TPM CV ≥ 0.5"
        else
          echo -e "${cbase}; mean ≤ 1% of samples are low expression outliers; mean TPM CV ≤ 0.3"
        fi
      fi
    done | paste -s \
  | sed 's/\t/\ \|\ /g' \
  | sed -e 's/^/\|\ /g' -e 's/$/\ \|/g'
  done
done \
>> genelist_table.html.txt
# HPO-associated genes
while read pheno hpo; do
  for wrapper in 1; do
    descrip=$( fgrep -w $hpo refs/HPOs_by_metacohort.table.tsv | cut -f2 )
    echo "$descrip (${hpo})-associated genes"
    cat $pheno.HPOdb.genes.list | wc -l | addcom
    echo "\`$pheno.HPOdb\`"
    echo "HPO database (accessed $( date | awk '{ print $2, $NF }' )) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478)"
    echo "Genes linked to $hpo"
  done | paste -s \
  | sed 's/\t/\ \|\ /g' \
  | sed -e 's/^/\|\ /g' -e 's/$/\ \|/g'
done < refs/test_phenotypes.list \
>> genelist_table.html.txt
# ClinGen dosage sensitive genes
for dos in haploinsufficient triplosensitive; do
  for conf in hc mc hmc lc all; do
    case $conf in
      "hc")
        clab="high"
        ;;
      "mc")
        clab="medium"
        ;;
      "hmc")
        clab="high or medium"
        ;;
      "lc")
        clab="low"
        ;;
      "all")
        clab="any"
        ;;
    esac
    for wrapper in 1; do
      echo "ClinGen dominant $dos genes (${clab} confidence)"
      cat ClinGen.${conf}_${dos}.genes.list | wc -l | addcom
      echo "\`ClinGen.${conf}_${dos}\`"
      echo "ClinGen gene curation map (accessed $( date | awk '{ print $2, $NF }' )) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198)"
      echo "Dominant $dos genes scored at $clab confidence by ClinGen"
    done | paste -s \
      | sed 's/\t/\ \|\ /g' \
      | sed -e 's/^/\|\ /g' -e 's/$/\ \|/g'
  done
done \
>> genelist_table.html.txt
# DECIPHER genes
for mech in lof gof other; do
  case $mech in
    "lof")
      mdescrip="loss-of-function"
      article="a"
      ;;
    "gof")
      mdescrip="activating/gain-of-function"
      article="an"
      ;;
    "other")
      mdescrip="other/unknown coding"
      article="an"
      ;;
  esac
  for conf in hc mc hmc lc all; do
    case $conf in
      "hc")
        clab="high"
        rating='"confirmed"'
        ;;
      "mc")
        clab="medium"
        rating='"probable"'
        ;;
      "hmc")
        clab="high or medium"
        rating='"confirmed" or "probable"'
        ;;
      "lc")
        clab="low"
        rating='"possible"'
        ;;
      "all")
        clab="any"
        rating='"confirmed", "probable", or "possible"'
        ;;
    esac
    for wrapper in 1; do
      echo "DECIPHER dominant dev. disorder genes (${mdescrip} mechanism; ${clab} confidence)"
      cat DDG2P.${conf}_${mech}.genes.list | wc -l | addcom
      echo "\`DDG2P.${conf}_${mech}\`"
      echo "DECIPHER/DDG2P (accessed $( date | awk '{ print $2, $NF }' )) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582)"
      echo "Dominant dev. disorder genes with $article $mdescrip mechanism scored as $rating in DECIPHER"
    done | paste -s \
      | sed 's/\t/\ \|\ /g' \
      | sed -e 's/^/\|\ /g' -e 's/$/\ \|/g'
  done
done \
>> genelist_table.html.txt
# COSMIC oncogenes
for wrapper in 1; do
  echo "High-confidence oncogenes"
  wc -l COSMIC.hc_oncogenes.genes.list \
  | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
  | sed 's/\.genes\.list//g' | addcom
  echo "COSMIC v91 [Sondka _et al._, _Nat. Rev. Cancer_, 2018](https://www.nature.com/articles/s41568-018-0060-1)"
  echo "Tier 1 dominant oncogenes due to amplification, missense, or 'other' mutations"
done | paste -s \
| sed 's/\t/\ \|\ /g' \
| sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
>> genelist_table.html.txt
cat genelist_table.html.txt

