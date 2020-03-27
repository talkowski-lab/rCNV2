# Dosage Sensitivity Scoring  

We developed a model to predict haploinsufficiency and triplosensitivity for each protein-coding gene across all 22 autosomes. This process is described below.  

The code to reproduce these analyses is contained in `gene_scoring_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `gene_scoring_analysis.wdl`.  

---  

## Dosage sensitivity scoring procedure

We trained & applied a model to predict dosage sensitivity scores for all [autosomal, canonical, protein-coding genes](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene/).  

The steps of this model are described below:  

### 1. Burden test of rare CNVs per gene between cases and controls  

We counted CNVs per gene between all [cases and controls](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/) for each metacohort, and then [meta-analyzed](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/genes/#3-combine-association-statistics-across-metacohorts) those CNV counts to derive an odds ratio for each gene.  

This process was executed identically to gene-level disease association analyses, described in detail [here](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/genes/#gene-based-burden-test-procedure).  

### 2. Empirical Bayes estimation of priors  

We next empirically estimated the average effect size and variance expected for true dosage sensitive genes based on the rCNV data in this study.  

We first defined a set of high-confidence haploinsufficient genes based on the following criteria:
1. 

### 3. Bayes factor calculation for each gene  


### 4. Prediction of dosage sensitivity scores


---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  

