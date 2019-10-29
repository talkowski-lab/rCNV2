# Gene-Based Analyses  

We compared CNV counts between cases and controls for each of 20,243 canonical, protein-coding genes across all 22 autosomes. This process is described below.  

The code to reproduce these analyses is contained in `gene_burden_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `gene_burden_analysis.wdl`.  

## Gene-based burden test procedure

We executed a standardized procedure to conduct gene-based burden tests for each [phenotype group](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/).  

The steps of this procedure are described below:  

### 1. Counting weighted CNVs per gene  

For each [metacohort](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/), we intersected rare CNVs against [all canonical, autosomal, protein-coding genes](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene/) separately for cases and controls.  

We conducted this procedure a total of three times per phenotype group & metacohort: once each for deletions, duplications, and all CNVs (deletions + duplications).  

The code to perform this step is contained in `count_cnvs_per_gene.py`.  

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  
