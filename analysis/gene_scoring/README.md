# Dosage Sensitivity Scoring  

We developed a model to predict haploinsufficiency and triplosensitivity for each protein-coding gene across all 22 autosomes. This process is described in the methods of Collins _et al._, 2022.  

The code to reproduce these analyses is contained in `gene_scoring_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `gene_scoring_analysis.wdl`.  

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  

