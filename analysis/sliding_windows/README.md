# Sliding Window Association Analyses  

We compared rare CNV (rCNV) counts between cases and controls across all 22 autosomes in regularized sliding windows. This process is described in the methods of Collins _et al._, 2022.  

The code to reproduce these analyses is contained in `sliding_window_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `sliding_window_analysis.wdl`.  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  

