# Sliding Window Analyses  

We compared CNV counts between cases and controls across all 22 autosomes in regularized sliding windows. This process is described below.  

The code to reproduce these analyses is contained in `sliding_window_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `sliding_window_analysis.wdl`.  

## Sliding window analysis procedure

We performed a standard procedure to conduct sliding window analyses for each [phenotype group](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/).  

The steps in this procedure are described below:  

### 1. Counting CNVs per window  

For each [metacohort](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/), we intersected rare CNVs against [genome-wide sliding windows](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/binned_genome/) separately for cases and controls.  

The code to perform this step is contained in `count_cnvs_per_window.py`.  

We required CNVs to overlap at least 50% of each window to be considered.  

We selected a cutoff of 50% based on the following two facts: 
1. the smallest CNV size considered in this study is 100kb; and  
2. the window size used for this analysis was 200kb.  

Thus, a 100kb CNV wholly contained within a 200kb window would count towards that window at 50% minimum overlap.  

### 2. Calculate burden statistics between cases & controls  

Once CNVs were tallied in cases and controls for each window, we next compared the ratios of CNV carriers between cases and controls using a one-sided Fisher's exact test.  

The code to perform this step is contained in `window_burden_test.R`.  

### 3. Combine statistics across metacohorts (_where applicable_)  

For phenotype groups with at least 100 cases represented in two or more metacohorts (see [the last table on this page](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/)), we XYZ...

### 4. Collapse overlapping significant windows  

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets or references.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in the [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  
