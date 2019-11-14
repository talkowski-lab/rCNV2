# Sliding Window Analyses  

We compared rare CNV counts between cases and controls across all 22 autosomes in regularized sliding windows. This process is described below.  

The code to reproduce these analyses is contained in `sliding_window_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `sliding_window_analysis.wdl`.  

## Sliding window analysis procedure

We executed a standardized procedure to conduct sliding window analyses for each [phenotype group](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/).  

The steps of this procedure are described below:  

### 1. Counting rare CNVs per window  

For each [metacohort](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/), we intersected [ultra-rare CNVs](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#curation-steps-rare-cnvs) against [genome-wide sliding windows](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/binned_genome/) separately for cases and controls.  

We conducted this procedure a total of three times per phenotype group & metacohort: once each for deletions, duplications, and all CNVs (deletions + duplications).  

The code to perform this step is contained in `count_cnvs_per_window.py`.  

We required CNVs to overlap at least 50% of each window to be considered.  

We selected a cutoff of 50% based on the following two facts: 
1. the smallest CNV size considered in this study is 100kb; and  
2. the window size used for this analysis was 200kb.  

Thus, a 100kb CNV wholly contained within a 200kb window would count towards that window at 50% minimum overlap.  

Finally, given that CNV breakpoint precision varies by locus, platform, and CNV calling algorithm, we systematically extended the breakpoints of each control CNV by +50kb. We did this to conservatively protect against loci where control CNVs breakpoints were underestimated compared to cases.  

### 2. Calculate burden statistics between cases & controls  

Once CNVs were tallied in cases and controls for each window, we next compared the ratios of CNV carriers between cases and controls using a one-sided Fisher's exact test.  

The code to perform this step is contained in `window_burden_test.R`.  

#### Output files  

For each combination of phenotype group, metacohort, and CNV type, the following files are generated:  

1. `$metacohort.$hpo.rCNV.$CNV_type.sliding_window.stats.bed.gz`: a bgzipped BED file containing CNV counts and association statistics for each window  
2. `$metacohort.$hpo.rCNV.$CNV_type.sliding_window.manhattan.png`: a Manhattan plot of association statistics for each window
3. `$metacohort.$hpo.rCNV.$CNV_type.sliding_window.qq.png`: a QQ plot of observed P-values compared to expected P-values under a uniform null  
4. `$metacohort.$hpo.rCNV.$CNV_type.sliding_window.manhattan_with_qq.png`: a two-panel composite plot combining plots `2` and `3`, above  

Furthermore, for each pair of phenotype group & metacohort, two additional files are generated:  
5. `$metacohort.$hpo.rCNV.sliding_window.miami.png`: a Miami plot of association statistics for each window, with duplications above and deletions below the x-axis  
6. `$metacohort.$hpo.rCNV.sliding_window.miami_with_qq.png`: a multi-panel composite plot of association statistics, combining a Miami plot with QQ plots for deletions and duplications separately  

All `.png` plots are annotated with 54 known genomic disorder loci (adapted from [Owen _et al._, _BMC Genomics_, 2018](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5292-7)), for reference.  

### 3. Combine association statistics across metacohorts  

We combined CNV association statistics across metacohorts for each sliding window using the Mantel-Haenszel meta-analysis method for fourfold-tables of count data.  

This model was implemented using `rma.mh()` from [the `metafor` package in R](https://cran.r-project.org/web/packages/metafor/metafor.pdf).  

Given that rare CNV counts per window are (_i_) sparse and (_ii_) zero-inflated, and furthermore that (_iii_) the case & control sample sizes are unbalanced for most phenotype groups (_e.g._, frequently >10- to 100-fold more controls than cases), we implemented an empirical continuity correction as proposed by [Sweeting _et al._, _Stat. Med._, 2004.](https://onlinelibrary.wiley.com/doi/10.1002/sim.1761)  

Each phenotype & CNV were meta-analyzed separately for a total of three meta-analyses per phenotype.  

Following meta-analysis, individual windows were labeled as genome-wide significant if the met the following three criteria:
1. _P_<sub>meta</sub> ≤ 10<sup>-6</sup> (_note: threshold of 10<sup>-6</sup> corresponds to a study-wide FDR of 1% via phenotype permutation analysis_);  
2. Lower bound of 95% confidence interval of meta-analysis odds ratio ≥ 2; and  
3. Nominally significant in ≥ 2 metacohorts.  

### 4. Refine associations to minimal credible regions  

TBD

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  
