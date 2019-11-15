# Sliding Window Analyses  

We compared rare CNV counts between cases and controls across all 22 autosomes in regularized sliding windows. This process is described below.  

The code to reproduce these analyses is contained in `sliding_window_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `sliding_window_analysis.wdl`.  

## Sliding window analysis procedure

We executed a standardized procedure to conduct sliding window analyses for each [phenotype group](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/).  

The four steps of this procedure are described below:  

### 1. Count rare CNVs per window  

For each [metacohort](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/), we intersected [rare CNVs](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#curation-steps-rare-cnvs) against [genome-wide sliding windows](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/binned_genome/) separately for cases and controls.  

We conducted this procedure a total of three times per phenotype group & metacohort: once each for deletions, duplications, and all CNVs (deletions + duplications).  

The code to perform this step is contained in `count_cnvs_per_window.py`.  

We required CNVs to overlap at least 50% of each window to be counted.  

Finally, given that CNV breakpoint precision varies by locus, platform, and CNV calling algorithm, we systematically extended the breakpoints of each control CNV by +50kb. Note that CNV breakpoints in cases were **not** extended. We did this to conservatively protect against spurrious associations arising from control CNVs breakpoints being systematically underestimated compared to case breakpoints.  

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

These files are stored in a protected Google Cloud bucket with one subdirectory per HPO group, here:  
```
$ gsutil ls gs://rcnv_project/analysis/sliding_windows/

gs://rcnv_project/analysis/sliding_windows/HP0000118/
gs://rcnv_project/analysis/sliding_windows/HP0000152/
gs://rcnv_project/analysis/sliding_windows/HP0000707/
gs://rcnv_project/analysis/sliding_windows/HP0000708/
gs://rcnv_project/analysis/sliding_windows/HP0000717/
gs://rcnv_project/analysis/sliding_windows/HP0000729/
gs://rcnv_project/analysis/sliding_windows/HP0000752/
gs://rcnv_project/analysis/sliding_windows/HP0000924/
gs://rcnv_project/analysis/sliding_windows/HP0001197/
gs://rcnv_project/analysis/sliding_windows/HP0001250/
gs://rcnv_project/analysis/sliding_windows/HP0001507/
gs://rcnv_project/analysis/sliding_windows/HP0001626/
gs://rcnv_project/analysis/sliding_windows/HP0001627/
gs://rcnv_project/analysis/sliding_windows/HP0002011/
gs://rcnv_project/analysis/sliding_windows/HP0002597/
gs://rcnv_project/analysis/sliding_windows/HP0002715/
gs://rcnv_project/analysis/sliding_windows/HP0002960/
gs://rcnv_project/analysis/sliding_windows/HP0003011/
gs://rcnv_project/analysis/sliding_windows/HP0011446/
gs://rcnv_project/analysis/sliding_windows/HP0012443/
gs://rcnv_project/analysis/sliding_windows/HP0012638/
gs://rcnv_project/analysis/sliding_windows/HP0012639/
gs://rcnv_project/analysis/sliding_windows/HP0012759/
gs://rcnv_project/analysis/sliding_windows/HP0025031/
gs://rcnv_project/analysis/sliding_windows/HP0031466/
gs://rcnv_project/analysis/sliding_windows/HP0100022/
gs://rcnv_project/analysis/sliding_windows/HP0100545/
gs://rcnv_project/analysis/sliding_windows/HP0100753/
gs://rcnv_project/analysis/sliding_windows/HP0100852/
gs://rcnv_project/analysis/sliding_windows/UNKNOWN/
```

### 3. Combine association statistics across metacohorts  

We combined CNV association statistics across metacohorts for each sliding window using the Mantel-Haenszel meta-analysis method for 2x2 contingency tables of count data.  

The code to perform this step is contained in `window_meta_analysis.R`.  

Given that rare CNV counts per window are (a) sparse and (b) zero-inflated, and furthermore that (c) the case & control sample sizes are unbalanced for most phenotype groups (_e.g._, frequently >10- to 100-fold more controls than cases), we implemented an empirical continuity correction as proposed by [Sweeting _et al._, _Stat. Med._, 2004.](https://onlinelibrary.wiley.com/doi/10.1002/sim.1761)  

Each phenotype & CNV type were meta-analyzed separately for a total of three meta-analyses per phenotype.  

Following meta-analysis, individual windows were labeled as genome-wide significant if the met the following three criteria:
1. _P_<sub>meta</sub> ≤ 10<sup>-6</sup>;  
2. Lower bound of 95% confidence interval of odds ratio ≥ 2; and  
3. Nominally significant in ≥ 2 metacohorts.  

_Note: a P<sub>meta</sub> threshold of 10<sup>-6</sup> was determined to correspond to a study-wide FDR of 1% via phenotype permutation analysis (not described here)_

#### Output files  

[As described above for Step 2](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/sliding_windows#output-files), we generated the same combination of plots and statistics files for the meta-analyses results of each phenotype group.  

These files are stored in the same location as the per-metacohort analysis results.  

### 4. Collapse associations across phenotypes & refine loci to minimal credible regions  

Lastly, we collapsed all significant windows across phenotypes and refined them to discrete intervals.  

For each locus with significant evidence of rCNVs being associated with one or more phenotypes, we aimed to identify:  
1. all statistically independent associations that were _individually_ genome-wide significant; and
3. refine these associations to the minimal interval that contained the causal element(s) with 90% confidence (_i.e._, analogous to 90% credible sets from common variant GWAS).  

This procedure is described below, and was performed separately for deletions and duplications using `refine_significant_regions.py`.  

First, all windows significant in at least one phenotype were collapsed into a list of larger nonredundant `query regions` to be refined. These query regions were created by merging overlapping significant windows, and clustering all merged windows within ±1Mb. Finally, query regions were extended by ±1Mb on each side, such that no significant window was closer than ±1Mb to the edge of each query region.  

Next, for each query region, we performed the following steps:
1. Designated the significant window with strongest P-value from any phenotype as the `sentinel window`.  
2. Designated all phenotypes with at least one significant association in the region as the `sentinel phenotypes`.  
3. Gathered all rCNVs overlapping the sentinel window in sentinel phenotypes & controls.  
4. Ordered all case CNVs based on the smallest max breakpoint distance to middle of sentinel window.  
5. Added case CNVs one at a time (in order) and ran a Mantel-Haenszel meta-analysis with Sweeting correction against all control CNVs until the test reaches P ≤ 10-6 and at least two metacohorts have Fisher’s exact P ≤ 0.05.  
  * _Note: Step 5 was attempted first with CNVs covering at least half of the sentinel window, and restricted to case CNVs ≤ 3Mb in size. If genome-wide significance was not achieved for this subset, the test was expanded to include all CNVs overlapping the sentinel window irrespective of size or window overlap._  
6. When the M-H test reached significance, we defined the `minimal credible region` as the middle 90% of the CNV density for the minimal set of case CNVs required to achieve genome-wide significance for the sentinel window.  
7. Exclude all sentinel CNVs and other case CNVs at least 50% covered by the minimal credible region from step 6.  
8-N. Recompute all association statistics for each phenotype at every window within the query region after conditioning on the removal of case CNVs from step 7. Repeat steps 1-7 until no more genome-wide significant windows exist within the query region.  

The output of this process was a set of minimal credible regions, where each credible region contained at least one genome-wide significant association between rCNVs and a phenotype group.  

Each minimal credible region was:
1. 90% confident to contain the causal element(s) for that association, and 
2. Statistically independent from all other minimal credible regions.  

#### Output files  

For each CNV type, we generated three files:  
1. `rCNV.$CNV.final_regions.loci.bed.gz`: a bgzipped BED file containing one row for each minimal credible region, including with summary association statistics pooled across all associated phenotype groups.  
2. `rCNV.$CNV.final_regions.associations.bed.gz`: a bgzipped BED file containing one row for each phenotype-region pair, including association statistics for that specific phenotype versus (i) the entire minimal credible region and (ii) the sentinel (peak) window within the credible region.  
3. `rCNV.$CNV.region_refinement.log`: a flat text log file containing details about the refinement process for each query region.  

These files are stored in a protected Google Cloud bucket, here:
```
$ gsutil ls gs://rcnv_project/results/sliding_windows/

gs://rcnv_project/results/sliding_windows/rCNV.DEL.final_regions.loci.bed.gz
gs://rcnv_project/results/sliding_windows/rCNV.DEL.final_regions.associations.bed.gz
gs://rcnv_project/results/sliding_windows/rCNV.DEL.region_refinement.log
gs://rcnv_project/results/sliding_windows/rCNV.DUP.final_regions.loci.bed.gz
gs://rcnv_project/results/sliding_windows/rCNV.DUP.final_regions.associations.bed.gz
gs://rcnv_project/results/sliding_windows/rCNV.DUP.region_refinement.log
```

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  
