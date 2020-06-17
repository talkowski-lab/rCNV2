# Sliding Window Analyses  

We compared rare CNV (rCNV) counts between cases and controls across all 22 autosomes in regularized sliding windows. This process is described below.  

The code to reproduce these analyses is contained in `sliding_window_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `sliding_window_analysis.wdl`.  

---  

## TL;DR: Summary of sliding window results  

Skipping ahead to the final outcome of this analysis, we identified a set of regions where [rCNVs](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#curation-steps-rare-cnvs) are statistically associated with at least one of [30 disease phenotypes](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype#tldr-final-outcome-of-phenotype-curation) at genome-wide significance.  

The details of this procedure are described below, but summary plots for the final results can be generated using `regions_summary.plot.R`, which includes the following multi-panel figure:  

![rCNV Region Summary Stats](https://storage.googleapis.com/rcnv_project/public/rCNV.final_segments.multipanel_summary.jpg)  

---  

## Sliding window analysis procedure

We executed a standardized procedure to conduct sliding window analyses for each [phenotype group](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/).  

The four steps of this procedure are described below:  

### 1. Count rare CNVs per window  

For each [metacohort](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/), we intersected [rCNVs](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#curation-steps-rare-cnvs) against [genome-wide sliding windows](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/binned_genome/) separately for cases and controls.  

We conducted this procedure twice per phenotype group & metacohort: once each for deletions and duplications.  

The code to perform this step is contained in `count_cnvs_per_window.py`.  

We required CNVs to overlap at least 50% of each window to be counted.  

Finally, given that CNV breakpoint precision varies by locus, platform, and CNV calling algorithm, we systematically extended the breakpoints of each control CNV by +50kb. CNV breakpoints in cases were **not** extended. We did this to conservatively protect against spurrious associations arising from situations where control CNVs breakpoints might be underestimated compared to case breakpoints.  

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

All plots also feature a horizontal dashed line indicating Bonferroni-corrected significance threshold (P ≤ 3.72x10<sup>-6</sup>) for all non-overlapping 200kb windows tested.  

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

We combined CNV association statistics across metacohorts for each sliding window using a fixed-effects meta-analysis. The P-value from this model was designated as the "`primary`" P-value.  

We also computed a "`secondary`" P-value, which was calculated from an identical fixed-effects meta-analysis model after excluding the single most significant metacohort per window.  

The code to perform this step is contained in `window_meta_analysis.R`.  

Given that rare CNV counts per window are (a) sparse and (b) zero-inflated, and furthermore that (c) the case & control sample sizes are unbalanced for most phenotype groups (_e.g._, frequently >10- to 100-fold more controls than cases), we implemented an empirical continuity correction as proposed by [Sweeting _et al._, _Stat. Med._, 2004.](https://onlinelibrary.wiley.com/doi/10.1002/sim.1761)  

Furthermore, as sample size imbalance between cases and controls has been shown to distort test statistics in genetic association studies, we applied a saddlepoint approximation correction to the test statistics for each phenotype per CNV type to recalibrate our primary and secondary P-values. This procedure is based on [Dey et al., _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5501775/).  

Each phenotype & CNV type were meta-analyzed separately for a total of two meta-analyses per phenotype.  

#### Determining & calibrating genome-wide significance threshold  

We next controlled false discovery rate (FDR) across all sliding window meta-analyses.  

Estimating the number of independent tests performed across all windows, phenotypes, and CNV types is difficult due to numerous necessary assumptions, such as the independence of samples between phenotypes, or the local correlation structure of CNV counts between neighboring windows, among others.    

Instead, we used a "genome-wide" significance threshold of P ≤ 3.72x10<sup>-6</sup>, which corresonds to a Bonferroni correction if applied to the number of non-overlapping 200kb windows tested in our analysis.  

We empirically assessed the calibration of this primary P-value threshold using a permutation-based approach, as follows:  

1. permute phenotype labels for all CNVs while matching on size (split by quantile) and CNV type (DEL/DUP);  
2. rerun all association tests (described above), including meta-analysis, for each phenotype & CNV combination;   
3. compute the fraction of significant windows for a broad range of primary P-value thresholds (_e.g._, -log<sub>10</sub>(_P_) ~ [0, 20] ); and
4. report the least significant primary P-value threshold that results in an empirical FDR ≤ 3.72x10<sup>-6</sup>.  

Steps 1-3 were repeated 50 times for each CNV type and phenotype, and the median value per phenotype was computed.

#### Output files  

[As described above for Step 2](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/sliding_windows#output-files), we generated the same combination of plots and statistics files for the meta-analyses results of each phenotype group.  

Unlike the plots from Step 2, the dashed lines from Step 3 correspond to empirically-derived genome-wide significance thresholds [as determined via permutation](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/sliding_windows#determining-genome-wide-significance-threshold).  

These files are stored in the same location as the per-metacohort analysis results.  

### 4. Refine associations to minimal credible regions & merge across phenotypes  

Lastly, we refined each association to the minimum interval(s) predicted to contain the causal factor(s).  

In practice, we considered a window to be genome-wide significant if its primary P-value P exceeded the genome-wide significance threshold, and it satisifed at least one of the following two criteria:
1. Secondary P-value (as [described above](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/sliding_windows/#3-combine-association-statistics-across-metacohorts)) was also nominally significant (P < 0.05); and/or
2. At least two metacohorts were nominally significant (P < 0.05) per Fisher's exact test.  

Absent a true replication sample, these _post hoc_ filters were required to protect against Winner's Curse.  

Next, for each locus with a significant association between rCNVs and one or more phenotypes, we aimed to identify the minimal interval(s) that contained the causal element(s) with 99% confidence.  

This procedure is described below, and was performed separately for deletions and duplications using `refine_significant_regions.py`.  

First, per phenotype, all significant windows were collapsed into nonredundant `blocks` by merging all windows within ±200kb of any significant window.  

Next, for each block, we performed the following steps:
1. Compute an approximate Bayes factor (ABF) for each window following the procedure specified by [Wakefield _et al._, _Genet. Epi._, 2009](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20359); 
2. Identify the set of windows that captured at least 99% of the total ABF for the block (_i.e._, define the 99% credible set); and
3. Merge all windows in the 99% credible set (with ±200kb padding) to define the set of credible intervals comprising the 99% credible set.  

For each refined association, we reported the following:  
* Case and control CNV frequencies were computed as the mean CNV carrier rate across all windows in the 99% credible set; and
* Pooled effect sizes were computed as the inverse variance-weighted mean across all windows in the 99% credible set.

After credible intervals were defined for each association, we collapsed associations across phenotypes to derive a final set of nonredundant disease-associated segments for downstream analyses.  

To reduce associations to loci, we clustered any overlapping credible intervals of the same CNV type (deletion or duplication) from different phenotypes. For each locus, we reported the following:  
* Case and control CNV frequencies were computed as the weighted mean CNV carrier rate per 99% credible set, where the weights corresponded to the square root of the sample size for each phenotype; and  
* Pooled effect sizes were computed as the inverse variance-weighted mean across all significant phenotypes.  

#### Output files  

From this analysis, we generated two files:  
1. `rCNV.final_segments.loci.bed.gz`: a bgzipped BED file containing one row for each disease-associated locus, including summary association statistics pooled across all significant phenotypes.  
2. `rCNV.final_segments.associations.bed.gz`: a bgzipped BED file containing one row for each phenotype-locus pair, including association statistics for that specific phenotype.  

These files are stored in a protected Google Cloud bucket, here:
```
$ gsutil ls gs://rcnv_project/results/segment_association/

gs://rcnv_project/results/segment_association/rCNV.final_segments.associations.bed.gz
gs://rcnv_project/results/segment_association/rCNV.final_segments.loci.bed.gz
```

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  

