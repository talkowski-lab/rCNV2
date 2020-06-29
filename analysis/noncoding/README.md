# Noncoding Association Analyses  

We compared rare CNV (rCNV) counts between cases and controls for clusters of _cis_-regulatory elements (_cis_-regulatory blocks, or CRBs) across all 22 autosomes. This process is described below.  

The code to reproduce these analyses is contained in `noncoding_association_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `noncoding_association_analysis.wdl`.  

---  

## Preprocessing  

There were two preprocessing steps used exclusively for these noncoding association analyses:

1. Noncoding CNV filtering ([described here](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#noncoding-subsets))  
2. CRB definition ([described here](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/genome_annotations))  

---  

## Noncoding association test procedure  

We executed a standardized procedure to conduct CRB-based burden tests for each [phenotype group](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/).  

The steps of this procedure are described below:  

### 1. Count rCNVs per CRB  

For each [metacohort](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/), we intersected [rCNVs](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#curation-steps-rare-cnvs) against [CRBs](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/genome_annotations) separately for cases and controls.  

As [described here](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#noncoding-subsets), rCNVs were restricted to a loose definition of "noncoding," where each rCNV was allowed to intersect only exons from a subset of [likely unconstrained genes](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene#gene-set-definitions) (or no genes at all; _i.e._, completely noncoding), but no other genes.  

Finally, given that CNV breakpoint precision varies by locus, platform, and CNV calling algorithm, we systematically extended the breakpoints of each control CNV by +50kb. CNV breakpoints in cases were **not** extended. We did this to conservatively protect against spurrious associations arising from situations where control CNVs breakpoints might be underestimated compared to case breakpoints.  

We conducted this procedure a total of two times per phenotype group & metacohort: once each for deletions and duplications.  

The code to perform this step is contained in `count_cnvs_per_crb.py`.  

We intersected noncoding rCNVs versus the individual elements from each CRB, and only counted rCNV-element pairs where the rCNV completely (100%) covered the element.  

After intersecting rCNVs with individual elements, we only tallied rCNV-CRB pairs where the rCNV hit at least 5% of all elements from that CRB.  

### 2. Calculate burden statistics between cases & controls  

Once CNVs were tallied in cases and controls for each CRB, we next compared the ratios of CNV carriers between cases and controls using a one-sided Fisher's exact test.  

The code to perform this step is contained in `crb_burden_test.R`.  

#### Output files  

For each combination of phenotype group, metacohort, and CNV type, the following files are generated:  

1. `$metacohort.$hpo.rCNV.loose_noncoding.$CNV_type.crb_burden.stats.bed.gz`: a bgzipped BED file containing CNV counts and association statistics for each CRB  
2. `$metacohort.$hpo.rCNV.loose_noncoding.$CNV_type.crb_burden.manhattan.png`: a Manhattan plot of association statistics for each CRB
3. `$metacohort.$hpo.rCNV.loose_noncoding.$CNV_type.crb_burden.qq.png`: a QQ plot of observed P-values compared to expected P-values under a uniform null  
4. `$metacohort.$hpo.rCNV.loose_noncoding.$CNV_type.crb_burden.manhattan_with_qq.png`: a two-panel composite plot combining plots `2` and `3`, above  

Furthermore, for each pair of phenotype group & metacohort, two additional files are generated:  

5. `$metacohort.$hpo.rCNV.crb_burden.miami.png`: a Miami plot of association statistics for each CRB, with duplications above and deletions below the x-axis  
6. `$metacohort.$hpo.rCNV.crb_burden.miami_with_qq.png`: a multi-panel composite plot of association statistics, combining a Miami plot with QQ plots for deletions and duplications separately  

All plots also feature a horizontal dashed line indicating Bonferroni-corrected significance threshold (P ≤ 3.76x10<sup>-6</sup>) for all 13,315 CRBs tested.  

These files are stored in a protected Google Cloud bucket with one subdirectory per HPO group, here:  
```
$ gsutil ls gs://rcnv_project/analysis/crb_burden/

gs://rcnv_project/analysis/crb_burden/HP0000118/
gs://rcnv_project/analysis/crb_burden/HP0000152/
gs://rcnv_project/analysis/crb_burden/HP0000707/
gs://rcnv_project/analysis/crb_burden/HP0000708/
gs://rcnv_project/analysis/crb_burden/HP0000717/
gs://rcnv_project/analysis/crb_burden/HP0000729/
gs://rcnv_project/analysis/crb_burden/HP0000752/
gs://rcnv_project/analysis/crb_burden/HP0000924/
gs://rcnv_project/analysis/crb_burden/HP0001197/
gs://rcnv_project/analysis/crb_burden/HP0001250/
gs://rcnv_project/analysis/crb_burden/HP0001507/
gs://rcnv_project/analysis/crb_burden/HP0001626/
gs://rcnv_project/analysis/crb_burden/HP0001627/
gs://rcnv_project/analysis/crb_burden/HP0002011/
gs://rcnv_project/analysis/crb_burden/HP0002597/
gs://rcnv_project/analysis/crb_burden/HP0002715/
gs://rcnv_project/analysis/crb_burden/HP0002960/
gs://rcnv_project/analysis/crb_burden/HP0003011/
gs://rcnv_project/analysis/crb_burden/HP0011446/
gs://rcnv_project/analysis/crb_burden/HP0012443/
gs://rcnv_project/analysis/crb_burden/HP0012638/
gs://rcnv_project/analysis/crb_burden/HP0012639/
gs://rcnv_project/analysis/crb_burden/HP0012759/
gs://rcnv_project/analysis/crb_burden/HP0025031/
gs://rcnv_project/analysis/crb_burden/HP0031466/
gs://rcnv_project/analysis/crb_burden/HP0100022/
gs://rcnv_project/analysis/crb_burden/HP0100545/
gs://rcnv_project/analysis/crb_burden/HP0100753/
gs://rcnv_project/analysis/crb_burden/HP0100852/
gs://rcnv_project/analysis/crb_burden/UNKNOWN/
```

### 3. Combine association statistics across metacohorts  

We combined CNV association statistics across metacohorts for each CRB using a fixed-effects meta-analysis. The P-value from this model was designated as the "`primary`" P-value.  

We also computed a "`secondary`" P-value, which was calculated from an identical fixed-effects meta-analysis model after excluding the single most significant metacohort per gene.  

The code to perform this step is contained in `crb_meta_analysis.R`.  

Given that rare CNV counts per gene are (a) sparse and (b) zero-inflated, and furthermore that (c) the case & control sample sizes are unbalanced for most phenotype groups (_e.g._, frequently >10- to 100-fold more controls than cases), we implemented an empirical continuity correction as proposed by [Sweeting _et al._, _Stat. Med._, 2004.](https://onlinelibrary.wiley.com/doi/10.1002/sim.1761)  

Furthermore, as sample size imbalance between cases and controls has been shown to distort test statistics in genetic association studies, we applied a saddlepoint approximation correction to the test statistics for each phenotype per CNV type to recalibrate our primary and secondary P-values. This procedure is based on [Dey et al., _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5501775/).  

Each phenotype & CNV type were meta-analyzed separately for a total of two meta-analyses per phenotype.  

#### Determining & calibrating genome-wide significance threshold  

We next controlled false discovery rate (FDR) across all gene meta-analyses.  

Estimating the number of independent tests performed across all CRBs, phenotypes, and CNV types is difficult due to numerous necessary assumptions, such as the independence of samples between phenotypes, or the local correlation structure of CNV counts between neighboring CRBs, among others.    

Instead, we used a genome-wide significance threshold of P ≤ 3.76x10<sup>-6</sup>, which corresonds to a Bonferroni correction if applied to the number of individual CRBs tested in our analysis (N=13,315).  

We empirically assessed the calibration of this primary P-value threshold using a permutation-based approach, as follows:  

1. permute phenotype labels for all loose noncoding rCNVs while matching on size (split by quantile) and CNV type (DEL/DUP);  
2. rerun all association tests (described above), including meta-analysis, for each phenotype & CNV combination;   
3. compute the fraction of significant CRBs for a broad range of primary P-value thresholds (_e.g._, -log<sub>10</sub>(_P_) ~ [0, 20] ); and
4. report the least significant primary P-value threshold that results in an empirical FDR ≤ 3.76x10<sup>-6</sup>.  

Steps 1-3 were repeated 50 times for each CNV type and phenotype, and the median value per phenotype was computed.

#### Output files  

[As described above for Step 2](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/noncoding#output-files), we generated the same combination of plots and statistics files for the meta-analyses results of each phenotype group.  

Unlike the plots from Step 2, the dashed lines from Step 3 correspond to empirically-derived exome-wide significance thresholds [as determined via permutation](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/noncoding#determining-and-calibrating-genome-wide-significance-threshold).  

These files are stored in the same location as the per-metacohort analysis results.  

#### 4. Refine correlated associations  

#### 5. Reporting of final association statistics  

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  
