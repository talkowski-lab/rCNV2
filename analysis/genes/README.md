# Gene-Based Analyses  

We compared rare CNV (rCNV) counts between cases and controls for each protein-coding gene across all 22 autosomes. This process is described below.  

The code to reproduce these analyses is contained in `gene_burden_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `gene_burden_analysis.wdl`.  

---  

## Gene-based burden test procedure

We executed a standardized procedure to conduct gene-based burden tests for each [phenotype group](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/).  

The steps of this procedure are described below:  

### 1. Count rare CNVs per gene  

For each [metacohort](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/), we intersected [rCNVs](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#curation-steps-rare-cnvs) against [canonical, autosomal, protein-coding genes](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene/) separately for cases and controls.  

We further excluded genes if their canonical transcript met any of the following criteria:
1. Less than 30% covered by somatically hypermutable sites (as applied in [Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674));  
2. Less than 30% covered by segmental duplications and/or simple, low-complexity, or satellite repeats; and  
3. Less than 30% covered by N-masked regions of the hg19 reference genome assembly.  

We also excluded any exons expressed in <20% of transcripts across all tissues in GTEx (see [Cummings _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/10.1101/554444v1) for details).

After all filtering, we retained 170,422 exons from 17,263 genes for these analyses.  

We conducted this procedure a total of three times per phenotype group & metacohort: once each for deletions, duplications, and all CNVs (deletions + duplications).  

The code to perform this step is contained in `count_cnvs_per_gene.py`.  

For each CNV-gene pair, we computed the total fraction of exonic bases (CDS) from the gene overlapped by the CNV.  

We considered a CNV to overlap a gene based on its CDS overlap and CNV type, as follows:  

| CNV type | Min. CDS overlap | Rationale |  
| :--- | ---: | :--- |  
| DEL | ≥10% | Most coding deletions should result in loss-of-function, but we wanted to impose a liberal minimum overlap to protect against spurious annotations given the coarse resolution of most microarrays. |  
| DUP | ≥50% | Functional annotation of duplications is more challenging than deletions, so we imposed a stricter CDS overlap to isolate CNVs predicted to duplicate most (if not all) of a gene.  |  


### 2. Calculate burden statistics between cases & controls  

Once CNVs were tallied in cases and controls for each gene, we next compared the ratios of CNV carriers between cases and controls using a one-sided Fisher's exact test.  

The code to perform this step is contained in `gene_burden_test.R`.  

#### Output files  

For each combination of phenotype group, metacohort, and CNV type, the following files are generated:  

1. `$metacohort.$hpo.rCNV.$CNV_type.gene_burden.stats.bed.gz`: a bgzipped BED file containing CNV counts and association statistics for each gene  
2. `$metacohort.$hpo.rCNV.$CNV_type.gene_burden.manhattan.png`: a Manhattan plot of association statistics for each gene
3. `$metacohort.$hpo.rCNV.$CNV_type.gene_burden.qq.png`: a QQ plot of observed P-values compared to expected P-values under a uniform null  
4. `$metacohort.$hpo.rCNV.$CNV_type.gene_burden.manhattan_with_qq.png`: a two-panel composite plot combining plots `2` and `3`, above  

Furthermore, for each pair of phenotype group & metacohort, two additional files are generated:  

5. `$metacohort.$hpo.rCNV.gene_burden.miami.png`: a Miami plot of association statistics for each gene, with duplications above and deletions below the x-axis  
6. `$metacohort.$hpo.rCNV.gene_burden.miami_with_qq.png`: a multi-panel composite plot of association statistics, combining a Miami plot with QQ plots for deletions and duplications separately  

All `.png` plots are annotated with a set of loss-of-function constrained genes associated with each HPO code. See [`gene curation/`](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene#gene-set-definitions) for more details about these gene sets.  

All plots also feature a horizontal dashed line indicating Bonferroni-corrected significance threshold (P ≤ 2.90x10<sup>-6</sup>) for all 17,263 windows tested.  

These files are stored in a protected Google Cloud bucket with one subdirectory per HPO group, here:  
```
$ gsutil ls gs://rcnv_project/analysis/gene_burden/

gs://rcnv_project/analysis/gene_burden/HP0000118/
gs://rcnv_project/analysis/gene_burden/HP0000152/
gs://rcnv_project/analysis/gene_burden/HP0000707/
gs://rcnv_project/analysis/gene_burden/HP0000708/
gs://rcnv_project/analysis/gene_burden/HP0000717/
gs://rcnv_project/analysis/gene_burden/HP0000729/
gs://rcnv_project/analysis/gene_burden/HP0000752/
gs://rcnv_project/analysis/gene_burden/HP0000924/
gs://rcnv_project/analysis/gene_burden/HP0001197/
gs://rcnv_project/analysis/gene_burden/HP0001250/
gs://rcnv_project/analysis/gene_burden/HP0001507/
gs://rcnv_project/analysis/gene_burden/HP0001626/
gs://rcnv_project/analysis/gene_burden/HP0001627/
gs://rcnv_project/analysis/gene_burden/HP0002011/
gs://rcnv_project/analysis/gene_burden/HP0002597/
gs://rcnv_project/analysis/gene_burden/HP0002715/
gs://rcnv_project/analysis/gene_burden/HP0002960/
gs://rcnv_project/analysis/gene_burden/HP0003011/
gs://rcnv_project/analysis/gene_burden/HP0011446/
gs://rcnv_project/analysis/gene_burden/HP0012443/
gs://rcnv_project/analysis/gene_burden/HP0012638/
gs://rcnv_project/analysis/gene_burden/HP0012639/
gs://rcnv_project/analysis/gene_burden/HP0012759/
gs://rcnv_project/analysis/gene_burden/HP0025031/
gs://rcnv_project/analysis/gene_burden/HP0031466/
gs://rcnv_project/analysis/gene_burden/HP0100022/
gs://rcnv_project/analysis/gene_burden/HP0100545/
gs://rcnv_project/analysis/gene_burden/HP0100753/
gs://rcnv_project/analysis/gene_burden/HP0100852/
gs://rcnv_project/analysis/gene_burden/UNKNOWN/
```

### 3. Combine association statistics across metacohorts  

We combined CNV association statistics across metacohorts for each gene using a random-effects meta-analysis. The P-value from this model was designated as the "`primary`" P-value.  

We also computed a "`secondary`" P-value, which was calculated from an identical random-effects meta-analysis model after excluding the single most significant metacohort per gene.  

The code to perform this step is contained in `gene_meta_analysis.R`.  

Given that rare CNV counts per gene are (a) sparse and (b) zero-inflated, and furthermore that (c) the case & control sample sizes are unbalanced for most phenotype groups (_e.g._, frequently >10- to 100-fold more controls than cases), we implemented an empirical continuity correction as proposed by [Sweeting _et al._, _Stat. Med._, 2004.](https://onlinelibrary.wiley.com/doi/10.1002/sim.1761)  

Each phenotype & CNV type were meta-analyzed separately for a total of two meta-analyses per phenotype.  

#### Determining exome-wide significance threshold  

We next controlled false discovery rate (FDR) across all gene meta-analyses using a permutation-based approach.  

Estimating the number of independent tests performed across all genes, phenotypes, and CNV types is difficult due to numerous necessary assumptions, such as the independence of samples between phenotypes, or the local correlation structure of CNV counts between neighboring genes, among others.    

Instead, we targeted an "exome-wide" significance threshold of P ≤ 2.90x10<sup>-6</sup>, which corresonds to a Bonferroni correction if applied to the number of individual genes tested in our analysis (N=17,263).  

We empirically determined the primary P-value threshold for each CNV type and phenotype that matched our desired exome-wide FDR as follows:  

1. permute phenotype labels for all CNVs while matching on size (split by quantile) and CNV type (DEL/DUP);  
2. rerun all association tests (described above), including meta-analysis, for each phenotype & CNV combination;   
3. compute the fraction of significant genes for a broad range of primary P-value thresholds (_e.g._, -log<sub>10</sub>(_P_) ~ [0, 20] ); and
4. report the least significant primary P-value threshold that results in an empirical FDR ≤ 2.90x10<sup>-6</sup>.  

Steps 1-3 were repeated 50 times for each CNV type and phenotype.  

Following permutation, we computed the mean primary P-value threshold per phenotype and CNV type. We used this value to assess exome-wide significance of primary P-values for gene discovery.    

#### Output files  

[As described above for Step 2](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/genes#output-files), we generated the same combination of plots and statistics files for the meta-analyses results of each phenotype group.  

Unlike the plots from Step 2, the dashed lines from Step 3 correspond to empirically-derived exome-wide significance thresholds [as determined via permutation](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/genes#determining-exome-wide-significance-threshold).  

These files are stored in the same location as the per-metacohort analysis results.  

### 4. Fine-map causal genes  

### Clustering significant genes into "gene blocks"  

We clustered all significant genes per phenotype into "significant gene blocks" by merging all genes within ±500kb that were significantly associated with the same phenotype.  

In practice, we considered a gene to be exome-wide significant if its primary P-value exceeded the exome-wide significance threshold for that phenotype and CNV type, and it satisifed at least one of the following two criteria:
1. Secondary P-value (as [described above](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/genes/#3-combine-association-statistics-across-metacohorts)) was also nominally significant (P < 0.05); and/or
2. At least two metacohorts were nominally significant (P < 0.05) per Fisher's exact test.  

Absent a true replication sample, these _post hoc_ filters were applied to protect against Winner's Curse.  

### Functional fine-mapping

Next, we aimed to define the minimal set of genes per block most likely to be causal for each phenotypic association.  

To accomplish this, we adapted several GWAS fine-mapping algorithms, as follows:

1. We transformed all association statistics for all genes per block (including non-significant genes) into approximate Bayes factors (ABFs) following the procedure specified by [Wakefield _et al._, _Genet. Epi._, 2009](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20359) and [Wakefield _et al._, _Am. J. Hum. Genet._, 2007](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950810/). For these calculations, we used an empirical Bayes approach to estimate null variance based on the observed variances for all exome-wide significant genes pooled across all phenotypes.  

2. We calculated posterior inclusion probabilities (PIPs) for each gene per block based on the ABFs from [1] while assuming a flat prior (_i.e._, each gene per block is equally likely to be causal).  

3. We developed an adaptation of E-M algorithms described in [Kichaev _et al._, _PLOS Genet._, 2014](https://doi.org/10.1371/journal.pgen.1004722) and [Wen _et al._, _PLOS Genet._, 2017](https://doi.org/10.1371/journal.pgen.1006646) to iteratively update gene priors for each block based on functional enrichments in a logistic regression framework. For this purpose, we used all gene features [described elsewhere in this repository](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene#gene-features).

4. Finally, we computed a 95% credible set of genes per block based on the cumulative sum of PIPs for genes ranked by causal likelihood.  

For downstream analyses, we considered any gene with PIP ≥ 0.1 to be a "candidate" causal gene, and any gene with PIP ≥ 0.9 to be a "strong candidate" causal gene.  

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  

