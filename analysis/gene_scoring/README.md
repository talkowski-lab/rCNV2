# Dosage Sensitivity Scoring  

We developed a model to predict haploinsufficiency and triplosensitivity for each protein-coding gene across all 22 autosomes. This process is described below.  

The code to reproduce these analyses is contained in `gene_scoring_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `gene_scoring_analysis.wdl`.  

---  

## Dosage sensitivity scoring procedure

We trained & applied a model to predict dosage sensitivity scores for all [autosomal, canonical, protein-coding genes](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene/).  

The steps of this model are described below:  

### 0. Curating likely dosage sensitive gene sets for model training  

To train our model, we first curated sets of likely dosage sensitive and insensitive genes based on existing evidence, as follows:

#### Curated dosage sensitive training gene sets  

We defined 291 likely haploinsufficient genes as those which met at least two of the three following criteria:  
1. LoF-constrained (N=3,036 genes); that is, either (i) pLI ≥ 0.9 and/or (ii) in the top LOEUF sextile in gnomAD v2.1; 
2. Known dominant disease gene via LoF mechanism (N=313); that is, either (i) ClinGen high-confidence dominant haploinsufficient genes (N=208) or (ii) DDG2P high-confidence dominant LoF genes (N=300);
3. Genes intolerant of low-expressing samples in GTEx (N=617); specifically, genes that satisfied all of the following: (i) tight expression distributions (mean expression coefficient of variation ≤ 0.3 across tissues), (ii) expressed (mean ≥ 5 TPM & first quartile ≥ 0) in at least 3 tissues, (iii) low rate of low expression outliers (mean < 0.1% samples per tissue are low expression outliers), and (iv) no samples with TPM < 1 in any tissue where gene is expressed.  

We also defined 195 likely triplosensitive genes as those which met at least two of the three following criteria:  
1. Missense-constrained genes (N=3,019 genes); that is, either (i) missense Z ≥ 3 and/or (ii) in the top MOEUF sextile in gnomAD v2.1;
2. Known dominant disease gene via triplosensitive, gain-of-function, or other (non-LoF) mechanisms (N=167); that is, either (i) ClinGen dominant triplosensitive genes at any confidence (N=15) and/or (ii) DDG2P high-confidence dominant gain-of-function or “other” mechanism genes (N=167); 
3. Genes intolerant of overexpression in GTEx (N=279); specifically, genes that satisfied all of the following: (i) tight expression distributions (mean expression coefficient of variation ≤ 0.3 across tissues), (ii) expressed (mean ≥ 5 TPM & first quartile ≥ 0) in at least 3 tissues, and (iii) low rate of high expression outliers (mean < 1% samples per tissue are high expression outliers).  

#### Curated dosage insensitive training gene sets  

We defined a set of 3,580 likely haplosufficient genes as those which met at least two of the three following criteria:
1. Mutation-tolerant (N=2,007 genes);
2. No known disease associations (N=14,659 genes);
3. Genes tolerant of low-expressing samples in GTEx (N=2,682); specifically, genes that satisfied all of the following: (i) variable expression levels (mean expression coefficient of variation ≥ 0.5 across tissues), (ii) expressed (mean ≥ 5 TPM & first quartile ≥ 0) in at least 3 tissues, and (iii) high rate of expression outliers (mean ≥ 0.5% samples per tissue are low expression outliers)  

We also defined a set of 3,065 likely triploinsensitive genes as those which met at least three of the four following criteria:
1. Mutation-tolerant (N=2,007 genes);
2. No known disease associations (N=14,659 genes);
3. Genes tolerant of overexpression in GTEx (N=1,924); specifically, genes that satisfied all of the following: (i) variable expression levels (mean expression coefficient of variation ≥ 0.5 across tissues), (ii) expressed (mean ≥ 5 TPM & first quartile ≥ 0) in at least 3 tissues, and (iii) high rate of expression outliers (mean ≥ 5% samples per tissue are low expression outliers)

For the above definitions:  
`Mutationally tolerant` genes refer to those meeting **both** of the following criteria in [gnomAD v2.1 (Karczewski _et al._, _bioRxiv_, 2019)](https://doi.org/10.1101/531210):  
1. pLI ≤ 0.01; **and**  
2. in the last third of LOEUF scores; **and** 
3. missense Z-score ≤ 0; **and**  
4. missense observed/expected upper fraction ≥ 1; **and**  
5. synonymous Z-score ~ (-3, 3)

`No known disease association` refers to genes not present in any of the following lists:  
1. Online Mendelian Inheritance in Man (OMIM) database; or
2. DECIPHER/DDG2P [DECIPHER/DDG2P (Wright _et al._, _Lancet_, 2015)](https://www.ncbi.nlm.nih.gov/pubmed/25529582); or
3. ClinGen dosage sensitivity map [ClinGen (Strande _et al._, _Am. J. Hum. Genet._, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/28552198).  

### 1. Burden test of rare CNVs per gene between cases and controls  

We counted CNVs per gene between all [cases and controls](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/) for each metacohort, and then [meta-analyzed](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/genes/#3-combine-association-statistics-across-metacohorts) those CNV counts to derive an odds ratio for each gene.  

This process was executed identically to gene-level disease association analyses, described in detail [here](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/genes/#gene-based-burden-test-procedure), with the following exceptions:  
1. Deletion and duplication CDS overlap increased to 80% to enrich for true loss-of-function and copy-gain effects, while still allowing for (i) CNV breakpoint imprecision and/or (ii) exon/gene body misannotation;  
2. CNVs restricted to ≤24 genes to better isolate gene-specific effect sizes;  
3. Phenotypes restricted to a subset of nine with nominally significant (P<0.05) enrichments of whole-gene deletions and duplications of curated dosage sensitive vs. insensitive genes;  
4. Genes restricted to those with at least ≥10 qualifying CNVs each (_note: only for model training_)

### 2a. Empirical Bayes estimation of prior effect sizes  

We next empirically estimated the average effect size expected for true dosage-sensitive genes and true dosage-insensitive genes based on the rCNV data in this study.  

We computed the median odds ratio and variance across these gold-standard genes (described below), and used those values as priors for subsequent steps.

#### Gene blacklist for priors & training  

To protect against certain regions of the genome with recurrent CNVs contributing noise to our effect size estimates for individual genes, we did not consider a subset of genes for any step of prior estimation or training.  

Specifically, we blacklisted genes which:  
1. Overlapped 114 known genomic disorder loci curated from the literature (described [here](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/other#genomic-disorders)); or
2. Were within 1Mb of any centromere or telomere; or
3. Had fewer than 10 CNVs observed (cases + controls across all cohorts).  

Merging our original blacklist with this extra criteria resulted in excluding:
* 14,400 genes while training the deletion model  
* 13,634 genes while training the duplication model  

However, when training the dosage sensitivity classifier, we re-introduced any genes from our curated training gene sets, even if they had insufficient CNV evidence, were near a centromere or telomere, or overlapped a known GD locus.  

### 2b. Estimation of prior fraction of dosage-sensitive genes  

We used Bayesian model averaging to estimate the fraction of all genes in the genome expected to be truly dosage sensitive as an unweighted combination of the following fractions:

* 16.3% of all genes are pLoF-constrained in gnomAD v2.1
* 7.5% of gold-standard genes (curated here) were haploinsufficient
* 20.5% of all genes have a disease association per OMIM
* 3.9% of all genes have some evidence for dominant haploinsufficiency per ClinGen and/or DDG2P (when including unknown/unconfirmed mechanism)

We averaged these values (12.0%), which was used as the prior probability of any gene being haploinsufficient.

Given that less is known (in general) about triplosensitivity, we used the same prior probability for both deletion and duplication models.  

### 3. Bayes factor calculation for each gene  

We next computed the likelihood ratio that each gene was CNV-tolerant versus being dosage sensitive.  

Given a gene with an estimated log odds ratio `T_hat`  and standard error `V_hat`, we compute the Bayes Factor (`BF`) as:  

```
H0 : T_hat ≤ T_null
H1 : T_hat ≥ T_alt

BF = ( 1 - N(T_hat - T_null, 1) ) / N(T_hat - T_alt, V_hat)
```

Where:  
*  `T_null` was the estimated null effect size computed from [gold-standard dosage-insensitive genes](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/gene_scoring#gold-standard-dosage-insensitive-genes), and  
*  `T_alt` was the estimated alternative effect size computed from [gold-standard haploinsufficient genes](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/gene_scoring#gold-standard-dosage-sensitive-genes).  

Per equation (2) from [Wakefield, _Am. J. Hum. Genet._, 2007](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950810/), we can further define a Bayesian false discovery probability (BFDP) for each gene as:  

```
BFDP = ( PO x BF ) / ( 1 + ( PO x BF ) )
```

Where `PO` is the prior odds of the null; _i.e._, the odds of any gene _not_ being dosage sensitive, which was estimated [as described above](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/gene_scoring#2b-estimation-of-prior-fraction-of-dosage-sensitive-genes).  

### 4. Prediction of dosage sensitivity scores

Using the BFDP values calculated for each gene [as described above](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/gene_scoring#3-bayes-factor-calculation-for-each-gene), we next trained and applied a battery of regression and machine learning models to predict the probability of a gene being haploinsufficient or triplosensitive based on [gene-level functional metadata](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene#gene-features) alone.  

Specifically, we evaluated each of the seven following models using the same strategy (outlined below):  
*  Logistic (logit) generalized linear model  
*  Linear support vector machine (SVM)  
*  Random forest  
*  Linear discriminant analysis (LDA)  
*  Naive Bayes  
*  Logistic stochastic gradient descent (SGD)  
*  Neural network (MLP) with logistic activation  

For this analysis, we additionally blacklisted any gene with fewer than three observed CNVs (counted separately for deletions and duplications). Five was determined as the minimum number of CNVs required to obtain a nominally significant (_i.e._, informative) enrichment in cases from a Fisher's exact test.  

For each model, we first split each autosome into its respective p and q arms, then partitioned all chromsome arms into 11 groups while balancing the total number of non-blacklisted genes per group, as follows:  

**Chromosome arm pairings for deletion model** 

| Chromosome arms | Genes used for training | Genes excluded from training |  
| :--- | ---: | ---: |  
| 2q, 16q, 12p | 459 | 876 |  
| 1p, 6p, 13q | 459 | 1,243 |  
| 4q, 7q, 9p | 460 | 734 |  
| 1q, 16p, 20p | 488 | 939 |  
| 17q, 2p, 20q, 18p | 490 | 1,139 |  
| 19q, 15q, 8p, 5p | 493 | 1,144 |  
| 11q, 6q, 22q | 497 | 1,047 |  
| 5q, 11p, 7p, 10p | 506 | 987 |  
| 3p, 14q, 9q, 4p | 509 | 1,291 |  
| 3q, 19p, 10q, 17p | 518 | 1,448 |  
| 12q, 8q, 21q, 18q | 520 | 1,016 |  

**Chromosome arm pairings for duplication model**  

| Chromosome arms | Genes used for training | Genes excluded from training |  
| :--- | ---: | ---: |  
| 14q, 9q, 4p | 459 | 857 |  
| 2q, 7q, 20p | 468 | 923 |  
| 12q, 11p, 3p | 472 | 1,122 |  
| 3q, 2p, 6p | 481 | 939 |  
| 17q, 15q, 10q | 481 | 1,345 |  
| 19q, 7p, 13q, 18p | 500 | 882 |  
| 5q, 16p, 21q, 10p | 518 | 907 |  
| 1q, 16q, 8p, 5p | 519 | 1,062 |  
| 4q, 19p, 12p, 22q | 520 | 1,199 |  
| 11q, 8q, 9p, 18q | 520 | 1,036 |  
| 1p, 6q, 17p, 20q | 537 | 1,516 |  

To predict gene-level scores for each group of chromosome arms, we used the other 10 groups to train & cross-validate a classification model.  

For each of 10 cross-validation iterations, we randomly allocated 90% of genes for training and held out 10% for testing. All genes considered at this stage were only sampled from the 10 chromosome arm groups _not_ including the held-out 11<sup>th</sup> group to be used for prediction.    

This process can be summarized as follows:

1. Fit classification model to the 90% of genes allocated for training;  
2. Compute root mean-squared error (RMSE) of predicted and actual BFDPs for the 10% of genes held out for testing; 
3. Repeat steps [1-2] once for each of the 10 folds; 
4. Select the most predictive CV fit based on lowest RMSE; and 
5. Apply the most predictive model from [4] to predict BFDPs for the genes from the 11<sup>th</sup> group of (completely held-out) prediction chromosome arms.  

After computing predicted BFDPs for all genes using this strategy described above, we subsequently standard-normalized the predicted BFDPs, and defined final gene scores as `1 – normal_cdf( normalized BFDP )`.  

Lastly, we compared each of the seven models (logit, SVM, random forest, LDA, naive Bayes, SGD, neural net) based on AUROC and AUPRC, and computed an ensemble probability per gene as the weighted average of all seven models, which was weighted by the accuracy of each model at their corresponding ROC-optimal cutoff.  

In practice we dubbed the final scores from our ensemble model as:  
*  **pHI**: probability of haploinsufficiency; and
*  **pTS**: probability of triplosensitivity.  

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  

