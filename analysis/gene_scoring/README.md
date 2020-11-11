# Dosage Sensitivity Scoring  

We developed a model to predict haploinsufficiency and triplosensitivity for each protein-coding gene across all 22 autosomes. This process is described below.  

The code to reproduce these analyses is contained in `gene_scoring_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `gene_scoring_analysis.wdl`.  

---  

## Dosage sensitivity scoring procedure

We trained & applied a model to predict dosage sensitivity scores for all [autosomal, canonical, protein-coding genes](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene/).  

The steps of this model are described below:  

### 1. Burden test of rare CNVs per gene between cases and controls  

We counted CNVs per gene between all [cases and controls](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/) for each metacohort, and then [meta-analyzed](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/genes/#3-combine-association-statistics-across-metacohorts) those CNV counts to derive an odds ratio for each gene.  

This process was executed identically to gene-level disease association analyses, described in detail [here](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/genes/#gene-based-burden-test-procedure), with the following exceptions:  
1. CNVs restricted to ≤5 genes to better isolate gene-specific effect sizes
2. Deletion and duplication CDS overlap increased to 80% to enrich for true loss-of-function and copy-gain effects, while still allowing for (i) CNV breakpoint imprecision and/or (ii) exon/gene body misannotation  

### 2a. Empirical Bayes estimation of prior effect sizes  

We next empirically estimated the average effect size expected for true dosage-sensitive genes and true dosage-insensitive genes based on the rCNV data in this study.  

We computed the median odds ratio and variance across these gold-standard genes (described below), and used those values as priors for subsequent steps.

We created these gene sets as follows:  

#### Gold-standard dosage-sensitive genes  

We defined 294 gold-standard haploinsufficient genes as those which met at least two of the three following criteria:  
1. LoF-constrained (N=3,036 genes); that is, either (i) pLI ≥ 0.9 and/or (ii) in the top LOEUF sextile in gnomAD v2.1; 
2. Known dominant disease gene via LoF mechanism (N=313); that is, either (i) ClinGen high-confidence dominant haploinsufficient genes (N=208) or (ii) DDG2P high-confidence dominant LoF genes (N=300);
3. Genes intolerant of low-expressing samples in GTEx (N=617); specifically, genes that satisfied all of the following: (i) tight expression distributions (mean expression coefficient of variation ≤ 0.3 across tissues), (ii) expressed (mean ≥ 5 TPM & first quartile ≥ 0) in at least 3 tissues, (iii) low rate of low expression outliers (mean < 0.1% samples per tissue are low expression outliers), and (iv) no samples with TPM < 1 in any tissue where gene is expressed.  

We also defined 151 gold-standard triplosensitive genes as those which met at least two of the three following criteria:  
1. Missense-constrained genes (N=3,019 genes); that is, either (i) missense Z ≥ 3 and/or (ii) in the top MOEUF sextile in gnomAD v2.1;
2. Known dominant disease gene via triplosensitive, gain-of-function, or other (non-LoF) mechanisms (N=167); that is, either (i) ClinGen dominant triplosensitive genes at any confidence (N=15) and/or (ii) DDG2P high-confidence dominant gain-of-function or “other” mechanism genes (N=158); 
3. Genes intolerant of overexpression in GTEx (N=279); specifically, genes that satisfied all of the following: (i) tight expression distributions (mean expression coefficient of variation ≤ 0.3 across tissues), (ii) expressed (mean ≥ 5 TPM & first quartile ≥ 0) in at least 3 tissues, and (iii) low rate of high expression outliers (mean < 1% samples per tissue are high expression outliers).  


#### Gold-standard dosage-insensitive genes  

We defined a set of 3,602 gold-standard haplosufficient genes as those which met at least two of the three following criteria:
1. Mutation-tolerant (N=2,007 genes);
2. No known disease associations (N=14,659 genes);
3. Genes tolerant of low-expressing samples in GTEx (N=2,682); specifically, genes that satisfied all of the following: (i) variable expression levels (mean expression coefficient of variation ≥ 0.5 across tissues), (ii) expressed (mean ≥ 5 TPM & first quartile ≥ 0) in at least 3 tissues, and (iii) high rate of expression outliers (mean ≥ 0.5% samples per tissue are low expression outliers)  

We also defined a set of 3,806 bronze-standard triploinsensitive genes as those which met at least three of the four following criteria:
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

#### Gene blacklist for priors & training  

To protect against certain regions of the genome with recurrent CNVs contributing noise to our effect size estimates for individual genes, we did not consider a subset of genes for any step of prior estimation or training.  

Specifically, we blacklisted genes which:  
1. Overlapped 114 known genomic disorder loci curated from the literature (described [here](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/other#genomic-disorders)); or
2. Were within 1Mb of any centromere or telomere.  

This resulted in blacklisting a total of 1,848 genes for all prior estimation and model training (described below).  

### 2b. Estimation of prior fraction of dosage-sensitive genes  

To estimate the fraction of all protein-coding genes expected to be truly dosage sensitive, we averaged the following four fractions:  

*  16.3% of genes are pLoF-constrained in gnomAD v2.1  
*  5.2% of all gold-standard genes (curated above) were haploinsufficient  
*  20.5% of all genes have a disease association per OMIM  
*  3.9% of all genes have some evidence for dominant haploinsufficiency per ClinGen and/or DDG2P (when including unknown/unconfirmed mechanism)  

For modeling, we assumed each gene had equal prior likelihood of being haploinsufficient (0.115).  

### 3. Bayes factor calculation for each gene  

We next computed the likelihood ratio that each gene was CNV-tolerant versus being dosage sensitive.  

Given a gene with an estimated log odds ratio `T_hat`  and standard error `V_hat`, we compute the Bayes Factor (`BF`) as:  

```
H0 : T_hat ≤ T_null
H1 : T_hat ≥ T_alt

BF = ( 1 - N(T_hat - T_null, 1) ) / N(T_hat - T_alt)
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

Specifically, we evaluated each of the following models using the same strategy (outlined below):  
*  Logistic (logit) generalized linear model  
*  Linear support vector machine (SVM)  
*  Random forest  
*  Linear discriminant analysis (LDA)  
*  Naive Bayes  
*  Logistic stochastic gradient descent (SGD)  
*  Neural network (MLP) with logistic activation  

For this analysis, we additionally blacklisted any gene with fewer than three observed CNVs (counted separately for deletions and duplications). Five was determined as the minimum number of CNVs required to obtain a nominally significant (_i.e._, informative) enrichment in cases from a Fisher's exact test.  

We also excluded genes with less informative intermediate BFDPs ( ~ [0.2, 0.8]) from training.  

For each model, we first split each autosome into its respective p and q arms, then partitioned all chromsome arms into 11 groups while balancing the total number of non-blacklisted genes per group, as follows:  

**Chromosome arm pairings for deletion model** 

| Chromosome arms | Genes used in training | Blacklisted genes |  
| ---: | ---: | ---: |  
| 2q, 3p, 7p | 1,400 | 58 |  
| 19q, 15q, 12p | 1,404 | 141 |  
| 12q, 6p, 22q | 1,416 | 105 |  
| 11q, 6q, 13q | 1,425 | 38 |  
| 5q, 7q, 17p, 18p | 1,452 | 123 |  
| 3q, 9q, 11p, 20p | 1,505 | 107 |  
| 1p, 16q, 4p | 1,509 | 50 |  
| 19p, 4q, 20q, 9p | 1,514 | 89 |  
| 1q, 8q, 18q, 5p | 1,519 | 96 |  
| 17q, 2p, 8p, 21q | 1,519 | 105 |  
| 14q, 10q, 16p, 10p | 1,533 | 155 |  

**Chromosome arm pairings for duplication model**  

| Chromosome arms | Genes used in training | Blacklisted genes |  
| ---: | ---: | ---: |  
| 2q, 3p, 7p | 1,400 | 58 |  
| 19q, 15q, 12p | 1,404 | 141 |  
| 12q, 6p, 22q | 1,416 | 105 |  
| 11q, 6q, 13q | 1,425 | 38 |  
| 5q, 7q, 17p, 18p | 1,452 | 123 |  
| 3q, 9q, 11p, 20p | 1,505 | 107 |  
| 1p, 16q, 4p | 1,509 | 50 |  
| 19p, 4q, 20q, 9p | 1,514 | 89 |  
| 1q, 8q, 18q, 5p | 1,519 | 96 |  
| 17q, 2p, 8p, 21q | 1,519 | 105 |  
| 14q, 10q, 16p, 10p | 1,533 | 155 |  

To predict gene-level scores for each group of chromosome arms, we used the other 10 groups to train & cross-validate a classification model.  

For each of 10 cross-validation iterations, we randomly allocated 90% of genes for training and held out 10% for testing. All genes considered at this stage were only sampled from the 10 chromosome arm groups _not_ including the held-out 11<sup>th</sup> group to be used for prediction.    

This process can be summarized as follows:

1. Fit classification model to the 90% of genes allocated for training;  
2. Compute root mean-squared error (RMSE) of predicted and actual BFDPs for the 10% of genes held out for testing; 
3. Repeat steps [1-2] once for each of the 10 folds; 
4. Select the most predictive CV fit based on lowest RMSE; and 
5. Apply the most predictive model from [4] to predict BFDPs for the genes from the 11<sup>th</sup> group of (completely held-out) prediction chromosome arms.  

After computing predicted BFDPs for all genes using this strategy described above, we subsequently standard-normalized the predicted BFDPs, and defined final gene scores as `1 – normal_cdf( normalized BFDP )`.  

Lastly, we compared each of the six models (logit, SVM, random forest, LDA, naive Bayes, SGD, neural net) based on AUROC and AUPRC, and selected the overall best-performing model based on the highest harmonic mean AUC across deletions and duplications for our final analysis.  

In practice we dubbed the final scores from our best-performing model as:  
*  **pHI**: probability of haploinsufficiency; and
*  **pTS**: probability of triplosensitivity.  

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  

