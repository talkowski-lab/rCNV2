# Dosage Sensitivity Scoring  

We developed a model to predict haploinsufficiency and triplosensitivity for each protein-coding gene across all 22 autosomes. This process is described below.  

The code to reproduce these analyses is contained in `gene_scoring_analysis.sh`.  

---  

## Dosage sensitivity scoring procedure

We trained & applied a model to predict dosage sensitivity scores for all [autosomal, canonical, protein-coding genes](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene/).  

The steps of this model are described below:  

### 1. Burden test of rare CNVs per gene between cases and controls  

We counted CNVs per gene between all [cases and controls](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/) for each metacohort, and then [meta-analyzed](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/genes/#3-combine-association-statistics-across-metacohorts) those CNV counts to derive an odds ratio for each gene.  

This process was executed identically to gene-level disease association analyses, described in detail [here](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/genes/#gene-based-burden-test-procedure).  

### 2a. Empirical Bayes estimation of prior effect sizes  

We next empirically estimated the average effect size and variance expected for true dosage-sensitive genes and true dosage-insensitive genes based on the rCNV data in this study.  

We created these gene sets as follows:  

#### Gold-standard dosage-sensitive genes  

We defined a set of 95 gold-standard haploinsufficient genes, which were **both**:  
1. mutationally constrained, **and**
2. known causes of dominant disease via a loss-of-function/haploinsufficient mechanism.  

We defined `mutationally constrained` genes as those meeting **either** of the following two criteria in [gnomAD v2.1 (Karczewski _et al._, _bioRxiv_, 2019)](https://doi.org/10.1101/531210):  
1. probability of loss-of-function intolerance (pLI) ≥ 0.9 ; **or**  
2. in the first sextile of loss-of-function observed/expected upper fraction (LOEUF) scores.  

And we defined `known causes of dominant disease via a loss-of-function/haploinsufficient mechanism` genes as those meeting **both** of the following criteria: 
1. dominant haploinsufficient disease genes scored at high confidence by [ClinGen (Strande _et al._, _Am. J. Hum. Genet._, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/28552198); **and**  
2. dominant developmental disorder genes with a loss-of-function mechanism scored as "confirmed" in [DECIPHER/DDG2P (Wright _et al._, _Lancet_, 2015)](https://www.ncbi.nlm.nih.gov/pubmed/25529582).  

#### Gold-standard dosage-insensitive genes  

We also defined a set of 1,728 gold-standard haplo*sufficient* genes, which were **both**:
1. mutationally tolerant, **and**
2. had no known disease associations. 

We defined `mutationally tolerant` genes as those meeting **both** of the following criteria in [gnomAD v2.1 (Karczewski _et al._, _bioRxiv_, 2019)](https://doi.org/10.1101/531210):  
1. pLI ≤ 0.01; **and**  
2. in the last third of LOEUF scores; **and** 
3. missense Z-score ≤ 0; **and**  
4. missense observed/expected upper fraction ≥ 1; **and**  
5. synonymous Z-score ~ (-3, 3)

We defined `no known disease association` genes as those present in none of the following lists:  
1. Online Mendelian Inheritance in Man (OMIM) database; or
2. DECIPHER/DDG2P [DECIPHER/DDG2P (Wright _et al._, _Lancet_, 2015)](https://www.ncbi.nlm.nih.gov/pubmed/25529582); or
3. ClinGen dosage sensitivity map [ClinGen (Strande _et al._, _Am. J. Hum. Genet._, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/28552198).  

We computed the median odds ratio and variance across these gold-standard genes, and used those values as priors for subsequent steps.

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

BF = ( 1 - N(T_hat - T_null, W) ) / N(T_hat - T_alt)
```

Where:  
*  `T_null` was the estimated null effect size computed from [gold-standard dosage-insensitive genes](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/gene_scoring#gold-standard-dosage-insensitive-genes),  
*  `T_alt` was the estimated alternative effect size computed from [gold-standard haploinsufficient genes](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/gene_scoring#gold-standard-dosage-sensitive-genes), and  
*  `W` was the null standard error computed from gold-standard haploinsufficient genes.  

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

For each model, we first divided all 22 autosomes into 11 pairs while balancing the total number of genes per chromosome pair, as follows:  

| Chromosome #1 | Chromosome #2 | Genes |  
| ---: | ---: | ---: |  
| 1 | 21 | 2,036 |  
| 19 | 18 | 1,586 |  
| 2 | 13 | 1,473 |  
| 11 | 22 | 1,520 |  
| 17 | 20 | 1,579 |  
| 3 | 15 | 1,566 |  
| 12 | 14 | 1,593 |  
| 5 | 8 | 1,464 |  
| 6 | 10 | 1,507 |  
| 7 | 9 | 1,476 |  
| 16 | 4 | 1,463 |  

To predict gene-level scores for each pair of chromosomes, we used the other 10 pairs to train & cross-validate a classification model.  

For each of 10 cross-validation iterations, we used 9/10 chromosome pairs for training and held out one pair for testing.  

This process can be summarized as follows:

1. Fit classification model for all genes from the 9/10 training pairs;  
2. Compute root mean-squared error (RMSE) of predicted and actual BFDPs for all genes on the 10<sup>th</sup> chromosome pair (held out from training); 
3. Repeat steps [1-2] once for each of the 10 training pairs as the held-out test pair; 
4. Select the most predictive CV fit based on lowest RMSE; and 
5. Apply the most predictive model from [4] to predict BFDPs for the genes from the 11<sup>th</sup> pair of (completely held-out) prediction chromosomes.  

After computing predicted BFDPs for all genes using this strategy described above, we subsequently standard-normalized the predicted BFDPs, and defined final gene scores as `1 – normal_cdf( normalized BFDP )`.  

Lastly, we compared each of the six models (logit, SVM, random forest, LDA, naive Bayes, and SGD) based on RMSE, AUROC, and AUPRC, and selected the overall best-performing model across deletions and duplications for our final analysis.  

In practice we dubbed the final scores from our best-performing model as:  
*  **pHI**: probability of haploinsufficiency; and
*  **pTS**: probability of triplosensitivity.  

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  

