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



### 4. Prediction of dosage sensitivity scores


---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  

