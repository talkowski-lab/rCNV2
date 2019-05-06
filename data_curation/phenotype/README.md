# Phenotype Data Curation

We collapsed cohorts into sets of case:control cohort pairs matched on phenotype and technical CNV properties. This process is described below.  

### Phenotype abbreviations  

Several abbreviations and acronyms are used below, as follows:  

| Abbreviation | Definition |
| --- | :--- |
| ASD | Autism Spectrum Disorders |
| CHD | Congenital Heart Disease |  
| CTRL | Control (Unaffected Sample) |
| DD | Developmental Disorders |
| NPD | Other neuropsychiatric disorders |
| SCZ | Schizophrenia |

### Phenotype curation  

Given the variability of phenotypic detail available between different cohorts, we applied a standardized hierarchical phenotype consolidation scheme uniformly across all samples.  

We used the [Human Phenotype Ontology (HPO)](https://hpo.jax.org/app/) to harmonize phenotype data across all [cohorts included in the overall CNV dataset](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV/).  

For each sample, we recursively matched all phenotype information available against HPO keywords, and assigned the corresponding HPO terms to each sample when a matche was found.  

The code to perform this HPO conversion is contained in `convert_phenotypes_to_hpo.sh`.  

### Phenotype hierarchy  

After assigning HPO terms to each sample, we next subgrouped samples based on concordant HPO terms into a reduced hierarchy of "phenotype groups" considered in this study.  

After tabulating the number of samples matching each HPO code, we retained HPO codes with at least 1,000 samples. In some situations, we collapsed redundant HPO codes up the hierarchy.  

This process yielded a hierarchical phenotype classification system with XX distinct HPO-based sample groupings, each which has ≥ 1,000 samples.  

The hierarchy of phenotypes considered is listed below. HPO numbers for each term are listed in brackets.  

```
Updated hierarchy image will go here when ready
```  

--- 

## Phenotype Case-Control Cohort Pairings  

For each phenotype, we broke down available case and control cohorts into pairs for burden testing. 

These case-control pairings are described below per phenotype.  

### ASD: Autism spectrum disorders  

**Case phenotypes:** autism spectrum disorders  

| Group ID | Case Source(s) | Number of Cases | Control Source(s) | Number of Controls | Stage |  
| :--- | :--- | ---: | :--- | ---: | :--- |  
| ASD_1 | Coe, BCH, GDX | TBD | Coe, Cooper | 19,585 | Discovery |    
| ASD_2 | SSC | TBD | TCGA, PGC, TSAICG | 33,040 | Discovery |  
| ASD_3 | CHOP | TBD | CHOP | 28,070 | Replication |  
| **Total** | - | **TBD** | - | **TBD** | - |  

### DD: Developmental disorders  

**Case phenotypes:** developmental disorders, developmental delay, or intellectual disability  

| Group ID | Case Source(s) | Number of Cases | Control Source(s) | Number of Controls | Stage |  
| :--- | :--- | ---: | :--- | ---: | :--- |  
| DD_1 | Coe, BCH, GDX | TBD | Coe, Cooper | 19,585 | Discovery |    
| DD_2 | SSC | TBD | TCGA, PGC, TSAICG | 33,040 | Discovery |  
| DD_3 | CHOP | TBD | CHOP | 28,070 | Replication |  
| **Total** | - | **TBD** | - | **TBD** | - |  


### CHD: Congenital heart disease  

**Case phenotypes:** congenital heart disease  

| Group ID | Case Source(s) | Number of Cases | Control Source(s) | Number of Controls | Stage |  
| :--- | :--- | ---: | :--- | ---: | :--- |  
| ASD_1 | Coe, BCH, GDX | TBD | Coe, Cooper | 19,585 | Discovery |    
| ASD_2 | SSC | TBD | TCGA, PGC, TSAICG | 33,040 | Discovery |  
| ASD_3 | CHOP | TBD | CHOP | 28,070 | Replication |  
| **Total** | - | **TBD** | - | **TBD** | - |  



### SCZ: Schizophrenia  

**Case phenotypes:** schizophrenia  


### NPD: Other neuropsychiatric disorders  

**Case phenotypes:** ADHD, Tourette Syndrome, anxiety, anorexia, depression, mood disorders, bipolar   



