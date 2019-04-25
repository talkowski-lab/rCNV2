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

We harmonized phenotype data across [all cohorts included in the overall CNV dataset](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV/).  

Given the variability of phenotypic detail available between different cohorts, we applied a standardized hierarchical phenotype consolidation scheme to each cohort.  

### Phenotype hierarchy  

The hierarchy of phenotypes considered is listed below. HPO numbers for each term are listed in brackets.  

```
Germline disease
|
|--Neurological defect
|   |  
|   |  
|   |  
|
|--Non-neurological defect
```  

### Phenotype assignment per cohort  

We processed each cohort as follows:

#### CHOP  
 * 

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



