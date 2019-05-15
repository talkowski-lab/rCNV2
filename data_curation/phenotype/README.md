# Phenotype Data Curation

We grouped samples into phenotype-matched subsets using a standardized hierarchical ontology.  

This process is described below.  

### Phenotype curation  

Given the variability of phenotypic detail available between different cohorts, we applied a standardized hierarchical phenotype consolidation scheme uniformly across all samples.  

We used the [Human Phenotype Ontology (HPO)](https://hpo.jax.org/app/) to harmonize phenotype data across all [cohorts included in the overall CNV dataset](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV/).  

For each sample, we recursively matched all phenotype information available against HPO keywords, and assigned the corresponding HPO terms to each sample when a match was found.  

The code to perform this HPO conversion is contained in `convert_phenotypes_to_hpo.sh`.  

### Phenotype hierarchy  

After assigning HPO terms to each sample, we next subgrouped samples based on concordant HPO terms into a reduced hierarchy of "phenotype groups" considered in this study.  

After tabulating the number of samples matching each HPO code, we retained all HPO terms with at least 1,000 samples. 

We next compared shared sample memberships between all pairs of related HPO terms:
 * If fewer than 1,000 samples differed between a pair of HPO terms, we retained the more general (i.e., higher-level) term, and excluded the more specific (i.e., lower-level) term.
 * If both terms were equally high-level, we retained the term with the larger sample size.

This process yielded a hierarchical phenotype classification system with 30 distinct HPO-based sample groups, each which has ≥ 1,000 samples and differs from another phenotype group by at least 1,000 distinct samples.  

The code to perform the phenotype group consolidation process described above is contained in `convert_phenotypes_to_hpo.sh`    

The final set of phenotypes considered is detailed below:  

| HPO Term | Description | Samples | HPO Tier | Parent Terms | Child Terms |  
| :--- | :--- | ---: | ---: | :--- | :--- |  
| HP:0000118 | Phenotypic abnormality | 223,116 | 1 | NA | HP:0007367, HP:0031466, HP:0001627, HP:0012759, HP:0000736, HP:0000752, HP:0012638, HP:0001626, HP:0001249, HP:0000152, HP:0004305, HP:0001510, HP:0012639, HP:0012443, HP:0000924, HP:0001507, HP:0002715, HP:0012758, HP:0002597, HP:0002197, HP:0100022, HP:0100753, HP:0000729, HP:0001250, HP:0025031, HP:0000708, HP:0003011, HP:0000707 |  
| HP:0000707 | Abnormality of the nervous system | 116,756 | 2 | HP:0000118 | HP:0012639, HP:0012443, HP:0001250, HP:0007367, HP:0012758, HP:0001249, HP:0000708, HP:0012759, HP:0000736, HP:0004305, HP:0000752, HP:0012638, HP:0002197, HP:0100022, HP:0100753, HP:0000729, HP:0031466 |  
| HP:0012638 | Abnormality of nervous system physiology | 91,042 | 3 | HP:0000707, HP:0000118 | HP:0000708, HP:0012759, HP:0000736, HP:0004305, HP:0000752, HP:0001250, HP:0012758, HP:0002197, HP:0100022, HP:0001249, HP:0100753, HP:0000729, HP:0031466 |  
| HP:0000708 | Behavioral abnormality | 54,033 | 4 | HP:0000707, HP:0012638, HP:0000118 | HP:0100753, HP:0000729, HP:0000736, HP:0031466 |  
| UNKNOWN | NA | 49,175 | 2 | HP:0000118 | NA |  
| HP:0012759 | Neurodevelopmental abnormality | 36,359 | 4 | HP:0000707, HP:0012638, HP:0000118 | HP:0001249, HP:0012758 |  
| HP:0012758 | Neurodevelopmental delay | 34,615 | 5 | HP:0000707, HP:0012638, HP:0012759, HP:0000118 | NA |  
| HP:0002715 | Abnormality of the immune system | 32,930 | 2 | HP:0000118 | NA |  
| HP:0012639 | Abnormality of nervous system morphology | 26,271 | 3 | HP:0000707, HP:0000118 | HP:0012443, HP:0007367 |  
| HP:0007367 | Atrophy/Degeneration affecting the central nervous system | 24,257 | 5 | HP:0000707, HP:0012639, HP:0000118 | NA |  
| HP:0001626 | Abnormality of the cardiovascular system | 22,367 | 2 | HP:0000118 | HP:0001627, HP:0002597 |  
| HP:0100753 | Schizophrenia | 22,301 | 5 | HP:0000707, HP:0012638, HP:0000708, HP:0000118 | NA |  
| HP:0000729 | Autistic behavior | 16,919 | 5 | HP:0000707, HP:0012638, HP:0000708, HP:0000118 | NA |  
| HP:0100022 | Abnormality of movement | 15,552 | 4 | HP:0000707, HP:0012638, HP:0000118 | HP:0004305, HP:0000752 |  
| HP:0002597 | Abnormality of the vasculature | 15,070 | 3 | HP:0000118, HP:0001626 | NA |  
| HP:0000752 | Hyperactivity | 12,940 | 5 | HP:0000707, HP:0012638, HP:0100022, HP:0000118 | NA |  
| HP:0000736 | Short attention span | 12,907 | 5 | HP:0000707, HP:0012638, HP:0000708, HP:0000118 | NA |  
| HP:0000152 | Abnormality of head or neck | 7,037 | 2 | HP:0000118 | NA |  
| HP:0001627 | Abnormal heart morphology | 6,749 | 4 | HP:0000118, HP:0001626 | NA |  
| HP:0001250 | Seizures | 4,961 | 4 | HP:0000707, HP:0012638, HP:0000118 | HP:0002197 |  
| HP:0001507 | Growth abnormality | 3,254 | 2 | HP:0000118 | HP:0001510 |  
| HP:0004305 | Involuntary movements | 2,538 | 5 | HP:0000707, HP:0012638, HP:0100022, HP:0000118 | NA |  
| HP:0031466 | Impairment in personality functioning | 2,171 | 5 | HP:0000707, HP:0012638, HP:0000708, HP:0000118 | NA |  
| HP:0001249 | Intellectual disability | 2,069 | 5 | HP:0000707, HP:0012638, HP:0012759, HP:0000118 | NA |  
| HP:0000924 | Abnormality of the skeletal system | 1,976 | 2 | HP:0000118 | NA |  
| HP:0012443 | Abnormality of brain morphology | 1,944 | 5 | HP:0000707, HP:0012639, HP:0000118 | NA |  
| HP:0002197 | Generalized-onset seizure | 1,781 | 5 | HP:0000707, HP:0012638, HP:0001250, HP:0000118 | NA |  
| HP:0025031 | Abnormality of the digestive system | 1,586 | 2 | HP:0000118 | NA |  
| HP:0003011 | Abnormality of the musculature | 1,373 | 2 | HP:0000118 | NA |  
| HP:0001510 | Growth delay | 1,030 | 3 | HP:0000118, HP:0001507 | NA | 

### HPO terms per cohort

The number of samples per phenotype group per cohort is outlined in the table, below.  

The code to generate this table is provided in `gather_hpo_per_cohort_table.py`.  

| HPO | description | **Total** | PGC | Cooper | Coe | SSC | UKBB | CHOP | GDX | TSAICG | BCH | TCGA |  
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| HP:0000118 | Phenotypic abnormality | 223116 | 21094 | 0 | 29104 | 2795 | 0 | 154140 | 9958 | 2434 | 3591 | 0 |  
| HP:0000707 | Abnormality of the nervous system | 116756 | 21094 | 0 | 29104 | 2795 | 0 | 54406 | 4009 | 2434 | 2914 | 0 |  
| HP:0012638 | Abnormality of nervous system physiology | 91042 | 21094 | 0 | 29104 | 2795 | 0 | 28948 | 3821 | 2434 | 2846 | 0 |  
| HP:0000708 | Behavioral abnormality | 54033 | 21094 | 0 | 1834 | 2795 | 0 | 26364 | 722 | 0 | 1224 | 0 |  
| UNKNOWN | NA | 49175 | 0 | 0 | 0 | 0 | 0 | 43873 | 5009 | 0 | 293 | 0 |  
| HP:0012759 | Neurodevelopmental abnormality | 36359 | 0 | 0 | 29104 | 814 | 0 | 671 | 3404 | 0 | 2366 | 0 |  
| HP:0012758 | Neurodevelopmental delay | 34615 | 0 | 0 | 29104 | 0 | 0 | 0 | 3302 | 0 | 2209 | 0 |  
| HP:0002715 | Abnormality of the immune system | 32930 | 0 | 0 | 19 | 0 | 0 | 32843 | 34 | 0 | 34 | 0 |  
| HP:0012639 | Abnormality of nervous system morphology | 26271 | 0 | 0 | 352 | 0 | 0 | 24936 | 706 | 0 | 277 | 0 |  
| HP:0007367 | Atrophy/Degeneration affecting the central nervous system | 24257 | 0 | 0 | 4 | 0 | 0 | 24247 | 4 | 0 | 2 | 0 |  
| HEALTHY_CONTROL | Unaffected control sample | 24161 | 20277 | 8329 | 11256 | 0 | 480501 | 28070 | 0 | 4093 | 0 | 8670 |  
| HP:0001626 | Abnormality of the cardiovascular system | 22367 | 0 | 0 | 748 | 0 | 0 | 21003 | 404 | 0 | 212 | 0 |  
| HP:0100753 | Schizophrenia | 22301 | 21094 | 0 | 2 | 0 | 0 | 1200 | 5 | 0 | 0 | 0 |  
| HP:0000729 | Autistic behavior | 16919 | 0 | 0 | 1569 | 2795 | 0 | 11071 | 631 | 0 | 853 | 0 |  
| HP:0100022 | Abnormality of movement | 15552 | 0 | 0 | 317 | 0 | 0 | 12394 | 82 | 2434 | 325 | 0 |  
| HP:0002597 | Abnormality of the vasculature | 15070 | 0 | 0 | 50 | 0 | 0 | 14948 | 56 | 0 | 16 | 0 |  
| HP:0000752 | Hyperactivity | 12940 | 0 | 0 | 281 | 0 | 0 | 12359 | 62 | 0 | 238 | 0 |  
| HP:0000736 | Short attention span | 12907 | 0 | 0 | 264 | 0 | 0 | 12359 | 57 | 0 | 227 | 0 |  
| HP:0000152 | Abnormality of head or neck | 7037 | 0 | 0 | 4535 | 0 | 0 | 0 | 1654 | 0 | 848 | 0 |  
| HP:0001627 | Abnormal heart morphology | 6749 | 0 | 0 | 166 | 0 | 0 | 6055 | 332 | 0 | 196 | 0 |  
| HP:0001250 | Seizures | 4961 | 0 | 0 | 1807 | 318 | 0 | 1762 | 535 | 0 | 539 | 0 |  
| HP:0001507 | Growth abnormality | 3254 | 0 | 0 | 1221 | 0 | 0 | 757 | 1058 | 0 | 218 | 0 |  
| HP:0004305 | Involuntary movements | 2538 | 0 | 0 | 10 | 0 | 0 | 35 | 5 | 2434 | 54 | 0 |  
| HP:0031466 | Impairment in personality functioning | 2171 | 0 | 0 | 18 | 0 | 0 | 1730 | 14 | 0 | 409 | 0 |  
| HP:0001249 | Intellectual disability | 2069 | 0 | 0 | 129 | 814 | 0 | 671 | 155 | 0 | 300 | 0 |  
| HP:0000924 | Abnormality of the skeletal system | 1976 | 0 | 0 | 338 | 0 | 0 | 0 | 1282 | 0 | 356 | 0 |  
| HP:0012443 | Abnormality of brain morphology | 1944 | 0 | 0 | 319 | 0 | 0 | 689 | 680 | 0 | 256 | 0 |  
| HP:0002197 | Generalized-onset seizure | 1781 | 0 | 0 | 1776 | 0 | 0 | 0 | 2 | 0 | 3 | 0 |  
| HP:0025031 | Abnormality of the digestive system | 1586 | 0 | 0 | 91 | 0 | 0 | 1258 | 144 | 0 | 93 | 0 |  
| HP:0003011 | Abnormality of the musculature | 1373 | 0 | 0 | 242 | 0 | 0 | 0 | 702 | 0 | 429 | 0 |  
| HP:0001510 | Growth delay | 1030 | 0 | 0 | 223 | 0 | 0 | 0 | 720 | 0 | 87 | 0 |  



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



