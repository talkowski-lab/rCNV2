# Phenotype Data Curation

We grouped samples into phenotype-matched subsets using a standardized hierarchical ontology.  

This process is described below, and [a summary graphic is provided at the bottom of this page](#summary-samples-per-hpo-term-per-metacohort).  

### Phenotype curation  

Given the variability of phenotypic detail available between different cohorts, we applied a standardized hierarchical phenotype consolidation scheme uniformly across all samples.  

We used the [Human Phenotype Ontology (HPO)](https://hpo.jax.org/app/) to harmonize phenotype data across all [cohorts included in the overall CNV dataset](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV/).  

For each sample, we recursively matched all phenotype information available against HPO keywords, and assigned the corresponding HPO terms to each sample when a match was found.  

We abided by all HPO phenotype definitions with a single exception: we matched all congenital anomalies to `HP:0001197 [Abnormality of prenatal development or birth]`, despite the description for that term excluding fetal structural anomalies. We made this exception because congenital anomalies were a phenotype of particular interest in this study, and there was no clear corresponding existing HPO code for these phenotypes.

The code to perform this HPO conversion is contained in `convert_phenotypes_to_hpo.sh`.  

#### A note on ICD-10 to HPO conversion  

For the UK BioBank dataset, all phenotypes were encoded in [ICD-10 format](https://www.cms.gov/medicare/coding/icd10/) according to [the UK BioBank-sanctioned version of ICD-10](https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=19&nl=1).  

Given the scope of analyses in this study, we reduced the 19,153 unique ICD-10 codes to a smaller subset relevant to this study. This was accomplished through:
1. Automated filtering (`prune_ukbb_icd10_dictionary.py`) to isolate ICD-10 codes with a cohort prevalance of at least 0.01% but no greater than 5%
2. Manual review by a physician to further identify terms unlikely to have a strong genetic component (e.g., infectious .  

Afterwards, phenotypes for each sample were converted from ICD-10 to plain-text indications with `icd10_to_indication.py`.  

Once converted to plain-text indications, UK BioBank samples were subjected to indication-to-HPO conversion as described above.  

### Phenotype hierarchy  

After assigning HPO terms to each sample, we next reduced the number of HPO terms used in this study to those with appreciable sample size and those that were partially non-overlapping with other related HPO terms.

After tabulating the number of samples matching each HPO code, we retained all HPO terms with at least 1,000 samples. 

We next compared shared sample memberships between all pairs of hierarchically related HPO terms:
 * If fewer than 2,000 samples differed between a pair of HPO terms, we retained the more general (i.e., higher-level) term, and excluded the more specific (i.e., lower-level) term.
 * If both terms were equally high-level and siblings (defined as reciprocally sharing at least 50% of their parent terms), we retained the term with the larger sample size.

 Finally, to help control for batch-specific technical artifacts, we required that all HPO terms be represented with at least 500 samples in at least two different [metacohorts (described here)](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV/#case-control-metacohorts).

This process yielded a condensed hierarchical phenotype classification system with 30 distinct HPO terms (including healthy controls), each which had ≥ 1,000 samples in total, ≥ 500 samples from at least two different metacohorts and differed from any related HPO terms by at least 2,000 distinct samples.  

The code to perform this HPO consolidation process described above is also contained in `convert_phenotypes_to_hpo.sh`    

The final set of HPO terms is detailed below:  

| HPO Term | Description | Samples | HPO Tier | Parent Terms | Child Terms |  
| :--- | :--- | ---: | ---: | :--- | :--- |  
| HEALTHY_CONTROL | Unaffected control sample | 399,961 | 1 | NA | NA |  
| HP:0000118 | Phenotypic abnormality | 293,207 | 1 | NA | HP:0031466, HP:0002960, HP:0003011, HP:0100852, HP:0100022, HP:0001627, HP:0100545, HP:0000152, HP:0100753, HP:0012638, HP:0002715, HP:0002597, HP:0025031, HP:0000717, HP:0012443, HP:0000708, HP:0001250, HP:0011446, HP:0012759, HP:0002011, HP:0000729, HP:0012639, HP:0001626, HP:0000752, HP:0000707, HP:0000924, HP:0001507, HP:0001197 |  
| HP:0000707 | Abnormality of the nervous system | 161,862 | 2 | HP:0000118 | HP:0000717, HP:0012759, HP:0031466, HP:0000729, HP:0012639, HP:0100753, HP:0012443, HP:0012638, HP:0000752, HP:0100852, HP:0000708, HP:0100022, HP:0001250, HP:0011446, HP:0002011 |  
| HP:0012638 | Abnormality of nervous system physiology | 122,887 | 3 | HP:0000707, HP:0000118 | HP:0000717, HP:0012759, HP:0031466, HP:0000729, HP:0000752, HP:0100852, HP:0000708, HP:0100022, HP:0001250, HP:0011446, HP:0100753 |  
| HP:0000708 | Behavioral abnormality | 68,954 | 4 | HP:0012638, HP:0000707, HP:0000118 | HP:0000717, HP:0100852, HP:0031466, HP:0100753, HP:0000729 |  
| UNKNOWN | NA | 49,819 | 2 | HP:0000118 | NA |  
| HP:0012639 | Abnormality of nervous system morphology | 39,815 | 3 | HP:0000707, HP:0000118 | HP:0012443, HP:0002011 |  
| HP:0002715 | Abnormality of the immune system | 39,239 | 2 | HP:0000118 | HP:0002960 |  
| HP:0012759 | Neurodevelopmental abnormality | 36,578 | 4 | HP:0012638, HP:0000707, HP:0000118 | NA |  
| HP:0002960 | Autoimmunity | 34,153 | 4 | HP:0002715, HP:0000118 | NA |  
| HP:0002011 | Morphological abnormality of the central nervous system | 28,733 | 4 | HP:0012639, HP:0000707, HP:0000118 | HP:0012443 |  
| HP:0001626 | Abnormality of the cardiovascular system | 26,509 | 2 | HP:0000118 | HP:0001627, HP:0002597, HP:0100545 |  
| HP:0100753 | Schizophrenia | 23,165 | 5 | HP:0000707, HP:0000118, HP:0012638, HP:0000708 | NA |  
| HP:0000729 | Autistic behavior | 19,029 | 5 | HP:0000707, HP:0000118, HP:0012638, HP:0000708 | HP:0000717 |  
| HP:0002597 | Abnormality of the vasculature | 18,663 | 3 | HP:0001626, HP:0000118 | HP:0100545 |  
| HP:0100022 | Abnormality of movement | 17,645 | 4 | HP:0012638, HP:0000707, HP:0000118 | HP:0000752 |  
| HP:0001250 | Seizures | 17,500 | 4 | HP:0012638, HP:0000707, HP:0000118 | NA |  
| HP:0100545 | Arterial stenosis | 16,492 | 7 | HP:0001626, HP:0000118, HP:0002597 | NA |  
| HP:0000717 | Autism | 16,036 | 6 | HP:0000707, HP:0000118, HP:0012638, HP:0000708, HP:0000729 | NA |  
| HP:0000752 | Hyperactivity | 13,386 | 5 | HP:0000707, HP:0100022, HP:0000118, HP:0012638 | NA |  
| HP:0001197 | Abnormality of prenatal development or birth | 10,784 | 2 | HP:0000118 | NA |  
| HP:0000924 | Abnormality of the skeletal system | 10,102 | 2 | HP:0000118 | NA |  
| HP:0031466 | Impairment in personality functioning | 8,472 | 5 | HP:0000707, HP:0000118, HP:0012638, HP:0000708 | HP:0100852 |  
| HP:0000152 | Abnormality of head or neck | 7,902 | 2 | HP:0000118 | NA |  
| HP:0001627 | Abnormal heart morphology | 7,347 | 4 | HP:0001626, HP:0000118 | NA |  
| HP:0025031 | Abnormality of the digestive system | 6,036 | 2 | HP:0000118 | NA |  
| HP:0001507 | Growth abnormality | 3,736 | 2 | HP:0000118 | NA |  
| HP:0100852 | Abnormal fear/anxiety-related behavior | 2,976 | 6 | HP:0000707, HP:0000118, HP:0012638, HP:0000708, HP:0031466 | NA |  
| HP:0012443 | Abnormality of brain morphology | 2,634 | 5 | HP:0000707, HP:0000118, HP:0012639, HP:0002011 | NA |  
| HP:0003011 | Abnormality of the musculature | 2,539 | 2 | HP:0000118 | NA |  
| HP:0011446 | Abnormality of higher mental function | 2,393 | 4 | HP:0012638, HP:0000707, HP:0000118 | NA |  

### HPO terms per cohort

The number of samples per phenotype group per cohort is outlined in the table, below.  

The code to generate this table is provided in `gather_hpo_per_cohort_table.py`.  

| HPO | description | Total | PGC | Cooper | Coe | SSC | UKBB | CHOP | GDX | TSAICG | BCH | TCGA | Epi25k | SickKids | IU |  
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| HEALTHY_CONTROL | Unaffected control sample | 460,750 | 20,277 | 8,329 | 11,256 | 0 | 375,800 | 24,161 | 0 | 4,093 | 0 | 8,670 | 8,164 | 0 | 0 |  
| HP:0000118 | Phenotypic abnormality | 293,207 | 21,094 | 0 | 29,104 | 2,795 | 54,071 | 153,870 | 9,958 | 2,434 | 3,591 | 0 | 12,023 | 2,691 | 1,576 |  
| HP:0000707 | Abnormality of the nervous system | 161,862 | 21,094 | 0 | 29,104 | 2,795 | 29,861 | 54,406 | 4,039 | 2,434 | 2,915 | 0 | 12,023 | 2,691 | 500 |  
| HP:0012638 | Abnormality of nervous system physiology | 122,887 | 21,094 | 0 | 29,104 | 2,795 | 16,679 | 28,948 | 3,847 | 2,434 | 2,846 | 0 | 12,023 | 2,691 | 426 |  
| HP:0000708 | Behavioral abnormality | 68,954 | 21,094 | 0 | 1,851 | 2,795 | 11,895 | 26,364 | 734 | 0 | 1,225 | 0 | 0 | 2,691 | 305 |  
| UNKNOWN | NA | 49,819 | 0 | 0 | 0 | 0 | 768 | 43,603 | 4,962 | 0 | 293 | 0 | 0 | 0 | 193 |  
| HP:0012639 | Abnormality of nervous system morphology | 39,815 | 0 | 0 | 364 | 0 | 13,484 | 24,936 | 714 | 0 | 279 | 0 | 0 | 0 | 38 |  
| HP:0002715 | Abnormality of the immune system | 39,239 | 0 | 0 | 19 | 0 | 6,307 | 32,843 | 36 | 0 | 34 | 0 | 0 | 0 | 0 |  
| HP:0012759 | Neurodevelopmental abnormality | 36,578 | 0 | 0 | 29,104 | 814 | 132 | 671 | 3,423 | 0 | 2,366 | 0 | 0 | 0 | 68 |  
| HP:0002960 | Autoimmunity | 34,153 | 0 | 0 | 1 | 0 | 1,305 | 32,843 | 2 | 0 | 2 | 0 | 0 | 0 | 0 |  
| HP:0002011 | Morphological abnormality of the central nervous system | 28,733 | 0 | 0 | 349 | 0 | 2,428 | 24,936 | 709 | 0 | 273 | 0 | 0 | 0 | 38 |  
| HP:0001626 | Abnormality of the cardiovascular system | 26,509 | 0 | 0 | 757 | 0 | 3,935 | 21,003 | 426 | 0 | 213 | 0 | 0 | 0 | 175 |  
| HP:0100753 | Schizophrenia | 23,165 | 21,094 | 0 | 2 | 0 | 659 | 1,200 | 6 | 0 | 0 | 0 | 0 | 204 | 0 |  
| HP:0000729 | Autistic behavior | 19,029 | 0 | 0 | 1,569 | 2,795 | 16 | 11,071 | 639 | 0 | 853 | 0 | 0 | 1,838 | 248 |  
| HP:0002597 | Abnormality of the vasculature | 18,663 | 0 | 0 | 55 | 0 | 3,581 | 14,948 | 59 | 0 | 17 | 0 | 0 | 0 | 3 |  
| HP:0100022 | Abnormality of movement | 17,645 | 0 | 0 | 317 | 0 | 1,645 | 12,394 | 83 | 2,434 | 325 | 0 | 0 | 427 | 20 |  
| HP:0001250 | Seizures | 17,500 | 0 | 0 | 1,807 | 318 | 1,814 | 1,762 | 539 | 0 | 539 | 0 | 10,673 | 0 | 48 |  
| HP:0100545 | Arterial stenosis | 16,492 | 0 | 0 | 1 | 0 | 1,543 | 14,948 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  
| HP:0000717 | Autism | 16,036 | 0 | 0 | 1,414 | 2,795 | 16 | 11,071 | 574 | 0 | 166 | 0 | 0 | 0 | 0 |  
| HP:0000752 | Hyperactivity | 13,386 | 0 | 0 | 281 | 0 | 0 | 12,359 | 63 | 0 | 238 | 0 | 0 | 427 | 18 |  
| HP:0001197 | Abnormality of prenatal development or birth | 10,784 | 0 | 0 | 939 | 0 | 2,654 | 6,055 | 408 | 0 | 387 | 0 | 0 | 0 | 341 |  
| HP:0000924 | Abnormality of the skeletal system | 10,102 | 0 | 0 | 344 | 0 | 8,012 | 0 | 1,300 | 0 | 357 | 0 | 0 | 0 | 89 |  
| HP:0031466 | Impairment in personality functioning | 8,472 | 0 | 0 | 20 | 0 | 6,297 | 1,730 | 15 | 0 | 410 | 0 | 0 | 0 | 0 |  
| HP:0000152 | Abnormality of head or neck | 7,902 | 0 | 0 | 4,535 | 0 | 809 | 0 | 1,662 | 0 | 848 | 0 | 0 | 0 | 48 |  
| HP:0001627 | Abnormal heart morphology | 7,347 | 0 | 0 | 173 | 0 | 397 | 6,055 | 351 | 0 | 196 | 0 | 0 | 0 | 175 |  
| HP:0025031 | Abnormality of the digestive system | 6,036 | 0 | 0 | 91 | 0 | 4,438 | 1,258 | 146 | 0 | 93 | 0 | 0 | 0 | 10 |  
| HP:0001507 | Growth abnormality | 3,736 | 0 | 0 | 1,221 | 0 | 6 | 757 | 1,066 | 0 | 218 | 0 | 0 | 0 | 468 |  
| HP:0100852 | Abnormal fear/anxiety-related behavior | 2,976 | 0 | 0 | 6 | 0 | 1,552 | 1,304 | 7 | 0 | 107 | 0 | 0 | 0 | 0 |  
| HP:0012443 | Abnormality of brain morphology | 2,634 | 0 | 0 | 325 | 0 | 646 | 689 | 686 | 0 | 256 | 0 | 0 | 0 | 32 |  
| HP:0003011 | Abnormality of the musculature | 2,539 | 0 | 0 | 242 | 0 | 1,130 | 0 | 705 | 0 | 429 | 0 | 0 | 0 | 33 |  
| HP:0011446 | Abnormality of higher mental function | 2,393 | 0 | 0 | 156 | 814 | 176 | 671 | 189 | 0 | 365 | 0 | 0 | 0 | 22 |  

### HPO terms per metacohort

The number of samples per phenotype group per cohort is outlined in the table, below.  

The code to generate this table is provided in `gather_hpo_per_cohort_table.py`.  

| HPO | description | Total | meta1 | meta2 | meta3 | meta4 | mega |  
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: |  
| HEALTHY_CONTROL | Unaffected control sample | 460,750 | 19,585 | 41,204 | 24,161 | 375,800 | 460,750 |  
| HP:0000118 | Phenotypic abnormality | 293,207 | 44,229 | 41,037 | 153,870 | 54,071 | 293,207 |  
| HP:0000707 | Abnormality of the nervous system | 161,862 | 36,558 | 41,037 | 54,406 | 29,861 | 161,862 |  
| HP:0012638 | Abnormality of nervous system physiology | 122,887 | 36,223 | 41,037 | 28,948 | 16,679 | 122,887 |  
| HP:0000708 | Behavioral abnormality | 68,954 | 4,115 | 26,580 | 26,364 | 11,895 | 68,954 |  
| UNKNOWN | NA | 49,819 | 5,448 | 0 | 43,603 | 768 | 49,819 |  
| HP:0012639 | Abnormality of nervous system morphology | 39,815 | 1,395 | 0 | 24,936 | 13,484 | 39,815 |  
| HP:0002715 | Abnormality of the immune system | 39,239 | 89 | 0 | 32,843 | 6,307 | 39,239 |  
| HP:0012759 | Neurodevelopmental abnormality | 36,578 | 34,961 | 814 | 671 | 132 | 36,578 |  
| HP:0002960 | Autoimmunity | 34,153 | 5 | 0 | 32,843 | 1,305 | 34,153 |  
| HP:0002011 | Morphological abnormality of the central nervous system | 28,733 | 1,369 | 0 | 24,936 | 2,428 | 28,733 |  
| HP:0001626 | Abnormality of the cardiovascular system | 26,509 | 1,571 | 0 | 21,003 | 3,935 | 26,509 |  
| HP:0100753 | Schizophrenia | 23,165 | 8 | 21,298 | 1,200 | 659 | 23,165 |  
| HP:0000729 | Autistic behavior | 19,029 | 3,309 | 4,633 | 11,071 | 16 | 19,029 |  
| HP:0002597 | Abnormality of the vasculature | 18,663 | 134 | 0 | 14,948 | 3,581 | 18,663 |  
| HP:0100022 | Abnormality of movement | 17,645 | 745 | 2,861 | 12,394 | 1,645 | 17,645 |  
| HP:0001250 | Seizures | 17,500 | 2,933 | 10,991 | 1,762 | 1,814 | 17,500 |  
| HP:0100545 | Arterial stenosis | 16,492 | 1 | 0 | 14,948 | 1,543 | 16,492 |  
| HP:0000717 | Autism | 16,036 | 2,154 | 2,795 | 11,071 | 16 | 16,036 |  
| HP:0000752 | Hyperactivity | 13,386 | 600 | 427 | 12,359 | 0 | 13,386 |  
| HP:0001197 | Abnormality of prenatal development or birth | 10,784 | 2,075 | 0 | 6,055 | 2,654 | 10,784 |  
| HP:0000924 | Abnormality of the skeletal system | 10,102 | 2,090 | 0 | 0 | 8,012 | 10,102 |  
| HP:0031466 | Impairment in personality functioning | 8,472 | 445 | 0 | 1,730 | 6,297 | 8,472 |  
| HP:0000152 | Abnormality of head or neck | 7,902 | 7,093 | 0 | 0 | 809 | 7,902 |  
| HP:0001627 | Abnormal heart morphology | 7,347 | 895 | 0 | 6,055 | 397 | 7,347 |  
| HP:0025031 | Abnormality of the digestive system | 6,036 | 340 | 0 | 1,258 | 4,438 | 6,036 |  
| HP:0001507 | Growth abnormality | 3,736 | 2,973 | 0 | 757 | 6 | 3,736 |  
| HP:0100852 | Abnormal fear/anxiety-related behavior | 2,976 | 120 | 0 | 1,304 | 1,552 | 2,976 |  
| HP:0012443 | Abnormality of brain morphology | 2,634 | 1,299 | 0 | 689 | 646 | 2,634 |  
| HP:0003011 | Abnormality of the musculature | 2,539 | 1,409 | 0 | 0 | 1,130 | 2,539 |  
| HP:0011446 | Abnormality of higher mental function | 2,393 | 732 | 814 | 671 | 176 | 2,393 |  

A graphical representation of the above table is provided below. This plot is generated by `plot_hpo_per_cohort.R`.  

#### Summary: samples per HPO term per metacohort
![Samples per HPO per metacohort](https://storage.googleapis.com/rcnv_project/public/HPOs_by_metacohort.barplot.jpg)



