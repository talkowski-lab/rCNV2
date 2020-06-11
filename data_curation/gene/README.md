# Gene Curation  

We curated and annotated all protein-coding genes in the genome for association testing. This process is described below.  

All commands executed to filter the CNV data as descrbed below are contained in `format_genes.sh`.    

## Gene definitions  

We restricted all gene-based analyses to canonical transcripts from autosomal protein-coding genes as defined by [Gencode v19](https://www.gencodegenes.org/human/release_19.html).  

The Gencode definition of "canonical transcript" is [provided here](http://www.ensembl.org/Help/Glossary?id=346).  

We extracted all canonical transcripts from protein-coding genes using `get_canonical_transcripts.py`.   

## Exon-level expression filtering  

We also restricted all gene-based analyses to exons expressed in at least 20% of transcripts in at least one human tissue catalogued by the [Genotype-Tissue Expression (GTEx) Project](https://gtexportal.org/).  

For this filter, we computed the per-exon max `pext` score as defined in [Cummings _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/10.1101/554444v1) across all available tissues.  

Bases from exons with missing `pext` scores were ignored when calculating mean `pext` scores, and exons entirely missing `pext` scores were considered as passing.  

Genes with no exons remaining after exon-level expression filtering were removed outright from all subsequent analyses.  

## Gene features  

For certain analyses, we used per-gene features across a variety of categories.  

These "gene features" are described below.  

For each feature set, we controlled for inter-feature correlation structure by decomposing all features with principal components analysis and retaining the top "eigenfeatures" that explained at least 95% of inter-gene variance.  

In practice, the feature collection and decomposition process was parallelized in [FireCloud/Terra](https://portal.firecloud.org) with `get_gene_metadata.wdl`.  

#### Gene features: genomic  

We collected the following genome-based gene features:  

| Feature name | Abbreviation | Description | Source |  
| :--- | :--- | :--- | :--- |  
| Gene length | `gene_length` | Length of gene body corresponding to canonical transcript (log<sub>10</sub>-scaled) | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) |  
| Dist. to nearest gene | `nearest_gene` | Distance to nearest other gene (log<sub>10</sub>-scaled). Defaults to zero for overlapping genes. | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Genes within 1Mb | `genes_within_1mb` | Count of all other genes within ±1Mb of gene body | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) |  
| Length of coding sequence | `cds_length` | Length of coding sequence (log<sub>10</sub>-scaled) | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Number of exons | `n_exons` | Number of non-overlapping exons | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Min. exon size | `min_exon_size` | Size of smallest exon (log<sub>10</sub>-scaled) | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Median exon size | `med_exon_size` | Median exon size (log<sub>10</sub>-scaled) | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Max. exon size | `max_exon_size` | Size of largest exon (log<sub>10</sub>-scaled) | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Min. intron size | `min_intron_size` | Size of smallest intron (log<sub>10</sub>-scaled). Defauts to zero for single-exon genes. | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Median intron size | `med_intron_size` | Median intron size (log<sub>10</sub>-scaled). Defauts to zero for single-exon genes. | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Max. intron size | `max_intron_size` | Size of largest intron (log<sub>10</sub>-scaled). Defauts to zero for single-exon genes. | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Dist. to telomere | `telomere_dist` | Absolute distance to nearest telomere | [UCSC Genome Browser](http://genome.ucsc.edu) |    
| Dist. to centromere | `centromere_dist` | Absolute distance to nearest centromere | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Norm. dist. to telomere | `telomere_dist_norm` | Distance to nearest telomere, normalized to chromosome arm length | [UCSC Genome Browser](http://genome.ucsc.edu) |    
| Norm. dist. to centromere | `centromere_dist_norm` | Distance to nearest centromere, normalized to chromosome arm length | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| GC content | `gc_pct` | Fraction of G/C nucleotides within gene body | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Repetitive sequence content | `repeat_cov` | Fraction of gene body covered by segmental duplications and/or simple/low-complexity/satellite repeats | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Hypermutable region content | `repeat_cov` | Fraction of gene body covered by hypermutable sequences | Various sources as described in [_Collins\*, Brand\*, et al._ (2019)](https://broad.io/gnomad_sv) |  
| Sequence alignability | `align_100mer` | Mean 100mer alignability score over the gene body | [UCSC Genome Browser](http://genome.ucsc.edu) |  

#### Gene features: expression  

We collected the following expression-based gene features:  

| Feature name | Abbreviation | Description | Source |  
| :--- | :--- | :--- | :--- |  
| Number of expressing tissues | `n_tissues_expressed` | Number of tissues with median expression (TPM) > 1 in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Min. average expression | `median_expression_min` | Minimum of median expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| First quartile of average expression | `median_expression_q1` | First quartile of median expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Mean average expression | `median_expression_mean` | Mean of median expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Third quartile of average expression | `median_expression_q3` | Third quartile of median expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Max. average expression | `median_expression_max` | Maximum of median expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Min. expression MAD | `expression_mad_min` | Minimum of median absolute deviation (MAD) of expression values  across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| First quartile of expression MAD | `expression_mad_q1` | First quartile of median absolute deviation (MAD) of expression values  across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Mean expression MAD | `expression_mad_mean` | Mean of median absolute deviation (MAD) of expression values  across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Third quartile of expression MAD | `expression_mad_q3` | Third quartile of median absolute deviation (MAD) of expression values  across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Max. expression MAD | `expression_mad_max` | Maximum of median absolute deviation (MAD) of expression values  across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Expression components 1-20 | `expression_component_1` ... `expression_component_20` | Top 20 principal components of gene X tissue matrix from GTEx v7, using mean expression across individuals per gene per tissue | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  

All expression values were derived from GTEx v7 (rather than the final v8 release) because GTEx v7 was the last version native to hg19/Gencode v19, which was a direct match for the Gencode version used in this study.  

Note that all expression values from GTEx were transformed as log<sub>10</sub>(x + 1).  

#### Gene features: constraint

We collected the following gene features related to evolutionary constraint:  

| Feature name | Abbreviation | Description | Source |  
| :--- | :--- | :--- | :--- |  
| Probability of loss-of-function intolerance | `gnomad_pLI` | Probability of loss-of-function intolerance calculated in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _bioRxiv_, 2019)](https://www.biorxiv.org/content/10.1101/531210v2) |  
| Probability of complete haplosufficiency | `gnomad_pNull` | Probability of haplosufficiency calculated in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _bioRxiv_, 2019)](https://www.biorxiv.org/content/10.1101/531210v2) |  
| Probability of recessive lethality | `gnomad_pRec` | Probability of recessive lethality calculated in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _bioRxiv_, 2019)](https://www.biorxiv.org/content/10.1101/531210v2) |  
| gnomAD missense obs:exp | `gnomad_oe_mis` | Observed : expected ratio for missense SNVs in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _bioRxiv_, 2019)](https://www.biorxiv.org/content/10.1101/531210v2) |  
| gnomAD loss-of-function obs:exp | `gnomad_oe_lof` | Observed : expected ratio for loss-of-function SNVs in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _bioRxiv_, 2019)](https://www.biorxiv.org/content/10.1101/531210v2) |  
| gnomAD missense obs:exp (upper 90% CI) | `gnomad_oe_mis_upper` | Upper 90% confidence interval of obs:exp ratio for missense SNVs in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _bioRxiv_, 2019)](https://www.biorxiv.org/content/10.1101/531210v2) |  
| gnomAD loss-of-function obs:exp (upper 90% CI) | `gnomad_oe_lof_upper` | Upper 90% confidence interval of obs:exp ratio for loss-of-function SNVs in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _bioRxiv_, 2019)](https://www.biorxiv.org/content/10.1101/531210v2) |  
| gnomAD missense Z-score | `gnomad_mis_z` | Z-score for missense obs:exp ratio in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _bioRxiv_, 2019)](https://www.biorxiv.org/content/10.1101/531210v2) |  
| gnomAD loss-of-function Z-score | `gnomad_lof_z` | Z-score for loss-of-function obs:exp ratio in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _bioRxiv_, 2019)](https://www.biorxiv.org/content/10.1101/531210v2) |  
| ExAC CNV intolerance Z-score | `exac_cnv_z` | Z-score for obs:exp ratio of coding CNVs in ExAC | ExAC v1.0 [(Ruderfer _et al._, _Nat. Genet._, 2016)](https://www.ncbi.nlm.nih.gov/pubmed/27533299) |  
| Haploinsufficiency score | `hurles_hi` | Haploinsufficiency scores predicted based on haplosufficient genes | [Huang _et al._, _PLOS Genetics_, 2010](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1001154) |  
| RVIS | `rvis` | Residual variation intolerance score (RVIS) computed on gnomAD v2.0 | [Petrovski _et al._, _PLOS Genetics_, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23990802) |  
| RVIS percentile | `rvis_pct` | RVIS percentile computed on gnomAD v2.0 | [Petrovski _et al._, _PLOS Genetics_, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23990802) |  
| Promoter GC content | `promoter_gc_pct` | Percentage of C/G nucleotides in gene promoter | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Promoter CpG count | `promoter_cpg_count` | Count of CpG dinucleotides in gene promoter | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Promoter conservation | `promoter_phastcons` | Mean 100-way phastCons score in gene promoter | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Exon conservation | `exon_phastcons` | Mean 100-way hastCons score in exons (weighted by exon size) | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Gene body conservation | `gene_body_phastcons` | Mean 100-way phastCons score across entire gene body | [UCSC Genome Browser](http://genome.ucsc.edu) |  

#### Gene features: chromatin  

_TBD_  

## Gene set definitions  

Throughout our analyses, we reference various subsets of genes. These are defined below.  

Note that all gene sets were required to have a unique match to a gene name from one of the canonical, protein-coding genes curated for this study (as defined above).  

All gene lists in the table below are available from `gs://rcnv_project/cleaned_data/genes/gene_lists/`, which is accessible to users with appropriate GCP permissions.  

| Gene set | Genes | Filename prefix | Source | Description |
| :--- | ---: | :--- | :--- | :---- | 
| All genes | 18,641 | `gencode.v19.canonical.pext_filtered` | Gencode v19 [Harrow _et al._, _Genome Res._, 2012](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | Canonical transcripts from autosomal, protein-coding genes |
| Constrained genes | 3,036 | `gnomad.v2.1.1.lof_constrained` | gnomAD v2.1.1 [Karczewski _et al._, _bioRxiv_, 2019](https://www.biorxiv.org/content/10.1101/531210v3) | pLI ≥ 0.9 or in the first LOEUF sextile |  
| Mutation-tolerant genes | 2,013 | `gnomad.v2.1.1.mutation_tolerant` | gnomAD v2.1.1 [Karczewski _et al._, _bioRxiv_, 2019](https://www.biorxiv.org/content/10.1101/531210v3) | pLI ≥ 0.01, the last third of LOEUF, missense Z-score ≤ 0, missense OEUF ≥ 1, and synonymous Z-score ~ (-3, 3) |  
| Phenotypic abnormality (HP:0000118)-associated genes | 3,825 | `HP0000118.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000118 |
| Abnormality of the nervous system (HP:0000707)-associated genes | 2,875 | `HP0000707.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000707 |
| Abnormality of nervous system physiology (HP:0012638)-associated genes | 2,678 | `HP0012638.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012638 |
| Behavioral abnormality (HP:0000708)-associated genes | 1,227 | `HP0000708.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000708 |
| NA (UNKNOWN)-associated genes | 3,825 | `UNKNOWN.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to UNKNOWN |
| Abnormality of nervous system morphology (HP:0012639)-associated genes | 2,177 | `HP0012639.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012639 |
| Abnormality of the immune system (HP:0002715)-associated genes | 1,427 | `HP0002715.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002715 |
| Neurodevelopmental abnormality (HP:0012759)-associated genes | 1,976 | `HP0012759.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012759 |
| Autoimmunity (HP:0002960)-associated genes | 146 | `HP0002960.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002960 |
| Morphological abnormality of the central nervous system (HP:0002011)-associated genes | 1,998 | `HP0002011.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002011 |
| Abnormality of the cardiovascular system (HP:0001626)-associated genes | 1,983 | `HP0001626.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001626 |
| Schizophrenia (HP:0100753)-associated genes | 52 | `HP0100753.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0100753 |
| Autistic behavior (HP:0000729)-associated genes | 330 | `HP0000729.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000729 |
| Abnormality of the vasculature (HP:0002597)-associated genes | 1,330 | `HP0002597.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002597 |
| Abnormality of movement (HP:0100022)-associated genes | 1,519 | `HP0100022.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0100022 |
| Seizures (HP:0001250)-associated genes | 1,340 | `HP0001250.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001250 |
| Arterial stenosis (HP:0100545)-associated genes | 129 | `HP0100545.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0100545 |
| Autism (HP:0000717)-associated genes | 140 | `HP0000717.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000717 |
| Hyperactivity (HP:0000752)-associated genes | 287 | `HP0000752.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000752 |
| Abnormality of prenatal development or birth (HP:0001197)-associated genes | 674 | `HP0001197.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001197 |
| Abnormality of the skeletal system (HP:0000924)-associated genes | 2,326 | `HP0000924.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000924 |
| Impairment in personality functioning (HP:0031466)-associated genes | 443 | `HP0031466.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0031466 |
| Abnormality of head or neck (HP:0000152)-associated genes | 2,373 | `HP0000152.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000152 |
| Abnormal heart morphology (HP:0001627)-associated genes | 1,016 | `HP0001627.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001627 |
| Abnormality of the digestive system (HP:0025031)-associated genes | 2,085 | `HP0025031.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0025031 |
| Growth abnormality (HP:0001507)-associated genes | 1,894 | `HP0001507.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001507 |
| Abnormal fear/anxiety-related behavior (HP:0100852)-associated genes | 227 | `HP0100852.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0100852 |
| Abnormality of brain morphology (HP:0012443)-associated genes | 1,828 | `HP0012443.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012443 |
| Abnormality of the musculature (HP:0003011)-associated genes | 2,057 | `HP0003011.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0003011 |
| Abnormality of higher mental function (HP:0011446)-associated genes | 1,912 | `HP0011446.HPOdb` | HPO database (accessed Mar 2020) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0011446 |
| ClinGen dominant haploinsufficient genes (high confidence) | 208 | `ClinGen.hc_haploinsufficient` | ClinGen gene curation map (accessed Mar 2020) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant haploinsufficient genes scored at high confidence by ClinGen |
| ClinGen dominant haploinsufficient genes (medium confidence) | 43 | `ClinGen.mc_haploinsufficient` | ClinGen gene curation map (accessed Mar 2020) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant haploinsufficient genes scored at medium confidence by ClinGen |
| ClinGen dominant haploinsufficient genes (high or medium confidence) | 251 | `ClinGen.hmc_haploinsufficient` | ClinGen gene curation map (accessed Mar 2020) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant haploinsufficient genes scored at high or medium confidence by ClinGen |
| ClinGen dominant haploinsufficient genes (low confidence) | 78 | `ClinGen.lc_haploinsufficient` | ClinGen gene curation map (accessed Mar 2020) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant haploinsufficient genes scored at low confidence by ClinGen |
| ClinGen dominant haploinsufficient genes (any confidence) | 329 | `ClinGen.all_haploinsufficient` | ClinGen gene curation map (accessed Mar 2020) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant haploinsufficient genes scored at any confidence by ClinGen |
| ClinGen dominant triplosensitive genes (high confidence) | 1 | `ClinGen.hc_triplosensitive` | ClinGen gene curation map (accessed Mar 2020) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant triplosensitive genes scored at high confidence by ClinGen |
| ClinGen dominant triplosensitive genes (medium confidence) | 1 | `ClinGen.mc_triplosensitive` | ClinGen gene curation map (accessed Mar 2020) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant triplosensitive genes scored at medium confidence by ClinGen |
| ClinGen dominant triplosensitive genes (high or medium confidence) | 2 | `ClinGen.hmc_triplosensitive` | ClinGen gene curation map (accessed Mar 2020) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant triplosensitive genes scored at high or medium confidence by ClinGen |
| ClinGen dominant triplosensitive genes (low confidence) | 13 | `ClinGen.lc_triplosensitive` | ClinGen gene curation map (accessed Mar 2020) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant triplosensitive genes scored at low confidence by ClinGen |
| ClinGen dominant triplosensitive genes (any confidence) | 15 | `ClinGen.all_triplosensitive` | ClinGen gene curation map (accessed Mar 2020) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant triplosensitive genes scored at any confidence by ClinGen |
| DECIPHER dominant dev. disorder genes (loss-of-function mechanism; high confidence) | 300 | `DDG2P.hc_lof` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with a loss-of-function mechanism scored as "confirmed" in DECIPHER |
| DECIPHER dominant dev. disorder genes (loss-of-function mechanism; medium confidence) | 135 | `DDG2P.mc_lof` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with a loss-of-function mechanism scored as "probable" in DECIPHER |
| DECIPHER dominant dev. disorder genes (loss-of-function mechanism; high or medium confidence) | 435 | `DDG2P.hmc_lof` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with a loss-of-function mechanism scored as "confirmed" or "probable" in DECIPHER |
| DECIPHER dominant dev. disorder genes (loss-of-function mechanism; low confidence) | 63 | `DDG2P.lc_lof` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with a loss-of-function mechanism scored as "possible" in DECIPHER |
| DECIPHER dominant dev. disorder genes (loss-of-function mechanism; any confidence) | 498 | `DDG2P.all_lof` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with a loss-of-function mechanism scored as "confirmed", "probable", or "possible" in DECIPHER |
| DECIPHER dominant dev. disorder genes (activating/gain-of-function mechanism; high confidence) | 82 | `DDG2P.hc_gof` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an activating/gain-of-function mechanism scored as "confirmed" in DECIPHER |
| DECIPHER dominant dev. disorder genes (activating/gain-of-function mechanism; medium confidence) | 25 | `DDG2P.mc_gof` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an activating/gain-of-function mechanism scored as "probable" in DECIPHER |
| DECIPHER dominant dev. disorder genes (activating/gain-of-function mechanism; high or medium confidence) | 107 | `DDG2P.hmc_gof` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an activating/gain-of-function mechanism scored as "confirmed" or "probable" in DECIPHER |
| DECIPHER dominant dev. disorder genes (activating/gain-of-function mechanism; low confidence) | 9 | `DDG2P.lc_gof` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an activating/gain-of-function mechanism scored as "possible" in DECIPHER |
| DECIPHER dominant dev. disorder genes (activating/gain-of-function mechanism; any confidence) | 116 | `DDG2P.all_gof` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an activating/gain-of-function mechanism scored as "confirmed", "probable", or "possible" in DECIPHER |
| DECIPHER dominant dev. disorder genes (other/unknown coding mechanism; high confidence) | 145 | `DDG2P.hc_other` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an other/unknown coding mechanism scored as "confirmed" in DECIPHER |
| DECIPHER dominant dev. disorder genes (other/unknown coding mechanism; medium confidence) | 96 | `DDG2P.mc_other` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an other/unknown coding mechanism scored as "probable" in DECIPHER |
| DECIPHER dominant dev. disorder genes (other/unknown coding mechanism; high or medium confidence) | 241 | `DDG2P.hmc_other` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an other/unknown coding mechanism scored as "confirmed" or "probable" in DECIPHER |
| DECIPHER dominant dev. disorder genes (other/unknown coding mechanism; low confidence) | 36 | `DDG2P.lc_other` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an other/unknown coding mechanism scored as "possible" in DECIPHER |
| DECIPHER dominant dev. disorder genes (other/unknown coding mechanism; any confidence) | 277 | `DDG2P.all_other` | DECIPHER/DDG2P (accessed Mar 2020) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an other/unknown coding mechanism scored as "confirmed", "probable", or "possible" in DECIPHER |
