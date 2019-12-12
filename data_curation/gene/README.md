# Gene Curation  

## Gene curation  

We curated and annotated all protein-coding genes in the genome for association testing. This process is described below.  

All commands executed to filter the CNV data as descrbed below are contained in `format_genes.sh`.    

## Gene definitions  

We restricted all gene-based analyses to canonical transcripts from autosomal protein-coding genes as defined by [Gencode v19](https://www.gencodegenes.org/human/release_19.html).  

The Gencode definition of "canonical transcript" is [provided here](http://www.ensembl.org/Help/Glossary?id=346).  

We extracted all canonical transcripts from protein-coding genes using `get_canonical_transcripts.py`.   

## Gene features  

For certain analyses, we used per-gene features across a variety of categories.  

These "gene features" are described below.  

In practice, the feature collection process was parallelized parallelized in [FireCloud/Terra](https://portal.firecloud.org) with `get_gene_metadata.wdl`.  

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

TBD

#### Gene features: chromatin  

TBD

#### Gene features: variation  

TBD

#### Other gene features TBD  

## Gene set definitions  

Throughout our analyses, we reference various subsets of genes. These are defined below.  

Note that all gene sets were required to have a unique match to a gene name from one of the canonical, protein-coding genes curated for this study (as defined above).  

All gene lists in the table below are available from `gs://rcnv_project/cleaned_data/genes/gene_lists/`, which is accessible to users with appropriate GCP permissions.  

| Gene set | Genes | Filename prefix | Source | Description |
| :--- | ---: | :--- | :--- | :---- | 
| All genes | 19,346 | `gencode.v19.canonical` | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | Canonical transcripts from autosomal, protein-coding genes |
| Constrained genes | 3,501 | `gnomad.v2.1.1.lof_constrained` | gnomAD v2.1.1 [Karczewski _et al._, _bioRxiv_, 2019](https://www.biorxiv.org/content/10.1101/531210v3) | pLI ≥ 0.9 or in the first LOEUF sextile |
| Phenotypic abnormality (HP:0000118)-associated genes | 4,231 | `HP0000118.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000118 |
| Abnormality of the nervous system (HP:0000707)-associated genes | 2,822 | `HP0000707.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000707 |
| Abnormality of nervous system physiology (HP:0012638)-associated genes | 2,628 | `HP0012638.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012638 |
| Behavioral abnormality (HP:0000708)-associated genes | 1,196 | `HP0000708.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000708 |
| NA (UNKNOWN)-associated genes | 4,231 | `UNKNOWN.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to UNKNOWN |
| Abnormality of nervous system morphology (HP:0012639)-associated genes | 2,144 | `HP0012639.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012639 |
| Abnormality of the immune system (HP:0002715)-associated genes | 1,398 | `HP0002715.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002715 |
| Neurodevelopmental abnormality (HP:0012759)-associated genes | 1,931 | `HP0012759.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012759 |
| Autoimmunity (HP:0002960)-associated genes | 137 | `HP0002960.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002960 |
| Morphological abnormality of the central nervous system (HP:0002011)-associated genes | 1,961 | `HP0002011.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002011 |
| Abnormality of the cardiovascular system (HP:0001626)-associated genes | 1,938 | `HP0001626.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001626 |
| Schizophrenia (HP:0100753)-associated genes | 52 | `HP0100753.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0100753 |
| Autistic behavior (HP:0000729)-associated genes | 302 | `HP0000729.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000729 |
| Abnormality of the vasculature (HP:0002597)-associated genes | 1,286 | `HP0002597.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002597 |
| Abnormality of movement (HP:0100022)-associated genes | 1,478 | `HP0100022.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0100022 |
| Seizures (HP:0001250)-associated genes | 1,303 | `HP0001250.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001250 |
| Arterial stenosis (HP:0100545)-associated genes | 129 | `HP0100545.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0100545 |
| Autism (HP:0000717)-associated genes | 129 | `HP0000717.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000717 |
| Hyperactivity (HP:0000752)-associated genes | 274 | `HP0000752.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000752 |
| Abnormality of prenatal development or birth (HP:0001197)-associated genes | 645 | `HP0001197.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001197 |
| Abnormality of the skeletal system (HP:0000924)-associated genes | 2,280 | `HP0000924.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000924 |
| Impairment in personality functioning (HP:0031466)-associated genes | 431 | `HP0031466.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0031466 |
| Abnormality of head or neck (HP:0000152)-associated genes | 2,312 | `HP0000152.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000152 |
| Abnormal heart morphology (HP:0001627)-associated genes | 989 | `HP0001627.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001627 |
| Abnormality of the digestive system (HP:0025031)-associated genes | 2,040 | `HP0025031.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0025031 |
| Growth abnormality (HP:0001507)-associated genes | 1,838 | `HP0001507.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001507 |
| Abnormal fear/anxiety-related behavior (HP:0100852)-associated genes | 223 | `HP0100852.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0100852 |
| Abnormality of brain morphology (HP:0012443)-associated genes | 1,795 | `HP0012443.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012443 |
| Abnormality of the musculature (HP:0003011)-associated genes | 2,023 | `HP0003011.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0003011 |
| Abnormality of higher mental function (HP:0011446)-associated genes | 1,872 | `HP0011446.HPOdb` | HPO database (accessed Nov 2019) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0011446 |  

