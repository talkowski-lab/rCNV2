# Gene Curation  

## Gene curation  

We curated and annotated all protein-coding genes in the genome for association testing. This process is described below.  

All commands executed to filter the CNV data as descrbed below are contained in `format_genes.sh`.    

## Gene definitions  

We restricted all gene-based analyses to canonical transcripts from autosomal protein-coding genes as defined by [Gencode v19](https://www.gencodegenes.org/human/release_19.html).  

The Gencode definition of "canonical transcript" is [provided here](http://www.ensembl.org/Help/Glossary?id=346).  

We extracted all canonical transcripts from protein-coding genes using `get_canonical_transcripts.py`.  

After filtering, we retained 19,346 genes for further analysis.  

## Gene set definitions  

For certain analyses, we considered various subsets of genes. These are defined below.  

Note that all gene sets were required to have a unique match to a gene name from one of the 19,346 canonical, protein-coding genes curated for this study (as defined above).  

| Gene set | # genes | Source | Description |  
| :--- | ---: | :--- | :--- |  
| Constrained genes | 2,939 | [Karczewski et al., bioRxiv (2019)](https://www.biorxiv.org/content/10.1101/531210v3) | Genes with pLI ≥ 0.95 |  


