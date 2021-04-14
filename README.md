# Rare Copy Number Variation Project

Aggregation and analyses of rare CNVs across human disease phenotypes  

Copyright (c) 2019-Present, [Ryan L. Collins](mailto:rlcollins@g.harvard.edu) and the Talkowski Laboratory.  
Distributed under terms of the [MIT License](/LICENSE) (see `LICENSE`).  


## Table of Contents  

_Note: most subdirectories have their own documentation, which outlines specific details pertaining to that topic_

| Subdirectory | Description |
| --- | :--- |
| [`analysis/`](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/) | Code for CNV analyses |
| [`config/`](https://github.com/talkowski-lab/rCNV2/tree/master/config/) | Configuration files for default parameters, etc. |
| [`data_curation/`](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/) | Code for filtering & curation of CNV and annotation data |
| [`docker/`](https://github.com/talkowski-lab/rCNV2/tree/master/docker/) | Dockerfile to build rCNV Docker image |
| [`refs/`](https://github.com/talkowski-lab/rCNV2/tree/master/refs/) | Small reference files |
| [`source/`](https://github.com/talkowski-lab/rCNV2/tree/master/source/) | Source code |  
| [`utils/`](https://github.com/talkowski-lab/rCNV2/tree/master/utils/) | Miscellaneous small utilities and scripts |  


## Citation  

If you use the code in this repository, please cite:  
> Collins _et al._, _A cross-disorder dosage sensitivity map of the human genome_, _medRxiv_ (2021). DOI: [10.1101/2021.01.26.21250098](https://doi.org/10.1101/2021.01.26.21250098)  


### Data access  

All files pertaining to this project are stored in a protected Google Cloud bucket, here:  
`gs://rcnv_project/`

Note that permissions must be granted per user prior to data access.  


### Docker  

System dependencies for rCNV data processing & analyses are managed via Docker.  

The public Docker image for this project is `gcr.io/gnomad-wgs-v2-sv/rcnv`.  
