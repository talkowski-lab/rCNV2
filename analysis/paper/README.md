# rCNV Paper Analyses  
## Formalized secondary data processing, analysis, and visualization for study manuscript  

This subdirectory contains all code to conduct the formal secondary analyses presented in the rCNV2 manuscript.  

Note that the code in this subdirectory makes extensive use of [processed data](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation) and results from [primary (discovery) analyses](https://github.com/talkowski-lab/rCNV2/tree/master/analysis). Please refer to those directories for descriptions of those analyses.  

## Organization  

This directory is divided into multiple subdirectories, as follows:  

| Subdirectory | Description |
| :--- | :--- |
| [`plot/`](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/paper/plot/) | Code to generate figure panels |  
| [`scripts/`](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/paper/scripts/) | Miscellaneous helper scripts |  
| [`shell/`](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/paper/data_processing/) | Bash code blocks to perform data processing |  
| [`wdl/`](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/paper/wdl/) | WDL wrappers to execute entire secondary analysis workflows |  

