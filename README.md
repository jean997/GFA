## Overview

Genetic Factor Analysis (GFA) is an R package that uses factor analysis to identify common genetic factors shared by multiple traits. These factors are often interpretable as shared unobserved biological processes. Our manuscript describing GFA can be foun [here](url). 

The vignette "Running GFA with Simulated Data" provides an overview of GFA and demonstrates how to run it. This vignette also discusses some nuances like choice of estimation method for nuisance parameters and understanding the scale of the output estimates. 

For applications of GFA to GWAS data, most of the computational effort goes into merging and formatting data, which can be onerous, especially with a large number of traits from many sources. To ease this process, we have provided a [Snakemake pipeline  to automate these steps here](https://github.com/jean997/gfa_pipeline). Using this pipeline does not require any coding or familiarity with Snakemake. You will only need to create a spreadsheet describing your data files and edit a human readable configuration file. In most cases, the pipeline will take care of all necessary pre-processing and you will not need to edit, process, or "munge" your data at all before using it. Detailed pipeline instructions can be found on the pipeline github, linked above. 


## Installation Instructions


Install `GFA` from GitHub

```
devtools::install_github("jean997/GFA")
```


or to install with vignette 

```
devtools::install_github("jean997/GFA", build_vignettes = TRUE)
browseVignettes("GFA")
```

