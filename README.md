# Installation Instructions


Install `GFA` from GitHub

```
devtools::install_github("jean997/GFA")
```


or to install with vignette 

```
devtools::install_github("jean997/GFA", build_vignettes = TRUE)
browseVignettes("GFA")
```


## Pipeline

We have written a Snakemake pipeline that implements all of the steps of a GFA analysis and requires only input of configuration parameters and a spreadsheet of data files. This pipeline can be found at 
