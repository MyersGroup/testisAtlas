[![Travis-CI Build Status](https://travis-ci.org/MyersGroup/testisAtlas.svg?branch=master)](https://travis-ci.org/MyersGroup/MtestisAtlas)
[![DOI](https://zenodo.org/badge/140632831.svg)](https://zenodo.org/badge/latestdoi/140632831)

# Single Cell RNAseq Testis Atlas

This repository contains the code used to analyse the data for the paper "Unified single-cell analysis of testis gene regulation and pathology in 5 mouse strains".

The repository is structured as an R package. Generic functions that are used multiple times are in the R/ directory, other analysis scripts are in the analysis/ directory.

To install this package, run:
```
remotes::install_githuib("myersgroup/testisAtlas")
```

Data used in the study is kept in the data/ directory but is not uploaded to github due to the large size. However R objects that many of the functions assume are loaded in the environment (e.g. "datat" a data.table containing cellwise metadata, and "results" a list object containing the results of the SDA composition) can be downloaded from: [10.5281/zenodo.1341372](https://doi.org/10.5281/zenodo.1341372). The raw data is avaliable at [GEO: GSE113293](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113293). The functions are not particularly designed for use outside of this analysis but rather to document the steps taken.

The shinyApp folder contains code for the web application for interactively exploring gene expression and components: http://www.stats.ox.ac.uk/~wells/testisAtlas.html