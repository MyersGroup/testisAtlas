[![Travis-CI Build Status](https://travis-ci.org/MyersGroup/testisAtlas.svg?branch=master)](https://travis-ci.org/MyersGroup/testisAtlas)
[![DOI](https://img.shields.io/badge/DOI-10.7554%2FeLife.43966-brightgreen)](https://doi.org/10.7554/eLife.43966)
[![DOI](https://zenodo.org/badge/140632831.svg)](https://zenodo.org/badge/latestdoi/140632831)

# Single Cell RNAseq Testis Atlas

This repository contains the code used to analyse the data for the paper "Unified single-cell analysis of testis gene regulation and pathology in 5 mouse strains" doi:[10.7554/eLife.43966](https://doi.org/10.7554/eLife.43966).

For a quick start guide on using this package to explore the data check out the >>> [vignette](vignettes/vignette.md) <<< (with pictures!).

To install this package, run:
```{r}
remotes::install_github("myersgroup/testisAtlas")
```

The repository is structured as an R package. Generic functions that are used multiple times are in the R/ directory, other analysis scripts are in the analysis/ directory. Data used in the study is kept in the data/ directory but is not uploaded to github due to the large size. However R objects required to run the functions (e.g. "cell_data" a data.table containing cellwise metadata, and "SDAresults" a list object containing the results of the SDA composition) can be downloaded from: [10.5281/zenodo.3233870](https://doi.org/10.5281/zenodo.3233870). The raw data is avaliable at [GEO: GSE113293](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113293).

The shinyApp folder contains code for the web application for interactively exploring gene expression and components: http://www.stats.ox.ac.uk/~wells/testisAtlas.html
