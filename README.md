# Single Cell RNAseq Testis Atlas

This repository contains the code used to analyse the data for the paper "".

The repository is structured as an R package. Generic functions that are used multiple times are in the R/ directory, other analysis scripts are in the analysis/ directory.

To install this package, run:
```
remotes::install_githuib("myersgroup/testisAtlas")
```

Data used in the study is kept in the data/ directory but is not uploaded to github due to the large size. Many of the functions assume that a common set of object are loaded in the environment: "datat" a data.table containing cellwise metadata, and "results" a list object containing the results of the SDA composition. These objects can be downloaded from: https://zenodo.org/ The functions are not particularly designed for use outside of this analysis but rather to document the steps taken.

The shinyApp contains code for the interactive web application for interactively exploring gene expression: http://www.stats.ox.ac.uk/~wells/testisAtlas.html