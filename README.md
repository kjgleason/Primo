# Primo: Package in R for Integrative Multi-Omics analysis

The goal of `Primo` is to provide computationally efficient tools to integrate data across phenotypes, cell/tissue types, populations or sources when performing joint analysis of multi-omics data.

## Setup

Please note that this package uses functions from the `limma` package, which is downloadable from [Bioconductor](https://www.bioconductor.org). If you have not yet installed the `limma` package, please run the following commands prior to installing `Primo`:

  ```R
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
  ```

Once you have installed `limma`, you can install and load functions from `Primo`:

  ```R
  devtools::install_github("kjgleason/primo")
  library("primo")
  ```

## Citation

To cite `Primo` in publications, please use:

Kevin J. Gleason, Fan Yang, Brandon L. Pierce, Xin He, and Lin S. Chen. Integrating GWAS and omics QTL summary statistics in mapping multiple complex and omics trait associations. Manuscript in preparation.
