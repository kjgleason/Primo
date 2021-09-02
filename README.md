# Primo: Package in R for Integrative Multi-Omics association analysis

The goal of `Primo` is to provide computationally efficient tools to integrate data across phenotypes, cell/tissue types, populations or sources when performing joint analysis of multi-omics data.

## Setup

Please note that this package uses functions from the `limma` package, which is downloadable from [Bioconductor](https://www.bioconductor.org), and the `lcmix` package, which is downloadable from [R-Forge](https://r-forge.r-project.org). If installation of `Primo` fails and you have not yet installed the `limma` or `lcmix` packages, please try running the following commands:

  ```R
  if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
  BiocManager::install("limma")
  
  install.packages("MASS")
  install.packages("matrixStats")
  install.packages("nnls")
  install.packages("R.methodsS3")
  install.packages("lcmix",repos="http://r-forge.r-project.org")
  ```

Once you have installed `limma` and `lcmix`, you can install and load functions from `Primo`:

  ```R
  devtools::install_github("kjgleason/Primo")
  library("Primo")
  ```

## Citation

To cite `Primo` in publications, please use:

Kevin J. Gleason, Fan Yang, Brandon L. Pierce, Xin He, and Lin S. Chen. Primo: integration of multiple GWAS and omics QTL summary statistics for elucidation of molecular mechanisms of trait-associated SNPs and detection of pleiotropy in complex traits. Genome Biol (2020); 21(1):236. doi:10.1186/s13059-020-02125-w.

## Licenses
The code is released under [GNU General Public License (GPL)](https://opensource.org/licenses/gpl-license).
