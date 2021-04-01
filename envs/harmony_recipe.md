## How to setup a local env for CrossTalkeR
> Motivation: CrossTalkeR is in development. Thus it is not provided via conda channels. Here we show how to setup it

Create a conda environment with all deps so far from channels
```bash
mamba create -p envs/harmony r-base=4.0 r-devtools r-dplyr r-cowplot r-tidyr r-ggplot2 r-irlba r-Matrix r-tibble r-rlang bioconductor-SingleCellExperiment r-seurat r-rmarkdown r-knitr r-testthat bioconductor-BiocStyle r-Rcpp
```

Init R session within the env
```bash
conda activate envs/harmony
$CONDA_PREFIX/bin/R 
```

Run in R session:
```r
devtools::install_github("https://github.com/immunogenomics/harmony")
```

Done.
