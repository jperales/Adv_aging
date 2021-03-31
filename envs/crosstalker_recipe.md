## How to setup a local env for CrossTalkeR
> Motivation: CrossTalkeR is in development. Thus it is not provided via conda channels. Here we show how to setup it

Create a conda environment with all deps so far from channels
```bash
mamba create -p envs/crosstalker r-base=4.0 r-devtools r-units r-sf r-matrix r-munsell r-nlme r-scales r-mgcv r-reshape2 r-ggplot2 r-bit64 r-RSQLite bioconductor-IRanges bioconductor-biobase bioconductor-GO.db bioconductor-GO.db r-lambda.r r-RcppEigen r-tweenr r-gridExtra r-futile.logger r-tidyr r-rvcheck r-graphlayouts r-tidygraph r-viridis r-ggrepel r-ggforce r-dplyr bioconductor-qvalue bioconductor-fgsea bioconductor-DO.db bioconductor-BiocParallel bioconductor-GOSemSim r-shadowtext r-scatterpie r-igraph r-ggraph bioconductor-DOSE r-cowplot r-progress r-cellranger r-readr r-readxl r-haven r-conquer r-MatrixModels r-sp r-lme4 r-rio r-maptools r-quantreg r-pbkrtest r-car r-broom r-rstatix r-ggsignif r-ggsci r-coda r-statnet.common r-class r-E1071 r-classInt r-maps r-rttf2pt1 r-extrafontdb r-Cairo r-png r-clue r-getoptlong r-circlize bioconductor-enrichplot r-downloader r-gridGraphics r-ggpubr r-factominer r-dendextend r-networkdynamic r-sna r-gsw r-tinytex r-mapproj r-dichromat r-qpdf r-extrafont bioconductor-complexheatmap bioconductor-org.Hs.eg.db bioconductor-clusterProfiler r-factoextra r-ggalluvial r-oce r-rmarkdown r-patchwork r-pals
```

Init R session within the env
```bash
conda activate envs/crosstalker
$CONDA_PREFIX/bin/R 
```

Run in R session:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DO.db")


devtools::install_github("https://github.com/CostaLab/CrossTalkeR")
```

Done.
