# Create feature plot for query
# Javier Perales-Paton (C)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)

# args <- c(
# 	  "out/ann/Contrasts/old-young/AA/logNormSCT_harmony_barcodes.tsv",
# 	"out/norm/Contrasts/old-young/AA/logNormSCT_features.tsv",
# 	"out/norm/Contrasts/old-young/AA/logNormSCT_counts.mtx",
# 	"out/norm/Contrasts/old-young/AA/logNormSCT_data.mtx",
# 	"young-AA",
#  	  "PDGFRA",
#  	  "Annotation.Level.2",
#  	  "png",
#  	  "/tmp/test.png"
# 	)
# 

  
INPUT_colData <- args[1]
INPUT_rowData <- args[2]
INPUT_cnts <- args[3]
INPUT_data <- args[4]

ds <- args[5]
query <- args[6]
ann <- args[7]
figOUT <- args[8]
OUTPUT <- args[9]
OUTDIR <- dirname(OUTPUT)

#---- Functions
getFIGdim <- function(nrows, ncols) {
	return(c(nrows*5, ncols*(8/25)))
}

### Load libraries
library("Matrix")
library("Seurat")
library("cowplot")
library("ggplot2")
source("workflow/src/10Xmat.R")

## Load data
CNT <- readMatrix(INPUT_cnts, INPUT_rowData, INPUT_colData, header=TRUE)
NORM <- readMatrix(INPUT_data, INPUT_rowData, INPUT_colData, header=TRUE)

colDat <- readBarcodes(INPUT_colData, header=TRUE)
rownames(colDat) <- colDat$CellName


## Create a SeuratObject
S <- Seurat::CreateSeuratObject(counts=CNT,
				project = ds)
S <- SetAssayData(S, slot="data", new.data=NORM) 
S <- AddMetaData(S, metadata = colDat)

## Visualization
if(!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive=TRUE);
figDim <- getFIGdim(length(query), length(unique(S@meta.data[, ann])))
if(figOUT=="png") png(OUTPUT, width = figDim[1], height = figDim[2], units = "in", res = 300);
if(figOUT=="pdf") pdf(OUTPUT, width = figDim[1], height = figDim[2], compress=FALSE);
if(figOUT=="tiff") tiff(OUTPUT, width = figDim[1], height = figDim[2], units = "in", compression="none", res = 300);

# Make it
gg <- VlnPlot(S, feature=query, group.by=ann, split.by="Group", 
	      split.plot=FALSE, col=c("grey", "red")) + 
	theme(axis.title.x=element_blank())

print(gg)

# Save it
dev.off()

