# Pooled deconvolution normalization with cluster annotation
# Javier Perales-Paton (C)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)

# INPUT
INPUT_colData <- args[1]
INPUT_rowData <- args[2]
INPUT_mtx <- args[3]
# OUTPUT
OUTPUT_colData <- args[4]
OUTPUT_rowData <- args[5]
OUTPUT_mtx <- args[6]
OUTDIR <- dirname(OUTPUT_mtx)

cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("Matrix")
library("SingleCellExperiment")
library("scran")
library("scater")
source("workflow/src/10Xmat.R")

## Load data
colDat <- readBarcodes(INPUT_colData, header=TRUE)
rowDat <- readFeatures(INPUT_rowData, header=TRUE)
dat <- readMatrix(INPUT_mtx, INPUT_rowData, INPUT_colData, header=TRUE)

rownames(rowDat) <- rowDat$hgnc_symbol
rownames(colDat) <- colDat$CellName

# Summary of cell metadata
## Create a Single-Cell Experiment
sce <- SingleCellExperiment(assays=list("counts"=dat),
			    colData=colDat,
			    rowData=rowDat)

## Normalize data
#NOTE: Params defined by M.Ibrahim
sce = scran::computeSumFactors(sce, 
			       sizes = seq(10, 200, 20),  
			       clusters = sce$Annotation.Level.3, 
			       positive = TRUE)
sce <- logNormCounts(sce)

# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

writeMM(logcounts(sce),	file=OUTPUT_mtx)
write.table(colData(sce), sep=",", file=OUTPUT_colData, 
	    col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(rowData(sce), sep=",", file=OUTPUT_rowData, 
	    col.names=TRUE, row.names=FALSE, quote=FALSE)
