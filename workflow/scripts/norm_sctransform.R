# SCTransform on single-cell data
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
OUTPUT_mtx2 <- args[7]
OUTPUT_sel <- args[8]
ID <- args[9]
OUTDIR <- dirname(OUTPUT_mtx)

cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("Matrix")
library("Seurat")
library("sctransform")
library("glmGamPoi")
source("workflow/src/10Xmat.R")
# src: https://stackoverflow.com/questions/40536067/how-to-adjust-future-global-maxsize-in-r
# 6 globals that need to be exported for the future expression (‘FUN()’) is 2.54 GiB
# This exceeds the maximum allowed size of 500.00 MiB
options(future.globals.maxSize= (1024*3)*1024^2)

## Load data
colDat <- readBarcodes(INPUT_colData, header=TRUE)
rowDat <- readFeatures(INPUT_rowData, header=TRUE)
rownames(rowDat) <- rowDat[,1]
dat <- readMatrix(INPUT_mtx, INPUT_rowData, INPUT_colData, header=TRUE)

## Create a Single-Cell Experiment
S <- Seurat::CreateSeuratObject(counts = dat,
				project = ID)
if(ncol(colDat)>1) {
	colDat2 <- colDat
	if("barcode" %in% colnames(colDat2)) {
		rownames(colDat2) <- colDat2$barcode
		col_dups <- c("barcode", "nCounts_RNA", "nFeatures_RNA")
	} else if ("CellName" %in% colnames(colDat2)) {
		rownames(colDat2) <- colDat2$CellName
		col_dups <- c("barcode", "nCounts_RNA", "nFeatures_RNA")
	} else {
		rownames(colDat2) <- colDat2[,1]
		col_dups <- c(colnames(colDat2)[1], "nCounts_RNA", "nFeatures_RNA")
	}
	colDat2 <- colDat2[,which(!colnames(colDat2) %in% col_dups)]
	S <- AddMetaData(S, metadata = colDat2)
}

## Normalization by scTransform on each individual seq run/library tissue sample
#NOTE: replaces NormalizeData(), ScaleData() and FindVariableFeatures() at once
#NOTE: creates a new assay called 'SCT' with the normalized data in 'data' slot
S <- SCTransform(S, method = "glmGamPoi", verbose = TRUE)

# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

writeMM(GetAssayData(S, assay="SCT", slot="data"), file=OUTPUT_mtx)
writeMM(as(GetAssayData(S, assay="SCT", slot="scale.data"),"dgCMatrix"), file=OUTPUT_mtx2)
cat(VariableFeatures(S), sep="\n", file=OUTPUT_sel)
write.table(S@meta.data, sep="\t", file=OUTPUT_colData, 
	    col.names=NA, row.names=TRUE, quote=FALSE)
write.table(rowDat[rownames(S), ], sep="\t", file=OUTPUT_rowData, 
	    col.names=TRUE, row.names=FALSE, quote=FALSE)
