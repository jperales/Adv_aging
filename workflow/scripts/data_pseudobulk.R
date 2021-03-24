# Pseudobulking cell type population
# Javier Perales-Paton (C)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)

# args <- c(list.files("out/data2/Contrasts/old-young/AA/", full.names=TRUE),
# 	  "Annotation.Level.1",
# 	  "out/minibulk/Contrasts/old-young/AA/")

# INPUT
INPUT_colData <- args[1]
INPUT_rowData <- args[2]
INPUT_mtx <- args[3]

# PARAMS
ann <- args[4]
OUTDIR <- args[5]

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

## Create a Single-Cell Experiment
sce <- SingleCellExperiment(assays=list("counts"=dat),
			    colData=colDat,
			    rowData=rowDat)

## Pseudobulking
summed <- aggregateAcrossCells(sce, id=colData(sce)[, c(ann, "Sample")])
colnames(summed) <- paste0(colData(summed)$Sample, "_", colData(summed)[[ann]])

# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

clusters <- unique(colData(summed)[[ann]])
for(cluster in clusters) {
	idx <- which(colData(summed)[[ann]] == cluster)
	fls <- paste0(OUTDIR, cluster, "_" ,c("counts.tsv", "samples.tsv"))

	cnt <- counts(summed)[, idx]
	target <- colData(summed)[idx, ]

	cat(paste0("[INFO] : saving data for '", cluster,"' at ",OUTDIR,"\n"),
	    file=stdout())
	write.table(cnt, file=fls[1], sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
	write.table(target, file=fls[2], sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
}

