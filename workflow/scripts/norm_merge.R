# Merge normalized matrices for a given group/Contrast
# Javier Perales-Paton (C)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)

# args <- c("out/data2/Groups/young/AA/barcodes.tsv",
# 	  "out/data2/Groups/young/AA/features.tsv",
# 	  "out/data2/Groups/young/AA/matrix.mtx",
# 	  "young",
# 	  "CA",
# 	"out/norm/Groups/young/AA/barcodes.tsv",
# 	"out/norm/Groups/young/AA/features.tsv",
# 	"out/norm/Groups/young/AA/counts.mtx",
# 	"out/norm/Groups/young/AA/data.mtx",
# 	"out/norm/Groups/young/AA/scaled.tsv",
# 	"out/norm/Groups/young/AA/HGV.txt",
#	)
 

# INPUT
INPUT_colData <- args[1]
INPUT_rowData <- args[2]
INPUT_mtx <- args[3]

# PARAM
GID <- args[4]
REGION <- args[5]
# OUTPUT
OUTPUT_colData <- args[6]
OUTPUT_rowData <- args[7]
OUTPUT_mtx1 <- args[8]
OUTPUT_mtx2 <- args[9]
OUTPUT_mtx3 <- args[10]
OUTPUT_sel <- args[11]
OUTDIR <- dirname(OUTPUT_mtx1)

## Print internal variables
cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), c("args"))) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("Matrix")
library("Seurat")
library("sctransform")
library("glmGamPoi")
source("workflow/src/10Xmat.R")
options(future.globals.maxSize= (1024*3)*1024^2)

## Load data
colDat <- readBarcodes(INPUT_colData, header=TRUE)
rowDat <- readFeatures(INPUT_rowData, header=TRUE)
rownames(rowDat) <- rowDat[,1]
dat <- readMatrix(INPUT_mtx, INPUT_rowData, INPUT_colData, header=TRUE)

#NOTE: we first create a SeuratObject with the merged data. Then we split it

S <- Seurat::CreateSeuratObject(counts = dat,
				project = GID)
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


# SCT normalization without transforming data (scaling and centering),
#NOTE: SCT normalization is done again as a "Group". See source link
# source : https://github.com/immunogenomics/harmony/issues/41#issuecomment-617032339
SL <- SplitObject(S, split.by = "orig.ident")
SL <- lapply(SL, function(S) {
		S <- SCTransform(S, method = "glmGamPoi", 
				 verbose=FALSE, 
				 do.scale = FALSE, 
				 do.center = FALSE)
				})

#NOTE: Get common HGV for all samples of the given group
HGV <- SelectIntegrationFeatures(object.list = SL, nfeatures = 3000)
# Scale data for those common HGV
SL <- lapply(SL, function(S) {
		     S <- ScaleData(S, features=HGV, verbose=FALSE)
	})

## Merge data
S <- merge(SL[[1]], SL[-1])
VariableFeatures(S) <- HGV

# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

writeMM(GetAssayData(S, assay="SCT", slot="counts"), file=OUTPUT_mtx1)
writeMM(GetAssayData(S, assay="SCT", slot="data"), file=OUTPUT_mtx2)
write.table(GetAssayData(S, assay="SCT", slot="scale.data"), 
	    sep="\t", file=OUTPUT_mtx3, row.names=TRUE, col.names=TRUE, quote=FALSE)
cat(VariableFeatures(S), sep="\n", file=OUTPUT_sel)
write.table(S@meta.data, sep="\t", file=OUTPUT_colData, 
	    col.names=NA, row.names=TRUE, quote=FALSE)
write.table(rowDat[rownames(S), ], sep="\t", file=OUTPUT_rowData, 
	    col.names=TRUE, row.names=FALSE, quote=FALSE)
