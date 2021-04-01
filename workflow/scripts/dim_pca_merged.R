# Merge normalized matrices for a given group/Contrast
# Javier Perales-Paton (C)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)

# args <- c(
# 	  "out/data2/Groups/young/AA/barcodes.tsv",
# 	  "out/data2/Groups/young/AA/features.tsv",
# 	  "out/data2/Groups/young/AA/matrix.mtx",
# 	  "out/norm/Groups/young/AA/logNormSCT_barcodes.tsv",
# 	"out/norm/Groups/young/AA/logNormSCT_features.tsv",
# 	"out/norm/Groups/young/AA/logNormSCT_counts.mtx",
# 	"out/norm/Groups/young/AA/logNormSCT_data.mtx",
# 	"out/norm/Groups/young/AA/logNormSCT_scaled.tsv",
# 	"out/norm/Groups/young/AA/logNormSCT_HGV.txt",
# 	"young",
# 	"AA",
# 	"out/dim/Groups/young/AA/logNormSCT_PCA_embeddings.tsv",
# 	"out/dim/Groups/young/AA/logNormSCT_PCA_loadings.tsv",
# 	"out/dim/Groups/young/AA/logNormSCT_PCA_stdev.txt"
# 	)
#  

# INPUT
INPUT_umi_colData <- args[1]
INPUT_umi_rowData <- args[2]
INPUT_umi_mtx <- args[3]

INPUT_norm_colData <- args[4]
INPUT_norm_rowData <- args[5]
INPUT_norm_cnt <- args[6]
INPUT_norm_data<- args[7]
INPUT_norm_scaled<- args[8]
INPUT_norm_hgv <- args[9]

# PARAM
GID <- args[10]
REGION <- args[11]
# OUTPUT
OUTPUT_embeddings <- args[12]
OUTPUT_loadings <- args[13]
OUTPUT_projected<- args[14]
OUTPUT_stdev <- args[15]
OUTPUT_nPCs<- args[16]
OUTPUT_horn<- args[17]
OUTDIR <- dirname(OUTPUT_embeddings)
 
## Print internal variables
cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), c("args"))) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("Matrix")
library("Seurat")
library("PCAtools")
source("workflow/src/10Xmat.R")

## Load data
# umi_colDat <- readBarcodes(INPUT_umi_colData, header=TRUE)
# umi_rowDat <- readFeatures(INPUT_umi_rowData, header=TRUE)
# rownames(umi_rowDat) <- umi_rowDat[,1]
umi <- readMatrix(INPUT_umi_mtx, INPUT_umi_rowData, INPUT_umi_colData, header=TRUE)
sct_cnt <- readMatrix(INPUT_norm_cnt, INPUT_norm_rowData, INPUT_norm_colData, header=TRUE)
sct_dat <- readMatrix(INPUT_norm_data, INPUT_norm_rowData, INPUT_norm_colData, header=TRUE)
sct_scal <- as.matrix(read.table(INPUT_norm_scaled, sep="\t", header=TRUE, stringsAsFactors=FALSE))
sct_hgv <- scan(INPUT_norm_hgv, what="character") 

## Build Seurat Object
S <- Seurat::CreateSeuratObject(counts = umi,
				project = GID)

## Add assay SCT as originally built
SCT <- CreateAssayObject(data = sct_dat)
SCT@counts <- sct_cnt
SCT@var.features <- sct_hgv
SCT@scale.data <- sct_scal
Seurat::Key(SCT) <- "SCT_"
S[["SCT"]] <- SCT

DefaultAssay(S) <- "SCT" 

colDat <- readBarcodes(INPUT_norm_colData, header=TRUE)

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


S <- RunPCA(S, assay="SCT", npcs=50, verbose=FALSE) 
pca <- Reductions(S, "pca")

### Choose number of PCs for downstream analysis
# Horn’s parallel analysis is commonly used to pick the number of PCs
set.seed(100010)
horn <- PCAtools::parallelPCA(as.matrix(S@assays$SCT@data[VariableFeatures(S),]),
    BSPARAM=BiocSingular::IrlbaParam(), niters=10)


# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

write.table(pca@cell.embeddings, 
	    sep="\t", file=OUTPUT_embeddings, 
	    row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(pca@feature.loadings, 
	    sep="\t", file=OUTPUT_loadings, 
	    row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(pca@feature.loadings.projected, 
	    sep="\t", file=OUTPUT_projected, 
	    row.names=TRUE, col.names=TRUE, quote=FALSE)
cat(pca@stdev, sep="\n", file=OUTPUT_stdev)
cat(horn$n, sep="\n", file=OUTPUT_nPCs)

pdf(OUTPUT_horn,  width=6, height=4, compress=FALSE)
plot(horn$original$variance, type="b", log="y", pch=16)
permuted <- horn$permuted
for (i in seq_len(ncol(permuted))) {
    points(permuted[,i], col="grey80", pch=16)
    lines(permuted[,i], col="grey80", pch=16)
}
abline(v=horn$n, col="red")
dev.off()


