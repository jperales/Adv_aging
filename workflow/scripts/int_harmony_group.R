# Integrate group/contrast via harmony
# Javier Perales-Paton (C)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)


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

INPUT_embeddings <- args[10]
INPUT_loadings <- args[11]
INPUT_projected <- args[12]
INPUT_stdev <- args[13]
INPUT_nPCs<- args[14]

# PARAM
GID <- args[15]
REGION <- args[16]
# OUTPUT
OUTPUT_embeddings <- args[17]
OUTPUT_loadings <- args[18]
OUTPUT_projected <- args[19]
OUTPUT_stdev <- args[20]
OUTPUT_nPCs <- args[21]
OUTDIR <- dirname(OUTPUT_embeddings)
 
## Print internal variables
cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("Matrix")
library("Seurat")
library("harmony")
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

## Add PCA dimensionality reduction
emb <- as.matrix(read.table(INPUT_embeddings, sep="\t", header=TRUE))
loa <- as.matrix(read.table(INPUT_loadings, sep="\t", header=TRUE))
# projected <- as.matrix(read.table(INPUT_projected, sep="\t", header=TRUE))
# projected <- new(Class="matrix")
std <- scan(file=INPUT_stdev, what=numeric())
nPCs <- scan(file=INPUT_nPCs, what=numeric())

PCA <- CreateDimReducObject(embeddings=emb, loadings = loa, assay="SCT", stdev= std, key="pca_")
S[["pca"]] <- PCA


## Run Harmony
S <- RunHarmony(S, assay.use="SCT", reduction="pca",
		group.by.vars="Sample",
		dim.use=1:nPCs, plot_convergence=FALSE)

# S <- RunUMAP(S, assay="SCT", reduction="harmony", dims=1:nPCs)
# S <- FindNeighbors(S, assay="SCT", reduction="harmony", dims = 1:nPCs)
# S <- FindClusters(S, resolution=0.5)


harmony <- Reductions(S, "harmony")

# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

write.table(harmony@cell.embeddings, 
	    sep="\t", file=OUTPUT_embeddings, 
	    row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(harmony@feature.loadings, 
	    sep="\t", file=OUTPUT_loadings, 
	    row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(harmony@feature.loadings.projected, 
	    sep="\t", file=OUTPUT_projected, 
	    row.names=TRUE, col.names=TRUE, quote=FALSE)

cat(harmony@stdev, sep="\n", file=OUTPUT_stdev)
cat(nPCs, sep="\n", file=OUTPUT_nPCs)

