# Export SeuratObject 
# Javier Perales-Paton (C)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)


# args <- c(
# 	  "out/data2/Contrasts/old-young/CA/barcodes.tsv",
# 	  "out/data2/Contrasts/old-young/CA/features.tsv",
# 	  "out/data2/Contrasts/old-young/CA/matrix.mtx",
# 	# 	  "out/norm/Contrasts/old-young/CA/logNormSCT_barcodes.tsv",
#  	  "out/ann/Contrasts/old-young/CA/logNormSCT_harmony_barcodes.tsv",
# 	"out/norm/Contrasts/old-young/CA/logNormSCT_features.tsv",
# 	"out/norm/Contrasts/old-young/CA/logNormSCT_counts.mtx",
# 	"out/norm/Contrasts/old-young/CA/logNormSCT_data.mtx",
# 	"out/norm/Contrasts/old-young/CA/logNormSCT_scaled.tsv",
# 	"out/norm/Contrasts/old-young/CA/logNormSCT_HGV.txt",
# 	"out/dim/Contrasts/old-young/CA/logNormSCT_harmony_embeddings.tsv",
# 	"out/dim/Contrasts/old-young/CA/logNormSCT_harmony_loadings.tsv",
# 	"out/dim/Contrasts/old-young/CA/logNormSCT_harmony_projected.tsv",
# 	"out/dim/Contrasts/old-young/CA/logNormSCT_harmony_stdev.txt",
# 	"out/dim/Contrasts/old-young/CA/logNormSCT_harmony_nPCs.txt",
# 	"out/clust/Contrasts/old-young/CA/logNormSCT_harmony_GraphNN.rds",
# 	"out/clust/Contrasts/old-young/CA/logNormSCT_harmony_GraphSNN.rds",
# 	"old-young",
# 	"CA",
# 	"harmony",
# 	"out/export/SeuratObject_CA.rds"
# 	)
 
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

INPUT_gNN <- args[15]
INPUT_gSNN <- args[16]

# PARAM
GID <- args[17]
REGION <- args[18]
RED <- args[19]
# OUTPUT
#UNK
OUTPUT <- args[20]
OUTDIR <- dirname(OUTPUT)
 
## Print internal variables
cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("Matrix")
library("Seurat")
library("cowplot")
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

# Remove dummy metadata cols
S@meta.data <- S@meta.data[, grep("^X(\\.[1-9]+)?", colnames(S@meta.data), invert=TRUE)]

## Add PCA dimensionality reduction
emb <- as.matrix(read.table(INPUT_embeddings, sep="\t", header=TRUE))
loa <- as.matrix(read.table(INPUT_loadings, sep="\t", header=TRUE))
# projected <- as.matrix(read.table(INPUT_projected, sep="\t", header=TRUE))
# projected <- new(Class="matrix")
std <- scan(file=INPUT_stdev, what=numeric())
nPCs <- scan(file=INPUT_nPCs, what=numeric())

red <- CreateDimReducObject(embeddings=emb, loadings = loa, assay="SCT", stdev= std, key=paste0(RED,"_"))
S[[RED]] <- red

## Add Graph
NN <- readRDS(INPUT_gNN)
SNN <- readRDS(INPUT_gSNN)
S[[paste0(DefaultAssay(S),"_nn")]] <- NN
S[[paste0(DefaultAssay(S),"_snn")]] <- SNN

if(!all(Idents(S) == S$Annotation.Level.2)) Idents(S) <- S$Annotation.Level.2;

## Run UMAP 
S <- RunUMAP(S, assay="SCT", reduction="harmony", dims=1:nPCs)

# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

saveRDS(S, file=OUTPUT)

