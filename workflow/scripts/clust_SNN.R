# Unsupervised Clustering
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
# 	"out/dim/Groups/young/AA/logNormSCT_harmony_embeddings.tsv",
# 	"out/dim/Groups/young/AA/logNormSCT_harmony_loadings.tsv",
# 	"out/dim/Groups/young/AA/logNormSCT_harmony_projected.tsv",
# 	"out/dim/Groups/young/AA/logNormSCT_harmony_stdev.txt",
# 	"out/dim/Groups/young/AA/logNormSCT_harmony_nPCs.txt",
# 	"young",
# 	"AA",
# 	"0.5,1.0,1.5",
# 	"harmony",
# 	"out/clust/Groups/young/AA/logNormSCT_harmony_GraphNN.rds",
# 	"out/clust/Groups/young/AA/logNormSCT_harmony_GraphSNN.rds",
# 	"out/clust/Groups/young/AA/logNormSCT_harmony_Idents"
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

# PARAM
GID <- args[15]
REGION <- args[16]
RES <- as.numeric(strsplit(args[17], split=",")[[1]])
RED <- args[18]
# OUTPUT
#UNK
OUTPUT_nn <- args[19]
OUTPUT_snn <- args[20]
OUTPUT_idents <- args[21]
OUTDIR <- dirname(OUTPUT_nn)
 
## Print internal variables
cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("Matrix")
library("Seurat")
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

red <- CreateDimReducObject(embeddings=emb, loadings = loa, assay="SCT", stdev= std, key=paste0(RED,"_"))
S[[RED]] <- red


## Unsupervised clustering
S <- FindNeighbors(S, assay="SCT", reduction="harmony", dims = 1:nPCs)
S <- FindClusters(S, resolution=RES)
#S <- FindSubCluster(S, cluster=5, graph.name="SCT_snn")
# S <- RunUMAP(S, assay="SCT", reduction="harmony", dims=1:nPCs)

# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

saveRDS(Graphs(S, "SCT_nn"), file=OUTPUT_nn)
saveRDS(Graphs(S, "SCT_snn"), file=OUTPUT_snn)

tab <- S@meta.data[,grep("snn_res\\.", colnames(S@meta.data))]
colnames(tab) <- paste0(GID, "_Annotation.Level.",1:ncol(tab)) 
write.table(tab,
	file=OUTPUT_idents, sep="\t", row.names=TRUE, col.names=TRUE,
	quote=FALSE)	

