# Merge UMI matrices for a given group/Contrast
# Javier Perales-Paton (C)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)

# args <- c(paste(grep("Y.*CA",list.files("out/data2/Samples", pattern = "^barcodes.tsv", recursive = TRUE, full.names = TRUE),value=TRUE), collapse=","),
# 	paste(grep("Y.*CA", list.files("out/data2/Samples", pattern = "^features.tsv", recursive = TRUE, full.names = TRUE), value=TRUE), collapse=","),
#  	paste(grep("Y.*CA", list.files("out/data2/Samples", pattern = "^matrix.mtx", recursive = TRUE, full.names = TRUE), value=TRUE), collapse=","),
# 	"young",
# 	"CA",
#  	"out/data2/Groups/Y/barcodes.tsv",
#  	"out/data2/Groups/Y/features.tsv",
#  	"out/data2/Groups/Y/matrix.mtx"
#  	)
# 


# INPUT
INPUT_colData <- strsplit(args[1], split=",")[[1]]
INPUT_rowData <- strsplit(args[2], split=",")[[1]]
INPUT_mtx <- strsplit(args[3], split=",")[[1]]
# PARAM
GID <- args[4]
REGION <- args[5]
# OUTPUT
OUTPUT_colData <- args[6]
OUTPUT_rowData <- args[7]
OUTPUT_mtx <- args[8]
OUTDIR <- dirname(OUTPUT_mtx)

## Sanity check on inputs
# Same length
stopifnot(all(length(INPUT_rowData)==length(INPUT_colData) & length(INPUT_rowData)== length(INPUT_mtx)))

# Same order of every sample
for (i in 1:length(INPUT_rowData)) {
	sidx1 <- basename(dirname(dirname(INPUT_rowData[i])))	
	sidx2 <- basename(dirname(dirname(INPUT_colData[i])))	
	sidx3 <- basename(dirname(dirname(INPUT_mtx[i])))	

	stopifnot(all(sidx1==sidx2 && sidx1==sidx3))
}

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
colDat_list <- sapply(INPUT_colData, function(fl) readBarcodes(fl, header=TRUE), simplify=FALSE)
rowDat_list <- sapply(INPUT_rowData, function(fl) readFeatures(fl, header=TRUE), simplify=FALSE)

dat_list <- vector("list", length=length(INPUT_mtx))
names(dat_list) <- sapply(INPUT_mtx, function(fl) basename(dirname(dirname(fl))), USE.NAMES=FALSE, simplify=TRUE)
for(i in 1:length(INPUT_mtx)) {
	dat_list[[i]] <- readMatrix(INPUT_mtx[i], INPUT_rowData[i], INPUT_colData[i], header=TRUE)
}

## Merge data
rowDat <- unique(do.call("rbind", rowDat_list))
rownames(rowDat) <- rowDat[,1]

#NOTE: We use Seurat C++ implementation of merging data for a better performance
# Creating Seurat objects for each sample

SL <- setNames(vector("list", length=length(dat_list)),
	       names(dat_list))
for(i in 1:length(dat_list)) {
	ID <- names(dat_list)[i]
	# Create SeuratObject
	S <- Seurat::CreateSeuratObject(counts = dat_list[[i]],
					project = ID)
	# Add metadata
	colDat <- colDat_list[[i]]
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
	S$orig.ident <- ID
	if(!grepl("-",GID)) S$Group <- GID;
	if(!grepl("-",GID)) S$Sample <- ID;

	stopifnot(all(S$Region == REGION))

	# Add it
	SL[[ID]] <- S
	rm(S, colDat, colDat2)
}

S <- merge(SL[[1]], SL[-1])

# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

writeMM(GetAssayData(S, assay="RNA", slot="counts"), file=OUTPUT_mtx)
write.table(S@meta.data, sep="\t", file=OUTPUT_colData, 
	    col.names=NA, row.names=TRUE, quote=FALSE)
write.table(rowDat[rownames(S), ], sep="\t", file=OUTPUT_rowData, 
	    col.names=TRUE, row.names=FALSE, quote=FALSE)

