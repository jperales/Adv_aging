#' Prepare single sample data for cell communication analysis
#' Javier Perales-Paton (c)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)
# args <- c(
# 	  # Inputs
# 	  "./out/norm/Samples/O.F.1/CA/logNormSCT_barcodes.tsv",
# 	  "./out/norm/Samples/O.F.1/CA/logNormSCT_features.tsv",
# 	  "./out/norm/Samples/O.F.1/CA/logNormSCT_matrix.mtx",
#  	  # Params
# 	  "O.F.1",
#  	  "Annotation.Level.1",
# 	  # Outputs
# 	  "./out/comm/Samples/O.F.1/CA/cellphonedb_meta.txt",
# 	  "./out/comm/Samples/O.F.1/CA/cellphonedb_count.txt"
# 	  )
 
# INPUT
INPUT_colData <- args[1]
INPUT_rowData <- args[2]
INPUT_mtx <- args[3]
# PARAMS
ID <- args[4]
ANN <- args[5]
# OUTPUT
OUTPUT_meta <- args[6]
OUTPUT_mtx <- args[7]
OUTDIR <- dirname(OUTPUT_mtx)

cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}


### Load libraries
library("Matrix")
source("workflow/src/10Xmat.R")
#

## Load data
colDat <- readBarcodes(INPUT_colData, header=TRUE)
dat <- readMatrix(INPUT_mtx, INPUT_rowData, INPUT_colData, header=TRUE)
mat <- as.matrix(dat)

## Sanity check
meta_idx <- colnames(colDat)[which(colnames(colDat)==ANN)]
if(length(meta_idx)!=1)
	stop("[ERROR]: Cluster annotation level does not found\n")

cellID_idx <- colnames(colDat)[1]
cat(paste0("[INFO] : Showing columns used for CellPhoneDB for cell metadata:\n",
	   "\t Barcode: ",cellID_idx,"\n",
	   "\t Cluster: ", meta_idx,"\n",
	   "---\n"),
    file=stdout())

## Drop EPI
idx <- which(colDat[, meta_idx]!= "EPI")
colDat <- colDat[idx, ]
mat <- mat[, colDat[, cellID_idx]]

## Write new data
# Format, source: https://www.cellphonedb.org/faq-and-troubleshooting#extract-cellphonedb-files-from-seurat 
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

write.table(mat, OUTPUT_mtx, sep="\t", quote=FALSE)
write.table(colDat[, c(cellID_idx, meta_idx)], 
	OUTPUT_meta, sep="\t", quote=FALSE, row.names=FALSE)

