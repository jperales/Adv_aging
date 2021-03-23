#' Perform a simple data processing from original merge, remove empty droplets and bad quality cells
#' Javier Perales-Paton (c)

#Original QC criteria:
# To remove the cells with low quality, cells with gene number over 2000,
# and the ratio of mitochondria lower than 0.15 were maintained, and genes with at
# least one feature count in more than three cells were used for the following analysis

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)
# args <- c(
# 	  # Inputs
# 	  "./data/GSE117715/CellType_info.csv",
# 	  "./data/GSE117715/CellType_abbn.csv",
# 	  "./data/GSE117715/GSE117715_Cynomolgus_monkey_aging_artery_count.txt.gz",
#  	  # Params
# 	  "O.F.1",
# 	  "CA",
#  	  "0.12",
#  	  "3",
#  	  "2000",
#  	  "999999999",
# 	  # Outputs
# 	  "./out/data2/GSE117715/Samples/CA/O.F.1/barcodes.tsv",
# 	  "./out/data2/GSE117715/Samples/CA/O.F.1/features.tsv",
# 	  "./out/data2/GSE117715/Samples/CA/O.F.1/matrix.mtx",
# 	  "./out/QC/GSE117715/Samples/CA/O.F.1/barcodes_QC1.tsv",
# 	  "./out/QC/GSE117715/Samples/CA/O.F.1/features_QC1.tsv"
# 	  )
# 
# INPUT
INPUT_colData <- args[1]
INPUT_colData2 <- args[2]
INPUT_mtx <- args[3]
# PARAMS
ID <- args[4]
REGION <- args[5]
MAX_percMT <- as.numeric(args[6])
MIN_nCells <- as.numeric(args[7])
MIN_nFeatures <- as.numeric(args[8])
MAX_nFeatures <- as.numeric(args[9])
# OUTPUT
OUTPUT_colData <- args[10]
OUTPUT_rowData <- args[11]
OUTPUT_mtx <- args[12]
OUTPUT_colRpt <- args[13]
OUTPUT_rowRpt <- args[14]
OUTDIR2 <- dirname(OUTPUT_mtx)
OUTDIR3 <- dirname(OUTPUT_colRpt)

cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("Matrix")

## Load data

# Col Data
colDat <- read.table(INPUT_colData, sep=",", header=TRUE, stringsAsFactors=FALSE)
colDat$Annotation.Level.1 <- gsub("^(AA_|CA_)", "", colDat$cell_type)

# NOTE: BUG, missing vals
abbn <- read.table(INPUT_colData2, sep=",", header=FALSE, stringsAsFactors=FALSE)
abbn <- setNames(abbn$V1, abbn$V2)
colDat$DESC_Annotation.Level.1 <- abbn[colDat$Annotation.Level.1]
rownames(colDat) <- colDat$CellName

# Matrix
dat <- read.table(INPUT_mtx, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Row Data (retrieve from matrix)
rowDat <- data.frame("hgnc_symbol"=rownames(dat))

## Subset to presence
#NOTE: It seems not all cells in the original matrix are present with annotation,
# meaning that absence was due to filtering criteria in the original paper
# so we apply same criteria accordingly
cellCommon <- intersect(colDat$CellName, colnames(dat))
cat(paste0("[INFO] : Subsetting to common cells in MERGED mtx and colData, n=", length(cellCommon),"\n"),
    file=stdout())
dat <- dat[, cellCommon]
colDat <- colDat[cellCommon, ]

## Match
#NOTE: CellName were named in a pattern:  Region.library_MonkeysInfo_MonkeyID_barcode
cell_idx_bol <- grepl(paste0("^",REGION,"[pn]"), colnames(dat)) &
		 grepl(paste0("_",ID,"_"), colnames(dat))

cat(paste0("[INFO]"," For Sample (id-region)= ", ID,"-",REGION," ... found nCells= ",sum(cell_idx_bol),"\n"),
    file=stdout())
stopifnot(sum(cell_idx_bol)>0)

cellID <- colnames(dat)[cell_idx_bol]

# Subset by Sample-Region, i.e. from merged deposited data to individual libraries
dat <- dat[, cellID]
colDat <- colDat[cellID, ]


#NOTE: Somehow there are not MT expression in the matrix
# # Find MT features
# mt.bool <- grepl("^MT-", rowDat[,1])
# # Verbose
# cat("[INFO] : Show MT-genes found\n", file=stdout())
# print(rowDat[mt.bool,])
# 
## Calculate cell metadata
colDat$nCounts_RNA <- colSums(dat)
colDat$nFeatures_RNA <- colSums(dat!=0)
# colDat$percMT_RNA <- colSums(dat[mt.bool,])/colSums(dat)*100
# if all 0s, then nan. So we reannotate them
# colDat$percMT_RNA[is.nan(colDat$percMT_RNA)] <- 0

# Boolean
colDat$QC1 <- colDat$nFeatures_RNA >= MIN_nFeatures &
		colDat$nFeatures_RNA <= MAX_nFeatures # &
# 		colDat$percMT_RNA <= MAX_percMT

col_passing <- colDat$QC1
cat("[INFO] : Show droplets passing QC\n", file=stdout())
print(table(col_passing))

## Calculate row metadata for passing criteria
rowDat$QC1 <- rowSums(dat[, col_passing]!=0) >= MIN_nCells
row_passing <- rowDat$QC1
cat("[INFO] : Show features passing QC\n", file=stdout())
print(table(row_passing))

## Subset the data by QC passing criteria
colDat2 <- colDat[col_passing,,drop=FALSE]
rowDat2 <- rowDat[row_passing, ]
dat2 <- dat[row_passing, col_passing]

dat2 <- as(as.matrix(dat2), Class = "dgCMatrix")

## Write new data
if(!dir.exists(OUTDIR2)) 
	dir.create(OUTDIR2, recursive=TRUE)

writeMM(dat2, file=OUTPUT_mtx)
write.table(colDat2, sep="\t", file=OUTPUT_colData, 
	    col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(rowDat2, sep="\t", file=OUTPUT_rowData, 
	    col.names=TRUE, row.names=FALSE, quote=FALSE)

if(!dir.exists(OUTDIR3)) 
	dir.create(OUTDIR3, recursive=TRUE)

write.table(colDat, sep="\t", file=OUTPUT_colRpt, 
	    col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(rowDat, sep="\t", file=OUTPUT_rowRpt, 
	    col.names=TRUE, row.names=FALSE, quote=FALSE)

