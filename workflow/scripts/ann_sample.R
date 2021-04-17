# Propagate cluster Annotation from integrated to lower levels
# Javier Perales-Paton (C)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)

# args <- c("out/norm/Samples/O.F.1/AA/logNormSCT_barcodes.tsv",
# 	  "out/ann/Groups/old/AA/logNormSCT_harmony_barcodes.tsv",
# 	  "out/ann/Samples/O.F.1/AA/logNormSCT_barcodes.tsv")

# INPUT
INPUT_LOW_colData <- args[1]
INPUT_HIGH_colData <- args[2]
# OUTPUT
OUTPUT_LOW_colData <- args[3]

OUTDIR <- dirname(OUTPUT_LOW_colData)

cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("Matrix")
source("workflow/src/10Xmat.R")

## Load data
low_colDat <- readBarcodes(INPUT_LOW_colData, header=TRUE)
high_colDat <- readBarcodes(INPUT_HIGH_colData, header=TRUE)
stopifnot(all(low_colDat[,1] %in% high_colDat[,1]))

## Subset
rownames(high_colDat) <- high_colDat[,1]
new_low_colDat <- high_colDat[low_colDat[,1],]

# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

write.table(new_low_colDat, sep="\t", file=OUTPUT_LOW_colData, 
	    col.names=TRUE, row.names=FALSE, quote=FALSE)

