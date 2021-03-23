# Read 10x cellranger count matrix
# Javier Perales-Paton (C)
# Date: 08.03.2021
# Source: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices 
# Dep: library(Matrix)

readFeatures <- function(feature.path, header=FALSE) {
	feature.names <- read.delim(feature.path,
				   header = header,
				   stringsAsFactors = FALSE)
	if(!header)
		colnames(feature.names)[1:2] <- c("ensembl_gene_id", "hgnc_symbol")

	return(feature.names)
}

readBarcodes <- function(barcode.path, header=FALSE) {
	barcode.names = read.delim(barcode.path,
				   header = header,
				   stringsAsFactors = FALSE)

	if(!header)
		colnames(barcode.names)[1] <- c("barcode")

	return(barcode.names)

}

readMatrix <- function(matrix.path, feature.path, barcode.path, header=FALSE, feature.loc=1, barcode.loc=1) {
	# Read sparse matrix
	mat <- Matrix::readMM(file = matrix.path)
	# Read dimension of sparse matrix
	barcode.names <- readBarcodes(barcode.path, header)
	feature.names <- readFeatures(feature.path, header)
	# Transform to dgCMatrix, default was dgTMatrix which leads 
	# issues with high mem demands when conver to 'matrix'
	mat <- as(mat, Class = "dgCMatrix")
	# Rename dim
	stopifnot(ncol(mat)==nrow(barcode.names))
	colnames(mat) <- barcode.names[,barcode.loc]
	stopifnot(nrow(mat)==nrow(feature.names))
	rownames(mat) <- feature.names[,feature.loc]
	# Return
	return(mat)
}


