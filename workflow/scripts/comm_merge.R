#' Merge CellPhoneDB output from multiple-samples in a tidy format
#' Javier Perales-Paton (c)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)
# args <- c(paste(grep("Y.*CA",list.files("out/comm/Samples", pattern = "^pvalues.txt", recursive = TRUE, full.names = TRUE),value=TRUE), collapse=","),
#  	paste(grep("Y.*CA", list.files("out/comm/Samples", pattern = "^means.txt", recursive = TRUE, full.names = TRUE), value=TRUE), collapse=","),
#  	paste(grep("Y.*CA", list.files("out/comm/Samples", pattern = "^significant_means.txt", recursive = TRUE, full.names = TRUE), value=TRUE), collapse=","),
#  	"out/comm/Groups/Y/merged.tsv",
# 	"out/comm/Groups/Y/legend.tsv"
#  	)
#  

# INPUT
INPUT_pvalues<- strsplit(args[1], split=",")[[1]]
INPUT_means <- strsplit(args[2], split=",")[[1]]
INPUT_signif <- strsplit(args[3], split=",")[[1]]
# PARAMS
# ID <- args[4]
# OUTPUT
OUTPUT<- args[4]
LEGEND <- args[5]
OUTDIR <- dirname(OUTPUT)

## Sanity check on inputs
# Same length
stopifnot(all(length(INPUT_pvalues)==length(INPUT_means) & length(INPUT_pvalues)== length(INPUT_signif)))

# Same order of every sample
for (i in 1:length(INPUT_pvalues)) {
	sidx1 <- basename(dirname(dirname(INPUT_pvalues[i])))	
	sidx2 <- basename(dirname(dirname(INPUT_means[i])))	
	sidx3 <- basename(dirname(dirname(INPUT_signif[i])))	

	stopifnot(all(sidx1==sidx2 && sidx1==sidx3))
}

## Print internal variables
cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("reshape2")

###------Functions
RMduplicated.id_cp_interactions <- function(df) {
	if(anyDuplicated(df$id_cp_interaction)) {
		dups <- df$id_cp_interaction[duplicated(df$id_cp_interaction)]
		cat(paste0("WARN: Duplicated protein interaction pairs:", length(dups), ". Retreiving 1/each\n"),
		    file=stdout())
		RMdups_idx <- setNames(vector("list", length(dups)), dups)
		# Collect idx to discard
		for(d in dups) {
			tmp <- df[df$id_cp_interaction %in% d, ]
			tmp[is.na(tmp)] <- ""
			sel <- names(which.max(rowSums(tmp!=""))[1])
			RMdups_idx[[d]] <- base::setdiff(rownames(tmp), sel)
		}
		df <- df[!rownames(df) %in% unlist(RMdups_idx), ]
	}
	return(df)

}

### Process data for merging

# Get all Sample IDs (SIDX)
SIDX <- sapply(INPUT_pvalues, function(fl) basename(dirname(dirname(fl))), USE.NAMES=FALSE, simplify=TRUE)

# Rename names acoordingly
names(INPUT_pvalues) <- sapply(INPUT_pvalues, function(fl) basename(dirname(dirname(fl))), USE.NAMES=FALSE, simplify=TRUE)
names(INPUT_means) <- sapply(INPUT_means, function(fl) basename(dirname(dirname(fl))), USE.NAMES=FALSE, simplify=TRUE)
names(INPUT_signif) <- sapply(INPUT_signif, function(fl) basename(dirname(dirname(fl))), USE.NAMES=FALSE, simplify=TRUE)


cpdb_tidy <- setNames(vector("list", length(SIDX)), SIDX)
cpdb_leg <- setNames(vector("list", length(SIDX)), SIDX)

for(sidx in SIDX) {

	# Read cpdb output files for single sample
	pval <- read.table(INPUT_pvalues[sidx], sep="\t", header=TRUE,
			   stringsAsFactors = FALSE, check.names = FALSE)
	avg <- read.table(INPUT_means[sidx], sep="\t", header=TRUE, 
			  stringsAsFactors = FALSE, check.names = FALSE)
	pavg <- read.table(INPUT_signif[sidx], sep="\t", header=TRUE,
			   stringsAsFactors = FALSE, check.names = FALSE)

	#--- Reformat 1: rows
	#NOTE: issue , sometimes there are id_cp_interaction duplicated.
	# It seems the ypresent exact same statistics and means, but one of
	# the duplicated forms for each are a bit more completed annotated,
	# so we retrive those

	pval <-RMduplicated.id_cp_interactions(pval)
	avg <-RMduplicated.id_cp_interactions(avg)
	pavg <-RMduplicated.id_cp_interactions(pavg)

	rownames(pval) <- pval$id_cp_interaction
	rownames(avg) <- avg$id_cp_interaction
	rownames(pavg) <- pavg$id_cp_interaction

	# Sanity check1: same order of rows
	stopifnot(all(rownames(pval) == rownames(avg)))
	stopifnot(all(rownames(pval) %in% rownames(pavg)))
	if(!all(rownames(pavg)==rownames(pval))) {
		cat(paste0("WARN : Reordering 'significant means'\n"), file=stdout())
		pavg <- pavg[rownames(pval),]
	}


	# Sanity check2: same order of columns (cell pairs)
	stopifnot(all(rownames(pval) == rownames(pavg)))
	stopifnot(all(rownames(avg) == rownames(avg)))
	# NOTE: pavg contains an extra column with the rankings of specificity
	stopifnot(colnames(pval)[1:11]==colnames(pavg)[1:11])

	#--- Reformat 2 : cols with stats
	leg <- pavg[, 1:11]
	rnk <- pavg[, "rank", drop=FALSE]
	pval <- pval[, 12:ncol(pval)]
	avg <- avg[, 12:ncol(avg)]


	#-- Reformat 3: tidy
	pval2 <- reshape2::melt(as.matrix(pval))
	colnames(pval2) <- c("id_cp_interaction", "cell_pair", "Pvalue")
	#NOTE: Multiple testing correction BH
	pval2$AdjPvalue <- p.adjust(pval2$Pvalue, method="fdr")
	pval2$key <- paste0(pval2$id_cp_interaction, "::", pval2$cell_pair)

	avg2 <- reshape2::melt(as.matrix(avg))
	colnames(avg2) <- c("id_cp_interaction", "cell_pair", "Mean")
	avg2$key <- paste0(avg2$id_cp_interaction, "::", avg2$cell_pair)

	rnk2 <- rnk
	colnames(rnk2) <- "Rank"
	rnk2$id_cp_interaction <- rownames(rnk)

	res <- merge(pval2, avg2[, 3:4], by= c("key", "key"), all=TRUE)
	res <- merge(res, rnk2, by = c("id_cp_interaction", "id_cp_interaction"), all=TRUE)

	res <- res[, colnames(res)!="key"]

	out_cols <- colnames(res)
	res$Sample <- sidx
	res <- res[, c("Sample", out_cols)]

	cpdb_leg[[sidx]] <- leg
	cpdb_tidy[[sidx]] <- res

}


# Collection of cp_id_interaction legends
LEG <- do.call("rbind", cpdb_leg)
rownames(LEG) <- NULL
LEG <- unique(LEG)

# Tidy results
OUT <- do.call("rbind", cpdb_tidy)
rownames(OUT) <- NULL

# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

write.table(LEG, file=LEGEND, sep="\t", 
	    row.names=FALSE, col.names=TRUE, quote=FALSE)

write.table(OUT, file=OUTPUT, sep="\t", 
	    row.names=FALSE, col.names=TRUE, quote=FALSE)

