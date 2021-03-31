#' Generate a 'significant means' CPDB output file with combined Group statistics
#' Javier Perales-Paton (c)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)
# args <- c("out/comm/Groups/young/AA/combined.tsv",
# 	"out/comm/Groups/young/AA/legend.tsv",
# 	"out/comm/Groups/young/AA/combined_significant_means.tsv"
#  	)
 
# INPUT
INPUT_COMBINED <- args[1]
INPUT_LEGEND <- args[2]
# OUTPUT
OUTPUT<- args[3]
OUTPUT2 <- args[4]
OUTDIR <- dirname(OUTPUT)

## Print internal variables
cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("reshape2")
library("tidyverse")

### Load data
LEGEND <- read.table(INPUT_LEGEND, sep="\t", header=TRUE, stringsAsFactors=FALSE)
COMBINED <- read.table(INPUT_COMBINED, sep="\t", header=TRUE, stringsAsFactors=FALSE)

#-- Add rank to legend as in significant means
rnk_dic <- setNames(COMBINED$madRank,COMBINED$id_cp_interaction)
LEGEND$rank <- rnk_dic[LEGEND$id_cp_interaction]
stopifnot(ncol(LEGEND)==12)

#-- Gather
means <- COMBINED %>% 
	select(id_cp_interaction, cell_pair, madMean) %>% 
	acast( id_cp_interaction ~ cell_pair)

pvals <- COMBINED %>% 
	select(id_cp_interaction, cell_pair, combinedAdjPval) %>% 
	acast( id_cp_interaction ~ cell_pair)

stopifnot(rownames(means)==rownames(pvals))
stopifnot(colnames(means)==colnames(pvals))

# Add NA to non-significnat interactions
means[pvals < 0.05] <- NA

RES <- data.frame(means, check.names=FALSE)
RES$id_cp_interaction <- rownames(RES)
rownames(RES) <- NULL

# Reformat for CrossTalker
cls <- unique(unlist(sapply(colnames(RES)[-ncol(RES)], function(z) strsplit(z, split="\\|")[[1]], simplify=FALSE)))
cls <- cls[order(nchar(cls), decreasing=TRUE)]
idx2cluster <- data.frame(idx=as.integer(factor(cls)),
			  name=cls)
idx2cluster_dic <- setNames(idx2cluster$idx, idx2cluster$name)
#NOTE: CrossTalkeR does not handle '_' in annotation from a single cluster
idx2cluster$name <- gsub("_","-",idx2cluster$name)

for(cl in cls) {
	colnames(RES) <- gsub(cl, idx2cluster_dic[cl], colnames(RES))
}


OUT <- merge(LEGEND, RES, by.x="id_cp_interaction", by.y="id_cp_interaction", all.x=FALSE, all.y=TRUE)

# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

write.table(OUT, file=OUTPUT, sep="\t", 
	    na = "",
	    row.names=FALSE, col.names=TRUE, quote=FALSE)

write.table(idx2cluster, file=OUTPUT2, sep="\t",
	    row.names=FALSE, col.names=TRUE, quote=FALSE)
