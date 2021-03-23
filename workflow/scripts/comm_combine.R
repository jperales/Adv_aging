#' Combine stats from multiple samples CellPhoneDB
#' Javier Perales-Paton (c)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)

# args <- c("out/comm/Groups/young/AA/merged.tsv",
# 	  "out/comm/Groups/young/AA/combined.tsv")

# INPUT
INPUT<- args[1]

# OUTPUT
OUTPUT<- args[2]
OUTDIR <- dirname(OUTPUT)

## Sanity check on inputs
# Same length

## Print internal variables
cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("reshape2")
library("dplyr")
library("metap")

### Load data
IN <- read.table(INPUT, sep="\t", header=TRUE, stringsAsFactors=TRUE)

# Combination of stats
cPval <- function(p) {
	p[p==0] <- 1e-4
	if(length(p) > 1) {
		cp <- sumlog(p)$p
	} else if(length(p) == 1) {
		cp <- p
	}
	return(p)
}

OUT <- IN %>% select(id_cp_interaction, cell_pair, AdjPvalue, Mean, Rank) %>% 
 	group_by(id_cp_interaction, cell_pair) %>%
 	summarise(combinedAdjPval=cPval(AdjPvalue), 
		   madMean=mad(Mean), 
		   madRank=mad(Rank), 
		   nCellPairs=n())
 
# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

write.table(OUT, file=OUTPUT, sep="\t", 
	    row.names=FALSE, col.names=TRUE, quote=FALSE)
