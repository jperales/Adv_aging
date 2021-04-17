#' Get idx2cluster from 'significant means' CPDB output file for crostalker
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
INPUT_signif<- args[1]
# OUTPUT
OUTPUT<- args[2]
OUTPUT_idx <- args[3]
OUTDIR <- dirname(OUTPUT)

## Print internal variables
cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

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

## Read data
pavg <- read.table(INPUT_signif, sep="\t", header=TRUE,
		 stringsAsFactors = FALSE, check.names = FALSE)

pavg <-RMduplicated.id_cp_interactions(pavg)

# Reformat for CrossTalker
cls <- unique(unlist(sapply(colnames(pavg)[13:ncol(pavg)], function(z) strsplit(z, split="\\|")[[1]], simplify=FALSE)))
cls <- cls[order(nchar(cls), decreasing=TRUE)]
idx2cluster <- data.frame(idx=as.integer(factor(cls)),
			  name=cls)
idx2cluster_dic <- setNames(idx2cluster$idx, idx2cluster$name)
#NOTE: CrossTalkeR does not handle '_' in annotation from a single cluster
idx2cluster$name <- gsub("_","-",idx2cluster$name)

for(cl in cls) {
	colnames(pavg)[13:ncol(pavg)] <- gsub(cl, idx2cluster_dic[cl], colnames(pavg)[13:ncol(pavg)])
}


# Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

write.table(pavg, file=OUTPUT, sep="\t", 
	    na = "",
	    row.names=FALSE, col.names=TRUE, quote=FALSE)

write.table(idx2cluster, file=OUTPUT_idx, sep="\t",
	    row.names=FALSE, col.names=TRUE, quote=FALSE)

