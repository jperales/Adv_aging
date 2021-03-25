# Differential gene expression from cell-type pseudobulk profiles
# Javier Perales-Paton (C)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)

# args <- c("out/minibulk/Contrasts/old-young/AA/AF_counts.tsv",
# 	  "out/minibulk/Contrasts/old-young/AA/AF_samples.tsv",
#  	  "old-young",
#  	  "out/minidge/Contrasts/old-young/AA/AF_topTags.tsv")
 
# INPUT
INPUT_cnt <- args[1]
INPUT_samples <- args[2]
comparison <- args[3]
OUTPUT <- args[4]

# PARAMS
OUTDIR <- dirname(OUTPUT)

cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}

### Load libraries
library("edgeR")

### Load data
cnt <- read.table(INPUT_cnt, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cnt <- apply(as.matrix(cnt), c(1,2), as.integer)

samples <- read.table(INPUT_samples, sep="\t", header=TRUE, stringsAsFactors=FALSE)

grs <- strsplit(comparison, split="-")[[1]]
samples$Group <- factor(samples$Group, levels=rev(grs))

#NOTE: At this point, data collection for the comparison has been straightforward
# But it might not present enough data for a proper statistical comparison
# Thus it is important to add some sanity checks to avoid crashes during a 
# 2-group comparison

#NOTE: This script can only go out via a data.frame, so if no results are generated
# an empty data.frame to be written. If statistical test is performed, then a full
# data.frame with genome-wide statistics is written
topT <- c("NoTest")

# Sanity check
if(length(unique(samples$Group))==2) {
	### Creating DGEList
	y <- DGEList(cnt, samples=samples, group=samples$Group)

	### Preprocessing
	discarded <- samples$ncells < 10
	cat(paste0("[INFO] Samples discarded :\n"), file=stdout())
	table(discarded)
	y <- y[,!discarded]

	#NOTE: At least 3 biological replicates per condition
	if(all(table(y$samples$Group) >= 3)) {
		keep <- filterByExpr(y)
		y <- y[keep,]
		cat(paste0("[INFO] Expressed genes kept :\n"), file=stdout())
		table(keep)
		stopifnot(length(keep)>0)

		### Normalization
		y <- calcNormFactors(y)

		### Statistical modelling
		design <- model.matrix(~ group, y$samples)

		y <- estimateDisp(y, design)
		summary(y$trended.dispersion)

		fit <- glmQLFit(y, design, robust=TRUE)
		res <- glmQLFTest(fit, coef=ncol(design))
		cat(paste0("[INFO] Differentially expressed genes (FDR 5%):\n"), file=stdout())
		summary(decideTests(res))

		topT <- topTags(res, n=Inf)
	} else {
		cat("[INFO] : No statistical modelling due to lack of samples:\n", file=stdout())
		cat(paste0("\t N per Group=",table(y$samples$Group),"\n"), file=stdout())
	}
}

### Save data
if(!dir.exists(OUTDIR)) 
	dir.create(OUTDIR, recursive=TRUE)

if(is(topT)[1]=="TopTags")
	write.table(topT, file=OUTPUT, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

if(is.vector(topT))
	cat(topT,sep="\n", file=OUTPUT)

