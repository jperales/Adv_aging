#' Run CrossTalker
#' Javier Perales-Paton (c)

## Setting the environment
### Internal variables
set.seed(1234)
args <- commandArgs(TRUE)

# args <- c("out/comm/Groups/old/AA/filtered_corrected.csv,out/comm/Groups/young/AA/filtered_corrected.csv",
# 	  "out/comm/Constrasts/old-young/AA/crosstalker/",
# 	  "out/comm/Constrasts/old-young/AA/crosstalker/report.html")

# INPUT
INPUT<- strsplit(args[1], split=",")[[1]]
# PARAMS
genes <- c("TGFB1", "PF4", "SPP1")
GROUPS <- strsplit(args[2], split="-")[[1]]
SUFFIX <- args[3]

# OUTPUT
OUTDIR <- args[4]

# Get absolute paths
#NOTE: It seems that rendering R-markdowns by crosstalker requires absolute paths
INPUT <- paste0(getwd(), "/", INPUT)
OUTDIR <- paste0(getwd(),"/",OUTDIR)
if(!grepl("/$", OUTDIR)) OUTDIR <- paste0(OUTDIR,"/");

# Add group names
stopifnot(length(INPUT)==length(GROUPS))
names(INPUT) <- GROUPS

## Print internal variables
cat("[INFO] : Print internal variables: \n", file=stdout())
for(i in setdiff(ls(), "args")) {
	cat(paste0(i,"\t=\t", get(i), "\n"), file=stdout())
}


## Load libraries
library("CrossTalkeR")


### Generating the report and the object
# source: https://github.com/CostaLab/CrossTalkeR/blob/master/CellPhoneDB%20Tutorial.md

if(!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive=TRUE);

data <- generate_report(lrpaths=INPUT, # paths list
			genes=genes, # Selected list
			out_path=OUTDIR, # output path
			threshold = 0, # threshold of prune edges 0=keep all
			out_file=SUFFIX # Suffix for report filename)
	)

# Save data
saveRDS(data, paste0(OUTDIR, "data_", SUFFIX,".rds"))
