#!/usr/bin/Rscript

library(methylKit)
library(ggplot2)

# source the script holding the helper functions
source('scripts/DMR_source_functions.R')

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
cov_files <- unlist(strsplit(args[1], ","))
samp_names <- unlist(strsplit(args[2], ","))
my_treat <- as.numeric(unlist(strsplit(args[3], ",")))
norm <- unlist(strsplit(args[4], ","))
norm_method <- unlist(strsplit(args[5], ","))
comparison_samples <- unlist(strsplit(args[6], ","))
comparison_treat <- as.numeric(unlist(strsplit(args[7], ",")))
min_per_group <- as.integer(unlist(strsplit(args[8], ",")))
win_size <- as.numeric(unlist(strsplit(args[9], ",")))
step_size <- as.numeric(unlist(strsplit(args[10], ",")))
covar1 <- unlist(strsplit(args[11], ","))
covar2 <- unlist(strsplit(args[12], ","))
overdisp <- unlist(strsplit(args[13], ","))
test <- unlist(strsplit(args[14], ","))
comp_name <- unlist(strsplit(args[15], ","))
meth_diff <- as.numeric(unlist(strsplit(args[16], ",")))
qval <- as.numeric(unlist(strsplit(args[17], ",")))


# Run first function to create a methylbase object from the cov files
myObj <- create.myObj(cov_files)

print("cov_files read in successfully")

# create and setwd to DMR_analysis directory
dir.create(file.path(getwd(), 'data/DMR_analysis'), showWarnings = FALSE)
setwd(file.path(getwd(), 'data/DMR_analysis'))

# filter and normalize
# defaults are 10x 99.9, normalized by median
myObj.filtered <- filter(myObj, low=10, high=99.9, normalize=norm, norm.method=norm_method)

# reorganize myObj to get samples belonging to the desired 2 groups for a comparison
my.meth <- reorg(myObj.filtered,
        samples=comparison_samples,
        treat=comparison_treat,
	win.size=win_size,
	step.size=step_size,
	min.per.group=min_per_group
)

# construct the covariates dataframe
covariates <- make_covar(covar1, covar2)

# Calculate DMRs
# exports:
#   1. Sig DMR results txt file
#   2. Sig DMR results bed file (ensembl contigs)
#   3. Hyper vs Hypo DMR pie chart
#   4. Both GREAT bed file
#   5. Hyper GREAT bed file
#   6. Hypo GREAT bed file
res <- calc.DMRs(my.meth,
                   covariate=covariates,
                   overdispersion=overdisp,
                   test=test,
                   comparison=comp_name,
                   meth.diff=meth_diff,
                   qval=qval
)

# make sig DMR bed file track
makeBED(res, comparison=comp_name)

foo <- "DMR Analysis Complete"
write.table(foo, "complete.txt", sep="\t", col.names=F, row.names=F, quote=F)
