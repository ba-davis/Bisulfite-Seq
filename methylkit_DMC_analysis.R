# Differentially Methylated CpG Analysis with methylKit #

library(methylKit)

#---------------------------------------#
# Variables to change for each project: #
#---------------------------------------#

# path to coverage files
dir <- "/home/groups/hoolock/u1/bd/Projects/ECP28/meth_extract/relaxed/cov_files"

# set the desired order of the cov.files (group samples together)
# example: cov.order <- c(19,15,23,9,16,1,21,20,17,18,24,22,5,6,11,3,12,14,7,4,13,2,10,8)
cov.order <- c(1,4,5,6,7,8,9,10,11,2,3)

# list the sample names (same order as desired cov.order)
samp.names <- c("MOR1","MOR2","MOR3","MOR4","MOR5","FOR1","FOR2","FOR3","FOR4","EWE1","EWE2")

# create treatment vector for diff analysis (can have more than 2 groups, we will reorganize for comparisons)
my.treat <- c(1,1,1,1,1,2,2,2,2,3,3)

# create a covariate data frame for differential analysis (may need more than 1)
# ex: ovar1 <- data.frame(sex=c("M","M","F","M","F","F","M"))
covar1 <- data.frame()

# canonical contigs
# ex: sheep: cans <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","X","Y")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# FUNCTION to create the methyl object from the cov files
create.myObj <- function() {
  # define path to coverage files
  cov.files <- list.files(dir, 'cov.gz$', full=T)
  
  # reorder the cov.files
  cov.files <- cov.files[cov.order]
  
  # capture the cov.file paths
  my.cov.files <- lapply(cov.files, function(x) x)
  # capture the samp.names
  my.samp.names <- lapply(samp.names, function(x) x)
  
  # use methRead function from methylKit to read in the cov.files
  myObj <- methRead(my.cov.files,
                    sample.id=my.samp.names,
                    assembly="genome",
                    treatment=my.treat,
                    pipeline="bismarkCoverage",
                    context="CpG",
                    mincov=1)
  return(myObj)
}

# FUNCTION to filter and normalize myObj
filter <- function(myObj, low=10, high=99.9, norm.method="mean") {
  # filter myObj CpGs
  myObj.filtered <- filterByCoverage(myObj,
                                     lo.count=low,
                                     lo.perc=NULL,
                                     hi.count=NULL,
                                     hi.perc=high
  )
  
  # normalize by coverage
  myObj.filt.norm <- normalizeCoverage(myObj.filtered,method=norm.method)

  return(myObj.filt.norm)  
}

# FUNCTION to reorganize myObj and merge CpGs
reorg <- function(myObj, samples, treat, min.per.group) {
  # reorganize the desired object
  myObj.reorg <- reorganize(myObj, 
                            sample.ids=samples,
                            treatment=treat
  )
  
  # merge the reorganized object
  meth <- unite(myObj.reorg, min.per.group=min.per.group, mc.cores=8)
  
  return(meth)
}

# FUNCTION to calculate DMCs and export results
# will utilize methylKit's "MN" overdispersion and "Chisq" test
calc.DMCs <- function(my.meth, covariate=NULL, comparison, meth.diff=10, qval=0.1) {
  # Calculate DMCs: Overdispersion:YES, Test:Chisq
  myDiff <- calculateDiffMeth(my.meth,
                              covariates=covariate,
                              overdispersion="MN",
                              test="Chisq",
                              mc.cores=8
  )
  # convert results to a data frame
  myDiff.df <- getData(myDiff)
  
  # Plot p-value distribution
  png(paste0(comparison, ".pval_dist.png"))
  hist(myDiff$pvalue)
  dev.off()
  
  # Subset to significant DMCs
  myDiff.sig <- getMethylDiff(myDiff, difference=meth.diff, qvalue=qval)
  # convert to data frame
  myDiff.sig.df <- getData(myDiff.sig)

  # remove DMCs on noncanonical contigs
  myDiff.sig.df2 <- myDiff.sig.df[myDiff.sig.df$chr %in% cans, ]
  
  # export sig results
  write.table(myDiff.sig.df2,
              paste0(comparison, ".sigDMCs.txt"),
              sep="\t",
              col.names=TRUE,
              row.names=FALSE,
              quote=FALSE
  )
}
