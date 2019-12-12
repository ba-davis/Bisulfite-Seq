
# Using methylKit to perform Initial Data Exploration on BS_seq data #

# Obtain:
#  - CpG Numbers at various coverage cutoffs (default: 1X, 10X, 40X, 10X_99.9%) per sample
#  - Methylation Stats Plots per Sample
#  - Coverage Plots per sample
#  - PCA Plots
#  - anything else? Coverage Boxplots? % methylation boxplots? distributions of some sort?


############------------------------------------------------#
# FUNCTION # to create the methyl object from the cov files #
############------------------------------------------------#
create.myObj <- function(cov_files) {
  # capture the cov.file paths
  my_cov_files <- lapply(cov_files, function(x) x)
  # capture the samp.names
  my_samp_names <- lapply(samp_names, function(x) x)

  # use methRead function from methylKit to read in the cov.files
  myObj <- methRead(my_cov_files,
                    sample.id=my_samp_names,
                    assembly="genome",
                    treatment=my_treat,
                    pipeline="bismarkCoverage",
                    context="CpG",
                    mincov=1)
  return(myObj)
}

############-------------------------------#
# FUNCTION # to produce CpG Coverage Table #
############-------------------------------#

# This function can primarily be used to filter the CpGs to a defined coverage cutoff
#  default is cov4 and perc4 (set to 10X and 99.9% respectively)
# if table is set to true (default), a CpG coverage table will be written to stdout

cpg_cov <- function(cov1=1, perc1=NULL, cov2=10, perc2=NULL, cov3=40, perc3=NULL, cov4=10, perc4=99.9, table=TRUE) {

  # cov1
  myObj.filtered=filterByCoverage(myObj,lo.count=cov1,lo.perc=NULL,hi.count=NULL,hi.perc=perc1)
  my.cov1 <- sapply(myObj.filtered, function(x) nrow(x))

  # cov2
  myObj.filtered=filterByCoverage(myObj,lo.count=cov2,lo.perc=NULL,hi.count=NULL,hi.perc=perc2)
  my.cov2 <- sapply(myObj.filtered, function(x) nrow(x))

  # cov3
  myObj.filtered=filterByCoverage(myObj,lo.count=cov3,lo.perc=NULL,hi.count=NULL,hi.perc=perc3)
  my.cov3 <- sapply(myObj.filtered, function(x) nrow(x))

  # cov4
  myObj.filtered=filterByCoverage(myObj,lo.count=cov4,lo.perc=NULL,hi.count=NULL,hi.perc=perc4)
  my.cov4 <- sapply(myObj.filtered, function(x) nrow(x))

  if (table == TRUE) {
    # Create data frame of CpG coverages
    cpg_covs <- data.frame(sample=unlist(samp_names, use.names=FALSE),
  	                   V1=my.cov1,
	                   V2=my.cov2,
	                   V3=my.cov3,
	                   V4=my.cov4
    )
  
    colnames(cpg_covs) <- c('sample',
  		            paste('CpG', cov1, perc1, sep='.'),
			    paste('CpG', cov2, perc2, sep='.'),
			    paste('CpG', cov3, perc3, sep='.'),
			    paste('CpG', cov4, perc4, sep='.')
    )

    # Export CpG Coverage Table
    write.table(cpg_covs, "CpG.coverage.table.cov_files.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  }
  
  return(myObj.filtered)
}

############-----------------------------------#
# FUNCTION # to create methylation stats plots #
############-----------------------------------#
makeMethStats <- function(outdir=getwd()) {

  dir.create(file.path(getwd(), 'meth_stats_plots'), showWarnings = FALSE)
  setwd(file.path(getwd(), 'meth_stats_plots'))

  for (i in 1:length(myObj)) {
      png(paste0(samp_names[[i]],"_methstats.png"))
      getMethylationStats(myObj[[i]], plot=TRUE, both.strands=FALSE)
      dev.off()
  }
  setwd('..')
}

############--------------------------#
# FUNCTION # to create coverage plots #
############--------------------------#
makeCovPlots <- function(outdir=getwd()) {

  dir.create(file.path(getwd(), 'coverage_plots'), showWarnings = FALSE)
  setwd(file.path(getwd(), 'coverage_plots'))
  
  for (i in 1:length(myObj)) {
      png(paste0(samp_names[[i]],"_covplot.png"))
      getCoverageStats(myObj[[i]], plot=TRUE, both.strands=FALSE)
      dev.off()
  }
  setwd('..')
}

############----------------------------------------------------------------#
# FUNCTION # to create coldata dataframe from command line vector variables #
############----------------------------------------------------------------#
make_coldata <- function(samp_names, my_treatment, my_condition, my_rep, my_sex, my_genotype, my_batch) {
  # construct coldata from given input variables

  coldata <- data.frame(sample_name = samp_names,
  	     		treatment = my_treatment,
			condition = my_condition,
			replicate = my_rep,
			sex = my_sex,
			genotype = my_genotype,
			batch = my_batch
  )

  # remove any columns that contain 'NULL'
  cols <- colSums(mapply('==', 'NULL', coldata))
  new_coldata <- coldata[,which(cols == 0)]

  return(new_coldata)
}


############------------------------------------------------------------------------------#
# FUNCTION # to convert the meth object into a useable format for various downstream use: #
############------------------------------------------------------------------------------#
# library(dplyr)
# requires dplyr for producing percent meth and coverage df
# can produce either:
#              percent meth matrix
#              percent meth dataframe
#              percent meth and coverage dataframe
# default: df=FALSE will produce a matrix
#          df=TRUE will produce a dataframe with an ID column
#          coverage: default is FALSE and will not include coverage info, just percent meth per sample
#                    TRUE will include coverage info (for use with gene_plots)
convertToPercMethMatrix <- function(meth, df=FALSE, coverage=FALSE) {

  # convert the meth object to a dataframe
  mydf <- getData(meth)

  # determine number of samples
  num <- (ncol(mydf) - 4) / 3

  # iteratively create perc.meth columns to be numCs / coverage per sample
  for (i in 1:num) {
    num.name <- paste0("numCs", i)
    den.name <- paste0("coverage", i)
    new.name <- paste0("percent.meth", i)

    mydf[[new.name]] <- mydf[[num.name]] / mydf[[den.name]]
  }

  # create ID column to be chr_start (1-based)
  mydf$ID <- paste(mydf[,1], mydf[,2], sep="_")

  # create either the matrix or df with ID column depending on input param 'df'
  if (df == FALSE) {
    # create perc.meth df without ID
    new.matrix <- mydf[ ,c((num*3+6):ncol(mydf)-1)]

    # final.matrix
    final.matrix <- data.matrix(new.matrix)
  } else if (df == TRUE & coverage == FALSE) {
      # create perc.meth df with ID as first column
      final.matrix <- mydf[ ,c(ncol(mydf), (num*3+6):ncol(mydf)-1)]
  } else if (df == TRUE & coverage == TRUE) {
      # create perc.meth and coverage df with ID as first column (use "select" from dplyr)
      final.matrix <- mydf %>% select(grep("ID", colnames(mydf)), grep("percent.meth", colnames(mydf)), grep("coverage", colnames(mydf)))
  }
    return(final.matrix)
}

############------------------------------------#
# FUNCTION # to describe a prcomp result object #
############------------------------------------#
describe_pca <- function(mypca) {
  # Determine how many PC's were returned
  print(paste0("Number of PC's returned: ", ncol(mypca$x)))

  # Obtain the eigenvalues (can get values proportional to eigenvalues by taking sd^2)
  eigs <- mypca$sdev^2

  # Determine number of PC's with eigenvalue > 1 (considered important)
  print(paste0("Number of PC's with eigenvalue > 1: ", sum(eigs > 1)))

  # Determine how much of total variance is explained by first PC
  print(paste0("Percent of total variance explained by first PC: ", round((eigs[1]/sum(eigs))*100, digits=2), "%"))

  # How many PC's are needed to explain at least 80% of total variance
  my_sum = 0
  num_pc = 0
  for (i in 1:ncol(mypca$x)) {
    my_sum = my_sum + (eigs[i] / sum(eigs))
    num_pc = num_pc + 1
    if (my_sum >= 0.8) {
      break
    }
  }
  print(paste0("Number of PC's required to explain at least 80% of the variance: ", num_pc))

  return(num_pc)
}

############-------------------------------#
# FUNCTION # to create a biplot from a PCA #
############-------------------------------#
plot_pca <- function(df, pc_a="PC1", pc_b="PC2", color_var, shape_var, label_var, eigs=eigs) {
  
  # get the percent variance explained by the two PC's
  pc_a_var <- round(((eigs[as.numeric(gsub("PC", "", pc_a))] / sum(eigs)) * 100), 1)
  pc_b_var <- round(((eigs[as.numeric(gsub("PC", "", pc_b))] / sum(eigs)) * 100), 1)

  # subtitle
  #subtitle <- paste0(pc_a, " by ", pc_b)

  # check to see if any variables are "NULL" and change them to NULL
  if (color_var == "NULL") {
    color_var <- NULL
  }
  if (shape_var == "NULL") {
    shape_var <- NULL
  }
  if (label_var == "NULL") {
    label_var <- NULL
  }

  # filename
  if (is.null(shape_var)) {
    filename <- paste0(pc_a, "_", pc_b, "_", "color_by_", color_var, ".png")
  } else {
    filename <- paste0(pc_a, "_", pc_b, "_", "color_", color_var, "_shape_", shape_var, ".png")
  }
  
  # General Format of plotting PCA
  myplot <- ggplot(df, aes_string(pc_a, pc_b, color=color_var, shape=shape_var)) +
    geom_point(size=5) +
    xlab(paste0(pc_a, ": ", pc_a_var, "% variance")) +
    ylab(paste0(pc_b, ": ", pc_b_var, "% variance")) +
    ggtitle("PCA: CpG Methylation") +
    geom_text(aes_string(label=label_var),hjust=0.7, vjust=-1.1, show.legend = FALSE) +
    theme_classic() +
    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
  ggsave(filename=filename)
}


############-----------------------------------------------------------------------------------#
# FUNCTION # to perform NMF on a merged cpg methylation object and return coefficient matrix h #
############-----------------------------------------------------------------------------------#
# references and uses the convertToPercMethMatrix function defined above
# inputs:
#        meth: merged cpg methylation object from methylkit (often called "meth")
#        ntop: number of most variable CpGs to use, default is 5000
#        nmf_method: choice of NMF algortihm to use from CRAN's NMF package, default is 'brunet'
#        nmf_seed: choice if initialization seeding method for NMF, default is 'random'
#        nrun: if using a random seed, must run multiple runs, default is 100

# returns: h, the mixture coefficient matrix

# NOTES: to add:
#  RMT is not applicable with 10 or fewer samples, the permutation procedure from Buja-Eyuboglu (1992) is better
#  before applying EstDimRMT, data matrix should be centered (each row (genomic feature) has mean zero).

run_nmf <- function(meth, ntop=5000, nmf_method="brunet", nmf_seed="random", nrun=100) { 
  # convert the "meth" object to a percent methylation matrix (cpgs as rows X samples as columns)
  mymeth <-convertToPercMethMatrix(meth)

  # keep only unique rows (remove redundant rows)
  mymeth.unique <- unique(mymeth)

  # Select the "ntop" most variable CpGs (default is 5000)
  Pvars <- rowVars(mymeth.unique)
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))]
  mymeth.top <- mymeth.unique[select, ]

  # before using EstDimRMT, we need to scale our data so that each row has mean 0
  mymeth.top_scaled <- t(scale(t(mymeth.top), scale=FALSE))

  # use Random Matrix Theory (RMT) from isva package to estimate number of significant components of variation
  #  (factorization rank), k, estimated intrinsic dimensionality of the data
  #  (requires rows as features and columns as samples)
  png('rmt.plot.png')
  my.est<- EstDimRMT(mymeth.top_scaled, plot=TRUE)
  dev.off()
  k <- my.est$dim
  print(paste0("Estimated dimensionality is: ", k))

  # Run NMF
  # use the unscaled mymeth.top (all nonnegative values)
  # default is to use brunet with random seed and 100 runs
  # good value for nrun is 100-200
  # only returns the 'best' fit over all the runs (factorization that achieved the lowest approximation error)
  set.seed(42)
  my.nmf <- nmf(mymeth.top, rank=k, method=nmf_method, seed=nmf_seed, nrun=nrun, .opt=list(parallel=16))

  # get the basis matrix, w
  w <- basis(my.nmf)
  # get the mixture coefficient matrix, h
  h <- coef(my.nmf)

  return(h)
}


############----------------------------------------------------------------------------#
# FUNCTION # to create a t-SNE plot of the mixture coefficient matrix returned from NMF #
############----------------------------------------------------------------------------#
# inputs:
#         h: mixture coefficient matrix returned from NMF function
#         tsne_pca: whether or not to perform                    
#         theta: speed/accuracy tradeoff, increase for less accuracy
#                set to 0.0 for exact TSNE
#                default is 0.5
#         check_duplicates: default FALSE (we have unique data from the NMF method)
#         max_iter: number of iterations, default is 1000
#         num_threads: number of threads to use, default is 8
#         coldata: coldata
#         color_var: variable from coldata to color by, same as on PCA
#         shape_var: variable from coldata to shape the points by, same as on PCA
#         label_var: variable from coldata to label points by, same as on PCA
#         v1: first column of tsne output to plot, default "x"
#         v2: second columns of tsne output to plot, default "y"

# NOTE: perplexity is defined within the function based on number of samples
#         3 * perplexity should be less than number of samples - 1
#         so, if you have 16 samples, perplexity can't be greater than 5
#         typical values for perplexity range between 5 and 50
    
plot_tsne <- function(h, tsne_pca=FALSE, theta=0.5, check_duplicates=FALSE, max_iter=1000, num_threads=8, coldata=coldata, color_var, shape_var, label_var, v1="x", v2="y") {
  # transpose the h matrix to get samples as rows
  mydat <- as.matrix(as.data.frame(t(h)))

  # set perplexity (p) equal to the max value it can be considering the number of samples
  p <- as.integer((nrow(mydat) - 1) / 3)

  # create tsne plot and use perplexity defined above
  filename <- paste0("nmf_tsne_perp", p, ".png")
  tsne_out <- Rtsne(mydat, pca=tsne_pca, perplexity=p, theta=theta, check_duplicates=check_duplicates, max_iter=max_iter, num_threads=num_threads)
  tsnedf <- data.frame(x=tsne_out$Y[,1], y=tsne_out$Y[,2])
  plot_df <- cbind(tsnedf, coldata)
  myplot <- ggplot(plot_df, aes_string(x=v1, y=v2, color=color_var, shape=shape_var)) +
    geom_point(size=5) +
    xlab("tSNE 1") +
    ylab("tSNE 2") +
    ggtitle("NMF+tSNE")) +
    geom_text(aes_string(label=label_var), hjust=0.7, vjust=-1.1, show.legend = FALSE) +
    theme_classic() +
    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
  ggsave(filename=filename)
}


############----------------------------------------------------------------------#
# FUNCTION # to perform PCA and NMF+t-SNE dimension reduction on methylation data #
############----------------------------------------------------------------------#
# plans fot this function:
#  use PCASamples function from methylkit to run PCA (prcomp) on meth (cpg % meth matrix)
#  maybe just run on the 5000 most variable CpGs?
#  use the screeplot to determine if PCA works well for the data
#    - PC's with eigenvalues >= 1 are considered important to look at
#    - or, the selected PC's should describe at least 80% of the total variance
#    - if you end up with too many PC's (>3), then PCA may not be best for the data
#        - so, run t-sne or MDS
#          - need to determine number of dimensions in which to search
#              - use Random Matrix Theory (RMT), even though data is non-gaussian
#                - ex: estimated 14 dimensions, so set K=14
# PCA: When implementing dimension reduction with Gaussian assumptions, PCA is widely used
#      After taking eigenvalue decomposition on the covariance matrix of the observed data,
#       PCA keeps K eigenvectors that correspond to the K largest eigenvalues.
#
# BG-NMF + VBBMM method performs the best 
# It looks like epicluster package has BG-NMF and other related functions, but difficult to install
# For now, can try NMF (performed by Adey Lab on methylation data), followed by t-sne

# may want to create a pca function and a NMF+tsne function, then call these functions in the meth_dim_reduct function


# inputs: myObj: desired methylRawList object to merge and run pca on
#         min_pc=3: sets how many pc's is too many
#                    if more than this many PCs are required to explain 80% of variance, try t-SNE
#         coldata: dataframe of coldata
#         color_var: variable from coldata to color PCA by, default is "treatment"
#         shape_var: variable from coldata to add shape of points on PCA by, default is NULL
#         label_var: variable from coldata to label points on PCA by, default is "replicate"


meth_dim_reduct <- function(myObj, min_pc=3, coldata, color_var="treatment", shape_var=NULL, label_var="replicate") {

  # create 'dim_reduct' working directory
  dir.create(file.path(getwd(), 'dim_reduct'), showWarnings = FALSE)
  setwd(file.path(getwd(), 'dim_reduct'))
 
  # normalize by coverage
  myObj=normalizeCoverage(myObj)

  # merge (destrand param is not used because no strand info is provided from coverage input files)
  meth=unite(myObj, mc.cores=4)
  print(paste0("Merging complete, number of cpgs: ", nrow(meth)))

  #######
  # PCA #
  #######
  # create 'pca' working directory
  dir.create(file.path(getwd(), 'pca'), showWarnings = FALSE)
  setwd(file.path(getwd(), 'pca'))

  # use methylkit function PCASamples to run PCA on % methylation matrix, default settings
  #  uses prcomp, scale=TRUE, center=TRUE, transpose=TRUE
  png('screeplot.png')
  mypca <- PCASamples(meth, obj.return=TRUE, screeplot=TRUE)
  dev.off()

  # describe the pca results and capture number of PC's required to explain 80% variance
  num_pc <- describe_pca(mypca)

  # grab the principle components as a new df
  my_pca_df <- as.data.frame(mypca$x)

  # add the relevant variable columns from coldata to the pca df
  my_pca_df2 <- cbind(my_pca_df, coldata)

  # get eigen values (standard deviation squared)
  eigs <- mypca$sdev^2

  # make 6 PCA plots of the various pairwise combos of PC1-4
  plot_pca(my_pca_df2, "PC1", "PC2", color_var=color_var, shape_var=shape_var, label_var=label_var, eigs=eigs)
  plot_pca(my_pca_df2, "PC1", "PC3", color_var=color_var, shape_var=shape_var, label_var=label_var, eigs=eigs)
  plot_pca(my_pca_df2, "PC1", "PC4", color_var=color_var, shape_var=shape_var, label_var=label_var, eigs=eigs)
  plot_pca(my_pca_df2, "PC2", "PC3", color_var=color_var, shape_var=shape_var, label_var=label_var, eigs=eigs)
  plot_pca(my_pca_df2, "PC2", "PC4", color_var=color_var, shape_var=shape_var, label_var=label_var, eigs=eigs)
  plot_pca(my_pca_df2, "PC3", "PC4", color_var=color_var, shape_var=shape_var, label_var=label_var, eigs=eigs)

  setwd('..')

  if (num_pc > min_pc) {
    print(paste0("More than ", min_pc, " PCs are required to explain at least 80% of the variance."))
    print("Will now try the NMF + t-SNE dimension reduction technique.")

    ###############
    # NMF + t-SNE #
    ###############

    # create directory
    dir.create(file.path(getwd(), 'nmf_tsne'), showWarnings = FALSE)
    setwd(file.path(getwd(), 'nmf_tsne'))

    # run nmf to get mixture coefficient matrix
    h <- run_nmf(meth, ntop=5000, nmf_method="brunet", nmf_seed="random", nrun=100)

    # create the tsne plot
    plot_tsne(h, tsne_pca=FALSE, theta=0.5, check_duplicates=FALSE, max_iter=1000, num_threads=8, coldata=coldata, color_var=color_var, shape_var=shape_var, label_var=label_var, v1="x", v2="y")

    setwd('..')
  }
}


############----------------------------------------------#
# FUNCTION # to run t-SNE on a percent methylation matrix #
############----------------------------------------------#
# requires library(Rtsne)
# requires library(genefilter)
# input:
# perc.meth: percent methylation matrix (output of convertToPercMethMatrix)
# ntop: number of most variable CpGs to include (default 5000)
# norm: whether to normalize the matrix or not (default TRUE)
# perplexity: value of perplexity, 3*perplexity should be less than number of samples - 1
# pca: whether initial PCA step should be performed (default: TRUE)


tsneSample <- function(perc.meth, ntop=5000, norm=TRUE) {
  library(Rtsne)
  library(genefilter)
  
  # make sure only unique values are retained
  mymeth <- unique(perc.meth)

  # subset to the 5000 most variable CpGs across all samples
  Pvars <- rowVars(mymeth)

  # select the ntop rows by variance
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))]
  mymeth.top <- mymeth[select, ]

  # transpose
  mydat <- as.matrix(as.data.frame(t(mymeth.top)))

  # normalize the input matrix if norm=TRUE
  # Mean centers each column of an input data matrix so that it has a mean of zero
  # Scales the entire matrix so that the largest absolute of the centered matrix is equal to unity
  mymeth.norm <- normalize_input(mydat)

  
  # Run t-SNE
  tsne_out <- Rtsne(mymeth.norm, pca=FALSE, perplexity=2, theta=0.0, num_cores=8)

  
  # Show the objects in 2D tsne representation
  #png('tsne.plot.png')
  #plot(tsne_out$Y, col=coldata$group)
  #dev.off()

  library(ggplot2)
  tsnedf <- data.frame(x=tsne_out$Y[,1], y=tsne_out$Y[,2], col=coldata$group)
  png('tsne.ggplot.maxiter5000.perplexity7.PCAtrue.png')
  ggplot(tsnedf) +
    geom_point(aes(x=x, y=y, color=col))
  dev.off()
}

#-------------------#
# Execute Functions #
#-------------------#

# define variables
#dir <- "/home/groups/hoolock/u1/bd/Projects/ECP1/spring2019.redo/rrbs/meth_extract/cov_files"
#cov.order <- c(19,21,3,5,10,16,18,4,7,13,14,17,20,23,2,8,9,15,1,11,22,6,12)
#samp.names <- c('E3_sham_1','E3_sham_2','E3_sham_3','E3_sham_4','E3_sham_5','E3_sham_6',
#           'E4_sham_1','E4_sham_2','E4_sham_3','E4_sham_4','E4_sham_5','E4_sham_6',
#           'E3_surg_1','E3_surg_2','E3_surg_3','E3_surg_4','E3_surg_5','E3_surg_6',
#           'E4_surg_1','E4_surg_3','E4_surg_4','E4_surg_5','E4_surg_6'
#)
#my.treat <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4)
#coldata <- data.frame(name=samp.names,
#genotype=c(rep("E3", 6), rep("E4", 6), rep("E3", 6), rep("E4", 5)),
#treatment=c(rep("sham", 12), rep("surg", 11)),
#sex=c('F','M','F','M','M','F','M','M','M','F','F','F','F','M','F','M','F','M','M','M','M','F','F'),
#group=c(rep("E3.sham", 6), rep("E4.sham", 6), rep("E3.surg", 6), rep("E4.surg", 5))
#)

# Run first function to create a methylbase object from the cov files
#myObj <- create.myObj()

# Run function to produce CpG coverage table, use default values
#  by default, returns the last specified coverage filtering
#myObj.filtered <- cpg_cov()

# make methylation stats plots
#makeMethStats()

 # make coverage plots
#makeCovPlots()

# normalize
# merge
# PCA
