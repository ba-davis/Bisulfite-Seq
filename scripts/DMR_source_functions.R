# Differentially Methylated CpG Analysis with methylKit #

#library(methylKit)
#library(genomation)
#library(ggplot2)

options(scipen = 999) # disables scientific notation


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
# FUNCTION # to filter and normalize myObj #
############-------------------------------#
# norm.method can be "median" or "mean"
filter <- function(myObj, low=10, high=99.9, normalize=norm, norm.method=norm_method) {
  # filter myObj CpGs
  myObj.filtered <- filterByCoverage(myObj,
                                     lo.count=low,
                                     lo.perc=NULL,
                                     hi.count=NULL,
                                     hi.perc=high
  )

  if (normalize=="TRUE") {
    # normalize by coverage
    myObj.filtered <- normalizeCoverage(myObj.filtered,method=norm.method)
  }
  
  return(myObj.filtered)
}

############---------------------------------------------------------------#
# FUNCTION # to reorganize myObj, tile the genome, and merge tiled regions #
############---------------------------------------------------------------#
reorg <- function(myObj, samples, treat, win.size=1000, step.size=1000, min.per.group) {
  # reorganize the desired object
  myObj.reorg <- reorganize(myObj,
                            sample.ids=samples,
                            treatment=treat
  )

  cpg.meth <- unite(myObj.reorg, min.per.group=min.per.group, mc.cores=8)
  print(paste0("Number of CpGs included in tiles: ", nrow(cpg.meth)))
  
  # tile the reorganized object
  my.tiles <- tileMethylCounts(myObj.reorg,
                               win.size=win.size,
                               step.size=step.size,
                               mc.cores=8
  )

  # merge the tiles
  meth <- unite(my.tiles, min.per.group=min.per.group, mc.cores=8)

  return(meth)
}

############------------------------------------------------------------------------#
# FUNCTION # to construct a covariates dataframe for use with differential analysis #
############------------------------------------------------------------------------#
make_covar <- function(covar1, covar2) {
  # check to see if both covariates are 'NULL'
  if (covar1 == 'NULL' & covar2 == 'NULL') {
    new_covariates <- NULL
  } else {
      covariates <- data.frame(covar1=covar1, covar2=covar2)

      # remove any columns that contain 'NULL'
      cols <- colSums(mapply('==', 'NULL', covariates))
      new_covariates <- covariates[,which(cols == 0)]
  }
  return(new_covariates)
}

############--------------------------------------#
# FUNCTION # to calculate DMRs and export results #
############--------------------------------------#
# will utilize methylKit's "MN" overdispersion and "Chisq" test
calc.DMRs <- function(my.meth, covariate=NULL, overdispersion="MN", test="Chisq", comparison, meth.diff=10, qval=0.1) {
  # Calculate DMRs: Overdispersion:YES, Test:Chisq
  myDiff <- calculateDiffMeth(my.meth,
                              covariates=covariate,
                              overdispersion=overdispersion,
                              test=test,
                              mc.cores=8
  )
  # convert results to a data frame
  myDiff.df <- getData(myDiff)
  # Plot p-value distribution
  png(paste0(comparison, ".pval_dist.png"))
  hist(myDiff$pvalue)
  dev.off()
  # Subset to significant DMRs
  myDiff.sig <- getMethylDiff(myDiff, difference=meth.diff, qvalue=qval)
  # convert to data frame
  myDiff.sig.df <- getData(myDiff.sig)
 
  #----------------------------#
  # HYPO and HYPER SEPARATIONS #
  #----------------------------#
  
  # separate hyper and hypo sig DMRs
  myDiff.hyper <- myDiff.sig.df[myDiff.sig.df$meth.diff > 0, ]
  myDiff.hypo <- myDiff.sig.df[myDiff.sig.df$meth.diff < 0, ]

  # df for pie chart of hyper and hypo sig DMRs
  pie.df <- data.frame(class=c("hyper", "hypo"),
  	               percent=c(nrow(myDiff.hyper) / nrow(myDiff.sig.df),
		                 nrow(myDiff.hypo) / nrow(myDiff.sig.df))
  )
  pie.df$label=paste0(round(as.numeric(pie.df$percent)*100, digits=0), '%')
  
  # define blank theme for prettier pie chart
  blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

  # plot ratio of hyper vs hypo sig DMRs
  myplot <- ggplot(pie.df, aes(x="", y=percent, fill=class)) +
             geom_bar(width = 1, stat = "identity") +
	     coord_polar("y", start=0) +
	     blank_theme +
	     theme(axis.text.x=element_blank()) +
	     #geom_text(aes(x=1, y=pos, label=label), size = 6) +
	     geom_text(aes(x=1, y=cumsum(percent/sum(percent)) - percent/sum(percent)/2,
	       label=label[c(2,1)]), size=6) +  #label = paste(round(percent/sum(percent)*100),"%")), size = 6)
	     ggtitle(paste0("Ratio of Hyper and Hypo Methylated sig DMRs\ntotal: ", nrow(myDiff.sig.df)))
  ggsave(filename=paste0(comparison, ".hyper.hypo.sigDMR.pieChart.png"))

  # -------------------------------#
  # continue with both hyper and hypo export
  
  # add DMR_ID column
  myDiff.sig.df$DMR_ID <- paste('DMR', seq(1:nrow(myDiff.sig.df)), sep='_')

  # reorder columns
  myDiff.sig.df <- myDiff.sig.df[ ,c(8,1,2,3,5,6,7)]
  
  # export sig results
  write.table(myDiff.sig.df,
              paste0(comparison, ".sigDMRs.txt"),
              sep="\t",
              col.names=TRUE,
              row.names=FALSE,
              quote=FALSE
  )

  #--------------#
  # sig dmr beds #
  #--------------#

  # BED file containing all sig DMRs and ensembl contig format
  # subset sig DMR results to a bed file for annotation
  res <- myDiff.sig.df[ ,c(2,3,4,1)]
  res$start <- res$start - 1

  # export sig DMR bed for annotation
  write.table(res,
              paste0(comparison, ".sigDMRs.bed"),
	      sep="\t",
	      col.names=FALSE,
	      row.names=FALSE,
	      quote=FALSE
  )

  # 3 BED files for GREAT
  # UCSC contig format, 1 bed file of hyper only, 1 of hypo only, 1 of both
  res$chr <- paste0("chr", res$chr)
  write.table(res,
              paste0(comparison, ".sigDMRs.both.GREAT.bed"),
              sep="\t",
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE
  )
  
  my.hyper <- myDiff.sig.df[myDiff.sig.df$meth.diff > 0, c(2,3,4,1)]
  my.hyper$start <- my.hyper$start - 1
  my.hyper$chr <- paste0("chr", my.hyper$chr)
  write.table(my.hyper,
              paste0(comparison, ".sigDMRs.hyper.GREAT.bed"),
              sep="\t",
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE
  )

  my.hypo <- myDiff.sig.df[myDiff.sig.df$meth.diff < 0, c(2,3,4,1)]
  my.hypo$start <- my.hypo$start - 1
  my.hypo$chr <- paste0("chr", my.hypo$chr)
  write.table(my.hypo,
              paste0(comparison, ".sigDMRs.hypo.GREAT.bed"),
              sep="\t",
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE
  )


  return(myDiff.sig.df)
}

# FUNCTION to make a sig DMR bed file track from results dataframe
makeBED <- function(res, comparison) {
  res2 <- res[ ,c(2,3,4,1,7)]
  res2$start <- res2$start - 1
  res2$chr <- paste0("chr", res2$chr)

  res2[ ,6] <- "."
  res2[ ,7] <- res2$start
  res2[ ,8] <- res2$end
  res2[ ,9] <- ifelse(res2[ ,5] > 0, '255,0,0', ifelse(res2[ ,5] < 0, '0,0,255', '0,0,0'))
  
  # add track line to exported file
  cat(paste0("track type=bed name=", comparison, " ", "itemRgb=On"), file=paste0(comparison, ".sigDMRtrack.bed"))
  cat("\n", file=paste0(comparison, ".sigDMRtrack.bed"), append=TRUE)

  # export bed info
  write.table(res2,
              paste0(comparison, ".sigDMRtrack.bed"),
	      sep="\t",
	      col.names=FALSE,
	      row.names=FALSE,
	      quote=FALSE,
	      append=TRUE
  )
}
