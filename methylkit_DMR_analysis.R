# Differentially Methylated CpG Analysis with methylKit #

library(methylKit)
library(genomation)
library(ggplot2)

options(scipen = 999) # disables scientific notation

#---------------------------------------#
# Variables to change for each project: #
#---------------------------------------#

# path to coverage files
# ex: dir <- "/home/groups/hoolock/u1/bd/Projects/ECP28/meth_extract/relaxed/cov_files"

# set the desired order of the cov.files (group samples together)
# example: cov.order <- c(19,15,23,9,16,1,21,20,17,18,24,22,5,6,11,3,12,14,7,4,13,2,10,8)
# ex: cov.order <- c(1,4,5,6,7,8,9,10,11,2,3)

# list the sample names (same order as desired cov.order)
# ex: samp.names <- c("MOR1","MOR2","MOR3","MOR4","MOR5","FOR1","FOR2","FOR3","FOR4","EWE1","EWE2")

# create treatment vector for diff analysis (can have more than 2 groups, we will reorganize)
# ex: my.treat <- c(1,1,1,1,1,2,2,2,2,3,3)

# create a covariate data frame for differential analysis (may need more than 1)
# ex: covar1 <- data.frame(sex=c("M","M","F","M","F","F","M"))

# canonical contigs
# ex: sheep: cans <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","X","Y")

# bed12 file
# ex: my.bed <- "/home/groups/hoolock/u1/genomes/Epigenetics_Core/sheep/oar_v3.1/annotation/Ovis_aries.Oar_v3.1.95.sorted.bed"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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

# FUNCTION to reorganize myObj, tile the genome, and merge tiled regions
reorg <- function(myObj, samples, treat, win.size=1000, step.size=1000, min.per.group) {
  # reorganize the desired object
  myObj.reorg <- reorganize(myObj,
                            sample.ids=samples,
                            treatment=treat
  )

  cpg.meth <- unite(myObj.reorg,  min.per.group=min.per.group, mc.cores=8)
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

# FUNCTION to calculate DMRs and export results
# will utilize methylKit's "MN" overdispersion and "Chisq" test
calc.DMRs <- function(my.meth, covariate=NULL, comparison, meth.diff=10, qval=0.1, bed) {
  # Calculate DMRs: Overdispersion:YES, Test:Chisq
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
  # Subset to significant DMRs
  myDiff.sig <- getMethylDiff(myDiff, difference=meth.diff, qvalue=qval)
  # convert to data frame
  myDiff.sig.df <- getData(myDiff.sig)
  # remove DMRs on noncanonical contigs
  myDiff.sig.df2 <- myDiff.sig.df[myDiff.sig.df$chr %in% cans, ]

  # read in the bed12 file
  gene.obj <- readTranscriptFeatures(bed,
				     up.flank=3000,
				     down.flank=0,
				     unique.prom=FALSE
  )

  # get percent overlaps of sig DMRs with genic parts
  foo <- annotateWithGeneParts(as(myDiff.sig.df2, "GRanges"), gene.obj)

  # plot the percentage of sig DMR overlaps with promoters, exons, introns, intergenic,
  #  with precedence set as promoter > exon > intron
  png(paste0(comparison, ".sigDMR.anno.pieChart.png"))
  plotTargetAnnotation(foo, precedence=TRUE, main=paste(comparison, "sig DMR overlaps", sep=" "))
  dev.off()

  #----------------------------#
  # HYPO and HYPER SEPARATIONS #
  #----------------------------#
  
  # separate hyper and hypo sig DMRs
  myDiff.hyper <- myDiff.sig.df2[myDiff.sig.df2$meth.diff > 0, ]
  myDiff.hypo <- myDiff.sig.df2[myDiff.sig.df2$meth.diff < 0, ]

  # df for pie chart of hyper and hypo sig DMRs
  pie.df <- data.frame(class=c("hyper", "hypo"),
  	               percent=c(nrow(myDiff.hyper) / nrow(myDiff.sig.df2),
		                 nrow(myDiff.hypo) / nrow(myDiff.sig.df2))
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
	     ggtitle(paste0("Ratio of Hyper and Hypo Methylated sig DMRs\ntotal: ", nrow(myDiff.sig.df2)))
  ggsave(filename="hyper.hypo.sigDMR.pieChart.png")

  # Now, produce genic overlap pie charts for hyper and hypo DMRs separately
  # get percent overlaps of sig hyper DMRs with genic parts
  foo.hyper <- annotateWithGeneParts(as(myDiff.hyper, "GRanges"), gene.obj)

  # plot the percentage of sig hyper DMR overlaps with promoters, exons, introns, intergenic,
  #  with precedence set as promoter > exon > intron
  png(paste0(comparison, ".sigDMR.hyper.anno.pieChart.png"))
  plotTargetAnnotation(foo.hyper, precedence=TRUE, main=paste(comparison, "sig hyper DMR overlaps", sep=" "))
  dev.off()

  # get percent overlaps of sig hypo DMRs with genic parts
  foo.hypo <- annotateWithGeneParts(as(myDiff.hypo, "GRanges"), gene.obj)

  # plot the percentage of sig hypo DMR overlaps with promoters, exons, introns, intergenic,
  #  with precedence set as promoter > exon > intron
  png(paste0(comparison, ".sigDMR.hypo.anno.pieChart.png"))
  plotTargetAnnotation(foo.hypo, precedence=TRUE, main=paste(comparison, "sig hypo DMR overlaps", sep=" "))
  dev.off()

  # -------------------------------#
  # continue with both hyper and hypo export
  
  # add DMR_ID column
  myDiff.sig.df2$DMR_ID <- paste('DMR', seq(1:nrow(myDiff.sig.df2)), sep='_')

  # reorder columns
  myDiff.sig.df2 <- myDiff.sig.df2[ ,c(8,1,2,3,5,6,7)]
  
  # export sig results
  write.table(myDiff.sig.df2,
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
  res <- myDiff.sig.df2[ ,c(2,3,4,1)]
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
              paste0(comparison, ".sigDMRs.both.bed"),
              sep="\t",
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE
  )
  
  my.hyper <- myDiff.sig.df2[myDiff.sig.df2$meth.diff > 0, c(2,3,4,1)]
  my.hyper$start <- my.hyper$start - 1
  write.table(my.hyper,
              paste0(comparison, ".sigDMRs.hyper.bed"),
              sep="\t",
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE
  )

  my.hypo <- myDiff.sig.df2[myDiff.sig.df2$meth.diff < 0, c(2,3,4,1)]
  my.hypo$start <- my.hypo$start - 1
  write.table(my.hypo,
              paste0(comparison, ".sigDMRs.hypo.bed"),
              sep="\t",
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE
  )


  return(myDiff.sig.df2)
}

# FUNCTION to make a sig DMR bed file track from results dataframe
makeBED <- function(res, name) {
  res2 <- res[ ,c(2,3,4,1,7)]
  res2$start <- res2$start - 1
  res2$chr <- paste0("chr", res2$chr)

  res2[ ,6] <- "."
  res2[ ,7] <- res2$start
  res2[ ,8] <- res2$end
  res2[ ,9] <- ifelse(res2[ ,5] > 0, '255,0,0', ifelse(res2[ ,5] < 0, '0,0,255', '0,0,0'))
  
  # add track line to exported file
  cat(paste0("track type=bed name=", name), file=paste(name, "bed", sep="."))
  cat("\n", file=paste(name, "bed", sep="."), append=TRUE)

  # export bed info
  write.table(res2,
              paste(name, "bed", sep="."),
	      sep="\t",
	      col.names=FALSE,
	      row.names=FALSE,
	      quote=FALSE,
	      append=TRUE
  )
}
