#!/usr/bin/Rscript


# INPUT:
# cov_file: path to desired gzipped cov file
# organism: model organism, either "mouse", "rat", "human"

# canonical contigs include all numbered contigs, and X and Y


# FUNCTION to turn the cov file into a gzipped canonical cov file
makeCanCovFiles <- function(cov, chrs) {
  # remove the suffix to get name for outfile
  outname <- gsub('_trimmed_bismark_bt2.bismark.cov.gz', '', cov_file)

  # keep only canonical contigs
  cov2 <- cov[cov$V1 %in% chrs, ]

  # create a gzipped file for writing
  gz1 <- gzfile(paste(outname, "canonical.cov.gz", sep="."), "w")

  # export the canonical bedgraph to a gzipped file
  write.table(cov2, gz1, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  close(gz1)
}



args <- commandArgs(trailingOnly = TRUE)
cov_file <- args[1]
organism <- args[2]

cov <- read.table(gzfile(cov_file), header=F)

if (organism == "mouse") {
  chrs <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X','Y')
} else if (organism == "rat") {
  chrs <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','X','Y')
} else if (organism == "human") {
  chrs <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')
}
    
# execute the function to make a canonical cov file
makeCanCovFiles(cov, chrs)
