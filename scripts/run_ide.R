#!/usr/bin/Rscript

library(methylKit)
library(ggplot2)
library(NMF)
library(isva)
library(genefilter)
library(Rtsne)

# source the script holding the helper functions
source('scripts/ide_source_functions.R')

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
cov_files <- unlist(strsplit(args[1], ","))
samp_names <- unlist(strsplit(args[2], ","))
my_treat <- as.numeric(unlist(strsplit(args[3], ",")))
my_treatment <- unlist(strsplit(args[4], ","))
my_condition <- unlist(strsplit(args[5], ","))
my_rep <- unlist(strsplit(args[6], ","))
my_sex <- unlist(strsplit(args[7], ","))
my_genotype <- unlist(strsplit(args[8], ","))
my_batch <- unlist(strsplit(args[9], ","))
my_pca_color_var <- unlist(strsplit(args[10], ","))
my_pca_shape_var <- unlist(strsplit(args[11], ","))
my_pca_label_var <- unlist(strsplit(args[12], ","))



# Run first function to create a methylbase object from the cov files
myObj <- create.myObj(cov_files)

# create and setwd to ide directory
dir.create(file.path(getwd(), 'data/ide'), showWarnings = FALSE)
setwd(file.path(getwd(), 'data/ide'))


# Run function to produce CpG coverage table, use default values
#  by default, returns the last specified coverage filtering (10x 99.9)
#  default produces a CpG coverage table
myObj.filtered <- cpg_cov()

# make methylation stats plots
makeMethStats()

 # make coverage plots
makeCovPlots()

# make coldata
coldata <- make_coldata(samp_names, my_treatment, my_condition, my_rep, my_sex, my_genotype, my_batch)

# make pca plot(s)
meth_dim_reduct(myObj.filtered, min_pc=3, coldata=coldata, color_var=my_pca_color_var, shape_var=my_pca_shape_var, label_var=my_pca_label_var)






write.table(coldata, "coldata.txt", sep="\t", col.names=T, row.names=F, quote=F)














# NOTES

# For now, every variable option in the nmf and tsne plot are hardcoded other than those relevant to PCA:
#  coldata (won't change)
#  color_var
#  shape_var
#  label_var

# and perplexity which is calculated within the function
# may need to add additional variables if want to try different parameters and create multiple tsne plots







# old:

# PCA function
#  before you run the entire snakemake pipeline, what info could you provide to explain
#  the design and variables, and which variables to be used for pca plot for:
#   - color
#   - shape
#   - label
# already define 6 possible variables for coldata in snakefile
#  my scripts will need to check for the presence of these variables
# maybe in snakefile i define new variables:
#   color_by =
#   shape_by =
#   label_by =

# will also be useful to have an easily accessible script to manually construct PCA plots
#  if after the snakerun there is more fine tuning required
