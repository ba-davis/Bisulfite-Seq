# Bisulfite-Seq
Scripts related to processing of bisulfite sequencing data via [Snakemake](https://snakemake.readthedocs.io/en/stable/) and the SLURM workload manager for parallel computing.

Makes use of [methylKit](https://bioconductor.org/packages/release/bioc/html/methylKit.html) for differential analysis

Standard Steps include:
1. Check quality of raw data with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Trim with [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) (-rrbs option for RRBS data)
3. Check quality of trimmed reads with FastQC
4. Align trimmed reads with [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
5. Collect trimming and alignment statistics into tab delimited text files
6. Extract methylation data with bismark_methylation_extractor
7. Produce canonical coverage files for downstream use
8. Perform initial data exploration, including:
  - number of CpGs covered per sample at various coverage cutoffs
  - "coverage plots" from methylKit
  - "methylation stats plots" from methylKit
  - dimension reduction and visualization of samples, including PCA and Nonnegative Matrix Factorization (NMF) + tSNE
9. Perform DMR analysis via methylKit to find differentially methylated regions between two groups
10. Produce a DMR bed file track for a genome browser
11. Produce 3 bed files per comparison for use with [GREAT](http://great.stanford.edu/public/html/) pathway analysis
  - Both hyper and hypomethylated regions
  - Hypermethylated regions only
  - Hypomethylated regions only
