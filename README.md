# Bisulfite-Seq
Scripts related to processing of bisulfite sequencing data

Standard Steps include:
1. Check quality of sequences with FastQC
2. Trim with TrimGalore (-rrbs option for RRBS data)
3. Align with Bismark
4. Extract methylation info with bismark_methylation_extractor
5. Initial Data Exploration
  - Methylation Stats Plots
  - Coverage Plots
  - CpG Coverage Table
  - PCA Plot
6. Differential Analysis
  - DMC Analysis (limma)
  - DMR Analysis (comb-p)
7. Annotation of DMCs and DMRs
8. Export to Excel

Added a change