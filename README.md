# Single-Cell TCR-Seq Data Analysis

This repository provides a full pipeline for analyzing paired single-cell RNA sequencing (scRNA-seq) and T-cell receptor sequencing (scTCR-seq) data using the R language.

## 📦 Overview

- Load 10x Genomics scRNA-seq and scTCR-seq data.
- Integrate TCR clonotype data with transcriptome data.
- Visualize clonal expansion and gene usage.
- Compute diversity metrics and clone distributions.

## 📁 Input Files

Place the following in the root or adjust paths:
- `sample1/outs/filtered_feature_bc_matrix/` (from `cellranger count`)
- `sample1_TCR/outs/filtered_contig_annotations.csv` (from `cellranger vdj`)

## ▶️ How to Run

1. **Install Required R Packages**:

```r
install.packages("Seurat")
remotes::install_github("ncborcherding/scRepertoire")
