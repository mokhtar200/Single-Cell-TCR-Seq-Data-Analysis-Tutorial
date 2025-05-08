# scTCRseq_analysis.R
# Single-cell TCR-Seq + RNA-Seq Data Analysis

# Load libraries
library(Seurat)
library(scRepertoire)
library(ggplot2)

# Define paths
rna_path <- "sample1/outs/filtered_feature_bc_matrix/"
tcr_path <- "sample1_TCR/outs/filtered_contig_annotations.csv"

# Load and process scRNA-Seq data
rna_counts <- Read10X(data.dir = rna_path)
seurat_obj <- CreateSeuratObject(counts = rna_counts, project = "TCR_SingleCell")

# Basic scRNA preprocessing
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Load and clean TCR data
tcr <- read.csv(tcr_path)
tcr$barcode <- gsub("-1", "", tcr$barcode)

# Combine TCR into list for scRepertoire
combined_tcr <- combineTCR(list(tcr), samples = "Sample1", ID = "TCRdata")

# Merge TCR with Seurat object
seurat_obj <- combineExpression(combined_tcr, seurat_obj, cloneCall = "gene+nt", groupBy = "seurat_clusters")

# Visualizations
pdf("TCR_Clonotype_Plots.pdf")

DimPlot(seurat_obj, group.by = "CTgene", label = TRUE, pt.size = 1) + ggtitle("TCR Gene Usage")
VlnPlot(seurat_obj, features = "cloneType", group.by = "seurat_clusters") + ggtitle("Clone Type by Cluster")
highlightClonotypes(seurat_obj, cloneCall = "gene+nt", sequence = head(unique(tcr$cdr3), 3))

# Diversity and expansion
clonalSpaceHomeostasis(combined_tcr, cloneCall = "gene+nt")
clonalHomeostasis(combined_tcr)
clonalDiversity(combined_tcr, cloneCall = "gene+nt")
clonalOverlap(combined_tcr, cloneCall = "gene+nt", method = "overlap")

dev.off()

# Export results
write.csv(combined_tcr[[1]], file = "TCR_clonotypes_sample1.csv", row.names = FALSE)
saveRDS(seurat_obj, file = "seurat_with_TCR.rds")

message("✅ Analysis complete. Outputs: TCR_clonotypes_sample1.csv, seurat_with_TCR.rds, and TCR_Clonotype_Plots.pdf")
