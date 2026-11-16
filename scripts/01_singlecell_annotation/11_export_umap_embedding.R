## 11_export_umap_embedding.R
## Export the RNA UMAP coordinates and cell type of every post-QC discovery cell, so that the
## module-score and gene-set UMAPs can be drawn without reloading the Seurat object.
##
## Inputs   results/singlecell/atac_motif.rds                 (10_motif_annotation_chromvar.R)
##          results/singlecell/final_cell_barcodes.csv
## Outputs  results/motif_sets/umap_embedding_celltype.csv    cell, UMAP_1, UMAP_2, celltype
## Run from the repository root:  Rscript scripts/01_singlecell_annotation/11_export_umap_embedding.R

TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")

suppressMessages({ library(Seurat); library(Signac); library(data.table) })

setwd(TNBC_ROOT)
dir.create("results/motif_sets", recursive = TRUE, showWarnings = FALSE)

ct_levels <- c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "T", "B", "Myl", "Pc")
Final_cell_barcode <- fread("results/singlecell/final_cell_barcodes.csv", header = TRUE)
obj <- readRDS("results/singlecell/atac_motif.rds")
obj <- subset(obj, cells = Final_cell_barcode$x)

## manuscript cell-type labels
new.cluster.ids <- c("LASP-1", "Fb", "LHS", "LASP-2", "B-Myo", "Endo", "T", "Pc", "Myl", "B")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
Idents(obj) <- factor(Idents(obj), levels = ct_levels)

d <- data.table(Embeddings(obj, "umap.rna")[, 1:2]); setnames(d, c("UMAP_1", "UMAP_2"))
d[, `:=`(cell = colnames(obj), celltype = as.character(Idents(obj)))]
fwrite(d, "results/motif_sets/umap_embedding_celltype.csv")
cat("wrote results/motif_sets/umap_embedding_celltype.csv |", nrow(d), "cells\n")
