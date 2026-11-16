# ==============================================================================
# 11_cell_counts.R
# Number of cells per cell type.
# Produces: Fig 2D
# Outputs : results/coloc/cellcounts.csv
# Run     : Rscript scripts/01_singlecell_annotation/11_cell_counts.R   (from the repository root)
# ==============================================================================
suppressMessages(library(Seurat))
a <- readRDS('Projects/GSE162798_RAW/results/singlecell/atac_ATAC_RNAlabeled.rds')
tb <- table(a$celltype_ATAC)
write.csv(as.data.frame(tb), 'results/coloc/cellcounts.csv', row.names=FALSE)
print(tb)
cat('total', ncol(a), '\n')
