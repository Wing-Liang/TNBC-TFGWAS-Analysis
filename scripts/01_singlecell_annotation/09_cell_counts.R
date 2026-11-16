# ==============================================================================
# 09_cell_counts.R
# Number of cells per cell type.
# Produces: Fig 2D
# Outputs : results/singlecell/cellcounts.csv
# Run     : Rscript scripts/01_singlecell_annotation/09_cell_counts.R   (from the repository root)
# ==============================================================================
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressMessages(library(Seurat))
setwd(TNBC_ROOT)
a <- readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')
tb <- table(a$celltype_ATAC)
write.csv(as.data.frame(tb), 'results/singlecell/cellcounts.csv', row.names=FALSE)
print(tb)
cat('total', ncol(a), '\n')
