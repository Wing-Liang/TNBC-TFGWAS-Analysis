#!/usr/bin/env Rscript
# ==============================================================================
# 06_celltype_composition.R
# Cell-type composition of the two cohorts.
# Produces: Fig S9
# Inputs  : results/singlecell/atac_motif.rds
#           replication/celltypes/cellxgene_celltypes_replication.csv
# Run     : Rscript scripts/06_dissociation_ap1/06_celltype_composition.R   (from the repository root)
# ==============================================================================
# 06_celltype_composition.R -- cell-type composition of the two cohorts, side by side.
# AP-1/TRE is the dominant motif of mesenchymal and basal chromatin, so a shift in the
# stromal fraction is the leading remaining explanation for the TRE deficit.
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({ library(Seurat); library(Signac); library(data.table) })
setwd(TNBC_ROOT)
ts <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

ts("loading discovery object")
o <- readRDS("results/singlecell/atac_motif.rds")
md <- o@meta.data
ts("metadata columns:")
print(names(md))
for (k in names(md)) {
  v <- md[[k]]
  if (is.character(v) || is.factor(v)) {
    u <- unique(as.character(v))
    if (length(u) >= 2 && length(u) <= 30) {
      cat("\n---", k, "---\n"); print(sort(table(as.character(v)), decreasing = TRUE))
    }
  }
}
cat("\nIdents:\n"); print(sort(table(as.character(Idents(o))), decreasing = TRUE))
cat("\ntotal discovery cells:", ncol(o), "\n")

cat("\n================ REPLICATION ================\n")
R <- fread("replication/celltypes/cellxgene_celltypes_replication.csv")
tb <- sort(table(R$author_cell_type), decreasing = TRUE)
print(tb)
cat("total replication cells:", nrow(R), "\n")
cat("\nreplication percentages:\n")
print(round(100 * tb / sum(tb), 1))
