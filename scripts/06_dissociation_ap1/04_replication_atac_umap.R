#!/usr/bin/env Rscript
# ==============================================================================
# 04_replication_atac_umap.R
# ATAC UMAP of the replication cohort.
# Produces: Fig S9
# Inputs  : replication/qc/replication_postQC.rds
#           replication/data/cellxgene/embeddings_replication.csv
# Run     : Rscript scripts/06_dissociation_ap1/04_replication_atac_umap.R   (from the repository root)
# ==============================================================================
# 04_replication_atac_umap.R -- a self-contained ATAC UMAP for the replication cohort
#
# WHY THIS IS NEEDED
# The replication panel of AP1_TRE_UMAP was drawn on data/cellxgene/embeddings_replication.csv, which
# holds the AUTHORS' embedding computed on all 51,367 cells of the 16-pool atlas. Plotting
# only the 10,700 European ATAC cells shows ~21% of the cells that defined that layout, so
# clusters look sparse and interleaved -- the "blur". 03b_extract_cellxgene_embeddings.py warns about
# exactly this: "the layout reflects the full atlas ... not the same as embedding this subset
# independently."
#
# Nothing already on disk fixes it. replication_WNN_recomputed.rds looks like it should, but
# 05d_wnn.R recomputes LSI only `if (!"lsi" %in% Reductions(obj))`, and the CELLxGENE import
# had already supplied `lsi` (X_lsi, full atlas). So its `lsi`, `pca`, `umap.atac` and even
# `wnn.umap.recomputed` all rest on full-atlas components.
#
# SOURCE OBJECT: qc/replication_postQC.rds -- the post-QC object, on the author's
# instruction and because it is the cleanest starting point available:
#   * exactly the 10,700 cells of tfgwas/replication_chromvar.rds (verified 100% overlap),
#     so the embedding covers precisely the cells scored in the figure;
#   * the 230,195-peak ATAC assay the TF-GWAS itself ran on;
#   * NO reductions of any kind, so no imported full-atlas component can leak in, unlike
#     replication_WNN_recomputed.rds which carries the authors' lsi/pca.
# NOTE the replication cohort is 10,700 cells and always was. The 51,367 figure is the size
# of the FULL published 16-pool atlas that the authors' embedding was fitted on -- it is not
# a count of our cohort.
#
# This derives the embedding from the European ATAC data alone -- TFIDF, top features, SVD,
# UMAP -- which is what the discovery panel already is (its umap.atac was computed on its own
# 16,451 cells). Only then are the two panels like-for-like.
#
# NOTE the embedding is used for DISPLAY only. No quantitative claim in this workstream
# depends on it: the AP-1 result rests on peak-level counts and on cell-type composition,
# both of which are embedding-independent.
#
# INPUT   replication/qc/replication_postQC.rds   (10,700 cells, ATAC assay, celltype labels)
# OUTPUT  replication/tfgwas_grouped/replication_atac_umap.csv
#           cell, UMAP1, UMAP2, celltype

## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(data.table)
})
setwd(TNBC_ROOT)
ts <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
set.seed(42)
OUT <- "replication/tfgwas_grouped/replication_atac_umap.csv"

ts("loading qc/replication_postQC.rds")
o <- readRDS("replication/qc/replication_postQC.rds")
ts(sprintf("  %d cells | assays: %s", ncol(o), paste(names(o@assays), collapse = ", ")))
ts(sprintf("  reductions present: %s",
           if (length(names(o@reductions))) paste(names(o@reductions), collapse = ", ")
           else "none (as intended -- nothing imported to inherit)"))
stopifnot(length(names(o@reductions)) == 0)
ts(sprintf("  ATAC features: %d", nrow(o[["ATAC"]])))

DefaultAssay(o) <- "ATAC"
ts("RunTFIDF on the European cells only")
o <- RunTFIDF(o, verbose = FALSE)
ts("FindTopFeatures(q5)")
o <- FindTopFeatures(o, min.cutoff = "q5")
ts(sprintf("  variable features: %d", length(VariableFeatures(o))))
ts("RunSVD -> reduction 'lsi.eur' (independent of the imported full-atlas lsi)")
o <- RunSVD(o, reduction.name = "lsi.eur", reduction.key = "LSIeur_", verbose = FALSE)

## Component 1 of an ATAC LSI almost always tracks sequencing depth; drop it, as the
## discovery pipeline and Signac's own vignettes do.
dc <- DepthCor(o, reduction = "lsi.eur", n = 10)$data
ts("depth correlation by component:")
print(round(setNames(dc$counts, paste0("LSI", dc$Component)), 3))
DIMS <- 2:30

ts(sprintf("RunUMAP on lsi.eur dims %d:%d", min(DIMS), max(DIMS)))
o <- RunUMAP(o, reduction = "lsi.eur", dims = DIMS,
             reduction.name = "umap.atac.eur", reduction.key = "atacEUR_", verbose = FALSE)

em <- Embeddings(o, "umap.atac.eur")[, 1:2, drop = FALSE]
md <- o@meta.data
ct <- if ("celltype" %in% names(md)) as.character(md$celltype) else as.character(Idents(o))
D <- data.table(cell = colnames(o), UMAP1 = em[, 1], UMAP2 = em[, 2], celltype = ct)
fwrite(D, OUT)
ts("written ", OUT)

cat("\n---- cells per type in the new embedding ----\n")
print(D[, .N, by = celltype][order(-N)], row.names = FALSE)

## Sanity: a good embedding separates cell types. Report the share of each cell's 15 nearest
## neighbours that share its label -- a plain, interpretable measure of how crisp the layout
## is, comparable between the old and new embeddings.
suppressPackageStartupMessages(library(FNN))
knn_purity <- function(xy, lab, k = 15) {
  nn <- get.knn(as.matrix(xy), k = k)$nn.index
  mean(rowMeans(matrix(lab[nn], nrow = nrow(nn)) == lab))
}
old <- fread("replication/data/cellxgene/embeddings_replication.csv")
i <- match(D$cell, old$cell)
cat(sprintf("\nkNN label purity (k=15), higher = crisper separation:\n"))
cat(sprintf("  authors' full-atlas umap.atac, European subset : %.3f\n",
            knn_purity(data.frame(old$umap.atac_1[i], old$umap.atac_2[i]), D$celltype)))
cat(sprintf("  European-only ATAC UMAP computed here          : %.3f\n",
            knn_purity(D[, .(UMAP1, UMAP2)], D$celltype)))
