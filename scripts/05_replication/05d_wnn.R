#!/usr/bin/env Rscript
# ==============================================================================
# 05d_wnn.R
# Import the published weighted-nearest-neighbour embedding.
# Produces: replication object
# Outputs : wnn/WNN04_imported_vs_recomputed.png
#           wnn/WNN03.1_ATACweight.png
#           wnn/WNN01_DimPlot_RNA_ATAC_WNN.png
#           wnn/WNN02_DimPlot_WNN_bypool.png
#           wnn/WNN03_DimPlot_Doublets.png
# Run     : Rscript scripts/05_replication/05d_wnn.R   (from the repository root)
# ==============================================================================
# 05d_wnn.R — weighted nearest neighbour (RNA + ATAC) joint embedding figures
#             European replication cohort (GSE244585 pools 9/11/12)
#
# RUN ORDER: 05_QC -> 05b_RNA_annotation -> 05c_ATAC_annotation -> 05d_WNN -> 06_chromvar_matrix
#
# ============ WNN IS IMPORTED, NOT RE-DERIVED ============
# The CELLxGENE h5ad already contains the authors' joint embedding and modality weights:
#     obsm X_wnn.umap    -> attached by 05b as reduction `wnn.umap`
#     obs  ATAC.weight   -> attached by 05b as `ATAC.weight`
# So this script does NOT run SCTransform/PCA, TF-IDF/SVD, or FindMultiModalNeighbors. All
# three would recompute, at real cost, something already present and published -- and the
# imported version is the one behind the source publication's figures.
#
# WHAT REMAINS: the figures and the modality-weight summary. ATAC.weight is the quantity
# that says how much chromatin (vs expression) drives the joint embedding per cell; it is
# the main quantitative result here, and it is now read rather than recomputed.
#
# CELL TYPES / CLUSTERS ARE NOT RE-DERIVED. The discovery script clusters the WNN graph
# (FindClusters graph.name="wsnn", algorithm=3, resolution=0.06) and assigns celltype_WNN
# via a hard-coded RenameIdents map. Not done here: cell types come from Idents(), set in 04
# from CELLxGENE author_cell_type. Note the h5ad also carries the authors' own WNN clusters
# as cxg_wsnn_res.0.6 / cxg_wsnn_res.0.8 if a partition is ever wanted -- again, no need to
# recompute one.
#
# Set RECOMPUTE_WNN <- TRUE to derive the joint embedding here instead, with the discovery
# parameters (pca 1:20 + lsi 2:30). Stored as wnn.umap.recomputed and plotted against the
# imported one.
#
# Note: 01_singlecell_annotation/06_wnn_integration.R sets future.globals.maxSize with 1024^3,
# a much larger limit than the 1024^2 used elsewhere; it only raises a memory ceiling.
#
# ============================= DATA SOURCES (imports) =============================
#  figures/atac/replication_ATAC_RNAlabeled.rds  <- 05c_annotate_atac.R
#      (RNA + ATAC + ACTIVITY; umap.rna / umap.atac / wnn.umap / lsi; celltype;
#       ATAC.weight; cxg_* authors' obs metadata; scDblFinder_score)
# Outputs -> wnn/ : figures WNN01-WNN04, replication_WNN.rds
# ==================================================================================

## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({
  library(Signac); library(Seurat); library(GenomicRanges); library(GenomeInfoDb)
  library(ggplot2); library(patchwork); library(data.table); library(dplyr)
})
setwd(file.path(TNBC_ROOT, "replication"))
dir.create("wnn", showWarnings = FALSE)
set.seed(1234)
ts <- function(m) message("[", format(Sys.time(), "%H:%M:%S"), "] ", m)
options(future.globals.maxSize = 8000 * 1024^2)

## ---- configuration ----
## Set via the environment so the default run is untouched:
##   sbatch --export=ALL,RECOMPUTE_WNN=TRUE 05d_WNN_recompute.sbatch
## Default FALSE = plot the authors' imported joint embedding, as before.
RECOMPUTE_WNN <- isTRUE(as.logical(Sys.getenv("RECOMPUTE_WNN", "FALSE")))
RNA_DIMS      <- 1:20   # only used when RECOMPUTE_WNN is TRUE
LSI_DIMS      <- 2:30   # ditto

colors <- c("#E64B35CC","#4DBBD5CC","#00A087CC","#B09C85CC","#F39B7FCC",
            "#8491B4CC","#91D1C2CC","#f9e79f","#7E6148CC","#3C5488CC")
pal <- function(n) rep(colors, length.out = n)

# ─────────────────────────────────────────────
#  1. Load 05c output
# ─────────────────────────────────────────────
IN <- "figures/atac/replication_ATAC_RNAlabeled.rds"
if (!file.exists(IN)) stop("Not found: ", IN, "\n  Run 05c_annotate_atac.R first.")
ts(paste("loading", IN))
obj <- readRDS(IN)
ts(paste("loaded:", ncol(obj), "cells; reductions:",
         paste(Reductions(obj), collapse = ", ")))

stopifnot("celltype" %in% colnames(obj@meta.data))
Idents(obj) <- "celltype"
n_ct <- nlevels(factor(obj$celltype))
print(table(obj$celltype))

if (!"wnn.umap" %in% Reductions(obj))
  stop("wnn.umap absent -- 05b did not attach the imported embeddings.\n",
       "  Check data/cellxgene/embeddings_replication.csv and re-run 05b, or set RECOMPUTE_WNN.")

# ─────────────────────────────────────────────
#  2. Optional: derive the joint embedding here too
# ─────────────────────────────────────────────
if (RECOMPUTE_WNN) {
  ts("RECOMPUTE_WNN=TRUE - deriving the joint embedding locally")
  if (!"pca" %in% Reductions(obj)) {
    ts("  RNA: SCTransform -> RunPCA")
    DefaultAssay(obj) <- "RNA"
    obj <- SCTransform(obj, verbose = FALSE); obj <- RunPCA(obj, verbose = FALSE)
  } else ts("  reusing 'pca'")
  if (!"lsi" %in% Reductions(obj)) {
    ts("  ATAC: RunTFIDF -> FindTopFeatures(q5) -> RunSVD")
    DefaultAssay(obj) <- "ATAC"
    obj <- RunTFIDF(obj, verbose = FALSE)
    obj <- FindTopFeatures(obj, min.cutoff = "q5")
    obj <- RunSVD(obj)
  } else ts("  reusing 'lsi'")
  obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"),
                                 dims.list = list(RNA_DIMS, LSI_DIMS), verbose = FALSE)
  obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap.recomputed",
                 reduction.key = "wnnrcUMAP_", verbose = FALSE)
  png("wnn/WNN04_imported_vs_recomputed.png", width = 12, height = 5,
      units = "in", res = 200)
  print((DimPlot(obj, reduction = "wnn.umap", group.by = "celltype", label = TRUE, repel = TRUE) +
           ggtitle("imported (CELLxGENE authors)") + NoLegend() |
         DimPlot(obj, reduction = "wnn.umap.recomputed", group.by = "celltype",
                 label = TRUE, repel = TRUE) + ggtitle("recomputed here") + NoLegend()) +
        plot_annotation(title = "WNN UMAP: authors' embedding vs re-derived"))
  dev.off()
}

# ─────────────────────────────────────────────
#  3. Modality weights
# ─────────────────────────────────────────────
## ATAC.weight in [0,1]: the share of each cell's joint neighbourhood determined by
## chromatin rather than expression. Near 0.5 means both modalities contribute; near 0
## means the joint embedding is effectively the RNA embedding and WNN adds little.
wcol <- intersect(c("ATAC.weight", "cxg_ATAC.weight", "cxg_integrated.weight"),
                  colnames(obj@meta.data))
if (length(wcol)) {
  ts("modality weights by cell type:")
  wt <- sapply(wcol, function(w) tapply(obj@meta.data[[w]], obj$celltype, mean, na.rm = TRUE))
  print(round(wt, 3))
  fwrite(as.data.table(as.data.frame(wt), keep.rownames = "celltype"),
         "wnn/WNN_modality_weights.csv")
  ts("overall mean:")
  print(round(sapply(wcol, function(w) mean(obj@meta.data[[w]], na.rm = TRUE)), 3))

  png("wnn/WNN03.1_ATACweight.png", width = 12, height = 4.5, units = "in", res = 200)
  w1 <- FeaturePlot(obj, features = wcol[1], reduction = "wnn.umap") +
        ggtitle("ATAC modality weight")
  w2 <- VlnPlot(obj, features = wcol[1], group.by = "celltype", pt.size = 0) +
        ggtitle("ATAC weight by cell type") + NoLegend()
  print(w1 + w2); dev.off()
} else ts("no modality weight column found - skipped")

OUT_RDS <- if (RECOMPUTE_WNN) "wnn/replication_WNN_recomputed.rds" else "wnn/replication_WNN.rds"
ts(paste("writing", OUT_RDS))
saveRDS(obj, OUT_RDS)
ts("saved -> wnn/replication_WNN.rds")

# ─────────────────────────────────────────────
#  4. Three embeddings, one label set (WNN01)
# ─────────────────────────────────────────────
png("wnn/WNN01_DimPlot_RNA_ATAC_WNN.png", width = 14, height = 4.5,
    units = "in", res = 300)
p1 <- DimPlot(obj, reduction = "umap.rna",  group.by = "celltype",
              label = TRUE, repel = TRUE) + ggtitle("RNA")  +
      scale_color_manual(values = pal(n_ct)) + NoLegend()
p2 <- DimPlot(obj, reduction = "umap.atac", group.by = "celltype",
              label = TRUE, repel = TRUE) + ggtitle("ATAC") +
      scale_color_manual(values = pal(n_ct)) + NoLegend()
p3 <- DimPlot(obj, reduction = "wnn.umap",  group.by = "celltype",
              label = TRUE, repel = TRUE) + ggtitle("WNN")  +
      scale_color_manual(values = pal(n_ct))
print((p1 + p2 + p3) +
      plot_annotation(title = sprintf(
        "CELLxGENE cell types on the authors' three embeddings (%d cells)", ncol(obj))) &
      theme(plot.title = element_text(hjust = 0.5)))
dev.off()

## by pool: batch structure is easier to see on the joint embedding than on either alone
png("wnn/WNN02_DimPlot_WNN_bypool.png", width = 7, height = 5, units = "in", res = 300)
print(DimPlot(obj, reduction = "wnn.umap", group.by = "pool") +
      ggtitle("WNN - by pool (batch check)"))
dev.off()
ts("cell-type composition within each pool:")
print(round(prop.table(table(celltype = obj$celltype, pool = obj$pool), 2), 3))

# ─────────────────────────────────────────────
#  5. Doublet score on the WNN embedding (WNN03)
# ─────────────────────────────────────────────
if ("scDblFinder_score" %in% colnames(obj@meta.data) && any(!is.na(obj$scDblFinder_score))) {
  png("wnn/WNN03_DimPlot_Doublets.png", width = 11, height = 4.5, units = "in", res = 300)
  d1 <- DimPlot(obj, reduction = "wnn.umap", group.by = "celltype",
                label = TRUE, repel = TRUE) + ggtitle("WNN") +
        scale_color_manual(values = pal(n_ct)) + NoLegend()
  d2 <- FeaturePlot(obj, features = "scDblFinder_score", reduction = "wnn.umap") +
        ggtitle("scDblFinder score") + labs(color = "DoubletScore")
  print(d1 + d2); dev.off()
  ts("mean doublet score by cell type:")
  print(round(tapply(obj$scDblFinder_score, obj$celltype, mean, na.rm = TRUE), 4))
} else ts("no scDblFinder_score - WNN03 skipped")

ts("=== 05d_wnn.R complete ===")
