#!/usr/bin/env Rscript
# 05c_annotate_atac.R
# ATAC gene activity and ATAC figures for the replication cohort.
#
# The ATAC UMAP and the LSI (100 components) are imported from the CELLxGENE object, as for RNA
# (05b_annotate_rna.R); cell types are the authors' labels. Gene activity is computed here from
# this cohort's fragments with the hg38 annotation carried in the ChromatinAssay
# (EnsDb.Hsapiens.v86). Set RECOMPUTE_EMBEDDING <- TRUE to also derive TF-IDF -> SVD -> UMAP
# with the discovery parameters (FindTopFeatures "q5", dims 2:30) and plot both.
#
# Inputs : rna/replication_RNA_annotated.rds              (05b_annotate_rna.R)
#          data/cellxgene/lsi_replication.csv              (03b_extract_cellxgene_embeddings.py)
# Outputs: figures/atac/replication_ATAC_RNAlabeled.rds, figures ATAC01-ATAC10 in figures/atac/
# Run    : Rscript scripts/05_replication/05c_annotate_atac.R   (from the repository root)

## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({
  library(Signac); library(Seurat); library(GenomicRanges); library(GenomeInfoDb)
  library(ggplot2); library(patchwork); library(data.table); library(dplyr)
  library(grid)
})
setwd(file.path(TNBC_ROOT, "replication"))
dir.create("figures/atac", showWarnings = FALSE)
set.seed(1234)
ts <- function(m) message("[", format(Sys.time(), "%H:%M:%S"), "] ", m)
options(future.globals.maxSize = 8000 * 1024^2)

## ---- configuration ----
LSI_CSV             <- "data/cellxgene/lsi_replication.csv"
RECOMPUTE_EMBEDDING <- FALSE
LSI_DIMS            <- 2:30      # only used when RECOMPUTE_EMBEDDING is TRUE
TOPFEAT_CUT         <- "q5"      # ditto (discovery value)
GA_MAX_WIDTH        <- 5e+08     # discovery: GeneActivity(max.width = 5e+08, biotypes = NULL)
GA_BIOTYPES         <- NULL

colors <- c("#E64B35CC","#4DBBD5CC","#00A087CC","#B09C85CC","#F39B7FCC",
            "#8491B4CC","#91D1C2CC","#f9e79f","#7E6148CC","#3C5488CC")
pal <- function(n) rep(colors, length.out = n)

# ─────────────────────────────────────────────
#  1. Load 05b output
# ─────────────────────────────────────────────
IN <- "rna/replication_RNA_annotated.rds"
if (!file.exists(IN)) stop("Not found: ", IN, "\n  Run 05b_annotate_rna.R first.")
ts(paste("loading", IN))
# loading rna/replication_RNA_annotated.rds

obj <- readRDS(IN)
ts(paste("loaded:", ncol(obj), "cells;", nrow(obj[["ATAC"]]), "peaks"))
# loaded: 10700 cells; 230195 peaks

ts(paste("reductions present:", paste(Reductions(obj), collapse = ", ")))
# [10:37:48] reductions present: umap.rna, umap.atac, wnn.umap
#
# Adi-1       Adi-2          BM      Endo-1      Endo-2 Fibroblasts
# LASP         LHS Macrophages     T-cells

stopifnot("celltype" %in% colnames(obj@meta.data))
Idents(obj) <- "celltype"
n_ct <- nlevels(factor(obj$celltype))
print(table(obj$celltype))

if (!"umap.atac" %in% Reductions(obj))
  stop("umap.atac absent -- 05b did not attach the imported embeddings.\n",
       "  Check data/cellxgene/embeddings_replication.csv and re-run 05b.")

## celltype_ATAC == the RNA/CELLxGENE labels, as in the discovery cohort. Same nuclei
## (10x Multiome, one barcode per cell), so this is an assignment, not a join.
stopifnot(identical(colnames(obj[["ATAC"]]), colnames(obj[["RNA"]])))
obj$celltype_ATAC <- Idents(obj)

# ─────────────────────────────────────────────
#  2. Import LSI (for DepthCor and for 05d)
# ─────────────────────────────────────────────
if (file.exists(LSI_CSV)) {
  ts(paste("importing LSI from", LSI_CSV))
  # [10:37:48] importing LSI from data/cellxgene/lsi_replication.csv

  lsi <- fread(LSI_CSV, header = TRUE); setkey(lsi, cell)
  if (length(setdiff(colnames(obj), lsi$cell)))
    stop("cells missing from ", LSI_CSV, " -- refusing a partial reduction.")
  l <- lsi[colnames(obj)]
  stopifnot(identical(l$cell, colnames(obj)))
  m <- as.matrix(l[, -1])
  rownames(m) <- colnames(obj); colnames(m) <- paste0("LSI_", seq_len(ncol(m)))
  obj[["lsi"]] <- CreateDimReducObject(embeddings = m, key = "LSI_", assay = "ATAC")
  ts(paste("  lsi:", ncol(m), "components"))
  # [10:37:49]   lsi: 100 components

} else ts(paste("  ", LSI_CSV, "absent - DepthCor skipped"))

## DepthCor: does LSI-1 track sequencing depth? The authors dropped it (dims 2:30) on that
## assumption; this checks it holds for THIS cohort's depth distribution.
if ("lsi" %in% Reductions(obj)) {
  DefaultAssay(obj) <- "ATAC"
  png("figures/atac/ATAC01_DepthCor_ATAC.png", width = 7, height = 4, units = "in", res = 300)
  print(DepthCor(obj)); dev.off()
  dc <- DepthCor(obj)$data
  ts("LSI component correlation with depth (first 5):")

  # [10:37:50] LSI component correlation with depth (first 5):
    # LSI1   LSI2   LSI3   LSI4   LSI5

  print(round(setNames(dc$counts[1:5], paste0("LSI", dc$Component[1:5])), 3))
  if (abs(dc$counts[1]) < 0.5)
    ts("  !! LSI-1 depth correlation is weak here - the authors' dims 2:30 choice may not transfer")
}

## optional: derive the ATAC embedding here as well and compare
if (RECOMPUTE_EMBEDDING) {
  ts("RECOMPUTE_EMBEDDING=TRUE - RunTFIDF -> FindTopFeatures -> RunSVD -> RunUMAP")
  DefaultAssay(obj) <- "ATAC"
  obj <- RunTFIDF(obj, verbose = FALSE)
  obj <- FindTopFeatures(obj, min.cutoff = TOPFEAT_CUT)
  obj <- RunSVD(obj, reduction.name = "lsi.recomputed", reduction.key = "LSIrc_")
  obj <- RunUMAP(obj, reduction = "lsi.recomputed", dims = LSI_DIMS,
                 reduction.name = "umap.atac.recomputed",
                 reduction.key = "atacrcUMAP_", verbose = FALSE)
  png("figures/atac/ATAC06_imported_vs_recomputed.png", width = 12, height = 5,
      units = "in", res = 200)
  print((DimPlot(obj, reduction = "umap.atac", group.by = "celltype", label = TRUE, repel = TRUE) +
           ggtitle("imported (CELLxGENE authors)") + NoLegend() |
         DimPlot(obj, reduction = "umap.atac.recomputed", group.by = "celltype",
                 label = TRUE, repel = TRUE) + ggtitle("recomputed here") + NoLegend()) +
        plot_annotation(title = "ATAC UMAP: authors' embedding vs re-derived"))
  dev.off()
}

# ─────────────────────────────────────────────
#  3. Figures on the imported ATAC embedding
# ─────────────────────────────────────────────
png("figures/atac/ATAC02_DimPlot_ATAC_celltype.png", width = 7, height = 5, units = "in", res = 200)
print(DimPlot(obj, reduction = "umap.atac", group.by = "celltype", label = TRUE, repel = TRUE) +
      ggtitle(sprintf("ATAC - CELLxGENE embedding + labels (%d cells)", ncol(obj))) +
      scale_color_manual(values = pal(n_ct)) + NoLegend())
dev.off()

## batch check: a region of the embedding that is ~one pool is an artifact, not biology.
## Relevant because of pool11's truncated fragments (05a_qc.R section 0b).
png("figures/atac/ATAC02.1_DimPlot_ATAC_bypool.png", width = 7, height = 5, units = "in", res = 200)
print(DimPlot(obj, reduction = "umap.atac", group.by = "pool") +
      ggtitle("ATAC - by pool (batch check)"))
dev.off()

ts("cell-type composition within each pool:")
# [10:37:52] cell-type composition within each pool:
#   pool
# celltype      pool11 pool12 pool9
# Adi-1        0.017  0.021 0.015
# Adi-2        0.003  0.007 0.005
# BM           0.277  0.419 0.384
# Endo-1       0.027  0.033 0.047
# Endo-2       0.052  0.038 0.024
# Fibroblasts  0.052  0.107 0.086
# LASP         0.373  0.211 0.272
# LHS          0.135  0.120 0.126
# Macrophages  0.025  0.017 0.015
# T-cells      0.039  0.028 0.024

print(round(prop.table(table(celltype = obj$celltype, pool = obj$pool), 2), 3))
fwrite(as.data.table(as.data.frame.matrix(table(obj$celltype, obj$pool)),
                     keep.rownames = "celltype"),
       "figures/atac/ATAC_celltype_by_pool.csv")

## Discovery's ATAC03.0/03.1/03.2 differed only in which label field coloured them -- but
## celltype_ATAC and celltype_RNA hold the same values, so all three were identical. One
## comparison is drawn instead: same cells, same labels, RNA embedding vs ATAC embedding.
## Where the two panels disagree is the informative part.
png("figures/atac/ATAC03_DimPlot_RNA_vs_ATAC.png", width = 11, height = 4.5,
    units = "in", res = 200)
p_atac <- DimPlot(obj, reduction = "umap.atac", group.by = "celltype",
                  label = TRUE, label.size = 4, repel = TRUE) +
          ggtitle("ATAC UMAP") + scale_color_manual(values = pal(n_ct))
if ("umap.rna" %in% Reductions(obj)) {
  p_rna <- DimPlot(obj, reduction = "umap.rna", group.by = "celltype",
                   label = TRUE, label.size = 4, repel = TRUE) +
           ggtitle("RNA UMAP") + scale_color_manual(values = pal(n_ct))
  print((p_rna | p_atac) + plot_layout(guides = "collect") +
        plot_annotation(title = "Same cells, same CELLxGENE labels, two modalities"))
} else print(p_atac)
dev.off()

png("figures/atac/ATAC04_DimPlot_ATAC.png", width = 6.5, height = 5, units = "in", res = 300)
print(DimPlot(obj, reduction = "umap.atac", group.by = "celltype_ATAC",
              label = TRUE, pt.size = 0.5) + ggtitle("ATAC") +
      scale_color_manual(values = pal(n_ct)) + NoLegend())
dev.off()

# ─────────────────────────────────────────────
#  4. Gene activity -> ACTIVITY assay
# ─────────────────────────────────────────────
## The one step here that CANNOT be imported: computed from this cohort's fragments against
## the hg38 annotation in the ChromatinAssay. Slow.
DefaultAssay(obj) <- "ATAC"
ts(paste0("GeneActivity(max.width=", format(GA_MAX_WIDTH, scientific = TRUE),
          ", biotypes=NULL) - slow"))
# [10:37:54] GeneActivity(max.width=5e+08, biotypes=NULL) - slow
# Extracting gene coordinates
# Extracting reads overlapping genomic regions
#   In SingleFeatureMatrix(fragment = fragments[[x]], features = features,  :
#                            1135 features are on seqnames not present in the fragment file. These will be removed.

gene.activities <- GeneActivity(obj, max.width = GA_MAX_WIDTH,
                                biotypes = GA_BIOTYPES, assay = "ATAC")
ts(paste("gene activity matrix:", nrow(gene.activities), "genes x",
         ncol(gene.activities), "cells"))
# [11:02:30] gene activity matrix: 27388 genes x 10700 cells
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

obj[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(obj) <- "ACTIVITY"
obj <- NormalizeData(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
rm(gene.activities); invisible(gc())

saveRDS(obj, "figures/atac/replication_ATAC_RNAlabeled.rds")
ts("saved -> figures/atac/replication_ATAC_RNAlabeled.rds")

# ─────────────────────────────────────────────
#  5. Doublet score on the ATAC embedding (ATAC05)
# ─────────────────────────────────────────────
if ("scDblFinder_score" %in% colnames(obj@meta.data) && any(!is.na(obj$scDblFinder_score))) {
  png("figures/atac/ATAC05_DimPlot_ATAC_Doublet.png", width = 11, height = 4.5,
      units = "in", res = 300)
  p14 <- DimPlot(obj, reduction = "umap.atac", group.by = "celltype",
                 label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ATAC") +
         scale_color_manual(values = pal(n_ct)) + NoLegend()
  p15 <- FeaturePlot(obj, features = "scDblFinder_score", reduction = "umap.atac") +
         ggtitle("scDblFinder score") + labs(color = "DoubletScore")
  print(p14 + p15); dev.off()
  ts("mean doublet score by cell type:")
  print(round(tapply(obj$scDblFinder_score, obj$celltype, mean, na.rm = TRUE), 4))
} else ts("no scDblFinder_score - ATAC05 skipped (run 05a_qc.R section 1b first)")
# [11:10:12] mean doublet score by cell type:
#   Adi-1       Adi-2          BM      Endo-1      Endo-2 Fibroblasts
# LASP         LHS Macrophages     T-cells

# ---------------------------------------------
#  6. Figure 5D -- RNA expression vs Gene Activity, side by side
#      Ported from the discovery cohort's ATAC10_Dotplot_Fig5D.
#
#      Two DotPlots on the SAME genes and the SAME cells, faceted:
#                                from GeneActivity in section 4)
#      Where the two panels disagree is the point of the figure: accessible-but-silent
#      genes are poised, expressed-but-inaccessible ones argue the activity score is
#      missing distal regulation.
#
#      NOTE the discovery run warned that 11 genes were absent from its ACTIVITY assay
#      (LAMA2, HPSE2, PRKG1, RGS6, ZNF385D, MCTP1, PTPRT, NRG1, AFF3, INPP4B, ...).
#      GeneActivity here used biotypes=NULL and dropped 1,135 features on seqnames absent
#      from the autosome-only fragment file, so the same gaps are expected. Missing genes
#      are reported per assay rather than silently omitted.
# ---------------------------------------------
FIG5D_FEATURES <- rev(c("ANKRD30A","ERBB4","SYTL2","INPP4B","AFF3",
                        "IGF2BP2","ALDH1A3","COBL","PIGR","CCL28",
                        "NRG1","KRT14","ACTA2","PTPRT","MYLK",
                        "PTPRC","PARP8","FYN","STAT4","AOAH",
                        "ENAM","SSR4","MZB1","DERL3","HERPUD1",
                        "HLA-DRA","CD74","HLA-DPA1","HLA-DRB1","HLA-DPB1",
                        "MCTP1","ADAMTS9","ZNF385D","ADGRL4","SELE",
                        "C11orf96","MT1A","RGS6","IGFBP5","PRKG1",
                        "DCN","APOD","HPSE2","CFD","LAMA2"))

ts("Fig5D: RNA expression vs gene activity")
dot_data <- function(assay, label) {
  if (!assay %in% Assays(obj)) { ts(paste("  assay", assay, "absent - skipped")); return(NULL) }
  DefaultAssay(obj) <- assay
  present <- FIG5D_FEATURES[FIG5D_FEATURES %in% rownames(obj[[assay]])]
  miss    <- setdiff(FIG5D_FEATURES, present)
  if (length(miss)) ts(paste0("  ", assay, ": ", length(miss), " genes absent -- ",
                              paste(miss, collapse = ", ")))
  if (!length(present)) return(NULL)
  d <- DotPlot(obj, features = present, group.by = "celltype")$data
  d <- d[!is.na(d$features.plot), ]
  d$category <- label
  d
}
d_rna  <- dot_data("RNA",      "RNA expression")
d_act  <- dot_data("ACTIVITY", "Gene Activity (ATAC)")

if (!is.null(d_rna) || !is.null(d_act)) {
  combined <- rbind(d_rna, d_act)
  combined$category <- factor(combined$category,
                              levels = c("RNA expression", "Gene Activity (ATAC)"))
  ## keep the gene order of FIG5D_FEATURES on the y axis rather than alphabetical
  combined$features.plot <- factor(as.character(combined$features.plot),
                                   levels = FIG5D_FEATURES)

  p <- ggplot(combined, aes(x = id, y = features.plot)) +
       geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
       scale_color_gradient2(low = "plum3", mid = "grey90", high = "red", midpoint = 0) +
       scale_size(range = c(0.3, 6)) +
       facet_grid(. ~ category, scales = "free", space = "free") +
       theme_minimal() +
       theme(axis.ticks       = element_line(color = "black"),
             axis.text        = element_text(color = "black"),
             axis.text.x      = element_text(angle = 45, hjust = 1),
             axis.title       = element_text(color = "black"),
             strip.text       = element_text(size = 14, face = "bold", color = "black"),
             panel.grid       = element_blank(),
             panel.spacing    = unit(0.5, "lines"),
             panel.background = element_blank(),
             panel.border     = element_blank(),
             legend.position  = "bottom") +
       labs(x = NULL, y = NULL)

  p <- p +
    annotation_custom(grid::linesGrob(y = c(1, 1), x = c(0, 1),
                                      gp = grid::gpar(col = "black", lwd = 1))) +
    annotation_custom(grid::linesGrob(y = c(0, 0), x = c(0, 1),
                                      gp = grid::gpar(col = "black", lwd = 1))) +
    annotation_custom(grid::linesGrob(y = c(0, 1), x = c(0, 0),
                                      gp = grid::gpar(col = "black", lwd = 1)))

  png("figures/atac/ATAC10_DotPlot_Fig5D.png",
      width = 8, height = 12, units = "in", res = 200)
  print(p); dev.off()

  ## numerical source data (journal requirement)
  fwrite(data.table(gene            = as.character(combined$features.plot),
                    cell_type       = as.character(combined$id),
                    assay           = as.character(combined$category),
                    mean_expression = round(combined$avg.exp, 4),
                    mean_expression_scaled = round(combined$avg.exp.scaled, 4),
                    pct_expressing  = round(combined$pct.exp, 2)),
         "figures/atac/ATAC10_DotPlot_Fig5D_source_data.csv")
  ts(paste("  wrote figures/atac/ATAC10_DotPlot_Fig5D.png +", "_source_data.csv"))
  DefaultAssay(obj) <- "ATAC"
} else ts("  no data for either assay - Fig5D skipped")

ts("=== 05c_annotate_atac.R complete ===")
