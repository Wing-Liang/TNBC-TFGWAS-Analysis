#!/usr/bin/env Rscript
# 05b_annotate_rna.R
# RNA embedding, cell-type labels and marker genes for the replication cohort.
#
# The RNA UMAP is imported from the CELLxGENE object (the source study's own embedding, computed
# on all 51,367 cells of the atlas) rather than re-derived, so that plots are directly comparable
# with the source publication; post-QC cells are a subset of it and are not re-projected. Cell
# types are the authors' labels (author_cell_type, set in 04_build_combined_object.R); no
# re-clustering is done. celltype_RNA is set equal to celltype, the field name used downstream.
# Set RECOMPUTE_EMBEDDING <- TRUE to also derive SCTransform -> PCA -> UMAP (dims 1:20, as in
# 01_singlecell_annotation/04_annotate_cell_types.R) and plot both.
#
# Inputs : qc/replication_postQC.rds                     (05a_qc.R)
#          data/cellxgene/embeddings_replication.csv      (03b_extract_cellxgene_embeddings.py)
#          data/cellxgene/obs_replication.csv, data/cellxgene/var_replication.csv
# Outputs: rna/replication_RNA_annotated.rds, marker tables and figures RNA01-RNA09 in rna/
# Run    : Rscript scripts/05_replication/05b_annotate_rna.R   (from the repository root)

## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({
  library(Signac); library(Seurat); library(GenomicRanges); library(GenomeInfoDb)
  library(ggplot2); library(patchwork); library(data.table); library(dplyr)
})
setwd(file.path(TNBC_ROOT, "replication"))
dir.create("rna", showWarnings = FALSE)
set.seed(1234)
ts <- function(m) message("[", format(Sys.time(), "%H:%M:%S"), "] ", m)
options(future.globals.maxSize = 8000 * 1024^2)

## ---- configuration ----
EMB_CSV            <- "data/cellxgene/embeddings_replication.csv"
OBS_CSV            <- "data/cellxgene/obs_replication.csv"   # all 69 authors' obs columns
VAR_CSV            <- "data/cellxgene/var_replication.csv"   # gene biotypes (feature_type)
KEEP_BIOTYPES      <- c("protein_coding", "lncRNA")  # NULL = keep all 35,477 genes
RECOMPUTE_EMBEDDING <- FALSE   # TRUE = also derive SCT -> PCA -> UMAP here and compare
PCA_DIMS           <- 1:20     # only used when RECOMPUTE_EMBEDDING is TRUE

## discovery palette (01_singlecell_annotation/04_annotate_cell_types.R)
colors <- c("#E64B35CC","#4DBBD5CC","#00A087CC","#B09C85CC","#F39B7FCC",
            "#8491B4CC","#91D1C2CC","#f9e79f","#7E6148CC","#3C5488CC")
pal <- function(n) rep(colors, length.out = n)

## Shared FeaturePlot rendering settings, used by BOTH 05b and 05b.2 so the same gene
## is drawn identically in every figure. RNA09 and RNA11 previously diverged on all three
## of these (order FALSE/TRUE, pt.size 0.4/0.1, ncol 3/4), which made the same data look
## very different between figures.
##   pt.size : fixed, not NULL -- with NULL, AutoPointSize() derives size from cell count
##             and the same size is drawn onto panels of different pixel dimensions.
##   order   : TRUE draws expressing cells LAST, on top of the grey non-expressing ones.
##             With FALSE, grey cells overpaint real signal and sparse markers (FOXA1)
##             look weaker than they are.
##   ncol    : fixed so per-panel aspect ratio is constant across figures.
FEAT_PT_SIZE <- 0.1
FEAT_ORDER   <- TRUE
FEAT_NCOL    <- 4

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
#  1. Load post-QC object; cell types come from Idents
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
ts("loading qc/replication_postQC.rds")
obj <- readRDS("qc/replication_postQC.rds")
ts(paste("loaded:", ncol(obj), "cells;", nrow(obj[["RNA"]]), "genes"))

stopifnot("celltype" %in% colnames(obj@meta.data))
Idents(obj) <- "celltype"                 # defensive: 04 sets this, assert it here
obj$celltype_RNA <- Idents(obj)           # discovery's field name; same labels
if (any(is.na(obj$celltype))) stop("NA cell types present -- 05b assumes every cell is labelled.")
n_ct <- nlevels(factor(obj$celltype))
ts(paste("cell types (from Idents):", n_ct))
print(table(obj$celltype))
print(table(obj$celltype, obj$pool))

# ---------------------------------------------
#  1b. Restrict the RNA assay to protein_coding + lncRNA
# ---------------------------------------------
## APPLIED HERE, NOT IN 05a_qc.R (decision 2026-08-03). 05_QC computes nCount_RNA and
## nFeature_RNA over ALL 35,477 genes and derives the MAD bounds from those distributions,
## so filtering upstream would shift which cells pass QC and invalidate the 10,700-cell
## result already produced. Filtering after QC leaves every QC number intact and changes
## only marker detection and the expression figures.
##
## The h5ad's feature_type: protein_coding 19,267 + lncRNA 15,492 = 34,759 (98.0%).
## The 718 dropped are mostly IG_*/TR_* segments (~646), pseudogenes (54), 17 'artifact'
## and 1 TEC. NOTE this removes immunoglobulin and T-cell-receptor gene segments, so
## plasma-cell / clonality signal is gone from the marker tables by construction -- the
## lineage markers themselves (PTPRC, CD79A, CD3D) are protein_coding and unaffected.
if (!is.null(KEEP_BIOTYPES) && file.exists(VAR_CSV)) {
  vg <- fread(VAR_CSV, header = TRUE)
  if (!"feature_type" %in% names(vg)) {
    ts(paste("  ", VAR_CSV, "has no feature_type column - biotype filter SKIPPED"))
  } else {
    keep_genes <- vg$feature_name[vg$feature_type %in% KEEP_BIOTYPES]  # feature_name = SYMBOL; rownames are symbols, not ENSG
    keep_genes <- intersect(rownames(obj[["RNA"]]), keep_genes)
    n_before   <- nrow(obj[["RNA"]])
    ts(sprintf("biotype filter (%s): keeping %d of %d genes (%.1f%%)",
               paste(KEEP_BIOTYPES, collapse = "+"), length(keep_genes), n_before,
               100 * length(keep_genes) / n_before))
    print(table(vg$feature_type[!vg$feature_type %in% KEEP_BIOTYPES]))
    obj <- subset(obj, features = c(keep_genes, rownames(obj[["ATAC"]])))
    ts(sprintf("  RNA assay now %d genes; ATAC unchanged at %d peaks",
               nrow(obj[["RNA"]]), nrow(obj[["ATAC"]])))
  }
} else ts("biotype filter not applied (KEEP_BIOTYPES NULL or var CSV absent)")

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
#  2. Import the authors' embeddings
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
if (!file.exists(EMB_CSV))
  stop("Not found: ", EMB_CSV, "\n  Run: python3 03b_extract_cellxgene_embeddings.py")

ts(paste("importing embeddings from", EMB_CSV))
emb <- fread(EMB_CSV, header = TRUE)
setkey(emb, cell)

missing <- setdiff(colnames(obj), emb$cell)
if (length(missing))
  stop(length(missing), " cells have no imported embedding (e.g. ",
       paste(head(missing, 3), collapse = ", "), ").\n",
       "  The h5ad and the object disagree -- do not proceed.")

## subset to this object's cells IN ITS OWN ORDER; a positional join here would silently
## scramble coordinates across cells
e <- emb[colnames(obj)]
stopifnot(identical(e$cell, colnames(obj)))

add_reduc <- function(object, mat, name, key, assay = "RNA") {
  colnames(mat) <- paste0(key, seq_len(ncol(mat)))
  rownames(mat) <- colnames(object)
  object[[name]] <- CreateDimReducObject(embeddings = mat, key = key, assay = assay)
  object
}
obj <- add_reduc(obj, as.matrix(e[, .(`umap.rna_1`, `umap.rna_2`)]),  "umap.rna",  "rnaUMAP_")
obj <- add_reduc(obj, as.matrix(e[, .(`umap.atac_1`, `umap.atac_2`)]), "umap.atac", "atacUMAP_", "ATAC")
obj <- add_reduc(obj, as.matrix(e[, .(`wnn.umap_1`, `wnn.umap_2`)]),   "wnn.umap",  "wnnUMAP_")
if ("ATAC.weight" %in% names(e)) obj$ATAC.weight <- e$ATAC.weight
ts(paste("reductions now present:", paste(Reductions(obj), collapse = ", ")))

## ---- full authors' obs metadata ----
## Every column the authors recorded per cell, PREFIXED `cxg_`. The prefix is not cosmetic:
## ~9 of their names collide with metrics this pipeline computes itself (nCount_RNA,
## nCount_ATAC, nFeature_*, percent.mt, TSS.enrichment, nucleosome_signal, high.tss,
## nucleosome_group). Merging unprefixed would silently overwrite OUR QC values -- the ones
## the applied filter was derived from -- with the authors' equivalents computed on a
## different peak set. Prefixing keeps both, so they can be compared instead.
if (file.exists(OBS_CSV)) {
  ts(paste("importing authors' obs metadata from", OBS_CSV))
  ob <- fread(OBS_CSV, header = TRUE); setkey(ob, cell)
  if (length(setdiff(colnames(obj), ob$cell)))
    stop("cells missing from ", OBS_CSV, " -- refusing a partial metadata merge.")
  o <- ob[colnames(obj)]
  stopifnot(identical(o$cell, colnames(obj)))
  o[, cell := NULL]
  setnames(o, paste0("cxg_", names(o)))
  for (cn in names(o)) obj[[cn]] <- o[[cn]]
  ts(paste("  attached", ncol(o), "columns with cxg_ prefix"))

  ## Sanity: the authors' QC metrics vs ours. They were computed on the authors' peak set,
  ## not the 230,195 peaks 04/05 produced, so agreement is informative but not expected to
  ## be exact.
  for (m in c("nCount_RNA", "TSS.enrichment", "nucleosome_signal", "percent.mt")) {
    ours <- obj@meta.data[[m]]; theirs <- obj@meta.data[[paste0("cxg_", m)]]
    if (!is.null(ours) && !is.null(theirs) && is.numeric(theirs))
      ts(sprintf("  %-18s ours vs cxg_: Spearman rho = %.3f",
                 m, suppressWarnings(cor(ours, theirs, method = "spearman", use = "complete.obs"))))
  }

  ## DONOR FIELDS ARE POOL-LEVEL, NOT PER-CELL. In this deposit donor_id reads
  ## "pooled [K107181,K109128,...]" -- the list of 5 donors in that pool, identical for
  ## every cell in it (3 distinct values across all 12,469 cells). The released object was
  ## NOT genetically demultiplexed, so there is no per-cell donor assignment: donor-level
  ## random effects and per-donor aggregation are NOT possible from this data.
  if ("cxg_donor_id" %in% colnames(obj@meta.data)) {
    nd <- length(unique(obj$cxg_donor_id))
    ts(sprintf("  cxg_donor_id: %d distinct value(s) for %d cells%s",
               nd, ncol(obj),
               if (nd <= 3) "  -- POOL-LEVEL, no per-cell donor assignment" else ""))
  }
} else ts(paste("  ", OBS_CSV, "absent - authors' obs metadata not attached"))

## optional: derive the embedding here as well and show both
if (RECOMPUTE_EMBEDDING) {
  ts("RECOMPUTE_EMBEDDING=TRUE - SCTransform -> RunPCA -> RunUMAP (discovery recipe)")
  DefaultAssay(obj) <- "RNA"
  obj <- SCTransform(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  png("rna/RNA01_ElbowPlot_RNA.png", width = 8, height = 5, units = "in", res = 200)
  print(ElbowPlot(obj, ndims = 50)); dev.off()
  ## discovery's numeric elbow: co1 = first PC with cumulative var > 90% and < 5% own
  ## contribution; co2 = last PC where the drop to the next still exceeds 0.1%
  pct  <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1  <- which(cumu > 90 & pct < 5)[1]
  co2  <- sort(which((pct[1:(length(pct)-1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  ts(sprintf("elbow: co1=%s co2=%s -> pc=%s (USING dims 1:%d)",
             co1, co2, min(co1, co2, na.rm = TRUE), max(PCA_DIMS)))
  obj <- RunUMAP(obj, dims = PCA_DIMS, reduction.name = "umap.rna.recomputed",
                 reduction.key = "rnarcUMAP_", verbose = FALSE)
  png("rna/RNA03_imported_vs_recomputed.png", width = 12, height = 5, units = "in", res = 200)
  print((DimPlot(obj, reduction = "umap.rna", group.by = "celltype", label = TRUE, repel = TRUE) +
           ggtitle("imported (CELLxGENE authors)") + NoLegend() |
         DimPlot(obj, reduction = "umap.rna.recomputed", group.by = "celltype",
                 label = TRUE, repel = TRUE) + ggtitle("recomputed here") + NoLegend()) +
        plot_annotation(title = "RNA UMAP: authors' embedding vs re-derived"))
  dev.off()
}

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
#  3. UMAP figures
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
png("rna/RNA02_DimPlot_RNA_celltype.png", width = 7, height = 5, units = "in", res = 200)
print(DimPlot(obj, reduction = "umap.rna", group.by = "celltype", label = TRUE, repel = TRUE) +
      ggtitle(sprintf("RNA - CELLxGENE embedding + labels (%d cells)", ncol(obj))) +
      scale_color_manual(values = pal(n_ct)) + NoLegend())
dev.off()

png("rna/RNA02.1_DimPlot_RNA_bypool.png", width = 7, height = 5, units = "in", res = 200)
print(DimPlot(obj, reduction = "umap.rna", group.by = "pool") +
      ggtitle("RNA - by pool (batch check)"))
dev.off()

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
#  4. Markers per cell type
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose = FALSE)   # FindAllMarkers needs a data layer
Idents(obj) <- "celltype"

ts("FindAllMarkers (min.diff.pct = 0.5, discovery setting)")
mk_paper <- FindAllMarkers(obj, min.diff.pct = 0.5, verbose = FALSE)
mk_paper <- mk_paper %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1.5)
fwrite(mk_paper, "rna/Markers_RNA_DEGs_perpaper.csv")
ts(paste("  per-paper markers:", nrow(mk_paper)))

ts("FindAllMarkers (unfiltered)")
mk_all <- FindAllMarkers(obj, verbose = FALSE)
fwrite(mk_all, "rna/Markers_RNA_DEGs_all.csv")
ts(paste("  all markers:", nrow(mk_all)))

## top 5 per cell type by effect size.
## NOTE: the discovery script writes `dplyr::filter(adj_p_val, des(avg_log2FC))`, which is
## malformed (no such column, no such function) -- written here as the evidently intended
## arrange(desc(.)) rather than copying the bug.
top5 <- mk_paper %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        slice_head(n = 5) %>% ungroup()
fwrite(top5, "rna/Markers_RNA_top5_per_celltype.csv")
ts("top 5 markers per cell type:")
print(as.data.frame(top5[, c("cluster", "gene", "avg_log2FC", "p_val_adj")]))

## Save the object HERE, before the figure sections. 05c consumes this file; a failure in a
## downstream plot should not throw away the analysis that produced it. (It did once: a
## DoHeatmap error discarded a 30-minute run.)
saveRDS(obj, "rna/replication_RNA_annotated.rds")
ts(paste("saved -> rna/replication_RNA_annotated.rds  (", ncol(obj), "cells )"))

if (nrow(top5)) {
  ## DoHeatmap reads the scale.data layer. NormalizeData does NOT create it, so this used to
  ## die with "No requested features found in the scale.data layer for the RNA assay".
  ## Scale only the genes being plotted -- scaling all 35,477 would be pointless work.
  hm_genes <- unique(top5$gene)
  obj <- ScaleData(obj, features = hm_genes, verbose = FALSE)
  png("rna/RNA06_DoHeatmap_top5gene.png", width = 10, height = 14, units = "in", res = 150)
  print(DoHeatmap(obj, features = hm_genes) + NoLegend()); dev.off()
}

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
#  5. Marker feature plots (discovery panels)
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
feat_plot <- function(genes, file, w, h, label) {
  present <- genes[genes %in% rownames(obj)]
  if (!length(present)) { ts(paste("  none of the", label, "genes present - skipped")); return(invisible()) }
  if (length(present) < length(genes))
    ts(paste("  missing from RNA assay:", paste(setdiff(genes, present), collapse = ", ")))
  png(file, width = w, height = h, units = "in", res = 150)
  # pt.size fixed (not left NULL) so the same gene is rendered identically across
  # panels. With pt.size = NULL Seurat calls AutoPointSize(), which sets size from the
  # cell count; the SAME size is then drawn onto panels of different pixel dimensions,
  # so e.g. FOXA1 looked denser in RNA09 (14x12 in, 7 panels) than in RNA04 (9x4 in,
  # 2 panels) despite identical data and identical colour scales.
  print(FeaturePlot(obj, features = present, reduction = "umap.rna",
                    ncol = min(FEAT_NCOL, length(present)))); dev.off()
}
ts("marker feature plots")
feat_plot(c("FOXA1","TP63"), "rna/RNA04_FeaturePlot_FOXA1_TP63.png", 9, 4, "FOXA1/TP63")
feat_plot(c("ANKRD30A","ERBB4","PROM1","PTN","KRT14","ACTA2","RGS5","PTPRC","CD3E",
            "IL7R","CD79A","JCHAIN","HLA-DRA","VWF","COL1A2","LUM"),
          "rna/RNA05_FeaturePlot_figS5D.png", 18, 13, "figS5D")
feat_plot(c("FOXA1","TP63","EPCAM","KIT","ITGA6","ESR1","GATA3"),
          "rna/RNA09_FeaturePlot_specific_genes.png", 14, 7, "specific genes")

# ─────────────────────────────────────────────
#  5b. Marker dot plots (same marker sets as the feature plots above)
#      DotPlot adds what FeaturePlot cannot show: the FRACTION of cells expressing
#      each marker per cell type, alongside mean expression. Useful for markers that
#      are high but sparse (transcription factors especially).
# ─────────────────────────────────────────────
dot_plot <- function(genes, file, w, h, label) {
  present <- genes[genes %in% rownames(obj)]
  if (!length(present)) { ts(paste("  none of the", label, "genes present - skipped")); return(invisible()) }
  if (length(present) < length(genes))
    ts(paste("  missing from RNA assay:", paste(setdiff(genes, present), collapse = ", ")))
  p <- DotPlot(obj, features = present, group.by = "celltype", dot.scale = 6) +
       ggtitle(paste("Marker expression by cell type --", label)) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
       labs(x = NULL, y = NULL)
  png(file, width = w, height = h, units = "in", res = 200)
  print(p); dev.off()
  # numerical source data for the figure
  write.csv(
    data.frame(marker                 = p$data$features.plot,
               cell_type              = p$data$id,
               mean_expression        = round(p$data$avg.exp, 4),
               mean_expression_scaled = round(p$data$avg.exp.scaled, 4),
               pct_expressing         = round(p$data$pct.exp, 2)),
    sub("\\.png$", "_source_data.csv", file), row.names = FALSE)
}
ts("marker dot plots")
dot_plot(c("FOXA1","TP63","EPCAM","KIT","ITGA6","ESR1","GATA3"),
         "rna/RNA07_DotPlot_specific_genes.png", 8, 5, "specific genes")
dot_plot(c("ANKRD30A","ERBB4","PROM1","PTN","KRT14","ACTA2","RGS5","PTPRC","CD3E",
           "IL7R","CD79A","JCHAIN","HLA-DRA","VWF","COL1A2","LUM"),
         "rna/RNA07.1_DotPlot_figS5D.png", 11, 5, "figS5D")

# ---------------------------------------------
#  5c. Canonical cell-type marker panel
#      MERGED IN from 05b.2_RNA_marker_dimplots.R (2026-08-03); that script is deleted.
#
#      EPCAM   epithelial            ITGA6 (CD49f)  epithelial/basal
#      FOXA1   luminal TF            KIT            luminal progenitor / LASP
#      ESR1    hormone-sensing TF    TP63           basal TF
#      CD3D    T cells               CD79A          B cells
#      ITGAX   myeloid               VWF            endothelial
#      RGS5    pericyte              PDGFRA         fibroblast
#      Only FOXA1, ESR1 and TP63 are transcription factors; the rest are surface/lineage
#      markers.
#
#      NOTE ON RENDERING: 05b.2 called FeaturePlot directly with pt.size = FEAT_PT_SIZE and
#      order = FEAT_ORDER. Here it goes through feat_plot(), the same helper used by
#      RNA04/RNA05/RNA09, where those two arguments are currently commented out. That is
#      deliberate -- one helper means every feature panel in this script renders the same
#      way. Uncomment them in feat_plot() to apply the fixed point size to ALL panels at
#      once, which is what FEAT_PT_SIZE/FEAT_ORDER were introduced for.
# ---------------------------------------------
MARKERS <- c("EPCAM","ITGA6","FOXA1","KIT","ESR1","TP63",
             "CD3D","CD79A","ITGAX","VWF","RGS5","PDGFRA")

ts("canonical marker panel")
feat_plot(MARKERS, "rna/RNA11_FeaturePlot_marker_panel.png", 14, 11, "marker panel")
dot_plot(MARKERS,  "rna/RNA12_DotPlot_marker_panel.png",      9,  5, "marker panel")

## 05b.2 also wrote RNA10_DimPlot_celltype_labelled.png -- a DimPlot of celltype on
## umap.rna with the same palette and NoLegend(). That is the same figure as RNA02 above,
## differing only in the title, so it is NOT reproduced here. RNA02 is the one to cite.

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
#  6. percent.mt cross-check (discovery)
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
mt_gene <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
             "MT-CO1","MT-CO2","MT-CO3","MT-ATP6","MT-ATP8","MT-CYB")
mt_det  <- mt_gene[mt_gene %in% rownames(obj)]
ts(paste("MT genes detected:", length(mt_det), "of 13"))
if (length(mt_det)) {
  obj$percent.mt.byGene    <- PercentageFeatureSet(obj, features = mt_det, assay = "RNA")
  obj$percent.mt.byPattern <- PercentageFeatureSet(obj, pattern  = "^MT-", assay = "RNA")
  print(summary(obj$percent.mt.byGene))
  cv <- suppressWarnings(cor.test(obj$percent.mt.byGene, obj$percent.mt.byPattern,
                                  method = "spearman"))
  ts(sprintf("  %%MT byGene vs byPattern: rho = %.3f", cv$estimate))
  ts("  mean %MT by cell type:")
  print(round(tapply(obj$percent.mt.byGene, obj$celltype, mean, na.rm = TRUE), 3))
}

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
#  7. Doublet score on the RNA UMAP
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
if ("scDblFinder_score" %in% colnames(obj@meta.data) && any(!is.na(obj$scDblFinder_score))) {
  png("rna/RNA08_DimPlot_RNA_Doublet.png", width = 12, height = 4.5, units = "in", res = 200)
  p8 <- DimPlot(obj, reduction = "umap.rna", group.by = "celltype", label = TRUE, repel = TRUE) +
        ggtitle("RNA (cell types)") + scale_color_manual(values = pal(n_ct)) + NoLegend()
  p9 <- FeaturePlot(obj, features = "scDblFinder_score", reduction = "umap.rna") +
        ggtitle("scDblFinder score") + labs(color = "DoubletScore")
  print(p8 + p9); dev.off()
  ts("  mean doublet score by cell type:")
  print(round(tapply(obj$scDblFinder_score, obj$celltype, mean, na.rm = TRUE), 4))
} else ts("no scDblFinder_score in object - RNA08 skipped")

ts("=== 05b_annotate_rna.R complete ===")
