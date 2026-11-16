# ==============================================================================
# 06_module_score_umap.R
# UMAP of the TNBC-specific motif module score.
# Produces: Fig 3D
# Inputs  : results/singlecell/final_cell_barcodes.csv
#           results/singlecell/atac_motif.rds
#           results/motif_sets/AddModuleScore_per_cell_all_sets.csv
# Run     : Rscript scripts/04_celltype_enrichment/06_module_score_umap.R   (from the repository root)
# ==============================================================================
## 06_module_score_umap.R — redraw the TNBC-specific integrated UMAP with the module-score gradient ending in
## the TNBC colour of the Fig 3A UpSet plot (#E64B35), instead of orange.
## Same embedding, hull, cutoffs and floor-at-0 as 05_motif_sets_upset_heatmap.R umap_plot(); the per-cell TNBC-specific
## module score is READ from the saved table (results/motif_sets/AddModuleScore_per_cell_all_sets.csv, column
## TNBC_specific) -- chromVAR is NOT recomputed, so the values are identical to the published panel.
## Needs Signac (loads the Seurat/Signac object for the UMAP coordinates and hull).
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressMessages({library(Seurat); library(Signac); library(data.table); library(ggplot2); library(fmesher)})
setwd(TNBC_ROOT)
OUTF <- "results/motif_sets/figures"; UMAP_RED <- "umap.rna"
HI_COL <- "#E64B35"                                        # TNBC colour, Fig 3A UpSet set_cols[["TNBC"]]
ct_levels <- c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "T", "B", "Myl", "Pc")
epi <- c("LASP-1", "LASP-2", "LHS", "B-Myo")

## ---- object + embedding + hull (as in 05_motif_sets_upset_heatmap.R) ----
Final_cell_barcode <- fread("results/singlecell/final_cell_barcodes.csv", header = TRUE)
obj <- readRDS("results/singlecell/atac_motif.rds")
obj <- subset(obj, cells = Final_cell_barcode$x)
new.cluster.ids <- c("LASP-1", "Fb", "LHS", "LASP-2", "B-Myo", "Endo", "T", "Pc", "Myl", "B")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
Idents(obj) <- factor(Idents(obj), levels = ct_levels)
epi_cells <- WhichCells(obj, idents = epi)

emb  <- Embeddings(obj, UMAP_RED)[epi_cells, 1:2]
hull <- fmesher::fm_nonconvex_hull_inla(as.matrix(emb), convex = -0.0021)
hull_seg <- data.frame(x = hull$loc[hull$idx[, 1], 1], y = hull$loc[hull$idx[, 1], 2],
                       xend = hull$loc[hull$idx[, 2], 1], yend = hull$loc[hull$idx[, 2], 2])
umap_df <- data.table(Embeddings(obj, UMAP_RED)[, 1:2]); setnames(umap_df, c("UMAP_1", "UMAP_2"))
umap_df[, `:=`(cell = colnames(obj), celltype = as.character(Idents(obj)))]
lab_df <- umap_df[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = celltype]
blk <- theme(text = element_text(colour = "black"), axis.text = element_text(colour = "black"),
             axis.title = element_text(colour = "black"), legend.text = element_text(colour = "black"))

## ---- per-cell TNBC-specific module score (saved; not recomputed) ----
sc_dt <- fread("results/motif_sets/AddModuleScore_per_cell_all_sets.csv")
sc <- setNames(sc_dt$TNBC_specific, sc_dt$cell)
stopifnot(all(umap_df$cell %in% names(sc)))

## ---- umap_plot: identical to 05_motif_sets_upset_heatmap.R except cols default = c("grey90", HI_COL) ----
umap_plot <- function(sc, cells = NULL, title = "", hull = TRUE, legend = TRUE,
                      cols = c("grey90", HI_COL), floor0 = TRUE, min.cutoff = 0.05, max.cutoff = 0.95, pt.size = 0.45) {
  d <- copy(umap_df); d[, score := sc[cell]]
  if (!is.null(cells)) d <- d[cell %in% cells]
  q  <- quantile(d$score, c(min.cutoff, max.cutoff))
  lo <- if (floor0) 0 else q[[1]]
  d[, score := pmin(pmax(score, lo), q[[2]])]
  setorder(d, score)
  p <- ggplot(d, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = score), size = pt.size, stroke = 0) +
    scale_colour_gradientn(colours = cols, limits = c(lo, q[[2]]), name = "Module\nscore") +
    geom_text(data = lab_df[celltype %in% unique(d$celltype)], aes(label = celltype), size = 3.5, fontface = "bold") +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) + theme_classic(base_size = 11) + blk +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  if (hull) p <- p + geom_segment(data = hull_seg, aes(x = x, y = y, xend = xend, yend = yend),
                                  colour = "black", linetype = "dotted", linewidth = 0.5)
  if (!legend) p <- p + theme(legend.position = "none")
  p
}
## integrated UMAP, title = "" (no "TNBC_specific" title to crop)
f <- file.path(OUTF, "UMAP_TNBC_specific_AddModuleScore_all_TNBCcol.jpg")
ggsave(f, umap_plot(sc, title = ""), width = 6, height = 5, dpi = 300, bg = "white")
cat("wrote", f, "| high colour", HI_COL, "\n")
