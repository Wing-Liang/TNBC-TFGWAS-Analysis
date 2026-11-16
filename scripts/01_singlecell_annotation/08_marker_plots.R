#!/usr/bin/env Rscript
# ==============================================================================
# 08_marker_plots.R
# Canonical marker feature plots for the scRNA-seq cell types.
# Produces: Fig 2B
# Inputs  : results/singlecell/final_cell_barcodes.csv
#           results/singlecell/atac_ATAC_RNAlabeled.rds
# Run     : Rscript scripts/01_singlecell_annotation/08_marker_plots.R   (from the repository root)
# ==============================================================================
# 08_marker_plots.R
#
# scRNA UMAP plots for the canonical cell-type marker panel, discovery cohort
# (GSE162798). Uses the RNA UMAP (umap.rna) and the QC'd final cell barcodes.
#
# Panel: EPCAM, ITGA6 (CD49f), FOXA1, KIT, ESR1, TP63, CD3D, CD79A, ITGAX,
#        VWF, RGS5, PDGFRA
# Note: FOXA1, ESR1 and TP63 are transcription factors; the remainder are
# surface/lineage markers. All 12 are present in the RNA assay.
#
# Output (PNG only): results/coloc/figures/
#   Fig_scRNA_A_umap_celltype.png     DimPlot coloured by cell type
#   Fig_scRNA_B_markers_grid.png      FeaturePlot grid, all 12 markers
#   Fig_scRNA_C_dotplot_markers.png   DotPlot marker x cell type
#   Fig_scRNA_marker_summary.csv      mean expression + pct expressing per cell type

## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork)
})

setwd(TNBC_ROOT)
OUT <- file.path(TNBC_ROOT, "results/coloc/figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ts <- function(m) message("[", format(Sys.time(), "%H:%M:%S"), "] ", m)

MARKERS <- c("EPCAM", "ITGA6", "FOXA1", "KIT", "ESR1", "TP63",
             "CD3D", "CD79A", "ITGAX", "VWF", "RGS5", "PDGFRA")

ts("Loading object (~10 GB, takes a few minutes)")
Final_cell_barcode <- fread('results/singlecell/final_cell_barcodes.csv', header = TRUE)
o <- readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')
o <- subset(o, cells = Final_cell_barcode$x)
DefaultAssay(o) <- "RNA"
ts(paste("cells:", ncol(o)))

# short cell-type labels for readability
o$celltype <- as.character(o$celltype_RNA)
short <- c("luminal secretory-1 (LS-1)" = "LS-1", "luminal secretory-2 (LS-2)" = "LS-2",
           "luminal HR (LHr)" = "LHr", "basal epithelial (Bs)" = "Bs",
           "fibroblast (Fb)" = "Fb", "endothelial (Endo)" = "Endo",
           "Pericyte (Pc)" = "Pc", "T-cell (T)" = "T",
           "Myeloid (Myl)" = "Myl", "B-cell (B)" = "B")
o$celltype_short <- ifelse(o$celltype %in% names(short), short[o$celltype], o$celltype)
lev <- c("LS-1", "LS-2", "LHr", "Bs", "Fb", "Pc", "Endo", "T", "B", "Myl")
o$celltype_short <- factor(o$celltype_short, levels = intersect(lev, unique(o$celltype_short)))
Idents(o) <- "celltype_short"

RED <- if ("umap.rna" %in% Reductions(o)) "umap.rna" else "umap"
ts(paste("using reduction:", RED))

# Okabe-Ito colourblind-safe palette
COLS <- c("#0072B2", "#56B4E9", "#E69F00", "#009E73", "#D55E00",
          "#CC79A7", "#F0E442", "#999999", "#661100", "#332288")

# ---- A. DimPlot by cell type -------------------------------------------------
ts("A: DimPlot by cell type")
pA <- DimPlot(o, reduction = RED, group.by = "celltype_short",
              label = TRUE, label.size = 3.5, repel = TRUE, pt.size = 0.15) +
  scale_colour_manual(values = COLS[seq_len(nlevels(o$celltype_short))]) +
  labs(title = "Normal breast scRNA-seq: cell types", x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.grid = element_blank(), legend.position = "right")
ggsave(file.path(OUT, "Fig_scRNA_A_umap_celltype.png"), pA,
       width = 7.2, height = 5.6, dpi = 300)

# ---- B. FeaturePlot grid -----------------------------------------------------
ts("B: FeaturePlot grid for 12 markers")
present <- intersect(MARKERS, rownames(o[["RNA"]]))
ts(paste("markers plotted:", paste(present, collapse = ", ")))

plots <- lapply(present, function(g) {
  FeaturePlot(o, features = g, reduction = RED, pt.size = 0.12, order = TRUE) +
    scale_colour_gradientn(colours = c("grey88", "#D1E5F0", "#4393C3", "#2166AC", "#67001F")) +
    labs(title = g) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold.italic", size = 10, hjust = 0.5),
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.key.width = unit(0.25, "cm"), legend.title = element_blank())
})
pB <- wrap_plots(plots, ncol = 4) +
  plot_annotation(title = "Cell-type marker expression (scRNA-seq)",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(OUT, "Fig_scRNA_B_markers_grid.png"), pB,
       width = 12, height = 9, dpi = 300)

# ---- C. DotPlot --------------------------------------------------------------
ts("C: DotPlot marker x cell type")
pC <- DotPlot(o, features = present, group.by = "celltype_short", dot.scale = 6) +
  scale_colour_gradientn(colours = c("#F7F7F7", "#92C5DE", "#2166AC", "#67001F"),
                         name = "Mean\nexpression") +
  labs(x = NULL, y = NULL, title = "Marker expression by cell type") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        plot.title = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank())
ggsave(file.path(OUT, "Fig_scRNA_C_dotplot_markers.png"), pC,
       width = 8.5, height = 4.6, dpi = 300)

# ---- source data -------------------------------------------------------------
ts("writing source data")
dp <- pC$data
setDT(dp)
fwrite(dp[, .(marker = features.plot, cell_type = id,
              mean_expression = round(avg.exp, 4),
              mean_expression_scaled = round(avg.exp.scaled, 4),
              pct_expressing = round(pct.exp, 2))],
       file.path(OUT, "Fig_scRNA_marker_summary.csv"))

ts("DONE - wrote 3 PNGs + source data to results/coloc/figures/")
