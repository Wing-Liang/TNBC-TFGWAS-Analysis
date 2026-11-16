# ==============================================================================
# 10_mdm4_rps18_celltype_dotplot.R
# MDM4 and RPS18 expression by cell type, discovery cohort.
# Produces: Fig S11
# Inputs  : results/singlecell/final_cell_barcodes.csv
#           results/singlecell/atac_motif.rds
# Run     : Rscript scripts/09_motifbreakR_gtex_coloc/10_mdm4_rps18_celltype_dotplot.R   (from the repository root)
# ==============================================================================
## 10_mdm4_rps18_celltype_dotplot.R — dot plot of MDM4 and RPS18 RNA expression by cell type, main analysis set.
## Seurat DotPlot on the RNA assay (log-normalised): dot size = % of cells expressing; colour = average expression scaled
## across cell types within each gene (Seurat default, z-score clipped to [-2.5, 2.5]). Same object, cell set, labels and
## cell-type order as the other main-set figures. PNG only.
## Output: results/coloc/MDM4_RPS18_expression/main_MDM4_RPS18_DotPlot.png
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({library(Seurat); library(Signac); library(data.table); library(ggplot2)})
setwd(TNBC_ROOT)
OUT <- "results/coloc/MDM4_RPS18_expression"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ORD <- c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "Pc", "T", "B", "Myl")
fb  <- fread("results/singlecell/final_cell_barcodes.csv", header = TRUE)
obj <- subset(readRDS("results/singlecell/atac_motif.rds"), cells = fb$x)
ids <- c("LASP-1", "Fb", "LHS", "LASP-2", "B-Myo", "Endo", "T", "Pc", "Myl", "B"); names(ids) <- levels(obj)
obj <- RenameIdents(obj, ids)
Idents(obj) <- factor(Idents(obj), levels = rev(ORD))          # rev(): first cell type drawn at the top
DefaultAssay(obj) <- "RNA"; obj <- NormalizeData(obj, verbose = FALSE)
cat(ncol(obj), "cells\n")
p <- DotPlot(obj, features = c("MDM4", "RPS18"), assay = "RNA", cols = c("grey90", "#B40426"), dot.scale = 7) +
  labs(x = NULL, y = NULL, title = "Main analysis set (GSE162798)") +
  guides(size = guide_legend(title = "Cells expressing (%)", order = 1),
         colour = guide_colourbar(title = "Scaled average\nexpression", order = 2)) +
  theme(axis.text.x = element_text(face = "italic", size = 11), axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8))
f <- file.path(OUT, "main_MDM4_RPS18_DotPlot.png")
ggsave(f, p, width = 4.2, height = 4.6, dpi = 300, bg = "white")
print(as.data.table(p$data)[, .(gene = features.plot, celltype = id, pct.exp = round(pct.exp, 1),
                                avg.exp = round(avg.exp, 3), avg.exp.scaled = round(avg.exp.scaled, 2))])
cat("wrote", f, "\n")
