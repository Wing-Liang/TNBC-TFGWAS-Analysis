# ==============================================================================
# 23_marker_detection_by_cohort.R
# Canonical marker detection by cell type in both cohorts.
# Produces: Table S6
# Inputs  : results/singlecell/final_cell_barcodes.csv
#           results/singlecell/atac_motif.rds
#           replication/rna/replication_RNA_annotated.rds
# Outputs : results/motif_sets/marker_compare_main_vs_replication.csv
# Run     : Rscript scripts/05_replication/23_marker_detection_by_cohort.R   (from the repository root)
# ==============================================================================
## 23_marker_detection_by_cohort.R — canonical marker detection (% cells > 0) and mean log-normalised expression per cell
## type, main (GSE162798) vs replication (GSE244585), (Table S6)..
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({library(Seurat); library(Signac); library(data.table)})
setwd(TNBC_ROOT)
G <- c("EPCAM","ITGA6","ESR1","FOXA1","KIT","TP63","CD3D","CD79A","ITGAX","VWF","PDGFRA","RGS5","KRT18","KRT8","KRT5","KRT14")
summ <- function(obj, lab, cohort) {
  DefaultAssay(obj) <- "RNA"; obj <- NormalizeData(obj, verbose = FALSE)
  g <- intersect(G, rownames(obj)); X <- GetAssayData(obj, assay = "RNA", layer = "data")[g, , drop = FALSE]
  rbindlist(lapply(unique(lab), function(ct) { i <- which(lab == ct)
    data.table(cohort, celltype = ct, n = length(i), gene = g,
               pct = 100 * Matrix::rowMeans(X[, i, drop = FALSE] > 0), mean = Matrix::rowMeans(X[, i, drop = FALSE])) }))
}
fb <- fread("results/singlecell/final_cell_barcodes.csv", header = TRUE)
m <- subset(readRDS("results/singlecell/atac_motif.rds"), cells = fb$x)
ids <- c("LASP-1","Fb","LHS","LASP-2","B-Myo","Endo","T","Pc","Myl","B"); names(ids) <- levels(m); m <- RenameIdents(m, ids)
A <- summ(m, as.character(Idents(m)), "Main"); rm(m); gc(verbose = FALSE)
r <- readRDS("replication/rna/replication_RNA_annotated.rds")
B <- summ(r, as.character(r$celltype), "Replication")
D <- rbind(A, B); fwrite(D, "results/motif_sets/marker_compare_main_vs_replication.csv")
epi <- c("LASP-1","LASP-2","LASP","LHS","B-Myo","BM")
cat("\n% cells expressing, epithelial cell types:\n")
print(dcast(D[celltype %in% epi & gene %in% c("EPCAM","ITGA6","KIT","ESR1","FOXA1","TP63","KRT18","KRT8","KRT5","KRT14")],
            gene ~ cohort + celltype, value.var = "pct"), digits = 3)
cat("\nmean log-normalised expression, epithelial:\n")
print(dcast(D[celltype %in% epi & gene %in% c("EPCAM","KRT18","KRT8")], gene ~ cohort + celltype, value.var = "mean"), digits = 3)
## pooled epithelium
P <- D[celltype %in% epi, .(pct = weighted.mean(pct, n), mean = weighted.mean(mean, n), n = sum(n)), by = .(cohort, gene)]
cat("\npooled epithelium:\n"); print(dcast(P, gene ~ cohort, value.var = c("pct", "mean")), digits = 3)
## sequencing depth per cell (RNA)
