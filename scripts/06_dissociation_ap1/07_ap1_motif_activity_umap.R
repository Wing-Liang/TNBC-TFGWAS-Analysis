#!/usr/bin/env Rscript
# ==============================================================================
# 07_ap1_motif_activity_umap.R
# AP-1 motif activity on the UMAP, both cohorts.
# Produces: Fig S9
# Inputs  : replication/data/cellxgene/embeddings_replication.csv
#           replication/celltypes/cellxgene_celltypes_replication.csv
#           replication/tfgwas_grouped/replication_atac_umap.csv
# Run     : Rscript scripts/06_dissociation_ap1/07_ap1_motif_activity_umap.R   (from the repository root)
# ==============================================================================
# 07_ap1_motif_activity_umap.R -- AP-1 (TRE / CRE motif) accessibility by cell, both cohorts.
#
# Shows which cells carry AP-1 TRE accessibility in each cohort, how the cohorts differ in those
# cells, and how much of the difference cell-type composition alone accounts for.
#  * chromVAR deviation z-scores are standardised within a dataset against its own background
#    peaks, so both cohort means are ~0 by construction and absolute z is not comparable between
#    cohorts. The UMAP panels are drawn on per-cohort scales; the between-cohort comparison rests
#    on peak-level counts and on composition.
#  * With epithelial cells included, LHS has the highest TRE activity of any discovery cell type
#    (+2.43), ahead of fibroblasts (+1.24); basal/myoepithelial cells have the lowest (-2.68).
#  * Basal/myoepithelial cells rise from 10.5% of discovery cells to 37.7% of replication cells,
#    while fibroblasts and LHS shrink; reweighting discovery cells to the replication composition
#    moves mean TRE by -0.684 (panel f).
#
# INPUT   results/singlecell/atac_motif.rds                         discovery chromVAR
#         replication/tfgwas/replication_chromvar.rds               replication chromVAR
#         replication/data/cellxgene/embeddings_replication.csv     replication UMAP
#         reference/jaspar_tf_family.csv                            TF family of each motif
# AP-1 motifs are the Jun-, Fos- and B-ATF-related matrices (45). They are split by the expected
# number of CpG dinucleotides in the JASPAR2022 matrix, sum_i P(C at i) * P(G at i+1): TRE
# (TGA(C/G)TCA, CpG-free) < 0.5 <= CRE (TGACGTCA), giving 29 TRE and 16 CRE matrices.
# OUTPUT  replication/tfgwas_grouped/ap1_tre_umap_cells.csv
#         replication/tfgwas_grouped/ap1_tre_reweight.csv
#         replication/figures/atac/AP1_TRE_UMAP.pdf / .png

## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(data.table)
  library(ggplot2); library(patchwork)
  library(JASPAR2022); library(TFBSTools)
})
setwd(TNBC_ROOT)
ts <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
RES <- "replication/tfgwas_grouped"
FIG <- "replication/figures/atac"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
CACHE <- file.path(RES, "ap1_tre_umap_cells.csv")

GROUP_REF <- Sys.getenv("GROUP_REF", "reference")
AP1 <- c("Jun-related", "Fos-related", "B-ATF-related")
M <- fread(file.path(GROUP_REF, "jaspar_tf_family.csv"))[grepl(paste(AP1, collapse = "|"), family)]
cpg_expected <- function(id) {
  m <- Matrix(getMatrixByID(JASPAR2022, ID = id)); p <- sweep(m, 2, colSums(m), "/")
  sum(p["C", -ncol(p)] * p["G", -1])
}
M[, cpg_expected := vapply(motif_id, cpg_expected, numeric(1))]
TRE <- M[cpg_expected <  0.5, motif_id]
CRE <- M[cpg_expected >= 0.5, motif_id]
stopifnot(length(TRE) == 29, length(CRE) == 16)
ts(sprintf("TRE matrices %d | CRE matrices %d", length(TRE), length(CRE)))

HARM <- c(
  "LS-1" = "LASP", "LS-2" = "LASP", "LASP" = "LASP",
  "LHr" = "LHS", "LHS" = "LHS",
  "Bs" = "Basal/myoep", "BM" = "Basal/myoep",
  "Fb" = "Fibroblast", "Fibroblasts" = "Fibroblast",
  "Pc" = "Pericyte",
  "Endo" = "Endothelial", "Endo-1" = "Endothelial", "Endo-2" = "Endothelial",
  "T" = "T-cell", "T-cells" = "T-cell",
  "Myl" = "Myeloid", "Macrophages" = "Myeloid",
  "B" = "B-cell", "Adi-1" = "Adipocyte", "Adi-2" = "Adipocyte")
harm <- function(x) { y <- HARM[as.character(x)]; ifelse(is.na(y), as.character(x), y) }

pick_umap <- function(o) {
  r <- names(o@reductions)
  for (cand in c("umap.atac", "umap", "wnn.umap", "umap.rna")) if (cand %in% r) return(cand)
  r[grep("umap", r, ignore.case = TRUE)][1]
}
emb_from_file <- function(cells, lab) {
  E <- fread("replication/data/cellxgene/embeddings_replication.csv")
  i <- match(cells, E$cell)
  ts(lab, ": embedding joined from CELLxGENE for ", sum(!is.na(i)), " of ", length(cells))
  if (mean(!is.na(i)) < 0.5) stop(lab, ": barcodes do not match the embedding table")
  cbind(E$umap.atac_1[i], E$umap.atac_2[i])
}

extract <- function(path, lab, ctcol_candidates) {
  ts(lab, ": loading")
  o <- readRDS(path)
  stopifnot("chromvar" %in% names(o@assays))
  cv <- GetAssayData(o, assay = "chromvar", layer = "data")
  ts(lab, ": chromvar ", nrow(cv), " motifs x ", ncol(cv), " cells")
  tre <- intersect(TRE, rownames(cv)); cre <- intersect(CRE, rownames(cv))
  stopifnot(length(tre) > 10, length(cre) > 5)
  red <- pick_umap(o)
  if (length(red) && !is.na(red) && nzchar(red)) {
    ts(lab, ": using reduction '", red, "'"); em <- Embeddings(o, red)[, 1:2, drop = FALSE]
  } else { ts(lab, ": no UMAP reduction in object"); em <- emb_from_file(colnames(cv), lab) }
  md <- o@meta.data; ct <- NULL
  for (cc in ctcol_candidates) if (cc %in% names(md)) { ct <- as.character(md[[cc]]); break }
  if (is.null(ct)) ct <- as.character(Idents(o))
  D <- data.table(cohort = lab, cell = colnames(cv), UMAP1 = em[, 1], UMAP2 = em[, 2],
                  celltype_raw = ct,
                  TRE = Matrix::colMeans(cv[tre, , drop = FALSE]),
                  CRE = Matrix::colMeans(cv[cre, , drop = FALSE]))
  rm(o, cv); gc()
  n0 <- nrow(D); D <- D[!is.na(UMAP1) & !is.na(UMAP2) & !is.na(TRE)]
  if (nrow(D) < n0) ts(lab, ": dropped ", n0 - nrow(D), " cells lacking an embedding")
  D
}

CO <- c("Discovery (SNuBar, 4 donors)", "Replication (10x Multiome, 15 donors)")
if (file.exists(CACHE)) {
  ts("using cached per-cell table (delete it to re-extract)")
  X <- fread(CACHE)
} else {
  A <- extract("results/singlecell/atac_motif.rds", CO[1],
               c("celltype_ATAC", "celltype_RNA"))
  B <- extract("replication/tfgwas/replication_chromvar.rds", CO[2],
               c("author_cell_type", "celltype_ATAC", "celltype", "predicted.id"))
  if (uniqueN(B$celltype_raw) < 3) {
    ts("replication: joining CELLxGENE cell types")
    R <- fread("replication/celltypes/cellxgene_celltypes_replication.csv")
    B[, celltype_raw := R$author_cell_type[match(cell, R$cell)]]
  }
  X <- rbind(A, B); X[, celltype := harm(celltype_raw)]
  X <- X[!is.na(celltype) & celltype != "NA"]
  fwrite(X, CACHE)
}
## drop cells without a cell-type label (some replication chromVAR cells have no CELLxGENE label)
n0 <- nrow(X)
X <- X[!is.na(celltype) & !celltype %in% c("NA", "")]
if (nrow(X) < n0) ts("dropped ", n0 - nrow(X), " cells without a cell-type label")
X[, cohort := factor(cohort, levels = CO)]

## ---- replication panels: ATAC embedding of the replication cells alone ----------------------
## (06_replication_atac_umap.R), so each cohort is drawn on an embedding of its own cells, as for
## discovery. The source atlas embedding (data/cellxgene/embeddings_replication.csv) was computed
## over all 51,367 atlas cells, of which these 10,700 are ~21%. kNN cell-type purity is similar on
## either embedding (0.843 atlas vs 0.854 own; discovery 0.954).
NE <- fread("replication/tfgwas_grouped/replication_atac_umap.csv")
X[, newU1 := NE$UMAP1[match(cell, NE$cell)]]
X[, newU2 := NE$UMAP2[match(cell, NE$cell)]]
nrep <- X[cohort == CO[2] & !is.na(newU1), .N]
X[cohort == CO[2] & !is.na(newU1), `:=`(UMAP1 = newU1, UMAP2 = newU2)]
X <- X[!(cohort == CO[2] & is.na(newU1))]
X[, c("newU1", "newU2") := NULL]
ts("replication panels re-anchored to the European-only ATAC UMAP (", nrep, " cells)")

## ---- composition reweighting ------------------------------------------------------------
tre_by <- X[cohort == CO[1], .(tre = mean(TRE)), by = celltype]
wD <- X[cohort == CO[1], .N, by = celltype][, .(celltype, w_disc = N / sum(N))]
wR <- X[cohort == CO[2], .N, by = celltype][, .(celltype, w_repl = N / sum(N))]
RW <- merge(merge(tre_by, wD, by = "celltype"), wR, by = "celltype", all.x = TRUE)
RW[is.na(w_repl), w_repl := 0][, w_repl := w_repl / sum(w_repl)]
obs <- RW[, sum(tre * w_disc)]; rew <- RW[, sum(tre * w_repl)]
fwrite(RW[order(-tre)], file.path(RES, "ap1_tre_reweight.csv"))
cat("\n---- TRE by cell type (discovery scale) and cohort weights ----\n")
print(RW[order(-tre)], row.names = FALSE)
cat(sprintf("\nobserved discovery mean TRE               : %+.3f\n", obs))
cat(sprintf("discovery reweighted to repl composition  : %+.3f\n", rew))
cat(sprintf("shift attributable to composition alone   : %+.3f\n", rew - obs))

## ---- plotting ---------------------------------------------------------------------------
base <- theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
        legend.key.height = unit(0.34, "cm"), legend.key.width = unit(0.28, "cm"),
        legend.text = element_text(size = 6.5), legend.title = element_text(size = 7))

PAL <- c("LASP" = "#1F78B4", "LHS" = "#A6CEE3", "Basal/myoep" = "#33A02C",
         "Fibroblast" = "#E31A1C", "Pericyte" = "#FB9A99", "Endothelial" = "#FF7F00",
         "T-cell" = "#6A3D9A", "Myeloid" = "#CAB2D6", "B-cell" = "#B15928",
         "Adipocyte" = "#FDBF6F")

ct_panel <- function(dt, title) {
  ggplot(dt[sample(.N)], aes(UMAP1, UMAP2, colour = celltype)) +
    geom_point(size = 0.16, stroke = 0, alpha = 0.85) +
    scale_colour_manual(values = PAL, drop = TRUE) +
    ggtitle(title) + base +
    guides(colour = guide_legend(override.aes = list(size = 1.7), ncol = 1, title = NULL))
}
## per-cohort limits: chromVAR z is within-dataset normalised (see header)
score_panel <- function(dt, title) {
  ## 5-95% rather than 1-99%: per-cell chromVAR z has heavy tails and the wider
  ## limits flatten the between-cell-type structure the panel exists to show
  L <- quantile(dt$TRE, c(0.05, 0.95), na.rm = TRUE); m <- max(abs(L)); L <- c(-m, m)
  ggplot(dt[order(abs(TRE))],
         aes(UMAP1, UMAP2, colour = pmin(pmax(TRE, L[1]), L[2]))) +
    geom_point(size = 0.16, stroke = 0, alpha = 0.85) +
    scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B",
                           midpoint = 0, limits = L,
                           name = "AP-1 TRE\nz (within\ncohort)") +
    ggtitle(title) + base
}

p1 <- ct_panel(X[cohort == CO[1]], CO[1])
p2 <- ct_panel(X[cohort == CO[2]], CO[2])
p3 <- score_panel(X[cohort == CO[1]], "AP-1 TRE accessibility - discovery")
p4 <- score_panel(X[cohort == CO[2]], "AP-1 TRE accessibility - replication")

## `ord` ranks by discovery TRE and therefore covers discovery cell types only. Replication
## -only types (adipocyte) must be appended, or factor() silently maps them to NA and the
## boxplot grows a phantom "NA" category.
ord  <- RW[order(-tre), celltype]
ordA <- c(ord, setdiff(sort(unique(X$celltype)), ord))
Y <- copy(X)[, celltype := factor(celltype, levels = ordA)]
p5 <- ggplot(Y, aes(celltype, TRE, fill = cohort)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_boxplot(outlier.size = 0.12, outlier.alpha = 0.2, linewidth = 0.25) +
  scale_fill_manual(values = c("#4C94D8", "#E08214"), name = NULL) +
  labs(x = NULL, y = "AP-1 TRE chromVAR z",
       subtitle = "Same ranking in both cohorts: LHS and fibroblast highest, basal/myoepithelial lowest") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top",
        plot.subtitle = element_text(size = 7.5, colour = "grey30"),
        panel.grid.minor = element_blank())

## RW covers only cell types present in discovery, because the reweighting needs a
## discovery TRE value. The composition panel must show everything, including adipocytes,
## which exist only in the replication cohort (2.3%) and are excluded from the reweighting.
FULL <- merge(X[cohort == CO[1], .N, by = celltype][, .(celltype, Discovery = N / sum(N))],
              X[cohort == CO[2], .N, by = celltype][, .(celltype, Replication = N / sum(N))],
              by = "celltype", all = TRUE)
FULL[is.na(Discovery), Discovery := 0][is.na(Replication), Replication := 0]
CB <- melt(FULL, id.vars = "celltype", variable.name = "cohort", value.name = "frac")
CB[, celltype := factor(celltype, levels = ordA)]
p6 <- ggplot(CB, aes(celltype, 100 * frac, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c("#4C94D8", "#E08214"), name = NULL) +
  labs(x = NULL, y = "% of cells",
       subtitle = sprintf(paste("Composition shifts mass onto the lowest-TRE cell type.",
                                "Reweighting discovery cells to the replication",
                                "\ncomposition moves mean TRE %+.2f to %+.2f (shift %+.2f)."),
                          obs, rew, rew - obs)) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top",
        plot.subtitle = element_text(size = 7.5, colour = "grey30"),
        panel.grid.minor = element_blank())

fig <- ((p1 | p2) / (p3 | p4) / p5 / p6) +
  plot_layout(heights = c(1.05, 1.05, 0.85, 0.85)) +
  plot_annotation(
    title = "AP-1/TRE accessibility differs between cohorts through cell composition",
    subtitle = paste0(
      "Per-cell mean chromVAR deviation over ", length(TRE), " CpG-free AP-1 TRE matrices. ",
      "chromVAR z is standardised within each dataset, so panels c and d use per-cohort\n",
      "scales and show WHICH cells carry TRE accessibility, not the absolute difference ",
      "between cohorts; that difference is established at the peak level (2 of 83 AP-1\n",
      "associations replicate, and AP-1-bearing discovery peaks are recovered at 20.5% ",
      "against 40.5% for GC-matched controls). Each cohort is embedded on its own cells; the
",
      "replication cell types separate less sharply in ATAC space than the discovery ones
",
      "(kNN label purity 0.85 against 0.95), which is a property of the data."),
    tag_levels = "a",
    theme = theme(plot.title = element_text(size = 11.5, face = "bold"),
                  plot.subtitle = element_text(size = 7.8, colour = "grey30")))

ggsave(file.path(FIG, "AP1_TRE_UMAP.pdf"), fig, width = 10, height = 14.5)
ggsave(file.path(FIG, "AP1_TRE_UMAP.png"), fig, width = 10, height = 14.5, dpi = 200)
ts("figure written to ", FIG)
