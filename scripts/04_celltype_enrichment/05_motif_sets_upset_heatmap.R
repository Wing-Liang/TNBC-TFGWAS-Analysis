# ==============================================================================
# 05_motif_sets_upset_heatmap.R
# Motif sets per subtype, cell-type enrichment and module scores.
# Produces: Fig 3A-B, Figs S5-S6
# Inputs  : results/singlecell/final_cell_barcodes.csv
#           results/singlecell/atac_motif.rds
#           results/singlecell/combined_peaks_gene_all.csv
#           results/singlecell/jaspar2022_motifs_727.csv
#           results/singlecell/jaspar2022_motif_list.csv
# Run     : Rscript scripts/04_celltype_enrichment/05_motif_sets_upset_heatmap.R   (from the repository root)
# ==============================================================================
## 05_motif_sets_upset_heatmap.R — Figure 3 motif sets, cell-type enrichment and module scores.
##   (1) motif enrichment on ALL DAPs of each cell type (no "unique peak" restriction)
##   (2) FindMotifs run once per cell type, cached, joined to every motif set
##   (3) motif sets: Bonferroni-significant per subtype (5); ER+ shared = LumA ∩ LumB ∩ LumB/HER2-;
##       ER- shared = TNBC ∩ HER2; subtype-specific = sig(s) minus union of the other four subtypes;
##       overall BC = reference layer only (never subtracted)
##   (4) outputs per set: bubble plot (Fig S6), family counts (Fig S5), AddModuleScore UMAP (Fig 3D),
##       box plots with Kruskal-Wallis + BH-Dunn two-sided (Fig 3E)
##   (5) Fig 3A UpSet of the five subtype sets; Fig 3B heatmap of mean module score per cell type
##   (6) TF expression for bubble sizes computed in-script from the RNA assay
##   Manhattan panels are NOT produced here (see 03_subtype_specific_motifs.R, Manhattan section)
## Run from the repository root:  Rscript scripts/04_celltype_enrichment/05_motif_sets_upset_heatmap.R
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(patchwork)
  library(data.table); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggpubr); library(ggsignif); library(dunn.test)
  library(fmesher); #library(INLA)
  library(ComplexHeatmap); library(circlize)
})
setMethod("seqinfo", "ChromatinAssay", function(x) x@seqinfo)
setMethod("genome",  "ChromatinAssay", function(x) {
  si <- x@seqinfo
  if (is.null(si)) NULL else GenomeInfoDb::genome(si)
})
library(ComplexHeatmap)
packageVersion("ComplexHeatmap")
exists("make_comb_mat")

setwd(TNBC_ROOT)

## ---------------------------------------------------------------- parameters
BONF     <- 0.05 / 727
SEED     <- 4242025
UMAP_RED <- "umap.rna"                            # reduction used for the UMAP panels + hull
MAX_MOTIFS_PER_BUBBLE <- 60                       # large sets are split into parts by family
OVERALLBC_FILE <- "replication/tfgwas/Discovery_OverallBC_ALL.csv"   # overall BC model, discovery cohort (09e_combine_overall_bc.R)
subtypes <- c("TNBC", "LumA", "LumB", "LumB_HER2Neg", "HER2")
ERpos    <- c("LumA", "LumB", "LumB_HER2Neg")
ERneg    <- c("TNBC", "HER2")
ct_levels <- c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "T", "B", "Myl", "Pc")
epi       <- c("LASP-1", "LASP-2", "LHS", "B-Myo")
nonepi    <- c("Fb", "Endo", "T", "B", "Myl", "Pc")
rename_ct <- function(x) { x <- sub("^LS-1$", "LASP-1", x); x <- sub("^LS-2$", "LASP-2", x)
x <- sub("^LHr$",  "LHS",    x); sub("^Bs$", "B-Myo", x) }
OUTF <- "results/motif_sets/figures"; OUTR <- "results/motif_sets/"
dir.create(OUTF, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTR, recursive = TRUE, showWarnings = FALSE)
xlab_bubble <- "Fold enrichment\n(over-representation in differentially accessible peaks of each cell type)"

## ---------------------------------------------------------------- object
Final_cell_barcode <- fread("results/singlecell/final_cell_barcodes.csv", header = TRUE)
obj <- readRDS("results/singlecell/atac_motif.rds")
obj <- subset(obj, cells = Final_cell_barcode$x)
new.cluster.ids <- c("LASP-1", "Fb", "LHS", "LASP-2", "B-Myo", "Endo", "T", "Pc", "Myl", "B")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
Idents(obj) <- factor(Idents(obj), levels = ct_levels)
obj$celltype <- Idents(obj)
epi_cells    <- WhichCells(obj, idents = epi)
nonepi_cells <- WhichCells(obj, idents = nonepi)
cat("cells:", ncol(obj), "\n"); print(table(Idents(obj)))
# cells: 16451
#
# LASP-1 LASP-2    LHS  B-Myo     Fb   Endo      T      B    Myl     Pc
## ---------------------------------------------------------------- DAPs -> FindMotifs, once per cell type
dap <- fread("results/singlecell/combined_peaks_gene_all.csv", header = TRUE)
setnames(dap, 1, "peaks")
dap[, cluster := rename_ct(cluster)]
dap <- unique(dap[, .(peaks, cluster)])           # one row per peak per cell type (file has one row per transcript annotation)
stopifnot(all(unique(dap$cluster) %in% ct_levels))
cat("\nDAPs per cell type:\n"); print(dap[, .N, by = cluster][match(ct_levels, cluster)])

DefaultAssay(obj) <- "ATAC"
enr_file <- file.path(OUTR, "FindMotifs_allDAP_by_celltype.rds")
if (file.exists(enr_file)) {
  enrich_all <- readRDS(enr_file)
} else {
  enrich_all <- rbindlist(lapply(ct_levels, function(ct) {
    feats <- unique(dap[cluster == ct, peaks]); feats <- feats[feats %in% rownames(obj)]
    cat(sprintf("FindMotifs %-7s n DAPs = %d\n", ct, length(feats)))
    em <- as.data.table(FindMotifs(object = obj, features = feats)); em[, motif.name := NULL]
    em[, cluster := ct]
  }))
  saveRDS(enrich_all, enr_file)
}

fwrite(enrich_all, file.path(OUTR, "FindMotifs_allDAP_by_celltype.csv"))
enrich_all[, pvalue := pmin(pmax(pvalue, 1e-300), 1 - 1e-16)]
enrich_all[, neglog10pval := -log10(pvalue)]
## ---------------------------------------------------------------- motif reference + GWAS results
motif_ref <- fread("results/singlecell/jaspar2022_motifs_727.csv")[, .(motif, motif.name)]
fam_ref   <- unique(fread("results/singlecell/jaspar2022_motif_list.csv")[, .(motif, family)])[!duplicated(motif)]
motif_ref <- merge(motif_ref, fam_ref, by = "motif", all.x = TRUE)
motif_ref[is.na(family) | family == "", family := "Unassigned"]
motif_ref[, family := ifelse(family == "Thyroid hormone receptor-related factors (NR1);RXR-related receptors (NR2)",
                             "Thyroid hormone receptor/RXR-related factor", family)]
stopifnot(!anyDuplicated(motif_ref$motif))

read_model <- function(f) fread(f)[motif %in% motif_ref$motif, .(motif, TF_p_value)][!duplicated(motif)]
gwas <- lapply(subtypes, function(s) read_model(sprintf("results/tfgwas/tfgwas_results_%s_ALL.csv", s)))
names(gwas) <- subtypes
sig <- lapply(gwas, function(d) d[TF_p_value < BONF, motif])
cat("\nBonferroni-significant motifs per subtype:\n"); print(sapply(sig, length))
# Bonferroni-significant motifs per subtype:
#   TNBC         LumA         LumB LumB_HER2Neg         HER2

## overall breast cancer (overall BC) model: reference layer, NOT used in the subtype-specific subtraction
OverallBC <- fread(OVERALLBC_FILE)[motif %in% motif_ref$motif, .(motif, TF_p_value = OverallBC_p)][!duplicated(motif)]
sig_OverallBC <- OverallBC[TF_p_value < BONF, motif]
cat("overall BC Bonferroni-significant motifs:", length(sig_OverallBC), "\n")
# overall BC Bonferroni-significant motifs: 548

## overlap of every motif set with the overall BC Bonferroni-significant set

# set n_motifs in_OverallBC pct_in_OverallBC not_in_OverallBC
# <char>    <int>    <int>        <num>        <int>
#   1:          TNBC_BonfSig      260      236         90.8           24
# 2:          LumA_BonfSig      433      431         99.5            2
# 3:          LumB_BonfSig      214      208         97.2            6
# 4:  LumB_HER2Neg_BonfSig      217      213         98.2            4
# 5:          HER2_BonfSig       38       28         73.7           10
# 6:          ERpos_shared      137      137        100.0            0
# 7:          ERneg_shared       14       14        100.0            0
# 8:         TNBC_specific       57       35         61.4           22
# 9:         LumA_specific      104      102         98.1            2
# 10:         LumB_specific       13        9         69.2            4
# 11: LumB_HER2Neg_specific        9        5         55.6            4
# 12:         HER2_specific       14        4         28.6           10

## ---------------------------------------------------------------- motif sets
union_of <- function(ss) Reduce(union, sig[ss])
specific <- setNames(lapply(subtypes, function(s) setdiff(sig[[s]], union_of(setdiff(subtypes, s)))), subtypes)

motif_sets <- list()
for (s in subtypes) motif_sets[[paste0(s, "_BonfSig")]] <- sig[[s]]
motif_sets[["ERpos_shared"]] <- Reduce(intersect, sig[ERpos])
motif_sets[["ERneg_shared"]] <- Reduce(intersect, sig[ERneg])
for (s in subtypes) motif_sets[[paste0(s, "_specific")]] <- specific[[s]]
motif_sets[["OverallBC_BonfSig"]] <- sig_OverallBC

set_sizes <- data.table(set = names(motif_sets), n_motifs = sapply(motif_sets, length))
cat("\nmotif set sizes:\n"); print(set_sizes)
fwrite(set_sizes, file.path(OUTR, "motif_set_sizes.csv"))

# motif set sizes:
#   set n_motifs
# <char>    <int>
#   1:          TNBC_BonfSig      260
# 2:          LumA_BonfSig      433
# 3:          LumB_BonfSig      214
# 4:  LumB_HER2Neg_BonfSig      217
# 5:          HER2_BonfSig       38
# 6:          ERpos_shared      137
# 7:          ERneg_shared       14
# 8:         TNBC_specific       57
# 9:         LumA_specific      104
# 10:         LumB_specific       13
# 11: LumB_HER2Neg_specific        9
# 12:         HER2_specific       14
# 13:         OverallBC_BonfSig      548

## power diagnostic: are "specific" motifs nominally significant elsewhere / significant in overall BC?
cat("\nspecific sets: nominal significance in other subtypes and in overall BC\n")
diag <- rbindlist(lapply(subtypes, function(s) {
  m <- specific[[s]]; others <- setdiff(subtypes, s)
  pm <- sapply(others, function(o) gwas[[o]]$TF_p_value[match(m, gwas[[o]]$motif)])
  if (is.null(dim(pm))) pm <- matrix(pm, nrow = length(m))
  data.table(subtype = s, n_specific = length(m),
             nominal_p05_other = sum(apply(pm < 0.05, 1, any, na.rm = TRUE)),
             p5e4_other        = sum(apply(pm < 5e-4, 1, any, na.rm = TRUE)),
             OverallBC_sig         = sum(m %in% sig_OverallBC))
}))
print(diag); fwrite(diag, file.path(OUTR, "specific_sets_power_diagnostic.csv"))
# subtype n_specific nominal_p05_other p5e4_other OverallBC_sig
# <char>      <int>             <int>      <int>     <int>
#   1:         TNBC         57                49         20        35
# 2:         LumA        104                98         44       102
# 3:         LumB         13                12          8         9
# 4: LumB_HER2Neg          9                 8          2         5
# 5:         HER2         14                11          3         4

## supplementary membership table
## supplementary membership table (Fig 3A/3B source data)
OverallBC_full <- fread(OVERALLBC_FILE)[motif %in% motif_ref$motif][!duplicated(motif)]
memb <- copy(motif_ref)[, .(motif, motif.name, family)]
for (s in subtypes) {
  memb[, (paste0(s, "_TF_p_value")) := gwas[[s]]$TF_p_value[match(motif, gwas[[s]]$motif)]]
  memb[, (paste0(s, "_BonfSig"))    := ifelse(motif %in% sig[[s]], "Y", "")]
}
memb[, OverallBC_B          := OverallBC_full$OverallBC_B[match(motif, OverallBC_full$motif)]]
memb[, OverallBC_SE         := OverallBC_full$OverallBC_SE[match(motif, OverallBC_full$motif)]]
memb[, OverallBC_TF_p_value := OverallBC_full$OverallBC_p[match(motif, OverallBC_full$motif)]]
memb[, OverallBC_BonfSig    := ifelse(motif %in% sig_OverallBC, "Y", "")]
memb[, n_subtypes_BonfSig := rowSums(sapply(subtypes, function(s) motif %in% sig[[s]]))]
for (s in subtypes) memb[, (paste0(s, "_BonfSig_and_OverallBC")) := ifelse(motif %in% sig[[s]] & motif %in% sig_OverallBC, "Y", "")]
memb[, ERpos_shared := ifelse(motif %in% motif_sets$ERpos_shared, "Y", "")]
memb[, ERneg_shared := ifelse(motif %in% motif_sets$ERneg_shared, "Y", "")]
for (s in subtypes) memb[, (paste0(s, "_specific")) := ifelse(motif %in% specific[[s]], "Y", "")]
memb[, specific_in := sapply(motif, function(m) { s <- subtypes[sapply(subtypes, function(x) m %in% specific[[x]])]; if (length(s)) s else "" })]
memb[, OverallBC_only := ifelse(motif %in% sig_OverallBC & n_subtypes_BonfSig == 0, "Y", "")]
setorder(memb, -n_subtypes_BonfSig, family, motif.name)
fwrite(memb, file.path(OUTR, "tfgwas_motif_sets_membership_all_subtypes.csv"))
for (n in names(motif_sets)) fwrite(motif_ref[motif %in% motif_sets[[n]], .(motif, motif.name, family)],
                                    file.path(OUTR, sprintf("motifs_%s.csv", n)))

## ---------------------------------------------------------------- Fig 3A: UpSet of Bonferroni-significant motifs, 5 subtypes
sub_lab  <- c(TNBC = "TNBC",LumA = "Luminal A", LumB = "Luminal B", LumB_HER2Neg = "Luminal B/HER2-", HER2 = "HER2-enriched")
row_ord  <- c("TNBC", "LumA", "LumB", "LumB_HER2Neg","HER2")
set_cols <- c(TNBC = "#E64B35", HER2 = "#F39B7F", LumA = "#3C5488", LumB = "#4DBBD5", LumB_HER2Neg = "#00A087")
all_five <- Reduce(intersect, sig)
combos <- list(
  list(id = "ER+ shared",          sets = ERpos,          motifs = motif_sets$ERpos_shared, col = "#8491B4"),
  list(id = "LumA-specific",       sets = "LumA",         motifs = specific$LumA,           col = set_cols[["LumA"]]),
  list(id = "LumB-specific",       sets = "LumB",         motifs = specific$LumB,           col = set_cols[["LumB"]]),
  list(id = "LumB/HER2- specific", sets = "LumB_HER2Neg", motifs = specific$LumB_HER2Neg,   col = set_cols[["LumB_HER2Neg"]]),
  list(id = "ER- shared",          sets = ERneg,          motifs = motif_sets$ERneg_shared, col = "#B0703C"),
  list(id = "TNBC-specific",       sets = "TNBC",         motifs = specific$TNBC,           col = set_cols[["TNBC"]]),
  list(id = "HER2-specific",       sets = "HER2",         motifs = specific$HER2,           col = set_cols[["HER2"]]),
  list(id = "All five",            sets = subtypes,       motifs = all_five,                col = "grey20"))
ids  <- sapply(combos, `[[`, "id")
cols <- setNames(sapply(combos, `[[`, "col"), ids)

bars <- rbindlist(lapply(combos, function(x)
  data.table(id = x$id, n = length(x$motifs),
             OverallBC = c("overall BC significant", "Not overall BC significant"),
             k = c(sum(x$motifs %in% sig_OverallBC), sum(!x$motifs %in% sig_OverallBC)))))
bars[, id := factor(id, levels = ids)]
bars[, OverallBC := factor(OverallBC, levels = c("Not overall BC significant", "overall BC significant"))]
dots <- rbindlist(lapply(combos, function(x) data.table(id = x$id, subtype = row_ord, on = row_ord %in% x$sets)))
dots[, id := factor(id, levels = ids)]
dots[, ypos := match(subtype, rev(row_ord))]
segs  <- dots[on == TRUE, .(ymin = min(ypos), ymax = max(ypos)), by = id]
sizes <- data.table(subtype = factor(row_ord, levels = rev(row_ord)), n = sapply(sig[row_ord], length))
blk <- theme(text = element_text(colour = "black", size = 12), axis.text = element_text(colour = "black", size = 11),
             axis.title = element_text(colour = "black", size = 12), legend.text = element_text(colour = "black", size = 10))
p_top <- ggplot(bars, aes(id, k, fill = id, alpha = OverallBC)) + geom_col(width = 0.6) +
  geom_text(data = unique(bars[, .(id, n)]), aes(id, n, label = n), inherit.aes = FALSE, vjust = -0.3, size = 3) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_alpha_manual(values = c(`Not overall BC significant` = 0.3, `overall BC significant` = 1), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "") + theme_classic(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        legend.position = c(0.82, 0.88), legend.key.size = unit(3.5, "mm"), legend.text = element_text(size = 8),
        legend.background = element_blank())+ blk
p_dots <- ggplot() +
  geom_segment(data = segs, aes(x = id, xend = id, y = ymin, yend = ymax), linewidth = 0.8, colour = "grey30") +
  geom_point(data = dots, aes(id, ypos, colour = on), size = 3.2) +
  scale_colour_manual(values = c(`TRUE` = "grey20", `FALSE` = "grey85"), guide = "none") +
  scale_y_continuous(breaks = seq_along(row_ord), labels = sub_lab[rev(row_ord)], expand = expansion(add = 0.6)) +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))+ blk
p_left <- ggplot(sizes, aes(subtype, n, fill = subtype)) + geom_col(width = 0.6) +
  geom_text(aes(label = n), hjust = 1.15, size = 3) + coord_flip() +
  scale_fill_manual(values = set_cols, guide = "none") +
  scale_y_reverse(expand = expansion(mult = c(0.25, 0))) +
  labs(x = NULL, y = "") + theme_classic(base_size = 10) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())+ blk
p_A <- plot_spacer() + p_top + p_left + p_dots + plot_layout(ncol = 2, widths = c(1.1, 4), heights = c(2.5, 1.5))
ggsave(file.path(OUTF, "Fig3A_UpSet_motif_sets.png"), p_A, width = 10, height = 5, dpi = 300, bg = "white")
fwrite(dcast(bars, id + n ~ OverallBC, value.var = "k")[, .(set = id, n_motifs = n,
                                                        in_OverallBC = `overall BC significant`, not_in_OverallBC = `Not overall BC significant`)],
       file.path(OUTR, "Fig3A_motif_set_counts.csv"))

## ---------------------------------------------------------------- TF expression per motif x cell type (RNA assay)
tf_genes <- function(nm) {
  nm <- sub("\\.[0-9]+$", "", nm)                           # make.unique suffix
  nm <- gsub("\\(var\\.[0-9]+\\)", "", nm)
  nm <- gsub("EWSR1-FLI1", "EWSR1::FLI1", nm, fixed = TRUE)
  unlist(strsplit(nm, "::", fixed = TRUE))
}
motif_ref[, genes := lapply(motif.name, tf_genes)]
all_genes <- unique(unlist(motif_ref$genes))
all_genes <- all_genes[all_genes %in% rownames(obj[["RNA"]])]
cat("\nTF genes found in RNA assay:", length(all_genes), "of", length(unique(unlist(motif_ref$genes))), "\n")
DefaultAssay(obj) <- "RNA"
dp <- as.data.table(DotPlot(obj, assay = "RNA", features = all_genes)$data)   # avg.exp, pct.exp, features.plot, id, avg.exp.scaled
motif_expr <- rbindlist(lapply(seq_len(nrow(motif_ref)), function(i) {
  g <- motif_ref$genes[[i]]; g <- g[g %in% all_genes]
  if (!length(g)) return(NULL)
  dp[features.plot %in% g, .(motif = motif_ref$motif[i], avg.exp.scaled = mean(avg.exp.scaled), pct.exp = mean(pct.exp)),
     by = .(cluster = as.character(id))]
}))
fwrite(motif_expr, file.path(OUTR, "TF_expression_by_motif_celltype.csv"))

## ---------------------------------------------------------------- hull around epithelial clusters (computed once)
emb  <- Embeddings(obj, UMAP_RED)[epi_cells, 1:2]
hull <- fmesher::fm_nonconvex_hull_inla(as.matrix(emb), convex = -0.0021)
hull_seg <- data.frame(x = hull$loc[hull$idx[, 1], 1], y = hull$loc[hull$idx[, 1], 2],
                       xend = hull$loc[hull$idx[, 2], 1], yend = hull$loc[hull$idx[, 2], 2])

## ---------------------------------------------------------------- module score (Seurat AddModuleScore algorithm, version-independent)
chromvar_mat <- as.matrix(LayerData(obj, assay = "chromvar", layer = "data"))   # motifs x cells, chromVAR deviation z-scores
module_score <- function(mat, features, ctrl = 15, nbin = 24, seed = SEED) {
  set.seed(seed)
  features <- intersect(features, rownames(mat))
  avg  <- rowMeans(mat)
  bins <- ggplot2::cut_number(avg + rnorm(length(avg)) / 1e30, n = nbin, labels = FALSE, right = FALSE)
  names(bins) <- names(avg)
  ctrl_use <- unique(unlist(lapply(features, function(f)
    sample(names(bins)[bins == bins[[f]]], size = min(ctrl, sum(bins == bins[[f]])), replace = FALSE))))
  colMeans(mat[features, , drop = FALSE]) - colMeans(mat[ctrl_use, , drop = FALSE])
}

## ---------------------------------------------------------------- UMAP panels drawn from the embedding (no Seurat plotting functions)
umap_df <- data.table(Embeddings(obj, UMAP_RED)[, 1:2]); setnames(umap_df, c("UMAP_1", "UMAP_2"))
umap_df[, `:=`(cell = colnames(obj), celltype = as.character(Idents(obj)))]
lab_df  <- umap_df[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = celltype]
blk <- theme(text = element_text(colour = "black"), axis.text = element_text(colour = "black"),
             axis.title = element_text(colour = "black"), legend.text = element_text(colour = "black"))

umap_plot <- function(sc, cells = NULL, title = "", hull = TRUE, legend = TRUE,
                      cols = c("grey90", "orange"), floor0 = TRUE, min.cutoff = 0.05,
                      max.cutoff = 0.95,
                      pt.size = 0.45) {
  d <- copy(umap_df); d[, score := sc[cell]]
  if (!is.null(cells)) d <- d[cell %in% cells]
  q  <- quantile(d$score, c(min.cutoff, max.cutoff))
  lo <- if (floor0) 0 else q[[1]]                       # floor0: everything <= 0 is grey; colour ramps over (0, q95]
  d[, score := pmin(pmax(score, lo), q[[2]])]
  setorder(d, score)
  p <- ggplot(d, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = score), size = pt.size, stroke = 0) +
    scale_colour_gradientn(colours = cols, limits = c(lo, q[[2]]), name = "Module\nscore") +
    geom_text(data = lab_df[celltype %in% unique(d$celltype)], aes(label = celltype), size = 3.5, fontface = "bold") +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) + theme_classic(base_size = 11) + blk +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  if (hull)    p <- p + geom_segment(data = hull_seg, aes(x = x, y = y, xend = xend, yend = yend),
                                     colour = "black", linetype = "dotted", linewidth = 0.5)
  if (!legend) p <- p + theme(legend.position = "none")
  p
}
## ---------------------------------------------------------------- helpers
sig_stars <- function(p) fifelse(is.na(p), "NA", fifelse(p <= 0.001, "***", fifelse(p <= 0.01, "**", fifelse(p <= 0.05, "*", "ns"))))

dunn_pairs <- function(df, comps) {
  dt  <- dunn.test(df$score, df$celltype, method = "bh", altp = TRUE, kw = FALSE, table = FALSE)
  res <- data.table(Comparison = dt$comparisons, Z = dt$Z, P_unadj = dt$altP, P_adj = dt$altP.adjusted)
  res[, c("g1", "g2") := tstrsplit(Comparison, " - ", fixed = TRUE)]
  keys <- sapply(comps, paste, collapse = "|")
  res[, key := fifelse(paste(g1, g2, sep = "|") %in% keys, paste(g1, g2, sep = "|"), paste(g2, g1, sep = "|"))]
  res <- res[key %in% keys][match(keys, key)]
  res[, signif := sig_stars(P_adj)]
  res[, .(Comparison, Z, P_unadj, P_adj, signif)]
}

box_panel <- function(df, comps, dunn, palette, ylab) {
  rng <- range(df$score); step <- 0.07 * diff(rng)
  ypos <- rng[2] + step * seq_along(comps)
  ggboxplot(as.data.frame(df), x = "celltype", y = "score", color = "celltype", palette = palette,
            legend = "none", outlier.size = 0.4) +
    stat_compare_means(method = "kruskal.test", label.y = max(ypos) + 1.5 * step, size = 3) +
    geom_signif(comparisons = comps, annotations = dunn$signif, y_position = ypos,
                textsize = 3, vjust = 0.5, tip_length = 0.01) +
    labs(x = "Cell type", y = ylab, title = "") + blk +
    theme(axis.title.y = element_text(size = 8))
}

bubble_plot <- function(dat) {
  ggplot(dat, aes(x = fold.enrichment, y = reorder(motif.name, fold.enrichment))) +
    geom_point(aes(size = avg.exp.scaled, color = neglog10pval)) +
    coord_cartesian(clip = "off") +
    scale_color_gradient(low = "#F88586", high = "#FC0D0F", name = "Motif enrichment\n-log10(P value)",
                         limits = c(-log10(0.05), NA), na.value = "grey85") +
    scale_x_continuous(breaks = 1) +
    scale_size_continuous(name = "Average gene expression\n(scaled)", range = c(0.2, 6)) +
    theme_bw() +
    labs(x = xlab_bubble, y = "TF motif (* also significant in overall BC model)", title = "") +
    facet_grid(family ~ cluster, scales = "free_y", space = "free_y", labeller = label_wrap_gen(width = 28)) +
    blk +
    theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 10), legend.position = "bottom",
          strip.text.y = element_text(size = 10, angle = 0), strip.text.x = element_text(size = 10),
          panel.spacing.y = unit(0.2, "lines"))
}

## ---------------------------------------------------------------- one motif set -> bubble, counts, UMAP, box plots
run_set <- function(name, motifs) {
  cat(sprintf("\n===== %s : %d motifs =====\n", name, length(motifs)))
  if (length(motifs) < 2) { cat("  skipped (<2 motifs)\n"); return(NULL) }

  ## enrichment in DAPs of each cell type (Fig 3F for TNBC_specific)
  dat <- merge(enrich_all[motif %in% motifs], motif_ref[, .(motif, motif.name, family)], by = "motif")
  dat[, motif.name := ifelse(motif %in% sig_OverallBC, paste0(motif.name, "*"), motif.name)]
  dat <- merge(dat, motif_expr, by = c("motif", "cluster"), all.x = TRUE)
  dat[is.na(avg.exp.scaled), avg.exp.scaled := 0]
  dat[, cluster := factor(cluster, levels = ct_levels)]
  fwrite(dat, file.path(OUTR, sprintf("enrichment_allDAP_%s.csv", name)))

  fam <- dat[, .(n = uniqueN(motif)), by = family][order(family)]
  fam[, part := cumsum(n) %/% MAX_MOTIFS_PER_BUBBLE + 1]
  for (pt in unique(fam$part)) {
    d <- dat[family %in% fam[part == pt, family]]
    suffix <- if (max(fam$part) > 1) sprintf("_part%d", pt) else ""
    ggsave(file.path(OUTF, sprintf("Bubble_%s%s.jpg", name, suffix)), bubble_plot(d),
           width = 15, height = min(40, 4 + 0.22 * uniqueN(d$motif)), dpi = 300, bg = "white", limitsize = FALSE)
  }

  ## family counts
  cnt <- unique(dat[, .(motif, family)])[, .(count = .N), by = family]
  p_cnt <- ggplot(cnt, aes(x = count, y = reorder(family, count))) +
    geom_bar(stat = "identity", fill = "#FF4500") +
    labs(x = "No. of TF motifs", y = "TF motif family", title = "") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    scale_x_continuous(breaks = scales::pretty_breaks()) + theme_classic() + blk
  ggsave(file.path(OUTF, sprintf("FamilyCount_%s.jpg", name)), p_cnt,
         width = 6, height = min(20, 2 + 0.25 * nrow(cnt)), dpi = 300, bg = "white")

  ## module score on chromVAR deviations + UMAPs (Fig 3E for TNBC_specific)
  sc <- module_score(chromvar_mat, motifs, ctrl = 15)
  ggsave(file.path(OUTF, sprintf("UMAP_%s_AddModuleScore_all.jpg", name)), umap_plot(sc, title = if (name == "OverallBC_BonfSig") "Overall BC" else name),
         width = 6, height = 5, dpi = 300, bg = "white")
  a <- umap_plot(sc, cells = epi_cells,    hull = FALSE, legend = FALSE)
  b <- umap_plot(sc, cells = nonepi_cells, hull = FALSE)
  ggsave(file.path(OUTF, sprintf("UMAP_%s_AddModuleScore_split.jpg", name)), a | b,
         width = 10, height = 5, dpi = 300, bg = "white")

  ## box plots, Kruskal-Wallis + BH-adjusted two-sided Dunn (Fig 3E for TNBC_specific)
  df <- data.table(celltype = factor(as.character(Idents(obj)), levels = ct_levels), score = sc)
  kw <- kruskal.test(score ~ celltype, data = df)
  d1 <- df[celltype %in% epi];    d1[, celltype := droplevels(celltype)]
  d2 <- df[celltype %in% nonepi]; d2[, celltype := droplevels(celltype)]
  comps1 <- combn(epi, 2, simplify = FALSE); comps2 <- combn(nonepi, 2, simplify = FALSE)
  dn1 <- dunn_pairs(d1, comps1); dn2 <- dunn_pairs(d2, comps2)
  ylab <- sprintf("AddModuleScore (%s TF motifs)", gsub("_", " ", sub("^OverallBC_BonfSig$", "Overall BC", name)))
  ggsave(file.path(OUTF, sprintf("Boxplot_%s_part1.jpg", name)), box_panel(d1, comps1, dn1, "jco", ylab),
         width = 3.5, height = 5, dpi = 300, bg = "white")
  ggsave(file.path(OUTF, sprintf("Boxplot_%s_part2.jpg", name)), box_panel(d2, comps2, dn2, "YlGnBu", ylab),
         width = 4.5, height = 6, dpi = 300, bg = "white")
  fwrite(rbind(dn1[, group := "epithelial"], dn2[, group := "non-epithelial"])[, set := name],
         file.path(OUTR, sprintf("dunn_%s.csv", name)))

  list(summary = data.table(set = name, n_motifs = length(motifs),
                            n_motifs_with_expr = uniqueN(dat[avg.exp.scaled != 0, motif]),
                            kruskal_chisq = unname(kw$statistic), kruskal_p = kw$p.value),
       score = sc)
}

res <- lapply(names(motif_sets), function(n) run_set(n, motif_sets[[n]]))
names(res) <- names(motif_sets)
res <- Filter(Negate(is.null), res)
summ <- rbindlist(lapply(res, `[[`, "summary"))
cat("\n"); print(summ)
fwrite(summ, file.path(OUTR, "motif_set_summary_kruskal.csv"))
score_mat <- data.table(cell = colnames(obj), celltype = as.character(Idents(obj)))
for (n in names(res)) score_mat[, (n) := res[[n]]$score]
fwrite(score_mat, file.path(OUTR, "AddModuleScore_per_cell_all_sets.csv"))

## ---------------------------------------------------------------- Fig 3B: mean module score per cell type x motif set (PNG)
rows <- data.table(
  set   = c("OverallBC_BonfSig", "ERpos_shared", "ERneg_shared",
            "TNBC_BonfSig", "LumA_BonfSig", "LumB_BonfSig", "LumB_HER2Neg_BonfSig","HER2_BonfSig",
            "TNBC_specific", "LumA_specific", "LumB_specific", "LumB_HER2Neg_specific", "HER2_specific"),
  label = c("Overall BC", "ER+ shared", "ER- shared",
            "TNBC", "Luminal A", "Luminal B", "Luminal B/HER2-","HER2-enriched",
            "TNBC-specific", "LumA-specific", "LumB-specific", "LumB/HER2- specific", "HER2-specific"),
  group = c(1, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4))
rows <- rows[set %in% names(score_mat)]
rows[, n := set_sizes$n_motifs[match(set, set_sizes$set)]]
rows[, label := sprintf("%s (n = %d)", label, n)]
ct_show <- ct_levels                                                    # all ten cell types
M <- score_mat[celltype %in% ct_show, lapply(.SD, mean), by = celltype, .SDcols = rows$set]
M <- as.matrix(M[match(ct_show, celltype), -1]); rownames(M) <- ct_show
Mz <- t(scale(t(M), center = TRUE, scale = TRUE))                       # z-score each set across the cell types shown
Mz <- t(Mz); rownames(Mz) <- rows$label                                  # rows = sets, cols = cell types
fwrite(data.table(set = rows$set, label = rows$label, round(Mz, 3)), file.path(OUTR, "Fig3B_mean_AddModuleScore_by_celltype_z.csv"))
fwrite(data.table(set = rows$set, label = rows$label, round(t(M), 4)), file.path(OUTR, "Fig3B_mean_AddModuleScore_by_celltype_raw.csv"))

lim <- (length(ct_show) - 1) / sqrt(length(ct_show))                    # max |z| possible for n columns
colB <- colorRamp2(c(-lim, -0.02, 0.02, lim), c("#483d8b", "#D8D2EC", "#FDDBC7", "red1"))   # purple < 0 < red, no grey
comp_col <- c(Epithelial = "#E7B800", `Non-epithelial` = "#8A8A8A")
commonB <- list(col = colB, cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
                column_names_rot = 45, column_names_gp = gpar(fontsize = 9, col = "black"),
                rect_gp = gpar(col = "white", lwd = 1), column_title = NULL)
htB_epi <- do.call(Heatmap, c(list(Mz[, epi, drop = FALSE], name = "epi", width = unit(11 * length(epi), "mm"),
                                   row_split = rows$group, row_gap = unit(2, "mm"), row_title = NULL,
                                   row_names_side = "left", row_names_gp = gpar(fontsize = 9, col = "black"),
                                   top_annotation = HeatmapAnnotation(Compartment = rep("Epithelial", length(epi)), col = list(Compartment = comp_col),
                                                                      show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(2.5, "mm"))),
                              commonB))
htB_non <- do.call(Heatmap, c(list(Mz[, nonepi, drop = FALSE], name = "non-epi", width = unit(4.5 * length(nonepi), "mm"),
                                   show_row_names = FALSE,
                                   top_annotation = HeatmapAnnotation(Compartment = rep("Non-epithelial", length(nonepi)), col = list(Compartment = comp_col),
                                                                      show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(2.5, "mm"))),
                              commonB))
htB <- htB_epi + htB_non

lgdB <- packLegend(
  Legend(col_fun = colB, title = "Mean module score (z)", direction = "horizontal", title_position = "topleft",
         at = c(-round(lim, 1), 0, round(lim, 1)), legend_width = unit(3.5, "cm"),
         title_gp = gpar(fontsize = 8, col = "black"), labels_gp = gpar(fontsize = 7, col = "black")),
  Legend(labels = names(comp_col), legend_gp = gpar(fill = comp_col), title = "Cell type", direction = "horizontal", title_position = "topleft", nrow = 1,
         title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7)),
  direction = "horizontal", gap = unit(8, "mm"))
png(file.path(OUTF, "Fig3B_heatmap_AddModuleScore_sets_by_celltype.png"), width = 6.2, height = 4, units = "in", res = 300)
draw(htB, ht_gap = unit(2.5, "mm"), padding = unit(c(4, 4, 4, 4), "mm"),
     annotation_legend_list = list(lgdB), annotation_legend_side = "bottom", align_annotation_legend = "global_center")
dev.off()

cat("\ndone. figures ->", OUTF, "| results ->", OUTR, "\n")

## ---------------------------------------------------------------- Fig 3D: UMAPs in three groups, bordered panels
darken   <- function(col, f = 0.45) rgb(t(col2rgb(col) * (1 - f)), maxColorValue = 255)
set_hi   <- c(set_cols, OverallBC_BonfSig = "purple", ERpos_shared = "green3", ERneg_shared = "orange")
ramp_for <- function(s) { key <- if (s %in% names(set_hi)) s else sub("_(BonfSig|specific)$", "", s)
c("grey90", set_hi[[key]], darken(set_hi[[key]])) }

sub_lab <- c(TNBC = "TNBC", HER2 = "HER2-enriched", LumA = "Luminal A", LumB = "Luminal B", LumB_HER2Neg = "Luminal B/HER2-")
umap_set <- function(s, lab) umap_plot(res[[s]]$score, legend = FALSE, pt.size = 0.35, cols = ramp_for(s),
                                       title = sprintf("%s (n = %d)", lab, set_sizes$n_motifs[set_sizes$set == s])) +
  theme(plot.title = element_text(size = 10),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.7),   # border
        axis.line = element_blank(),
        plot.margin = margin(4, 6, 4, 6))

## (1) all Bonferroni-significant motifs per subtype
p_all <- wrap_plots(lapply(subtypes, function(s) umap_set(paste0(s, "_BonfSig"), sub_lab[[s]])), ncol = length(subtypes))
ggsave(file.path(OUTF, "Fig3C1_UMAP_allBonfSig_by_subtype.jpg"), p_all, width = 20, height = 4.3, dpi = 300, bg = "white")

## (2) subtype-specific motifs: TNBC | LumA, LumB/HER2- | LumB, HER2
umap_set <- function(s, lab, pt.size = 0.30) umap_plot(res[[s]]$score, legend = FALSE, pt.size = pt.size, cols = ramp_for(s),
                                                       title = sprintf("%s (n = %d)", lab, set_sizes$n_motifs[set_sizes$set == s])) +
  theme(plot.title = element_text(size = 10),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.7),
        axis.line = element_blank(),
        plot.margin = margin(4, 6, 4, 6))
u <- setNames(lapply(subtypes, function(s) umap_set(paste0(s, "_specific"), paste(sub_lab[[s]], "specific"),
                                                    pt.size = if (s == "TNBC") 0.4 else 0.25)), subtypes)
design <- "
AABD
AACE
"
p_spec <- u$TNBC + u$LumA + u$LumB_HER2Neg + u$LumB + u$HER2 + plot_layout(design = design)
ggsave(file.path(OUTF, "Fig3C2_UMAP_specific_by_subtype.jpg"), p_spec, width = 13, height = 5, dpi = 300, bg = "white")

## (3) overall BC and the shared programs
p_ref <- wrap_plots(list(umap_set("OverallBC_BonfSig", "Overall BC"), umap_set("ERpos_shared", "ER+ shared"), umap_set("ERneg_shared", "ER- shared")), ncol = 3)
ggsave(file.path(OUTF, "Fig3C3_UMAP_OverallBC_ERpos_ERneg.jpg"), p_ref, width = 12, height = 4.3, dpi = 300, bg = "white")

## ---------------------------------------------------------------- per-motif enrichment heatmap, transposed: cell types (rows) x motifs (columns)
enrich_heatmap_t <- function(name, cap = 100, file = NULL, width = NULL, height =4.5,
                             row_h_epi = 6, row_h_non = 3.5, fam_ncol = 6) {       # row heights in mm; family key columns
  m <- motif_sets[[name]]
  d <- merge(enrich_all[motif %in% m, .(motif, cluster, neglog10pval)], motif_ref[, .(motif, motif.name, family)], by = "motif")
  d[, label := ifelse(motif %in% sig_OverallBC, paste0(motif.name, "*"), motif.name)]
  mat <- as.matrix(dcast(d, motif ~ cluster, value.var = "neglog10pval")[, c("motif", ct_levels), with = FALSE], rownames = "motif")
  mat[is.na(mat)] <- 0
  top <- ct_levels[apply(mat, 1, which.max)]
  ord <- order(match(top, ct_levels), -apply(mat, 1, max))
  mat <- mat[ord, , drop = FALSE]; top <- top[ord]
  info <- d[match(rownames(mat), motif), .(motif, label, family)]
  tm <- t(mat); colnames(tm) <- info$label

  ## colour: flat dark blue below P = 0.05; above it, stops placed where the data are (1.3 / 10 / 50 / cap)
  thr <- -log10(0.05)
  col_fun <- colorRamp2(c(0, thr - 1e-6, thr, 10, 50, cap),
                        c("#483d8b", "#7B68B8", "#FEE0D2", "#FC9272", "#EF3B2C", "red"))

  fam_lv  <- sort(unique(info$family))
  fam_col <- setNames(grDevices::hcl.colors(length(fam_lv), "Dark 3"), fam_lv)
  comp_col <- c(Epithelial = "#E7B800", `Non-epithelial` = "#8A8A8A")
  ta <- HeatmapAnnotation(Family = info$family, col = list(Family = fam_col), show_legend = FALSE,
                          annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "left", simple_anno_size = unit(3, "mm"))
  common <- list(col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, row_title = NULL, show_heatmap_legend = FALSE,
                 row_names_side = "left", row_names_gp = gpar(fontsize = 9, col = "black"), rect_gp = gpar(col = "grey99"))
  ht_epi <- do.call(Heatmap, c(list(tm[epi, , drop = FALSE], name = "epi", top_annotation = ta,
                                    height = unit(row_h_epi * length(epi), "mm"), show_column_names = FALSE,
                                    left_annotation = rowAnnotation(Compartment = rep("Epithelial", length(epi)), col = list(Compartment = comp_col),
                                                                    show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(2.5, "mm"))),
                               common))
  ht_non <- do.call(Heatmap, c(list(tm[nonepi, , drop = FALSE], name = "non-epi",
                                    height = unit(row_h_non * length(nonepi), "mm"),
                                    column_names_rot = 60, column_names_gp = gpar(fontsize = 7, col = "black"),
                                    left_annotation = rowAnnotation(Compartment = rep("Non-epithelial", length(nonepi)), col = list(Compartment = comp_col),
                                                                    show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(2.5, "mm"))),
                               common))
  ht <- ht_epi %v% ht_non

  ## legends built by hand: row 1 = -log10 P + compartment, row 2 = family
  lgd_p <- Legend(col_fun = col_fun, title = "-log10 P", direction = "horizontal", title_position = "lefttop",
                  at = c(0, thr, 10, 50, cap), labels = c("0", "1.3", "10", "50", paste0("\u2265", cap)),
                  legend_width = unit(4, "cm"), title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
  lgd_c <- Legend(labels = names(comp_col), legend_gp = gpar(fill = comp_col), title = "Cell type", direction = "horizontal", title_position = "lefttop",
                  nrow = 1,
                  title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
  lgd_f <- Legend(labels = fam_lv, legend_gp = gpar(fill = fam_col), title = "TF family", ncol = fam_ncol,
                  title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
  lgd <- packLegend(packLegend(lgd_p, lgd_c, direction = "horizontal", gap = unit(8, "mm")), lgd_f,
                    direction = "vertical", gap = unit(3, "mm"))

  if (is.null(width)) width <- 2.5 + 0.16 * ncol(tm)
  if (is.null(file))  file  <- file.path(OUTF, sprintf("Heatmap_enrichment_%s_transposed.png", name))
  png(file, width = width, height = height, units = "in", res = 300)
  draw(ht, ht_gap = unit(2, "mm"), padding = unit(c(4, 4, 4, 4), "mm"),
       annotation_legend_list = list(lgd), annotation_legend_side = "bottom", align_annotation_legend = "global_center")
  dev.off()
  fwrite(data.table(motif = info$motif, motif.name = info$label, family = info$family, top_celltype = top, round(mat, 2)),
         file.path(OUTR, sprintf("Heatmap_enrichment_%s_transposed.csv", name)))
  invisible(tm)
}
enrich_heatmap_t("TNBC_specific")                                                  # Fig 3F
for (s in setdiff(subtypes, "TNBC")) enrich_heatmap_t(paste0(s, "_specific"))      # Supp
