#!/usr/bin/env Rscript
# ==============================================================================
# 12_peaks_per_motif.R
# Count peaks carrying each motif in both cohorts.
# Produces: Table S8 (input)
# Inputs  : reference/st4_motif_ids.txt
#           reference/group_JASPAR122.csv
#           reference/group_TFfamily.csv
# Run     : Rscript scripts/05_replication/12_peaks_per_motif.R   (from the repository root)
# ==============================================================================
# Per-motif peak counts for both cohorts, joined to the per-motif variant counts
# from the two GWAS matrices.
#
#   peaks_replication / variants_replication = same for GSE244585
#
# Peak counts come from GetMotifData(object) -- a peaks x motifs binary matrix --
# on the same objects the two pipelines used:
#   discovery   results/singlecell/atac_motif.rds
#   replication replication/tfgwas/replication_chromvar.rds        (06b_motif_matrix.R)

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(Seurat); library(Signac)
})
setDTthreads(4)
ts <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
OUT <- "replication/tfgwas_grouped/peaks_per_motif.csv"

peaks_per_motif <- function(path, lab) {
  ts(lab, ": loading ", path)
  o <- readRDS(path)
  assay <- if ("ATAC" %in% names(o@assays)) "ATAC" else DefaultAssay(o)
  ts(lab, ": assay=", assay, "  features=", nrow(o[[assay]]))
  md <- Signac::GetMotifData(object = o, assay = assay, slot = "data")
  ts(lab, ": motif matrix ", nrow(md), " peaks x ", ncol(md), " motifs")
  v <- Matrix::colSums(md > 0)
  res <- data.table(motif = colnames(md), peaks = as.numeric(v))
  attr(res, "npeaks") <- nrow(md)
  rm(o, md); gc()
  res
}

D <- peaks_per_motif("results/singlecell/atac_motif.rds", "discovery")
nD <- attr(D, "npeaks"); setnames(D, "peaks", "peaks_discovery")
R <- peaks_per_motif("replication/tfgwas/replication_chromvar.rds", "replication")
nR <- attr(R, "npeaks"); setnames(R, "peaks", "peaks_replication")

# per-motif variant counts straight from the two GWAS matrices
varcount <- function(f, lab) {
  hdr <- names(fread(f, nrows = 0))
  mc  <- setdiff(grep("^MA", hdr, value = TRUE), "Chromatin")
  ts(lab, ": counting variants over ", length(mc), " motif columns")
  d <- fread(f, select = mc, showProgress = FALSE)
  data.table(motif = mc, v = vapply(mc, function(cc) sum(d[[cc]]), numeric(1)))
}
VD <- varcount("results/tfgwas/variant_peak_motif_subtypes_ALL.csv", "discovery")
setnames(VD, "v", "variants_discovery")
VR <- varcount("replication/tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz", "replication")
setnames(VR, "v", "variants_replication")

ids <- readLines("reference/st4_motif_ids.txt"); ids <- ids[nzchar(ids)]
map <- fread("reference/group_JASPAR122.csv")
famm <- fread("reference/group_TFfamily.csv")
tab <- data.table(motif = ids)
tab <- merge(tab, D,  by = "motif", all.x = TRUE)
tab <- merge(tab, VD, by = "motif", all.x = TRUE)
tab <- merge(tab, R,  by = "motif", all.x = TRUE)
tab <- merge(tab, VR, by = "motif", all.x = TRUE)
tab <- merge(tab, map[,  .(motif = motif_id, jaspar122 = group)], by = "motif", all.x = TRUE)
tab <- merge(tab, famm[, .(motif = motif_id, tf_family = group)], by = "motif", all.x = TRUE)
tab[, `:=`(peak_ratio = round(peaks_replication/peaks_discovery, 3),
           variant_ratio = round(variants_replication/variants_discovery, 3))]
setcolorder(tab, c("motif","jaspar122","tf_family","peaks_discovery","variants_discovery",
                   "peaks_replication","variants_replication","peak_ratio","variant_ratio"))
setorder(tab, -variants_discovery)
fwrite(tab, OUT)

cat("\n=== totals ===\n")
cat(sprintf("peaks in object      : discovery %s | replication %s\n",
            format(nD, big.mark=","), format(nR, big.mark=",")))
cat(sprintf("summed motif-peak hits: discovery %s | replication %s\n",
            format(sum(tab$peaks_discovery, na.rm=TRUE), big.mark=","),
            format(sum(tab$peaks_replication, na.rm=TRUE), big.mark=",")))
cat(sprintf("motifs per peak       : discovery %.1f | replication %.1f\n",
            sum(tab$peaks_discovery, na.rm=TRUE)/nD,
            sum(tab$peaks_replication, na.rm=TRUE)/nR))
cat(sprintf("\nmedian peak ratio %.3f | median variant ratio %.3f\n",
            median(tab$peak_ratio, na.rm=TRUE), median(tab$variant_ratio, na.rm=TRUE)))
cat("\n=== head ===\n"); print(head(tab, 12))
ts("written -> ", OUT)
