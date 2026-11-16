#!/usr/bin/env Rscript
# ==============================================================================
# 00_describe_data.R
# Describe the replication data: cells, fragments and pools.
# Produces: data check
# Run     : Rscript scripts/05_replication/00_describe_data.R   (from the repository root)
# ==============================================================================
# 00_describe_data.R — descriptive summary of downloaded data (both cohorts).
# Cheap metrics (cells, peaks, sizes) computed LIVE; fragment counts read from
# cache (they take minutes on multi-GB files).
# Usage: Rscript 00_describe_data.R [--count] [--verify]
#   --count   recompute fragment counts live (slow: reads all fragment files)
#   --verify  gzip -t integrity check on fragment files (slow)

## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
args     <- commandArgs(trailingOnly = TRUE)
do_count <- "--count"  %in% args
do_verify<- "--verify" %in% args

suppressPackageStartupMessages(library(GenomicRanges))  # so readRDS(GRanges) length() works

ORIG <- TNBC_ROOT
NEW  <- file.path(TNBC_ROOT, "replication")
setwd(NEW)

## tee all output to a script-named log (split=TRUE -> console AND file)
dir.create("logs", showWarnings = FALSE)
.logcon <- file(file.path(NEW, "logs", "00_describe_data.log"), open = "wt")
sink(.logcon, split = TRUE)
on.exit({ sink(); close(.logcon) }, add = TRUE)

## cached fragment counts (refresh with --count)
frag_total <- c(orig=528404395, pool9=362591080, pool11=326699715, pool12=439342349)
frag_valid <- c(pool9=263563523, pool11=265658436, pool12=246998362)
## pool11 GEO file is truncated (chrX/chrY/chr9-tail missing). TRUE totals from per_barcode_metrics:
frag_true_total <- c(pool11=339089275)   # vs recovered 326,699,715 -> 96.3%
frag_true_valid <- c(pool11=275224947)   # vs recovered 265,658,436 -> 96.5%

cm   <- function(x) if (is.na(suppressWarnings(as.numeric(x)))) as.character(x) else formatC(as.numeric(x), format="d", big.mark=",")
sz   <- function(f) if (file.exists(f)) sprintf("%.1f GB", file.info(f)$size/1e9) else "MISSING"
nlz  <- function(f) if (file.exists(f)) length(readLines(f, warn=FALSE)) else NA          # lines in (gz) file
cfg  <- function(key) { L <- readLines(file.path(ORIG,"config.yaml"), warn=FALSE); v <- grep(paste0("^",key,":"), L, value=TRUE); if(length(v)) trimws(sub(".*:","",v[1])) else "NA" }
countfrags <- function(f) if (file.exists(f)) as.integer(system(sprintf("zcat %s 2>/dev/null | grep -vc '^#'", shQuote(f)), intern=TRUE)) else NA
grn  <- function(rds) if (file.exists(rds)) length(readRDS(rds)) else NA                   # length of GRanges RDS

rule <- function() cat(strrep("=", 78), "\n", sep="")
rule()
cat(sprintf(" DATA SUMMARY  -  TNBC single-cell multiome cohorts        %s\n", Sys.Date()))
rule()

## ---------------- ORIGINAL ----------------
cat("\n## ORIGINAL COHORT  (GSE162798, 4 donors, hg19, SNuBar multiome)\n")
cat("   location:", ORIG, "\n")
ofrag <- file.path(ORIG, "data/GSE162798/fragments/fragments.tsv.gz")
ft    <- if (do_count) countfrags(ofrag) else frag_total["orig"]
cat(sprintf("   fragments file   : %s\n", sz(ofrag)))
cat(sprintf("   fragments        : %s\n", cm(ft)))
cat(sprintf("   ATAC barcodes    : %s\n", cm(nlz(file.path(ORIG,"data/GSE162798/scATAC/barcodes.tsv.gz")))))
cat(sprintf("   analysis cells   : %s  (QC'd)\n", cm(cfg("num_cells"))))
cat(sprintf("   peaks (MACS2)    : %s\n", cm(cfg("num_features"))))

## ---------------- NEW ----------------
cat("\n## NEW COHORT  (GSE244585, European non-Ashkenazi, 15 donors, hg38, 10x Multiome)\n")
cat("   location:", NEW, "\n")
pools <- c("pool9","pool11","pool12")
nm    <- c(pool9="GSM7821199_KTBpool9_atac_fragments.tsv.gz",
           pool11="data/fragments_hg38/pool11_atac_fragments.hg38.tsv.gz",
           pool12="GSM7821201_KTBpool12_atac_fragments.tsv.gz")
hdr <- sprintf("   %-7s %10s %17s %17s %9s %11s %-12s %s", "Pool","is_cell","TotalFrags","ValidCellFrags","Peaks","scRNA_src","GEO_h5","ATAC_integ")
cat(hdr, "\n"); cat("  ", strrep("-", nchar(hdr)-3), "\n", sep="")
tot_cell <- 0; tot_valid <- 0
intlog <- c(sprintf("# h5/ATAC integrity check  %s", Sys.time()),
            "pool\tRNA_actual\tRNA_stored_eof\tRNA_opens\tATAC_chr9tail_frags\tGEO_h5\tATAC_integ")
for (p in pools) {
  cells <- nlz(file.path("geo_cells", paste0(p,"_geo_cells.txt")))
  pk    <- grn(file.path("peaks", paste0("peaks_", p, ".rds")))
  tf    <- if (do_count) countfrags(nm[p]) else frag_total[p]
  vf    <- if (do_count) countfrags(file.path("cell_fragments", paste0(p,"_cells.fragments.tsv.gz"))) else frag_valid[p]
  # scRNA now sourced uniformly from CELLxGENE (data/cellxgene/rna_replication/); GEO h5 status shown live for reference
  rna_src <- "CELLxGENE"                                     # uniform scRNA source (all pools)
  geo_h5 <- c(pool9="GSM7821215_KTBpool9_filtered_feature_bc_matrix.h5",
              pool11="GSM7821216_KTBpool11_filtered_feature_bc_matrix.h5",
              pool12="GSM7821217_KTBpool12_filtered_feature_bc_matrix.h5")[p]
  hi <- system(paste("python3 00_check_h5_info.py", shQuote(geo_h5)), intern=TRUE)   # actual / stored_eof / opens(1)
  h5s <- if (hi[3]=="1") "GEO-ok" else sprintf("corrupt %dM/%dM", round(as.numeric(hi[1])/1e6), round(as.numeric(hi[2])/1e6))
  # ATAC integrity: autosome-relevant only (chr9-tail; chrX/chrY ignored). tabix range read = fast.
  atac_f <- c(pool9="GSM7821199_KTBpool9_atac_fragments.tsv.gz",
              pool11="data/fragments_hg38/pool11_atac_fragments.hg38.tsv.gz",
              pool12="GSM7821201_KTBpool12_atac_fragments.tsv.gz")[p]
  n9tail <- length(system(paste("tabix", shQuote(atac_f), "chr9:138000000-138394717"), intern=TRUE))
  intg <- if (n9tail > 0) "OK" else "chr9-trunc*"
  intlog <- c(intlog, sprintf("%s\t%s\t%s\t%s\t%d\t%s\t%s", p, hi[1], hi[2], hi[3], n9tail, h5s, intg))
  tf_s  <- if (p=="pool11") paste0(cm(tf),"*") else cm(tf)
  vf_s  <- if (p=="pool11") paste0(cm(vf),"*") else cm(vf)
  cat(sprintf("   %-7s %10s %17s %17s %9s %11s %-12s %s\n", p, cm(cells), tf_s, vf_s, cm(pk), rna_src, h5s, intg))
  tot_cell  <- tot_cell + cells; tot_valid <- tot_valid + vf
}
dir.create("logs", showWarnings=FALSE)
writeLines(intlog, "logs/h5_integrity.log")
cat(sprintf("   [integrity log -> logs/h5_integrity.log]\n"))
cat(sprintf("   %-7s %10s %16s %16s\n", "TOTAL", cm(tot_cell), "", cm(tot_valid)))

cat(sprintf("   * pool11 GEO file truncated (chrX/chrY/chr9>133.35Mb absent). Shown values are RECOVERED.
"))
cat(sprintf("     true total=%s (recovered %s = 96.3%%);  true valid-cell=%s (recovered %s = 96.5%%).
",
            cm(frag_true_total["pool11"]), cm(frag_total["pool11"]), cm(frag_true_valid["pool11"]), cm(frag_valid["pool11"])))
cat(sprintf("     missing fragments are entirely in chrX/chrY/chr9-tail, all excluded from the autosomal analysis.
"))
cat(sprintf("   scRNA + cell types: ALL pools from CELLxGENE (data/cellxgene/rna_replication/ 35,477 genes x 12,469 cells;
     obs$author_cell_type). GEO_h5 column = GEO file status only (pool11 corrupt) - NOT the source used.
"))

cb  <- "peaks/peaks_combined.rds"
cat(sprintf("\n   combined peak set (autosomes, hg38): %s\n", cm(grn(cb))))
obj <- "peaks/replication_combined.rds"
if (file.exists(obj)) {
  o <- readRDS(obj)
  cat(sprintf("   ChromatinAssay object: %s  [%d peaks x %d cells]  %s\n",
              sz(obj), nrow(o), ncol(o), obj))
  print(table(o$pool))
} else cat("   ChromatinAssay object: not yet built (04_build_combined_object.R)\n")

if (do_verify) {
  cat("\n## INTEGRITY (gzip -t)\n")
  for (p in pools) {
    ok <- system(sprintf("gzip -t %s 2>/dev/null", shQuote(nm[p]))) == 0
    cat(sprintf("   %-7s fragments: %s\n", p, if (ok) "OK" else "FAILED"))
  }
}
cat(sprintf("\n   (fragment counts %s; --count to recompute; --verify for integrity)\n",
            if (do_count) "computed live" else "from cache"))
rule()
