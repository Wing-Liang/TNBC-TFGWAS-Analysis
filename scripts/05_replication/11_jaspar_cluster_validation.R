#!/usr/bin/env Rscript
# ==============================================================================
# 11_jaspar_cluster_validation.R
# Validate cluster assignments made by nearest-matrix similarity.
# Produces: reference/fallback_validation.csv
# Run     : Rscript scripts/05_replication/11_jaspar_cluster_validation.R   (from the repository root)
# ==============================================================================
# 11_jaspar_cluster_validation.R — why version-based fallback was rejected
#
# For the 35 tested motifs whose exact version JASPAR's 2022 clustering does not
# contain, an obvious shortcut is to assign each to the cluster of ANOTHER VERSION
# of the same accession ("base-ID fallback"). This script tests that shortcut by
# measuring how similar the scanned matrix actually is to the version JASPAR
# clustered, using the same metric used to co-cluster motifs (offset-aligned PCC,
# reverse complement allowed, min overlap 6).
#
# RESULT: 12 of 35 fall below the r > 0.80 bar — NKX3-1 0.611,
# RUNX1 0.649, PPARG::RXRA 0.661, ATOH1 0.680, LHX1 0.698, NR1H4 0.702,
# ELF5 0.704, ZIC3 0.718, AR 0.729, MSGN1 0.782, PRDM4 0.792, NEUROD2 0.795.
# The shortcut is therefore NOT used; 10_jaspar_cluster_assignment.R assigns by actual
# similarity instead. This also vindicates matching the primary analysis on full
# versioned accessions against the CORE *redundant* set: base-ID matching would
# have fed a materially different matrix into the scan for those motifs.
#
# SOURCES — fetched from JASPAR, not from any derived local copy (same as 44).
#
# INPUT   $GROUP_REF/st4_motif_ids.txt
# OUTPUT  $GROUP_REF/fallback_validation.csv
#
# Needs universalmotif (installed in the analysis conda environment).

suppressPackageStartupMessages({ library(universalmotif); library(data.table) })

REFDIR  <- Sys.getenv("GROUP_REF", "reference")
REFRESH <- "--refresh" %in% commandArgs(trailingOnly = TRUE)
dir.create(REFDIR, showWarnings = FALSE, recursive = TRUE)
ts <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

PFM_URL <- "https://jaspar.elixir.no/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_redundant_pfms_jaspar.txt"
CLU_URL <- "https://jaspar.elixir.no/static/clustering/2022/vertebrates/CORE/interactive_trees/clusters.tab"
PFM <- file.path(REFDIR, basename(PFM_URL))
CLU <- file.path(REFDIR, "jaspar_clusters_2022.tab")

get <- function(url, dest) {
  if (REFRESH || !file.exists(dest) || file.size(dest) == 0) {
    ts("downloading ", basename(dest))
    ok <- tryCatch({ download.file(url, dest, quiet = TRUE, mode = "wb"); TRUE },
                   error = function(e) FALSE)
    if (!ok) stop("could not fetch ", url)
  } else ts("using cached ", basename(dest))
  dest
}
get(PFM_URL, PFM); get(CLU_URL, CLU)

ids <- readLines(file.path(REFDIR, "st4_motif_ids.txt")); ids <- ids[nzchar(ids)]

memb <- rbindlist(lapply(readLines(CLU), function(l) {
  p <- strsplit(l, "\t", fixed = TRUE)[[1]]
  if (length(p) < 2 || !grepl("^cluster", p[1])) return(NULL)
  g <- regmatches(p[2], gregexpr("MA[0-9]+_[0-9]+", p[2]))[[1]]
  if (!length(g)) return(NULL)
  data.table(motif = sub("_", ".", g, fixed = TRUE), cluster = p[1])
}))

need <- setdiff(ids, memb$motif)
base_of <- function(x) sub("\\.[0-9]+$", "", x)
partner <- vapply(need, function(m) {
  cand <- memb$motif[base_of(memb$motif) == base_of(m)]
  if (length(cand)) cand[1] else NA_character_
}, character(1))
ts(length(need), " motifs lack an exact-version entry; ",
   sum(!is.na(partner)), " have another version present")

um  <- read_jaspar(PFM)
acc <- vapply(um, function(m) as.character(m["name"]),    character(1))
tfn <- vapply(um, function(m) as.character(m["altname"]), character(1))
names(um) <- acc

res <- rbindlist(lapply(seq_along(need), function(i) {
  a <- need[i]; b <- partner[i]
  if (is.na(b) || !(a %in% acc) || !(b %in% acc))
    return(data.table(scanned = a, jaspar = b, tf = NA_character_,
                      w_scan = NA_integer_, w_jas = NA_integer_,
                      pcc = NA_real_, ok = NA))
  s <- compare_motifs(list(um[[a]], um[[b]]), method = "PCC", tryRC = TRUE,
                      min.overlap = 6, score.strat = "a.mean")
  data.table(scanned = a, jaspar = b, tf = tfn[match(a, acc)],
             w_scan = ncol(um[[a]]["motif"]), w_jas = ncol(um[[b]]["motif"]),
             pcc = round(s[1, 2], 3), ok = s[1, 2] > 0.80)
}))
setorder(res, pcc)
print(res, row.names = FALSE)
fwrite(res, file.path(REFDIR, "fallback_validation.csv"))

cat(sprintf("\npairs tested: %d | clear r > 0.80: %d | BELOW 0.80: %d\n",
            sum(!is.na(res$pcc)), sum(res$ok, na.rm = TRUE), sum(!res$ok, na.rm = TRUE)))
if (any(!res$ok, na.rm = TRUE)) {
  cat("\nfallback NOT justified for:\n")
  print(res[!is.na(ok) & !ok, .(scanned, jaspar, tf, w_scan, w_jas, pcc)], row.names = FALSE)
}
