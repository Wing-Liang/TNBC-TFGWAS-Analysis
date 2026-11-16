#!/usr/bin/env Rscript
# ==============================================================================
# 10_jaspar_cluster_assignment.R
# Assign the 727 motifs to JASPAR 2022 matrix clusters.
# Produces: reference/group_JASPAR122.csv
# Run     : Rscript scripts/05_replication/10_jaspar_cluster_assignment.R   (from the repository root)
# ==============================================================================
# 10_jaspar_cluster_assignment.R — build the JASPAR122 grouping for the 727 tested motifs
#
# Places every tested motif into JASPAR's PUBLISHED 2022 vertebrates CORE
# matrix-clustering (RSAT; Castro-Mondragon et al. 2017), giving 122 clusters.
#
#   692 of 727 match at exact versioned accession.
#    35 do not — JASPAR clustered a different version of the same accession.
#       Each is assigned to the cluster of its MOST SIMILAR clustered matrix,
#       using the same offset-aligned Pearson correlation (reverse complement
#       allowed, min overlap 6). Assigning by the clustered version of the same
#       accession was not used: for 12 of the 35 the versions differ substantially
#       (e.g. RUNX1 .1 vs .2, r = 0.649; ELF5 0.704; AR 0.729).
#
# SOURCES — fetched from JASPAR, not from any derived local copy:
#   PFMs     https://jaspar.elixir.no/download/data/2022/CORE/
#              JASPAR2022_CORE_vertebrates_redundant_pfms_jaspar.txt
#            The REDUNDANT set is required: it retains historical matrix versions,
#            and 494 of the 727 tested IDs are .1 versions (the scan used
#            JASPAR2020-era accessions). Versioned accessions are immutable, so
#            MA0002.1 here is the matrix that was scanned.
#   Clusters https://jaspar.elixir.no/static/clustering/2022/vertebrates/CORE/
#              interactive_trees/clusters.tab
#            (2022 is the closest published clustering; no 2020 clustering exists.
#             Token format is MA<base>_<version>, not MA<base>.<version>.)
#
# Downloads are cached in GROUP_REF so repeat runs are offline; --refresh re-fetches.
#
# INPUT   $GROUP_REF/st4_motif_ids.txt        the 727 tested versioned accessions
# OUTPUT  $GROUP_REF/group_JASPAR122.csv      motif_id, group
#         $GROUP_REF/ST4_motifs_to_JASPAR2022_clusters.csv   with assignment provenance
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
    if (!ok || !file.exists(dest) || file.size(dest) == 0)
      stop("could not fetch ", url)
  } else ts("using cached ", basename(dest))
  dest
}
get(PFM_URL, PFM); get(CLU_URL, CLU)

ids <- readLines(file.path(REFDIR, "st4_motif_ids.txt")); ids <- ids[nzchar(ids)]
ts("tested motifs: ", length(ids))

## ---- published clustering: token MA<base>_<version> -> MA<base>.<version> -----
memb <- rbindlist(lapply(readLines(CLU), function(l) {
  p <- strsplit(l, "\t", fixed = TRUE)[[1]]
  if (length(p) < 2 || !grepl("^cluster", p[1])) return(NULL)
  g <- regmatches(p[2], gregexpr("MA[0-9]+_[0-9]+", p[2]))[[1]]
  if (!length(g)) return(NULL)
  data.table(motif = sub("_", ".", g, fixed = TRUE), cluster = p[1])
}))
ts("published clustering: ", uniqueN(memb$cluster), " clusters over ",
   nrow(memb), " matrices")

exact <- ids[ids %in% memb$motif]
need  <- setdiff(ids, memb$motif)
ts("exact-version matches: ", length(exact), " | to assign by similarity: ", length(need))

## ---- nearest-neighbour assignment for the remainder --------------------------
um  <- read_jaspar(PFM)
acc <- vapply(um, function(m) as.character(m["name"]), character(1))
if (!all(grepl("^MA[0-9]", acc)))
  stop("PFM headers did not parse to accessions in 'name'")
names(um) <- acc
tfn <- vapply(um, function(m) as.character(m["altname"]), character(1))

ref <- memb$motif[memb$motif %in% acc]
miss <- setdiff(c(need, ref), acc)
if (length(miss)) stop(length(miss), " accessions absent from the PFM file: ",
                       paste(head(miss, 5), collapse = ", "))

ts("comparing ", length(need), " motifs against ", length(ref), " clustered matrices")
sim <- compare_motifs(um[c(need, ref)], method = "PCC", tryRC = TRUE,
                      min.overlap = 6, score.strat = "a.mean")
rownames(sim) <- colnames(sim) <- c(need, ref)
blk <- sim[need, ref, drop = FALSE]; blk[is.na(blk)] <- -1

asg <- rbindlist(lapply(need, function(m) {
  v <- blk[m, ]; j <- which.max(v)
  data.table(motif_id = m, jaspar2022_cluster = memb$cluster[match(ref[j], memb$motif)],
             assignment = "nearest_neighbour", nearest = ref[j],
             nearest_tf = tfn[match(ref[j], acc)], pcc = round(v[j], 3))
}))
ts("assignment similarity range: ", sprintf("%.3f - %.3f", min(asg$pcc), max(asg$pcc)))

out <- rbind(
  data.table(motif_id = exact,
             jaspar2022_cluster = memb$cluster[match(exact, memb$motif)],
             assignment = "exact_version", nearest = "", nearest_tf = "", pcc = NA_real_),
  asg)[match(ids, motif_id)]

stopifnot(nrow(out) == length(ids), !anyNA(out$jaspar2022_cluster))
fwrite(out, file.path(REFDIR, "ST4_motifs_to_JASPAR2022_clusters.csv"))
fwrite(out[, .(motif_id, group = jaspar2022_cluster)],
       file.path(REFDIR, "group_JASPAR122.csv"))
ts("FINAL: ", nrow(out), " motifs -> ", uniqueN(out$jaspar2022_cluster), " clusters")
print(out[, .N, by = assignment])
