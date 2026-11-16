#!/usr/bin/env Rscript
# ==============================================================================
# 11_tf_family_annotation.R
# Annotate motifs with TFClass transcription-factor families.
# Produces: reference/group_TFfamily.csv
# Run     : Rscript scripts/05_replication/11_tf_family_annotation.R   (from the repository root)
# ==============================================================================
# 11_tf_family_annotation.R — TFClass family / class for the 727 tested motifs
#
# R port of 46_tf_family_fetch.ps1, so the TF-family grouping can be regenerated
# on the HPC cluster rather than only on the Windows workstation.
#
# Produces the 94-family grouping used by the grouped TF-GWAS models.
#
# WHY PER-MATRIX CALLS AND NOT A BULK FILE
#   JASPAR publishes no annotation table for the 2022 vertebrates CORE clustering
#   (.../2022/vertebrates/CORE/radial_trees/annotation_table.txt returns 404), and
#   the list endpoint /api/v1/matrix/?collection=CORE&tax_group=vertebrates returns
#   only matrix_id, name, collection, base_id, version, sequence_logo, url — class
#   and family exist ONLY on the per-matrix endpoint. Hence 727 calls.
#   Runtime a few minutes. the compute host needs outbound access to the JASPAR API.
#
# INPUT   reference/st4_motif_ids.txt
# OUTPUT  reference/jaspar_tf_family.csv   motif_id, tf, class, family
#         reference/group_TFfamily.csv     motif_id, group  (group = family)
#
# Verify against the existing mapping rather than overwriting blind:
#   Rscript scripts/05_replication/11_tf_family_annotation.R --check
# writes to *_new.csv and diffs against the current files.

suppressPackageStartupMessages({ library(jsonlite); library(data.table) })

args   <- commandArgs(trailingOnly = TRUE)
CHECK  <- "--check" %in% args
REFDIR <- Sys.getenv("GROUP_REF", "reference")
IDFILE <- file.path(REFDIR, "st4_motif_ids.txt")
API    <- "https://jaspar.elixir.no/api/v1/matrix/%s/"
ts <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

ids <- readLines(IDFILE); ids <- ids[nzchar(ids)]
ts("motifs to fetch: ", length(ids))

fetch_one <- function(id, tries = 3L) {
  for (k in seq_len(tries)) {
    r <- tryCatch(jsonlite::fromJSON(sprintf(API, id)), error = function(e) NULL)
    if (!is.null(r)) {
      j <- function(x) if (is.null(x) || !length(x)) "" else paste(x, collapse = "; ")
      return(list(tf = j(r$name), class = j(r$class), family = j(r$family)))
    }
    Sys.sleep(1)                       # transient API failure: back off and retry
  }
  list(tf = "", class = "", family = "")
}

res <- vector("list", length(ids)); fail <- 0L
for (i in seq_along(ids)) {
  res[[i]] <- fetch_one(ids[i])
  if (!nzchar(res[[i]]$family)) fail <- fail + 1L
  if (i %% 100 == 0) ts("  ", i, "/", length(ids))
}
tab <- data.table(motif_id = ids,
                  tf     = vapply(res, `[[`, character(1), "tf"),
                  class  = vapply(res, `[[`, character(1), "class"),
                  family = vapply(res, `[[`, character(1), "family"))
ts("fetched ", nrow(tab), " | without family annotation: ", fail)

grp <- tab[, .(motif_id, group = fifelse(nzchar(family), family,
                                         paste0("unclassified_", motif_id)))]
ts("distinct TF families: ", uniqueN(grp$group))

suffix <- if (CHECK) "_new" else ""
fwrite(tab, file.path(REFDIR, sprintf("jaspar_tf_family%s.csv", suffix)))
fwrite(grp, file.path(REFDIR, sprintf("group_TFfamily%s.csv",  suffix)))
ts("written to ", REFDIR, if (CHECK) " (as *_new.csv)" else "")

if (CHECK) {
  old <- file.path(REFDIR, "group_TFfamily.csv")
  if (file.exists(old)) {
    o <- fread(old); m <- merge(o, grp, by = "motif_id", suffixes = c("_old","_new"))
    d <- m[group_old != group_new]
    cat(sprintf("\nCHECK: %d motifs compared, %d differ\n", nrow(m), nrow(d)))
    if (nrow(d)) print(head(d, 20)) else cat("mapping reproduces exactly\n")
  } else cat("\nCHECK: no existing group_TFfamily.csv to compare\n")
}
