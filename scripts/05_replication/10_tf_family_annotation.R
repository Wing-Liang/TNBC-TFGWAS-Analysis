#!/usr/bin/env Rscript
# ==============================================================================
# 10_tf_family_annotation.R
# Annotate motifs with TFClass transcription-factor families.
# Produces: reference/group_TFfamily.csv
# Run     : Rscript scripts/05_replication/10_tf_family_annotation.R   (from the repository root)
# ==============================================================================
# 10_tf_family_annotation.R — TFClass family / class for the 727 tested motifs
#
# Produces the 94-family grouping used by the TF-family models and Table S8.
#
# TF class and family are available only from the per-matrix endpoint of the JASPAR
# REST API, so the script makes one call per motif (727 calls, a few minutes; the
# host needs outbound access to the JASPAR API).
#
# INPUT   reference/st4_motif_ids.txt
# OUTPUT  reference/jaspar_tf_family.csv   motif_id, tf, class, family
#         reference/group_TFfamily.csv     motif_id, group  (group = family)
#
# Verify against the existing mapping rather than overwriting blind:
#   Rscript scripts/05_replication/10_tf_family_annotation.R --check
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
