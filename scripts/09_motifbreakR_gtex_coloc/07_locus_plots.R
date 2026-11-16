# ==============================================================================
# 07_locus_plots.R
# LocusCompare and regional association plots per locus.
# Produces: Fig 5C
# Inputs  : results/gwas/bcac_gwas_ldblock.csv
#           results/motifbreakR/motifbreakR_TNBC_all.csv
#           reference/promoter_snp_genes.txt
# Run     : Rscript scripts/09_motifbreakR_gtex_coloc/07_locus_plots.R   (from the repository root)
# ==============================================================================
## 07_locus_plots.R — one locuscompare + regional figure PER LOCUS, not per gene.
##
## One figure per gene would repeat the same window: MDM4 and PIK3C2B share 1q32; SPDYA/WDR43/SNORD53/SNORD92 all
## share 2p23; nine genes share 19p13.11. The GWAS panel and the gene track are identical
## within a locus -- only the eQTL panel differs -- so this script draws each locus once, uses
## the best-colocalising gene for the eQTL panel, and titles by locus.
##
## Nothing is discarded: the left column carries the full PP0-PP4 table for EVERY gene tested
## at that locus, so the per-gene evidence 32 spread over nine figures sits on the one figure.
##
## BUILD: hg19 throughout -- the native build of the BCAC GWAS, of motifbreakR, and of
## EnsDb.Hsapiens.v75. GTEx (GRCh38) joins to the GWAS BY rsID (100% rsID coverage in
## QTD000211), so there is no liftOver and no crosswalk loss.
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
GTEX_DIR <- Sys.getenv("GTEX_DIR", "data/GTEx")
suppressMessages({library(data.table); library(locuszoomr); library(EnsDb.Hsapiens.v75)
                  library(ggplot2); library(cowplot)})
setwd(TNBC_ROOT)
TOKEN <- Sys.getenv("LDLINK_TOKEN")
FLANK <- 3e5L; LABEL_N <- 5L; P12 <- 5e-5
## OUT holds input tables and the LD cache; FIG is where figures go. Keeping these apart stops
## PNGs accumulating in the results directory, which is what had happened (30 PNGs in
## results_celltype/ against 11 in figures/).
OUT   <- "results/coloc/"
FIG   <- "results/coloc/figures"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
GZ    <- file.path(GTEX_DIR, "eQTL_v8/QTD000211.all.tsv.gz")
cn    <- c("mtid","chr","pos","ref","alt","variant","ma_samples","maf","pvalue","beta","se",
           "type","ac","an","r2","mtoid","gene_id","median_tpm","rsid")
chrn  <- function(x) as.integer(sub("^chr", "", as.character(x)))

## EnsDb.Hsapiens.v75 is Ensembl 75 and predates some HGNC renames. Only the annotation
## LOOKUP is aliased; titles, tables and filenames keep the current symbol.
ALIAS    <- c(POGLUT3 = "KDELC2")   # ENSG00000178202, chr11:108,342,832
ens_name <- function(g) if (g %in% names(ALIAS)) ALIAS[[g]] else g

## ---- GWAS (hg19) + rsID ----
bc <- fread("results/gwas/bcac_gwas_ldblock.csv",
            select = c("Chr", "Pos", "TNBC_pval", "chrpos37"))
si <- fread("snp_index.txt", select = c("ID", "CHR", "POS"))
si <- si[, .(chrpos37 = paste0("chr", CHR, ":", POS), rsid = ID)]
gw <- merge(bc, unique(si), by = "chrpos37")
gw <- gw[!is.na(TNBC_pval) & TNBC_pval > 0][!duplicated(rsid)]
gw[, `:=`(chrom = chrn(Chr), pos = as.integer(Pos), p = TNBC_pval)]
cat("GWAS variants with rsID:", nrow(gw), "\n")

## ---- motifbreakR TNBC-significant SNPs, PROMOTER ONLY (hg19, rsID) ----
## Restricted to annotation_category == "Promoter". This makes the
## figures consistent with how the colocalisation gene set was defined in the first place:
## 04_coloc_gtex_v8.R already picked its anchor genes with grepl("Promoter", annotation_category).
## The restriction takes 209 TNBC-significant SNPs down to 57, and those 57 map to exactly the
## 20 genes listed in reference/promoter_snp_genes.txt.
## Categories in this table are unbinned ("Promoter", not "Promoter (<=1kb)"), so an equality
## test and a grepl give the same 57 SNPs.
mb  <- fread("results/motifbreakR/motifbreakR_TNBC_all.csv")
mb[, sigflag := gsub("\"", "", TNBC_sig)]
sig <- unique(mb[sigflag == "Y" & annotation_category == "Promoter",
                 .(rsid = SNP_id, chrom = chrn(seqnames), pos = as.integer(start))])
cat("motifbreakR TNBC-sig PROMOTER SNPs:", uniqueN(sig$rsid),
    "(of", uniqueN(mb[sigflag == "Y"]$SNP_id), "at all annotations)\n")

## ---- GTEx v8 DAP-G 95% credible sets (script 30) ----
## Members of fine-mapped breast eQTL signals with SPIP > 0.95. Overlaid on the eQTL panel and
## the locuscompare scatter so the reader can see whether the GWAS signal coincides with the
## fine-mapped causal set, not merely with a nominal eQTL. hg19 positions via the BCAC crosswalk.
CSF <- file.path(OUT, "GTExv8_DAPG_breast_CS95_membership_full.csv")
CS  <- if (file.exists(CSF)) fread(CSF)[!is.na(pos19)] else data.table()
if (nrow(CS)) cat("DAP-G CS95 members loaded:", nrow(CS), "| genes:", uniqueN(CS$gene), "
") else
  cat("DAP-G CS95 file absent - credible sets will not be drawn
")

## ---- loci, taken from the coloc run that matches the plotted eQTL resource ----
## coloc_GTEx_bulk_strandharmonized_300kb.csv is the strand-harmonised GTEx bulk +/-300 kb run,
## i.e. the same resource and window as the eQTL panel. Its `locus` column defines the grouping.
co  <- fread(file.path(OUT, "coloc_GTEx_bulk_strandharmonized_300kb.csv"))[p12 == P12]
## The 20 promoter-SNP genes (Ensembl ID, symbol).
ens <- fread(file.path(TNBC_ROOT, "reference/promoter_snp_genes.txt"),
             header = FALSE, col.names = c("ensg", "gene"))
co  <- merge(co, ens, by = "gene")
setorder(co, locus, -PP4)
## Representative gene per locus = highest PP4, i.e. the gene whose eQTL is most likely to share
## a causal variant with the GWAS. That makes it the informative eQTL panel for the locus.
reps <- co[, .(gene = gene[1], ensg = ensg[1], n_tested = .N), by = locus]

## Genes with no coloc row are not tested eGenes in GTEx breast (QTD000211) -- verified by a
## whole-chromosome scan of column 1, not just the plotting window. SNORD53 and SNORD92 fall
## inside the 2p23 window already drawn for SPDYA, so only TERT needs its own eQTL-free figure.
untested <- setdiff(ens$gene, co$gene)
cat("not tested in GTEx breast:", paste(untested, collapse = ", "), "\n")
reps <- rbind(reps, data.table(locus = "5p15.33", gene = "TERT",
                               ensg = ens[gene == "TERT"]$ensg, n_tested = 0L))

## `gene` is the LOCUS gene (window centre, title, motifbreakR table). `eg` is the gene whose
## eQTL is drawn. They coincide except where the anchor gene is not a tested eGene: there the
## same best-PP4 rule is applied to the locus-wide bulk coloc (script 34) so the eQTL panel still
## shows the association in the GWAS region. At 5p15.33 that is MRPL36 (PP4 0.60 at +/-300 kb,
## p12 = 5e-5 -- prior-sensitive, 0.23 at 1e-5), which the panel labels as such.
reps[, `:=`(eg = gene, eensg = ensg)]
for (i in seq_len(nrow(reps))) if (reps$n_tested[i] == 0L) {
  f <- file.path(OUT, sprintf("coloc_%s_GTEx_bulk_allgenes.csv", reps$locus[i]))
  if (file.exists(f)) {
    A <- fread(f)[p12 == P12 & flank_kb == 300][order(-PP4)]
    if (nrow(A)) { reps$eg[i] <- A$symbol[1]; reps$eensg[i] <- A$ensg[1]
      cat(sprintf("%s: anchor %s untestable -> eQTL panel uses locus-wide best %s (PP4 %.3f)\n",
                  reps$locus[i], reps$gene[i], A$symbol[1], A$PP4[1])) }
  }
}

## Chromosome per locus. Parsed from the cytoband label, with the one label that does not
## start with its chromosome number handled explicitly.
LOCUS_CHR <- function(lc) if (lc == "6p21_HLA") 6L else
  as.integer(sub("^([0-9]+).*$", "\\1", lc))

ld_bin <- function(r2) cut(r2, c(-Inf, .2, .4, .6, .8, Inf),
                           labels = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">0.8"))
LDCOL  <- c("<0.2" = "#4A6FE3", "0.2-0.4" = "#63C7E6", "0.4-0.6" = "#4CAF50",
            "0.6-0.8" = "#F5A623", ">0.8" = "#D0021B")

## Overlay marking EVERY motifbreakR SNP with a black triangle. Built from the PLOT's own data:
## gg_scatter() rescales pos from bp to Mb internally, so a layer built on loc$data lands
## off-scale and every point is dropped with only a "Removed N rows" warning.
## gg_scatter() REMOVES the index SNP from the plot's main data frame and draws it as its own
## one-row layer. Any overlay built from g$data therefore silently skips the index SNP -- which
## is exactly the variant that matters most here (at 6p21 the index IS the PIP 0.999 credible-set
## member and the only motifbreakR SNP). panel_data() puts it back by unioning g$data with the
## single-row point layers, all of which carry pos (already rescaled to Mb) and logP.
panel_data <- function(g) {
  d <- g$data
  keep <- intersect(c("rsid", "pos", "logP"), names(d))
  out <- d[, keep, drop = FALSE]
  for (ly in g$layers) {
    dd <- ly$data
    if (is.data.frame(dd) && nrow(dd) && all(keep %in% names(dd)))
      out <- rbind(out, dd[, keep, drop = FALSE])
  }
  out[!duplicated(out$rsid), , drop = FALSE]
}

## Credible-set overlay: open diamonds, ringed in dark red, sized by PIP. Drawn from the PLOT's
## own data because gg_scatter() rescales pos to Mb (same trap as add_mb).
add_cs <- function(g, cs_rsids, pip) {
  if (!length(cs_rsids)) return(g)
  P <- panel_data(g); D <- P[P$rsid %in% cs_rsids, , drop = FALSE]
  if (!nrow(D)) return(g)
  D$pip <- pip[match(D$rsid, names(pip))]
  g + geom_point(data = D, aes(x = pos, y = logP, size = pip), shape = 23,
                 fill = NA, colour = "#8B0000", stroke = 1.3, inherit.aes = FALSE) +
      scale_size_continuous(range = c(3, 7), limits = c(0, 1), name = "DAP-G PIP",
                            breaks = c(0.1, 0.5, 0.9))
}

add_mb <- function(g, rs) {
  P <- panel_data(g); D <- P[P$rsid %in% rs, , drop = FALSE]
  if (!nrow(D)) return(g)
  g + geom_point(data = D, aes(x = pos, y = logP), shape = 24, size = 2.4,
                 fill = NA, colour = "black", stroke = .7, inherit.aes = FALSE)
}

## Fixed-width PP0-PP4 table for every gene tested at the locus, rendered as a panel.
## Built as a plain ggplot rather than cowplot::ggdraw(): ggdraw() leaves plot.background
## transparent, and once nested inside two plot_grid() calls that region composites to solid
## black in the PNG, hiding the text entirely. An explicit white rect avoids that.
text_panel <- function(txt, note = NULL, size = 7.4, family = "mono", colour = "black") {
  g <- ggplot() +
    annotate("text", x = 0, y = 1, label = txt, hjust = 0, vjust = 1,
             size = size / .pt, family = family, colour = colour) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(6, 6, 6, 10))
  if (!is.null(note))
    g <- g + annotate("text", x = 0, y = 0.02, label = note, hjust = 0, vjust = 0,
                      size = 7 / .pt, colour = "grey35")
  g
}

## motifbreakR panel for loci with NO bulk eQTL (currently only 5p15.33 / TERT). Rather than
## leave the left column as a placeholder, show the evidence the locus does have: every
## TNBC-significant PROMOTER SNP assigned to the gene, with its GWAS strength, LD to the index,
## the TF motifs it disrupts and the direction. This is bulk-eQTL-independent by construction.
## Locus-level table of every TNBC-significant PROMOTER motifbreakR SNP in the plotting window
## (the triangles), with its assigned promoter gene, GWAS strength, LD to the index SNP, the TF
## motifs it disrupts and the direction. Shown under the locuscompare scatter at every locus
## Capped at MB_ROWS rows by GWAS P; 19p13.11 has 23 in the window.
MB_ROWS <- 12L
mb_panel <- function(lc, mbs, d) {
  t <- mb[SNP_id %in% mbs & sigflag == "Y" & annotation_category == "Promoter"]
  if (!nrow(t)) return(text_panel(sprintf("%s: no promoter motifbreakR SNPs in window", lc),
                                  size = 9, family = ""))
  t <- merge(t, d[, .(SNP_id = rsid, p, ld)], by = "SNP_id", all.x = TRUE)
  idx <- d$rsid[which.min(d$p)]
  s <- t[order(-abs(alleleDiff)),
         .(pos = start[1], gene = SYMBOL[1], gwas = -log10(p[1]), r2 = ld[1],
           dTSS = distanceToTSS[1], n_tf = .N, tfs = paste(unique(geneSymbol), collapse = ","),
           maxdiff = alleleDiff[which.max(abs(alleleDiff))]), by = SNP_id][order(-gwas)]
  genes <- s[, .N, by = gene][order(-N)]
  hdr  <- sprintf("%-11s %-9s %-8s %5s %5s %5s %-19s %6s", "SNP", "hg19 pos", "gene", "GWAS",
                  "r2", "dTSS", "TF motifs disrupted", "dAll")
  rows <- sprintf("%-11s %-9d %-8s %5.1f %5s %5d %-19s %+6.2f",
                  ifelse(s$SNP_id == idx, paste0(s$SNP_id, "*"), s$SNP_id), s$pos,
                  substr(s$gene, 1, 8), s$gwas, ifelse(is.na(s$r2), "NA", sprintf("%.2f", s$r2)),
                  s$dTSS, substr(s$tfs, 1, 19), s$maxdiff)
  shown <- head(rows, MB_ROWS)
  more  <- if (nrow(s) > MB_ROWS) sprintf("... %d more (lower GWAS P)", nrow(s) - MB_ROWS) else NULL
  txt  <- paste(c(sprintf("%s: motifbreakR promoter SNPs in window", lc),
                  sprintf("%d TNBC-significant promoter SNPs; promoter gene(s): %s", nrow(s),
                          paste(sprintf("%s (%d)", genes$gene, genes$N), collapse = ", ")), "",
                  hdr, strrep("-", nchar(hdr)), shown, more), collapse = "\n")
  note <- paste0("GWAS = -log10 P (TN-GWAS); r2 = LD to index SNP (*); dTSS = bp to the assigned\n",
                 "gene's TSS; dAll = largest allele score difference (alt - ref); ",
                 "all effects 'strong'.")
  text_panel(txt, note)
}

## Locus-wide bulk coloc (script 34) for loci whose anchor gene is not a tested eGene. Shows
## every GTEx eGene in the window ranked by PP4, so "no coloc" is demonstrated, not asserted.
pp_panel_locuswide <- function(lc) {
  f <- file.path(OUT, sprintf("coloc_%s_GTEx_bulk_allgenes.csv", lc))
  if (!file.exists(f)) return(text_panel(
    paste0("No gene at this locus is a tested eGene\n",
           "in GTEx breast mammary tissue (QTD000211)"), size = 9, family = "", colour = "grey25"))
  ## Compact note only -- the full locus-wide table is a supplementary table (script 35).
  A  <- fread(f)[p12 == P12 & flank_kb == 300][order(-PP4)]
  A1 <- fread(f)[p12 == 1e-5 & flank_kb == 300]
  n5 <- sum(A$PP4 >= 0.5); top <- head(A, 2); top1 <- A1[match(top$symbol, symbol)]
  txt <- paste0("Anchor gene not a tested eGene in GTEx breast (QTD000211).\n\n",
                "Locus-wide bulk coloc, all ", nrow(A), " GTEx eGenes within\n",
                "+/-300 kb of the index SNP (p12 = 5e-5):\n",
                "  best PP4 = ", sprintf("%.2f", A$PP4[1]), " (", A$symbol[1], "); ",
                n5, " gene(s) >= 0.5, none >= 0.8.\n",
                if (n5 > 0) paste0("  Prior-sensitive: ",
                    paste(sprintf("%s %.2f -> %.2f", top$symbol, top$PP4, top1$PP4), collapse = ", "),
                    " at p12 = 1e-5.\n") else "",
                "\nThe eQTL panel shows ", A$symbol[1], ", the locus-wide best.\n",
                "Full table: Supplementary Table (locus-wide coloc).")
  note <- "coloc.abf, GTEx v8 breast bulk, N = 394; priors and harmonisation as the 20-gene run (script 34)."
  text_panel(txt, note, size = 8.5, family = "")
}

pp_panel <- function(lc, rep_gene) {
  D <- co[locus == lc]
  if (!nrow(D)) return(pp_panel_locuswide(lc))
  hdr  <- sprintf("%-9s %5s %5s %5s %5s %5s %6s",
                  "gene", "PP0", "PP1", "PP2", "PP3", "PP4", "nsnps")
  rows <- sprintf("%-9s %5.3f %5.3f %5.3f %5.3f %5.3f %6d",
                  ifelse(D$gene == rep_gene, paste0(D$gene, "*"), D$gene),
                  D$PP0, D$PP1, D$PP2, D$PP3, D$PP4, D$nsnps)
  txt  <- paste(c(sprintf("Colocalisation at %s (%d gene%s tested)", lc, nrow(D),
                          if (nrow(D) == 1) "" else "s"), "",
                  hdr, strrep("-", nchar(hdr)), rows), collapse = "\n")
  note <- sprintf("coloc, GTEx breast bulk, +/-300 kb, p12 = %g\n* gene shown in the eQTL panel above",
                  P12)
  text_panel(txt, note)
}

## LD comes from the LDlink web service, which times out often enough to break a batch run
## halfway through. Retry a few times, then cache the rsid -> r2 map to disk so reruns are
## do not depend on the service. Delete the cache directory to refresh.
LDCACHE <- file.path(OUT, "ld_cache")
dir.create(LDCACHE, showWarnings = FALSE, recursive = TRUE)

attach_ld <- function(loc, lc, tries = 3L) {
  f <- file.path(LDCACHE, sprintf("ld_%s.rds", lc))
  if (file.exists(f)) {
    map <- readRDS(f)
    loc$data$ld <- map[match(loc$data$rsid, names(map))]
    cat("  LD from cache:", sum(!is.na(loc$data$ld)), "of", nrow(loc$data), "\n")
    return(loc)
  }
  for (k in seq_len(tries)) {
    L <- tryCatch(link_LD(loc, token = TOKEN, genome_build = "grch37"),
                  error = function(e) {
                    cat("  link_LD attempt ", k, " failed: ",
                        conditionMessage(e), "\n", sep = ""); NULL })
    if (!is.null(L) && "ld" %in% names(L$data) && sum(!is.na(L$data$ld)) > 0) {
      map <- setNames(L$data$ld, L$data$rsid)
      saveRDS(map[!is.na(map)], f)
      cat("  LD matched:", sum(!is.na(L$data$ld)), "of", nrow(L$data), "-> cached\n")
      return(L)
    }
  }
  ## Fall back to no LD rather than dying: the figure is still readable, every point is grey,
  ## and the log says so. Never silently pretend r2 was retrieved.
  cat("  WARNING: no LD available for", lc, "-- points will be grey (NA)\n")
  loc$data$ld <- NA_real_
  loc
}

## gg = locus gene (window centre); ensg = ENSG of the gene whose eQTL rows are pulled (eensg).
make <- function(lc, gg, ensg) {
  ic  <- LOCUS_CHR(lc)
  d   <- as.data.frame(gw[chrom == ic, .(chrom, pos, rsid, p)])
  loc <- tryCatch(locus(data = d, gene = ens_name(gg), flank = FLANK,
                        ens_db = "EnsDb.Hsapiens.v75", chrom = "chrom", pos = "pos",
                        p = "p", labs = "rsid"),
                  error = function(e) {
                    cat("  locus() failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(loc)) return(NULL)
  loc <- attach_ld(loc, lc)
  tf <- tempfile()
  system(paste("tabix", shQuote(GZ),
               sprintf("%d:%d-%d", ic, max(1, loc$xrange[1] - 2e6), loc$xrange[2] + 2e6),
               ">", shQuote(tf)))
  if (!file.exists(tf) || file.info(tf)$size == 0) return(list(loc = loc, eq = NULL, ic = ic))
  E  <- fread(tf, header = FALSE, col.names = cn)[mtid == ensg & !is.na(pvalue)][!duplicated(rsid)]
  ld <- as.data.table(loc$data)[, .(rsid, pos, ld)]   # attach_ld() guarantees the ld column
  list(loc = loc, eq = merge(ld, E[, .(rsid, ep = pvalue)], by = "rsid"), ic = ic)
}

## eg = gene whose eQTL is shown (defaults to the locus gene gg).
draw <- function(lc, gg, X, eg = gg) {
  loc  <- X$loc; d <- as.data.table(loc$data)
  elab <- if (eg == gg) eg else paste0(eg, " (locus-wide best)")
  idx  <- d$rsid[which.min(d$p)]
  mbs  <- intersect(sig$rsid, d$rsid)
  ## Label the index SNP plus the LABEL_N most GWAS-significant motifbreakR SNPs only.
  ## Labelling all is unreadable -- 19p13.11 has 63 in the window. The rest stay as triangles.
  labs <- unique(c(idx, head(d[rsid %in% mbs][order(p)]$rsid, LABEL_N)))
  ## credible set for the gene whose eQTL is plotted; map hg19 pos -> the panel's rsids
  csm <- if (nrow(CS)) CS[locus == lc & gene == eg] else data.table()
  if (nrow(csm)) {
    csm <- merge(csm, d[, .(rsid, pos19 = pos)], by = "pos19")
    cs_pip <- setNames(csm$pip, csm$rsid)
    for (sg in sort(unique(csm$sig_idx))) {
      z <- csm[sig_idx == sg]
      cat(sprintf("  DAP-G CS95 %s signal %d: %d members on the panel, top PIP %.3f (%s)
",
                  eg, sg, nrow(z), max(z$pip), z$rsid[which.max(z$pip)]))
    }
  } else { cs_pip <- setNames(numeric(0), character(0)) }
  cat("  locus", lc, "| eQTL gene", gg, "| index", idx,
      "| motifbreakR SNPs", length(mbs), "\n")
  p_g <- add_mb(gg_scatter(loc, showLD = TRUE, labels = labs, index_snp = idx,
                           ylab = "TN-GWAS  -log10 P", xticks = FALSE,
                           legend_pos = "topright"), mbs) + ggtitle(lc)

  p_e <- NULL; p_lc <- NULL
  if (is.null(X$eq) || nrow(X$eq) <= 20)
    cat("  NO eQTL PANEL: ", if (is.null(X$eq)) 0L else nrow(X$eq),
        " shared variants -- gene not tested in GTEx breast\n", sep = "")
  if (!is.null(X$eq) && nrow(X$eq) > 20) {
    e <- copy(X$eq); e[, `:=`(chrom = X$ic, p = ep)]
    loce <- tryCatch(locus(data = as.data.frame(e), gene = ens_name(gg), flank = FLANK,
                           ens_db = "EnsDb.Hsapiens.v75", chrom = "chrom", pos = "pos",
                           p = "p", labs = "rsid"), error = function(e) NULL)
    if (!is.null(loce)) {
      loce$data$ld <- e$ld[match(loce$data$rsid, e$rsid)]
      p_e <- add_mb(gg_scatter(loce, showLD = TRUE,
                               labels = intersect(labs, loce$data$rsid),
                               index_snp = idx, ylab = paste0(elab, " eQTL  -log10 P"),
                               xticks = FALSE, legend_pos = NULL), mbs)
      p_e <- add_cs(p_e, names(cs_pip), cs_pip)
    }
    m <- merge(d[, .(rsid, gp = p, ld)], e[, .(rsid, ep)], by = "rsid")
    m[, bin := ld_bin(ld)]
    p_lc <- ggplot(m, aes(-log10(ep), -log10(gp))) +
      geom_point(aes(fill = bin), shape = 21, size = 2, colour = "grey30",
                 stroke = .2, na.rm = TRUE) +
      geom_point(data = m[rsid %in% mbs], shape = 24, size = 2.4, fill = NA,
                 colour = "black", stroke = .7) +
      geom_point(data = m[rsid == idx], shape = 23, size = 4, fill = "purple",
                 colour = "black") +
      ## CS ring drawn last: the top-PIP member is often the index SNP itself, and an
      ## underneath ring would be hidden by the purple diamond.
      geom_point(data = m[rsid %in% names(cs_pip)], shape = 23, size = 5.2, fill = NA,
                 colour = "#8B0000", stroke = 1.3) +
      ggrepel::geom_text_repel(data = m[rsid %in% labs], aes(label = rsid), size = 2.8,
                               box.padding = .7, min.segment.length = 0,
                               max.overlaps = Inf, segment.size = .2) +
      scale_fill_manual(values = LDCOL, na.value = "grey80", name = expression(r^2)) +
      labs(x = paste0(elab, " eQTL  -log10 P"), y = "TN-GWAS  -log10 P",
           title = if (eg == gg) paste0(lc, " — ", gg) else
                   paste0(lc, " — ", gg, "  (eQTL shown: ", eg, ", locus-wide best PP4)"),
           subtitle = sprintf("%d shared variants; %d motifbreakR SNPs (triangles)%s",
                              nrow(m), sum(m$rsid %in% mbs),
                              if (length(cs_pip)) sprintf("
%d DAP-G CS95 members (dark red diamonds), top PIP %.3f",
                                                          sum(m$rsid %in% names(cs_pip)), max(cs_pip)) else "")) +
      theme_bw(base_size = 9) + theme(legend.position = "top")
  }
  ## Left column. Normal locus: scatter + anchor-gene PP table. Anchor gene untestable in
  ## GTEx (5p15.33/TERT): scatter for the locus-wide best gene (if any), the anchor gene's
  ## promoter motifbreakR table, and the locus-wide PP table -- so the eQTL association in the
  ## GWAS region is shown AND the reason the anchor has no coloc is on the figure.
  if (is.null(p_lc)) p_lc <- text_panel(sprintf("%s — %s\n\nno eQTL available for the scatter", lc, gg),
                                        size = 10, family = "")

  gt    <- gg_genetracks(loc, filter_gene_biotype = c("protein_coding", "lincRNA"),
                         cex.text = .7)
  right <- plot_grid(plotlist = Filter(Negate(is.null), list(p_g, p_e, gt)),
                     ncol = 1, align = "v",
                     rel_heights = if (is.null(p_e)) c(1, .95) else c(1, 1, .95))
  ## Left column: scatter on top, PP0-PP4 table for every gene at the locus below. Splitting
  ## the column keeps the scatter roughly square rather than stretching it to the tall canvas.
  ## Left column: locuscompare scatter over the anchor-gene PP table. The promoter motifbreakR
  ## SNP table and the locus-wide coloc table are NOT drawn (they are reported as tables).
  left <- plot_grid(p_lc, pp_panel(lc, gg), ncol = 1, rel_heights = c(1.5, 1))
  plot_grid(left, right, ncol = 2, rel_widths = c(1, 1.2))
}

sel <- commandArgs(TRUE)
R   <- if (length(sel)) reps[locus %in% sel] else reps
for (i in seq_len(nrow(R))) {
  cat("\n==", R$locus[i], "==\n")
  X <- tryCatch(make(R$locus[i], R$gene[i], R$eensg[i]),
                error = function(e) {
                  cat("  make() error:", conditionMessage(e), "\n"); NULL })
  if (is.null(X)) next
  p <- tryCatch(draw(R$locus[i], R$gene[i], X, eg = R$eg[i]),
                error = function(e) {
                  cat("  draw() error:", conditionMessage(e), "\n"); NULL })
  if (is.null(p)) next
  f <- file.path(FIG, sprintf("locuszoom_%s_%s.png", R$locus[i], R$gene[i]))
  ## bg="white": several panels are theme_void/annotation-only and would otherwise leave
  ## transparent regions that composite to black.
  ggsave(f, p, width = 12, height = 11, dpi = 200, bg = "white")
  cat("  wrote", f, "\n")
}
