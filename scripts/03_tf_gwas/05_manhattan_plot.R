# ==============================================================================
# 05_manhattan_plot.R
# TNBC Manhattan plot highlighting variants in peaks with a TNBC-specific motif.
# Produces: Fig 3C
# Inputs  : results/tfgwas/variant_peak_motif_subtypes_ALL.csv
#           results/tfgwas/tfgwas_significant_TNBC_unique.csv
# Run     : Rscript scripts/03_tf_gwas/05_manhattan_plot.R   (from the repository root)
# ==============================================================================
## 05_manhattan_plot.R — redraw the TNBC Manhattan plot with the highlighted SNPs in the TNBC colour
## of the Fig 3A UpSet plot (#E64B35) instead of orange.
## Data prep is identical to 03_subtype_specific_motifs.R (the TNBC block, ~l.3410-3486); only the
## highlight colour changes. Highlighted = SNPs in an ATAC peak carrying >= 1 TNBC-specific motif (any MA* == 1).
## Uses the patched manhattan_functions.R (manhattan_local gains a `highlight.col` argument; default
## "orange", so it is backward-compatible). Base graphics -- no qqman needed.
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressMessages({library(data.table); library(dplyr)})
setwd(TNBC_ROOT)
source("manhattan_functions.R")
HI_COL <- "#E64B35"

BCAC_all    <- fread("results/tfgwas/variant_peak_motif_subtypes_ALL.csv", header = TRUE)
TNBC_unique <- fread("results/tfgwas/tfgwas_significant_TNBC_unique.csv", header = TRUE)
BCAC_TNBC_unique <- BCAC_all %>% select(Chr, Pos, chrpos37, TNBC_pval, TNBC_unique$motif) %>% unique()
BCAC_TNBC_unique$TNBC_pval <- ifelse(BCAC_TNBC_unique$TNBC_pval == 0, 1e-300,
                              ifelse(BCAC_TNBC_unique$TNBC_pval == 1, 0.99999999999999999999999, BCAC_TNBC_unique$TNBC_pval))
snp_data <- fread("snp_index.txt", header = TRUE)
snp_data$chrpos37 <- paste0("chr", snp_data$CHR, ":", snp_data$POS)
names(snp_data)[1] <- "rsid"
snp_data <- snp_data %>% dplyr::select(chrpos37, rsid, REF, ALT)
BCAC_TNBC_unique <- left_join(BCAC_TNBC_unique, snp_data, by = "chrpos37"); setDT(BCAC_TNBC_unique)
BCAC_TNBC_unique$Chr <- as.numeric(gsub("chr", "", BCAC_TNBC_unique$Chr))
BCAC_TNBC_unique <- BCAC_TNBC_unique[!is.na(rsid) & !is.na(TNBC_pval), ]
ma_cols  <- grep("^MA[0-9]+\\.[0-9]+", names(BCAC_TNBC_unique), value = TRUE)
has_ma_1 <- BCAC_TNBC_unique[, rowSums(.SD == 1) > 0, .SDcols = ma_cols]
TNBC_sig_snp <- BCAC_TNBC_unique[has_ma_1]
cat("SNPs:", nrow(BCAC_TNBC_unique), "| highlighted (>= 1 TNBC-specific motif):", nrow(TNBC_sig_snp), "\n")

f <- "results/tfgwas/manhattan_TNBC_TNBCcol.png"
png(f, width = 14, height = 8, units = "in", res = 150)
par(mar = c(5, 5, 1, 1))
manhattan_local(BCAC_TNBC_unique, chr = "Chr", bp = "Pos", p = "TNBC_pval", snp = "rsid",
                col = c("black", "grey"), chrlabs = NULL, suggestiveline = FALSE,
                genomewideline = -log10(5e-08), highlight = TNBC_sig_snp$rsid, highlight.col = HI_COL,
                logp = TRUE, annotateTop = FALSE)
dev.off()
cat("wrote", f, "| highlight colour", HI_COL, "\n")
