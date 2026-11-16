## 01_filter_gwas.R
## Quality-filter the BCAC intrinsic-subtype GWAS summary statistics (Zhang et al. 2020) and
## compute a per-variant chi-square and p-value for each subtype.
##
## Filter: control allele frequency in the OncoArray component between 0.01 and 0.99
## (imputation r2 >= 0.3 is already applied in the released file). 9,819,442 variants remain.
## Chi-square = (log OR / SE)^2 from the iCOGS + OncoArray meta-analysis; p is two-sided.
##
## Inputs   $BCAC_DIR/icogs_onco_meta_intrinsic_subtypes_summary_stats_zhang_2020.zip
## Outputs  results/gwas/bcac_gwas_filtered.csv
## Run from the repository root:  Rscript scripts/02_gwas_ld_blocks/01_filter_gwas.R

TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
BCAC_DIR  <- Sys.getenv("BCAC_DIR", "data/BCAC")

library(data.table)

setwd(TNBC_ROOT)
dir.create("results/gwas", recursive = TRUE, showWarnings = FALSE)

BCAC <- fread(file.path(BCAC_DIR, "icogs_onco_meta_intrinsic_subtypes_summary_stats_zhang_2020.zip"))

## minor-allele filter on the OncoArray control frequency
BCAC_filt <- BCAC[BCAC$EAFcontrols.Onco >= 0.01 & BCAC$EAFcontrols.Onco <= 0.99, ]
nrow(BCAC_filt)

## GRCh37 chr:pos key
BCAC_filt$chrpos37 <- paste0(BCAC_filt$chr.Onco, ":", BCAC_filt$Position.Onco)

## chi-square and two-sided p per subtype (meta-analysis log OR and SE)
SUBTYPES <- list(                      # output prefix = c(log OR column, SE column)
  Triple_Neg        = c("Triple_Neg_log_or_meta",        "Triple_Neg_sd_meta"),
  Luminal_A         = c("luminal_A_log_or_meta",         "Luminal_A_sd_meta"),
  Luminal_B         = c("Luminal_B_log_or_meta",         "Luminal_B_sd_meta"),
  Luminal_B_HER2Neg = c("Luminal_B_HER2Neg_log_or_meta", "Luminal_B_HER2Neg_sd_meta"),
  HER2              = c("HER2_Enriched_log_or_meta",     "HER2_Enriched_sd_meta")
)
for (s in names(SUBTYPES)) {
  z <- BCAC_filt[[SUBTYPES[[s]][1]]] / BCAC_filt[[SUBTYPES[[s]][2]]]
  BCAC_filt[[paste0(s, "_X2_meta")]] <- z^2
  BCAC_filt[[paste0(s, "_P_meta")]]  <- 2 * pnorm(-abs(z))
}

fwrite(BCAC_filt, "results/gwas/bcac_gwas_filtered.csv", row.names = FALSE)
