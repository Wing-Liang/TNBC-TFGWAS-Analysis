## 01_prepare_gwas_inputs.R
## Per-subtype GWAS tables for MAGMA / scDRS (Fig S10) from the filtered BCAC summary statistics
## (GRCh37), with rsIDs added from the 1000 Genomes EUR reference VCF distributed with scPagwas.
## Each table has columns chrom, pos, p, se, maf, beta, chrposID, rsid; beta is oriented to the
## minor allele.
##
## Inputs   results/gwas/bcac_gwas_filtered.csv              (02_gwas_ld_blocks/01_filter_gwas.R)
##          $KG_DIR/1000G.EUR.QC.plink.simple.vcf             (scPagwas reference data)
## Outputs  results/gwas/bcac_gwas_filtered_<subtype>.csv
## Run from the repository root:  Rscript scripts/07_scdrs/01_prepare_gwas_inputs.R

TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
KG_DIR    <- Sys.getenv("KG_DIR", "data/1KG")

library(data.table)
library(dplyr)
library(VariantAnnotation)

setwd(TNBC_ROOT)

## ---- rsID index from the 1000 Genomes EUR VCF (GRCh37) ----
## The VCF's #CHROM line carries a FORMAT column without sample columns; keep the first 8 fields.
lines <- readLines(file.path(KG_DIR, "1000G.EUR.QC.plink.simple.vcf"))
lines[grep("^#CHROM", lines)] <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
body <- grep("^[^#]", lines)
lines[body] <- vapply(strsplit(lines[body], "\t"), function(x) paste(x[1:8], collapse = "\t"), "")
vcf_fixed <- tempfile(fileext = ".vcf"); writeLines(lines, vcf_fixed)
vcf  <- readVcf(vcf_fixed, "hg19")
snps <- rowRanges(vcf)
snp_data <- data.frame(rsid = names(snps), chrposID = paste0("chr", seqnames(snps), ":", start(snps)))

## ---- one table per subtype ----
BCAC_filt <- fread("results/gwas/bcac_gwas_filtered.csv", header = TRUE)
BCAC_COLS <- list(                      # subtype = c(p, se, log OR)
  TNBC         = c("Triple_Neg_P_meta",        "Triple_Neg_sd_meta",        "Triple_Neg_log_or_meta"),
  LumA         = c("Luminal_A_P_meta",         "Luminal_A_sd_meta",         "luminal_A_log_or_meta"),
  LumB         = c("Luminal_B_P_meta",         "Luminal_B_sd_meta",         "Luminal_B_log_or_meta"),
  LumB_HER2Neg = c("Luminal_B_HER2Neg_P_meta", "Luminal_B_HER2Neg_sd_meta", "Luminal_B_HER2Neg_log_or_meta"),
  HER2         = c("HER2_P_meta",              "HER2_Enriched_sd_meta",     "HER2_Enriched_log_or_meta")
)
for (s in names(BCAC_COLS)) {
  cc  <- BCAC_COLS[[s]]
  eaf <- BCAC_filt$EAFcontrols.Onco
  d <- data.frame(chrom = BCAC_filt$chr.Onco, pos = BCAC_filt$Position.Onco,
                  p = BCAC_filt[[cc[1]]], se = BCAC_filt[[cc[2]]],
                  maf  = ifelse(eaf < 0.5, eaf, 1 - eaf),
                  beta = ifelse(eaf < 0.5, BCAC_filt[[cc[3]]], -BCAC_filt[[cc[3]]]),
                  chrposID = paste0("chr", BCAC_filt$chr.Onco, ":", BCAC_filt$Position.Onco))
  d <- left_join(d, snp_data, by = "chrposID")
  fwrite(d, sprintf("results/gwas/bcac_gwas_filtered_%s.csv", s))
}

