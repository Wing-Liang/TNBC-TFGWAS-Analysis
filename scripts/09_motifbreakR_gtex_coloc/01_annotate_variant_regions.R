## 01_annotate_variant_regions.R
## Genomic-region annotation (ChIPseeker, GRCh37 UCSC knownGene, promoter = TSS +/- 3 kb) of every
## GWAS variant in the TF-GWAS matrix; used to select promoter variants for motifbreakR.
##
## Inputs   results/tfgwas/variant_peak_motif_subtypes_ALL.csv   (03_tf_gwas/01_build_variant_motif_matrix.R)
## Outputs  results/tfgwas/BCAC_SNP_annotation.csv
## Run from the repository root:  Rscript scripts/09_motifbreakR_gtex_coloc/01_annotate_variant_regions.R

TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")

library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(data.table)
library(dplyr)

setwd(TNBC_ROOT)

BCAC_all <- fread("results/tfgwas/variant_peak_motif_subtypes_ALL.csv", header = TRUE)
BCAC_all_gr <- GRanges(
  seqnames = BCAC_all$Chr,
  ranges   = IRanges(start = BCAC_all$Pos, end = BCAC_all$Pos),
  chrpos37 = BCAC_all$chrpos37,
  Pos37    = BCAC_all$Pos
)

BCAC_all_Anno <- annotatePeak(BCAC_all_gr, tssRegion = c(-3000, 3000),
                              TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db")
BCAC_all_Anno_df <- as.data.frame(BCAC_all_Anno) %>%
  mutate(annotation_category = case_when(
    grepl("Promoter", annotation) ~ "Promoter",
    grepl("5' UTR", annotation) ~ "5' UTR",
    grepl("3' UTR", annotation) ~ "3' UTR",
    grepl("exon", annotation, ignore.case = TRUE) ~ "Exon",
    grepl("intron", annotation, ignore.case = TRUE) ~ "Intron",
    grepl("Downstream \\(<=300bp\\)", annotation) ~ "Downstream (<=300bp)",
    grepl("Distal Intergenic", annotation) ~ "Distal Intergenic",
    TRUE ~ "Other"
  ))
table(BCAC_all_Anno_df$annotation_category)

fwrite(BCAC_all_Anno_df, "results/tfgwas/BCAC_SNP_annotation.csv", row.names = FALSE)
