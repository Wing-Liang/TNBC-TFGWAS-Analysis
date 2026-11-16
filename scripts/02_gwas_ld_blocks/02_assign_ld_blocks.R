## 02_assign_ld_blocks.R
## Assign every BCAC GWAS variant to a European LD block (Table S2).
##
## GWAS positions are GRCh37. They are lifted to GRCh38 with the UCSC chain file and overlapped
## with the 1,336 autosomal pyrho EUR LD blocks (GRCh38; MacDonald et al., LDblocks_GRCh38).
## Variants that fail the liftover, or fall between blocks, take the block of the nearest
## assigned variant (distance measured on GRCh37).
##
## Inputs   results/gwas/bcac_gwas_filtered.csv          (from 03_tf_gwas/01_build_variant_motif_matrix.R)
##          $LDBLOCK_DIR/pyrho_EUR_LD_blocks.bed
##          $LDBLOCK_DIR/hg19ToHg38.over.chain
## Outputs  results/ld_blocks/ld_blocks_hg38_chr<1-22>.csv
##          results/ld_blocks/ld_blocks_hg38_ALL.csv       (Table S2)
## Run from the repository root:  Rscript scripts/02_gwas_ld_blocks/02_assign_ld_blocks.R

TNBC_ROOT   <- Sys.getenv("TNBC_ROOT", ".")
LDBLOCK_DIR <- Sys.getenv("LDBLOCK_DIR", "data/LDblocks")

library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(dplyr)

setwd(TNBC_ROOT)
OUT <- "results/ld_blocks"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

BCAC_filt <- fread("results/gwas/bcac_gwas_filtered.csv", header = TRUE)

## ---- LD blocks, numbered within chromosome ----
pyrho_EUR <- fread(file.path(LDBLOCK_DIR, "pyrho_EUR_LD_blocks.bed"), header = TRUE)
pyrho_EUR$range <- pyrho_EUR$end - pyrho_EUR$start
pyrho_EUR$block <- NA
for (chrom in unique(pyrho_EUR$chr)) {
  idx <- which(pyrho_EUR$chr == chrom)
  pyrho_EUR$block[idx] <- paste0(chrom, "_", seq_along(idx))
}
pyrho_gr <- GRanges(
  seqnames = pyrho_EUR$chr,
  ranges   = IRanges(start = pyrho_EUR$start, end = pyrho_EUR$end),
  block    = pyrho_EUR$block
)

chain <- import.chain(file.path(LDBLOCK_DIR, "hg19ToHg38.over.chain"))

## ---- per chromosome: liftover, overlap, nearest-block fallback ----
for (i in 1:22) {
  chr <- paste0("chr", i)

  BCAC_filt_chr <- BCAC_filt[chr.iCOGs == as.character(i)]
  if (nrow(BCAC_filt_chr) == 0) {
    cat(sprintf("No data for chromosome %d, skipping\n", i))
    next
  }
  BCAC_filt_chr <- BCAC_filt_chr %>%
    dplyr::select(chr.iCOGs, Position.iCOGs) %>%
    setnames(c("Chr", "Pos37"))
  BCAC_filt_chr$Chr <- paste0("chr", BCAC_filt_chr$Chr)
  BCAC_filt_chr$chr_pos37 <- paste(BCAC_filt_chr$Chr, BCAC_filt_chr$Pos37, sep = ":")
  gwas <- unique(BCAC_filt_chr) %>% arrange(Pos37)

  gwas_gr <- GRanges(
    seqnames  = gwas$Chr,
    ranges    = IRanges(start = gwas$Pos37, end = gwas$Pos37),
    chr_pos37 = gwas$chr_pos37,
    Pos37     = gwas$Pos37
  )

  # GRCh37 -> GRCh38
  gwas_lifted_df <- as.data.frame(unlist(liftOver(gwas_gr, chain)))
  gwas_hg38 <- data.frame(
    chr          = gwas_lifted_df$seqnames,
    pos_hg38     = gwas_lifted_df$start,
    Pos37        = gwas_lifted_df$Pos37,
    chr_pos37    = gwas_lifted_df$chr_pos37,
    chr_pos_hg38 = paste(gwas_lifted_df$seqnames, gwas_lifted_df$start, sep = ":")
  )
  gwas_hg38 <- gwas_hg38[gwas_hg38$chr == chr, ]
  cat(sprintf("Chromosome %d: %d SNPs unmapped during liftOver\n", i, nrow(gwas_lifted_df) - nrow(gwas_hg38)))
  if (nrow(gwas_hg38) == 0) {
    cat(sprintf("No lifted data for chromosome %d, skipping\n", i))
    next
  }

  # overlap with LD blocks
  gwas_gr_hg38 <- GRanges(
    seqnames     = gwas_hg38$chr,
    ranges       = IRanges(start = gwas_hg38$pos_hg38, end = gwas_hg38$pos_hg38),
    chr_pos_hg38 = gwas_hg38$chr_pos_hg38,
    Pos37        = gwas_hg38$Pos37,
    chr_pos37    = gwas_hg38$chr_pos37
  )
  overlaps <- suppressWarnings(findOverlaps(gwas_gr_hg38, pyrho_gr))
  gwas_hg38$block <- NA
  if (length(overlaps) > 0) {
    gwas_hg38$block[queryHits(overlaps)] <- pyrho_gr$block[subjectHits(overlaps)]
  }
  gwas_hg38 <- gwas_hg38 %>% filter(!is.na(block))
  cat(sprintf("Chromosome %d: %d SNPs without block overlap\n", i, nrow(gwas) - nrow(gwas_hg38)))

  gwas_block_joint <- left_join(gwas, gwas_hg38, by = "chr_pos37") %>% unique()

  # unassigned variants take the block of the nearest assigned variant (GRCh37 positions)
  no_overlap_idx <- which(is.na(gwas_block_joint$block))
  if (length(no_overlap_idx) > 0) {
    gwas_hg19 <- GRanges(
      seqnames     = gwas_block_joint$Chr,
      ranges       = IRanges(start = gwas_block_joint$Pos37.x, end = gwas_block_joint$Pos37.x),
      chr_pos_hg38 = gwas_block_joint$chr_pos_hg38,
      block        = gwas_block_joint$block
    )
    gwas_gr_37block <- GRanges(
      seqnames     = gwas_hg38$chr,
      ranges       = IRanges(start = gwas_hg38$Pos37, end = gwas_hg38$Pos37),
      chr_pos_hg38 = gwas_hg38$chr_pos_hg38,
      block        = gwas_hg38$block
    )
    nearest_idx <- nearest(gwas_hg19[no_overlap_idx], gwas_gr_37block)
    if (length(nearest_idx) > 0 && all(nearest_idx <= length(gwas_gr_37block))) {
      gwas_block_joint$block[no_overlap_idx] <- gwas_gr_37block$block[nearest_idx]
    }
  }
  gwas_block_joint <- gwas_block_joint %>% select(Chr, Pos37.x, chr_pos37, pos_hg38, block)

  write.csv(gwas_block_joint, file.path(OUT, sprintf("ld_blocks_hg38_%s.csv", chr)), row.names = FALSE)
  cat(sprintf("Completed chromosome %d\n", i))
}

## ---- combine chromosomes (Table S2) ----
ld_list <- list()
for (i in 1:22) {
  f <- file.path(OUT, sprintf("ld_blocks_hg38_chr%d.csv", i))
  if (file.exists(f)) ld_list[[i]] <- fread(f, header = TRUE)
}
LDblock <- rbindlist(ld_list)
write.csv(LDblock, file.path(OUT, "ld_blocks_hg38_ALL.csv"), row.names = FALSE)
cat("Variants assigned:", nrow(LDblock), "| blocks:", uniqueN(LDblock$block), "\n")
