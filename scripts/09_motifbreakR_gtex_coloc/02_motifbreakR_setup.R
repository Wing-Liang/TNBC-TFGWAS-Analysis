# ==============================================================================
# 02_motifbreakR_setup.R
# Select genome-wide significant TNBC variants for motifbreakR.
# Produces: input to 03_motifbreakR_tnbc.R
# Inputs  : results/tfgwas/variant_peak_motif_subtypes_ALL.csv
#           results/tfgwas/TNBC_uniqueSNP_motifbreak.csv
#           results/singlecell/jaspar2022_motifs_727.csv
#           results/tfgwas/tfgwas_results_TNBC_ALL.csv
#           results/gwas/bcac_gwas_filtered.csv
# Outputs : results/tfgwas/bcac_TNBC_gwsig_snps.csv
#           results/tfgwas/TNBC_uniqueSNP_motifbreak.csv
#           results/motifbreakR/chr1_204497913_rs2169137.jpg
#           results/motifbreakR/chr2_29006516_rs60287247.jpg
#           results/motifbreakR/chr5_1279790_rs10069690.jpg
#           results/motifbreakR/chr6_33239869_rs17215231.jpg
#           results/motifbreakR/chr8_124744912_rs4418323.jpg
#           ...
# Run     : Rscript scripts/09_motifbreakR_gtex_coloc/02_motifbreakR_setup.R   (from the repository root)
# ==============================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
KG_DIR <- Sys.getenv("KG_DIR", "data/1KG")
library(motifbreakR)
setwd(TNBC_ROOT)

# Step 1 | Read in Single Nucleotide Variants
library(BSgenome)
library(data.table)
available.SNPs()

BCAC_motif <- fread("results/tfgwas/variant_peak_motif_subtypes_ALL.csv", header=TRUE)
BCAC_motif_TNBC_unqiue <- BCAC_motif %>%
  select(colnames(BCAC_motif)[1:13])

table(BCAC_motif_TNBC_unqiue$TNBC_pval<5e-8)

table(BCAC_motif_TNBC_unqiue$TNBC_pval<5e-8 &
        BCAC_motif_TNBC_unqiue$LuminalA_pval>=5e-8 &
        BCAC_motif_TNBC_unqiue$LuminalB_P_meta>=5e-8 &
        BCAC_motif_TNBC_unqiue$LuminalB_HER2Neg_pval>=5e-8&
        BCAC_motif_TNBC_unqiue$HER2_pval>=5e-8)

BCAC_motif_TNBC_unqiue <- BCAC_motif_TNBC_unqiue %>%
  filter(TNBC_pval<5e-8 & LuminalA_pval>=5e-8 & LuminalB_P_meta>=5e-8 &LuminalB_HER2Neg_pval>=5e-8& HER2_pval>=5e-8)
dim(BCAC_motif_TNBC_unqiue)
table(BCAC_motif_TNBC_unqiue$Chr)
# chr1 chr11 chr19  chr2  chr5  chr6  chr8

snp_data<-fread("snp_index.txt", header = TRUE)
snp_data$chrpos37 <- paste0('chr',snp_data$CHR,':',snp_data$POS)
names(snp_data)[1] <- 'rsid'
snp_data  <- snp_data  %>% dplyr::select(chrpos37, rsid, REF, ALT)

BCAC_motif_TNBC_unqiue <- left_join(BCAC_motif_TNBC_unqiue, snp_data, by='chrpos37')
table(is.na(BCAC_motif_TNBC_unqiue$rsid))
table(BCAC_motif_TNBC_unqiue$Chr,is.na(BCAC_motif_TNBC_unqiue$rsid))

# loopk for rsid in BSgenome
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(BSgenome)
library(tidyverse)
library(data.table)

### granges
list_missing_snp <- BCAC_motif_TNBC_unqiue %>%
  filter(is.na(rsid)) %>%
  dplyr::mutate(ranges=str_c(gsub("^chr", "", Chr),":",Pos,"-",Pos))
### search snps
snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37
snp.res <- snpsByOverlaps(snps, GRanges(list_missing_snp$ranges))
### combine
list_missing_snp <- as.data.frame(snp.res) %>%
  dplyr::mutate(ranges=str_c(seqnames,":",pos,"-",pos)) %>%
  inner_join(list_missing_snp,by="ranges") %>%
  select(RefSNP_id, chrpos37)
names(list_missing_snp)[1] <- 'rsid'
dim(list_missing_snp)

BCAC_motif_TNBC_unqiue$rsid <- ifelse(BCAC_motif_TNBC_unqiue$chrpos37 %in%list_missing_snp$chrpos37,
                                 list_missing_snp$rsid, BCAC_motif_TNBC_unqiue$rsid)
table(is.na(BCAC_motif_TNBC_unqiue$rsid))
table(BCAC_motif_TNBC_unqiue$Chr,is.na(BCAC_motif_TNBC_unqiue$rsid))
# chr1    160   10
# chr11    57    2
# chr19   176    8
# chr2    251    8
# chr5     29    2
# chr6     23    2
# chr8      1    0

fwrite(BCAC_motif_TNBC_unqiue,'results/tfgwas/bcac_TNBC_gwsig_snps.csv', row.names = FALSE)

library(TwoSampleMR)

ld_clumped<-ld_clump_local(
  dplyr::tibble(rsid = BCAC_motif_TNBC_unqiue$rsid,
                pval = BCAC_motif_TNBC_unqiue$TNBC_pval,
                id = BCAC_motif_TNBC_unqiue$chrpos37),
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = file.path(KG_DIR, "1kg.v3/EUR")
)
dim(ld_clumped)
ld_clumped

BCAC_motif_TNBC_unqiue_clumped<-BCAC_motif_TNBC_unqiue %>% filter(rsid %in% ld_clumped$rsid)
fwrite(BCAC_motif_TNBC_unqiue_clumped, "results/tfgwas/TNBC_uniqueSNP_motifbreak.csv", row.names = FALSE)
# Chr       Pos    TNBC_pval TNBC_chisquare LuminalA_pval LuminalA_chisquare Luminal_B_P_meta
# <char>     <int>        <num>          <num>         <num>              <num>            <num>
#   1:  chr19  17401404 2.759964e-49      217.77876  5.854089e-02           3.578247     0.0175049722
# 2:   chr5   1279790 8.645201e-47      206.33831  6.141506e-04          11.732590     0.1966090751
# 3:   chr2 121088182 1.060199e-12       50.72939  1.191449e-03          10.503558     0.0019976929
# 4:   chr6 151970639 1.809322e-14       58.72924  7.144975e-05          15.772064     0.0004290088
# 5:   chr6  33239869 2.365168e-10       40.13902  6.070439e-01           0.264502     0.8419807466
# 6:   chr1 204497913 4.813557e-22       93.16401  1.708859e-01           1.875155     0.3714040858
# 7:   chr8 124744912 3.143479e-08       30.61687  2.727902e-04          13.248554     0.0246684353
# 8:  chr11 108311655 1.738098e-11       45.24561  1.055832e-01           2.619117     0.0529068944
# 9:   chr2  29006516 3.412329e-11       43.92483  5.574948e-02           3.659526     0.9597552603
# LuminalB_chisquare    HER2_pval HER2_chisquare        chrpos37       rsid    REF    ALT
# <num>        <num>          <num>          <char>     <char> <char> <char>
#   1:        5.645024441 7.509715e-01    0.100716497  chr19:17401404 rs67397200      C      G
# 2:        1.667385622 4.257841e-01    0.634295267    chr5:1279790 rs10069690      C      G
# 3:        9.551654334 8.109523e-02    3.042814010  chr2:121088182 rs11903787      G      G
# 4:       12.401414328 6.983877e-06   20.197647570  chr6:151970639  rs1293960      G      G
# 5:        0.039744614 5.668735e-01    0.327941857   chr6:33239869 rs17215231      C      G
# 6:        0.798961365 1.523152e-02    5.889489408  chr1:204497913  rs2169137      C      G
# 7:        5.047012568 5.818056e-01    0.303324780  chr8:124744912  rs4418323      G      G
# 8:        3.746857412 9.331534e-01    0.007035523 chr11:108311655  rs4754322      T      G
# 9:        0.002546283 7.870721e-01    0.072961892   chr2:29006516 rs60287247      G      G

library(data.table)
library(dplyr)
BCAC_motif_TNBC_unqiue_clumped<-fread("results/tfgwas/TNBC_uniqueSNP_motifbreak.csv", header = TRUE)

# 3.1.1 SNPs from rsID:
library(SNPlocs.Hsapiens.dbSNP155.GRCh37) # dbSNP155 in hg19
library(BSgenome.Hsapiens.UCSC.hg19)

snps.mb <- snps.from.rsid(rsid = BCAC_motif_TNBC_unqiue_clumped$rsid,
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh37,
                          search.genome = BSgenome.Hsapiens.UCSC.hg19)
snps.mb
unique(snps.mb$SNP_id)

# Step 2 | Find Broken Motifs
library(MotifDb)
# MotifDb
table(mcols(MotifDb)$organism, mcols(MotifDb)$dataSource=='jaspar2022')

data(motifbreakR_motif)
motifbreakR_motif

subset<- subset(MotifDb,
                dataSource %in% c("jaspar2022") & mcols(MotifDb)$organism=='Hsapiens')
length(subset@elementMetadata$providerName)

results <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                       pwmList = subset,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::SerialParam())
head(results)
unique(results$SNP_id)
BCAC_motif_TNBC_unqiue_clumped[!rsid %in% unique(results$SNP_id)]

results$providerName
results$geneSymbol

motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
TNBC_sig<- left_join(motif_id_list_final, TNBC, by=c('motif','motif.name')) %>% dplyr::filter(TF_p_value <= 0.05/727)

TNBC_motifscan_sig <- TNBC_sig[motif %in% results$providerId,]
TNBC_motifscan_sig$motif.name
# 6 TF in significant TNBC motifs

TNBC_motifscan_uniq <- TNBC_unique[motif %in% results$providerId,]
TNBC_motifscan_uniq$motif.name
# "TBX5"  "ZNF85"

# check the alt, ref
BCAC_filt <- fread('results/gwas/bcac_gwas_filtered.csv', header = TRUE)
SNP_orig_info <- BCAC_filt %>%
  filter(paste0('chr',chrpos37)%in% BCAC_motif_TNBC_unqiue_clumped$chrpos37) %>%
  select(chrpos37, Effect.Meta, Baseline.Meta,Triple_Neg_log_or_meta, luminal_A_log_or_meta, Luminal_B_log_or_meta, HER2_Enriched_log_or_meta)
SNP_orig_info$chrpos37 <- paste0('chr',SNP_orig_info$chrpos37)
SNP_orig_info <- left_join(SNP_orig_info, select(BCAC_motif_TNBC_unqiue_clumped,Chr, rsid, chrpos37),by='chrpos37')
SNP_orig_info[order(Chr)]
# chrpos37 Effect.Meta Baseline.Meta    Chr       rsid
# <char>      <char>        <char> <char>     <char>
#   1:  chr1:204497913           C             G   chr1  rs2169137
# 2: chr11:108311655           C             T  chr11  rs4754322
# 3:  chr19:17401404           G             C  chr19 rs67397200
# 4:  chr2:121088182           A             G   chr2 rs11903787
# 5:   chr2:29006516           A             G   chr2 rs60287247
# 6:    chr5:1279790           T             C   chr5 rs10069690
# 7:  chr6:151970639           G             C   chr6  rs1293960
# 8:   chr6:33239869           T             C   chr6 rs17215231
# 9:  chr8:124744912           T             G   chr8  rs4418323

rs2169137 <- results[results$SNP_id == "rs2169137"]
rs2169137 <- calculatePvalue(rs2169137, granularity = 1e-6)
rs2169137
jpeg('results/motifbreakR/chr1_204497913_rs2169137.jpg',width=7, height = 4, units = 'in', res = 150)
plotMB(results = results, rsid = "rs2169137", effect = "strong", altAllele = "C")
dev.off()

rs11903787 <- results[results$SNP_id == "rs11903787"]
rs11903787 <- calculatePvalue(rs11903787, granularity = 1e-6)
rs11903787

rs60287247 <- results[results$SNP_id == "rs60287247"]
rs60287247 <- calculatePvalue(rs60287247, granularity = 1e-6)
rs60287247
jpeg('results/motifbreakR/chr2_29006516_rs60287247.jpg',width=7, height = 3, units = 'in', res = 150)
plotMB(results = results, rsid = "rs60287247", effect = "strong", altAllele = "A")
dev.off()

rs10069690 <- results[results$SNP_id == "rs10069690"]
rs10069690 <- calculatePvalue(rs10069690, granularity = 1e-6)
rs10069690
jpeg('results/motifbreakR/chr5_1279790_rs10069690.jpg',width=7, height = 3, units = 'in', res = 150)
plotMB(results = results, rsid = "rs10069690", effect = "strong", altAllele = "T")
dev.off()

rs17215231 <- results[results$SNP_id == "rs17215231"]
rs17215231 <- calculatePvalue(rs17215231, granularity = 1e-6)
rs17215231
jpeg('results/motifbreakR/chr6_33239869_rs17215231.jpg',width=7, height = 4, units = 'in', res = 150)
plotMB(results = results, rsid = "rs17215231", effect = "strong", altAllele = "T")
dev.off()

rs4418323 <- results[results$SNP_id == "rs4418323"]
rs4418323 <- calculatePvalue(rs4418323, granularity = 1e-6)
rs4418323
jpeg('results/motifbreakR/chr8_124744912_rs4418323.jpg',width=7, height = 4, units = 'in', res = 150)
plotMB(results = results, rsid = "rs4418323", effect = "strong", altAllele = "T")
dev.off()

rs4754322 <- results[results$SNP_id == "rs4754322"]
rs4754322 <- calculatePvalue(rs4754322, granularity = 1e-6)
rs4754322
jpeg('results/motifbreakR/chr11_108311655_rs4754322.jpg',width=8, height = 5, units = 'in', res = 150)
plotMB(results = results, rsid = "rs4754322", effect = "strong", altAllele = "C")
dev.off()

rs67397200 <- results[results$SNP_id == "rs67397200"]
rs67397200 <- calculatePvalue(rs67397200, granularity = 1e-6)
rs67397200
jpeg('results/motifbreakR/chr19_17401404_rs67397200.jpg',width=7, height = 4, units = 'in', res = 150)
plotMB(results = results, rsid = "rs67397200", effect = "strong", altAllele = "G")
dev.off()

# connect to 02_motif_enrichment_by_subtype.R
# for the motif score heatmap
motif_df_motifbreakR <- motif_df[row.names(motif_df) %in% TNBC_motifscan_sig$motif.name,]
motif_matrix <- as.matrix(motif_df_motifbreakR)

library(ComplexHeatmap)
library(dendextend)
set.seed(2025)
heatmap_obj_orig <-Heatmap(motif_matrix,
                           name = "-log10 pval (chromVAR)",          # Title of legend
                           row_names_gp = gpar(fontsize = 7),
                           row_names_side = "left",
                           column_names_centered = TRUE,
                           width = unit(50, "mm"),
                           col = circlize::colorRamp2(c(0,35), c("#483d8b","red1")),
                           km =4,                        # K-means clustering with 10 clusters
                           rect_gp = gpar(col = "grey99"), # Set cell border color (instead of border_color)
                           show_row_dend = FALSE,         # Hide row dendrogram to avoid km labels
                           show_column_dend = FALSE,
                           column_split = 4,                 # Ensure 3 clusters without labels
                           row_title = NULL ,
                           column_title = NULL
)
heatmap_obj_orig
pdf("results/motifbreakR/heatmap_enrichedmotif_TNBC_sig.pdf", width = 4, height =6)
draw(heatmap_obj_orig)
dev.off()

# diff motif JASPAR2022 CORE VS ALL JASPAR2022
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)

library(MotifDb)
subset<- subset(MotifDb,
                dataSource %in% c("jaspar2022") & mcols(MotifDb)$organism=='Hsapiens')
length(subset@elementMetadata$providerName)
jaspar2022_core<-subset@elementMetadata$providerName

motif_not_in_motifbreakR<-motif_id_list_final[!motif %in% jaspar2022_core, ]

fwrite(motif_not_in_motifbreakR,'results/singlecell/Motif_not_in_motifbreakR.csv', row.names=FALSE)
