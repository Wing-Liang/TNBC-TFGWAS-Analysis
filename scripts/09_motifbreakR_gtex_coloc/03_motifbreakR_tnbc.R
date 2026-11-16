# ==============================================================================
# 03_motifbreakR_tnbc.R
# Predict allele-specific TF-binding changes (motifbreakR, JASPAR 2022) for TNBC variants.
# Produces: Table S12, Fig 5A
# Inputs  : results/tfgwas/bcac_TNBC_gwsig_snps.csv
#           results/gwas/bcac_gwas_filtered.csv
#           results/motifbreakR/ALLresults_TNBCsig.csv
#           results/singlecell/jaspar2022_motifs_727.csv
#           results/tfgwas/tfgwas_results_TNBC_ALL.csv
#           results/tfgwas/BCAC_SNP_annotation.csv
#           results/tfgwas/tfgwas_significant_TNBC_unique.csv
#           ...
# Outputs : results/motifbreakR/motifbreakR_SNPinfo_TNBC.csv
#           results/motifbreakR/ALLresults_TNBCsig.csv
#           results/motifbreakR/MotifbreakR_TNBC_significant.csv
#           results/motifbreakR/ALLresults_TNBCsig_annotation.csv
#           results/motifbreakR/heatmap_enrichedmotif_TNBC_all_5e-4.pdf
#           results/motifbreakR/motifbreakR_manhattan_TNBC.png
#           results/motifbreakR/motifbreakR_TNBC_all.csv
#           ...
# Run     : Rscript scripts/09_motifbreakR_gtex_coloc/03_motifbreakR_tnbc.R   (from the repository root)
# ==============================================================================
# TNBC Unique SNPs without clumping
#                        "Biostrings", "GenomeInfoDb", "GenomicRanges", "Gviz", "S4Vectors",
#                        "rtracklayer", "IRanges", "MotifDb", "BSgenome.Hsapiens.UCSC.hg19",
#                        "SNPlocs.Hsapiens.dbSNP.20120608", "SNPlocs.Hsapiens.dbSNP155.GRCh37",
#                        "VariantAnnotation", "matrixStats", "BiocStyle"))
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
DEPMAP_DIR <- Sys.getenv("DEPMAP_DIR", "data/DepMap")
TCGA_DIR <- Sys.getenv("TCGA_DIR", "data/TCGA")
library(data.table)
library(dplyr)
library(ggplot2)
setwd(TNBC_ROOT)

BCAC_motif_TNBC_unqiue<-fread('results/tfgwas/bcac_TNBC_gwsig_snps.csv', header=TRUE)
BCAC_motif_TNBC_unqiue <- BCAC_motif_TNBC_unqiue[!is.na(rsid),] %>%unique()
BCAC_motif_TNBC_unqiue <- BCAC_motif_TNBC_unqiue %>%
  filter(LuminalA_pval >=5e-4 & LuminalB_P_meta>=5e-4 & LuminalB_HER2Neg_pval >=5e-4&HER2_pval>=5e-4)

# GWAS snp
BCAC_filt <- fread('results/gwas/bcac_gwas_filtered.csv', header = TRUE)
BCAC_filt_select<-BCAC_filt%>%select(Effect.iCOGs,Baseline.iCOGs, chrpos37)
BCAC_filt_select$chrpos37<-paste0('chr',BCAC_filt_select$chrpos37)

# COMBINE
BCAC_motif_TNBC_unqiue<-left_join(BCAC_motif_TNBC_unqiue, BCAC_filt_select, by='chrpos37')
BCAC_motif_TNBC_unqiue$REF<-NULL
BCAC_motif_TNBC_unqiue$ALT<-NULL
names(BCAC_motif_TNBC_unqiue)[15:16]<-c("EA","OA")
BCAC_motif_TNBC_unqiue<-BCAC_motif_TNBC_unqiue[rsid!="",]

length(unique(BCAC_motif_TNBC_unqiue$chrpos37))
unique_SNP<-BCAC_motif_TNBC_unqiue%>%select(Chr, Pos, chrpos37) %>%unique()
table(unique_SNP$Chr)
# chr1 chr11 chr19  chr2  chr5  chr6

# 3.1.1 SNPs from rsID:
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37) # dbSNP155 in hg19
library(BSgenome.Hsapiens.UCSC.hg19)
rsid <- BCAC_motif_TNBC_unqiue$rsid
rsid <- rsid[grepl("^rs", rsid)]

BCAC_filt <- fread('results/gwas/bcac_gwas_filtered.csv', header = TRUE)
SNP_orig_info <- BCAC_filt %>%
  filter(paste0('chr',chrpos37)%in% BCAC_motif_TNBC_unqiue$chrpos37) %>%
  select(chrpos37, chr.Onco, Effect.Meta, Baseline.Meta,Triple_Neg_log_or_meta,Triple_Neg_P_meta,
         luminal_A_log_or_meta, Luminal_A_P_meta,
         Luminal_B_log_or_meta, Luminal_B_P_meta,
         HER2_Enriched_log_or_meta, HER2_P_meta)%>%
  unique()
SNP_orig_info$chrpos37 <- paste0('chr',SNP_orig_info$chrpos37)
SNP_orig_info <- left_join(SNP_orig_info, select(BCAC_motif_TNBC_unqiue,rsid, chrpos37)%>%unique(),by='chrpos37')
SNP_orig_info<-SNP_orig_info[order(chr.Onco)]
head(SNP_orig_info)
fwrite(SNP_orig_info, 'results/motifbreakR/motifbreakR_SNPinfo_TNBC.csv', row.names = FALSE)

# Find Broken Motifs
library(MotifDb)
subset<- subset(MotifDb,
                dataSource %in% c("jaspar2022") & mcols(MotifDb)$organism=='Hsapiens')
length(unique(subset@elementMetadata$providerName))

snps.mb1 <- snps.from.rsid(rsid = rsid[1:20],
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh37,
                          search.genome = BSgenome.Hsapiens.UCSC.hg19)

# check the alt, ref
results1 <- motifbreakR(snpList = snps.mb1, filterp = TRUE,
                       pwmList = subset,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::SerialParam())
results1 <- results1[paste0(results1$SNP_id,":",results1$ALT) %in% paste0(SNP_orig_info$rsid,":",SNP_orig_info$Effect.Meta)]
results1<-calculatePvalue(results1, granularity = 1e-6)
results_df1 <- as.data.frame(results1, row.names = NULL)

all_result_df <- results_df1

# Initialize empty result container
all_result <- NULL
batch_size <- 10
total_snps <- length(rsid)

# Loop through rsid in batches
for (i in seq(11, total_snps, by = batch_size)) {
  # Calculate end index for current batch
  end_idx <- min(i + batch_size - 1, total_snps)

  print(i)

  # Load SNP data for current batch
  snps.mb1 <- snps.from.rsid(rsid = rsid[i:end_idx],
                             dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh37,
                             search.genome = BSgenome.Hsapiens.UCSC.hg19)

  # Run motifbreakR with error handling
  results1 <- tryCatch({
    motifbreakR(snpList = snps.mb1, filterp = TRUE,
                pwmList = subset, threshold = 1e-4,
                method = "ic", bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                BPPARAM = BiocParallel::SerialParam())
  }, error = function(e) {
    warning("motifbreakR failed for batch starting at index ", i, ": ", e$message)
    NULL
  })

  if (is.null(results1) || length(results1) == 0) {
    warning("No results from motifbreakR for batch starting at index ", i)
    next
  }

  # Filter results
  results1 <- results1[paste0(results1$SNP_id, ":", results1$ALT) %in% paste0(SNP_orig_info$rsid, ":", SNP_orig_info$Effect.Meta)]

  # Process only if results1 is non-empty
  if (length(results1) > 0) {
    results1 <- calculatePvalue(results1, granularity = 1e-6)
    results_df1 <- as.data.frame(results1, row.names = NULL)
    all_result_df <- rbind(all_result_df, results_df1)
  } else {
    warning("No matching SNPs found in batch starting at index ", i)
  }

}

fwrite(all_result_df, "results/motifbreakR/ALLresults_TNBCsig.csv", row.names = FALSE)

all_result_df<-fread("results/motifbreakR/ALLresults_TNBCsig.csv",header = TRUE)
all_result_df[SNP_id=='rs4245736',]
unique(all_result_df[SNP_id=='rs75146030',])

#

all_result_df <- all_result_df %>%
  filter(SNP_id %in% BCAC_motif_TNBC_unqiue$rsid)
length(unique(all_result_df$geneSymbol))

all_result_df <- all_result_df[effect=='strong',]
length(unique(all_result_df$geneSymbol))
length(unique(all_result_df$SNP_id))

motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
TNBC_sig<- left_join(motif_id_list_final, TNBC, by=c('motif','motif.name')) %>% dplyr::filter(TF_p_value <= 0.05/727)

TNBC_motifscan_sig <- TNBC_sig[motif %in% all_result_df$providerId,]
length(unique(TNBC_motifscan_sig$motif.name))
fwrite(TNBC_motifscan_sig,"results/motifbreakR/MotifbreakR_TNBC_significant.csv", row.names = FALSE)

BCAC_all_Anno_df<-fread("results/tfgwas/BCAC_SNP_annotation.csv", header = TRUE)
all_result_df$chrpos37<-paste0(all_result_df$seqnames,":",all_result_df$start)
all_result_df<-left_join(all_result_df,BCAC_all_Anno_df,by=c('chrpos37',"seqnames","start","end"))

TNBC_sig_motifid<-TNBC_sig$motif
all_result_df$TNBC_sig<-ifelse(all_result_df$providerId%in%TNBC_sig_motifid,"Y","")

TNBC_unique<-fread("results/tfgwas/tfgwas_significant_TNBC_unique.csv",header = TRUE)
TNBC_unique_motifid<-TNBC_unique$motif
all_result_df$TNBC_unique<-ifelse(all_result_df$providerId%in%TNBC_unique_motifid,"Y","")

all_result_df<-unique(all_result_df)
all_result_df<-left_join(all_result_df, BCAC_motif_TNBC_unqiue, by='chrpos37')
all_result_df <- all_result_df %>% filter((REF==EA & ALT==OA)|(REF==OA & ALT==EA))
all_result_df<-unique(all_result_df)
fwrite(all_result_df, "results/motifbreakR/ALLresults_TNBCsig_annotation.csv", row.names = FALSE)

MDM4 <- all_result_df%>% filter(SYMBOL=='MDM4')%>%unique()
# 67 SNPs

# # [1] 40
#
# TNBC_unique[!TNBC_unique$motif.name %in% TNBC_motifscan_uniq$motif.name,motif.name]
# # [1] "HLF"         "NFE2L2"      "NFIC"        "ERG"         "GLI2"        "CEBPG"       "TEF"         "E2F2"
# # [9] "PPARA::RXRA" "BHLHA15"     "HOXA4"       "HOXA6"       "HOXB6"       "TBX6"        "TCF21"       "TFAP4"
# # [17] "NKX2-4"
#
# heatmap
enriched_motifs_list<-fread("results/singlecell/enriched_motifs_list_uniqDAP_ALL.csv", header = TRUE)
# Get all unique motif names across all clusters
all_motifs <- unique(enriched_motifs_list$motif.name)

# Get all cluster names
all_clusters <- unique(enriched_motifs_list$cluster)

# Initialize an empty matrix
motif_matrix <- matrix(
  0,
  ncol = length(all_clusters),
  nrow = length(all_motifs),
  dimnames = list( all_motifs,all_clusters)
)

# Fill the matrix with p-values for enriched motifs in each cluster
for (cluster1 in all_clusters) {
  # Get the enriched motifs and their p-values for the current cluster
  enriched_data <- enriched_motifs_list[enriched_motifs_list$cluster==cluster1,]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif.name[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# Convert the matrix to a data frame for better readability
motif_df<- as.data.frame(motif_matrix)
motif_df <- motif_df[rownames(motif_df) %in% TNBC_motifscan_uniq$motif.name,]
motif_matrix <- as.matrix(motif_df)

library(ComplexHeatmap)
library(dendextend)
set.seed(2025)
heatmap_obj_orig <-Heatmap(motif_matrix,
                           name = "-log10 pval",          # Title of legend
                           row_names_gp = gpar(fontsize = 7),
                           row_names_side = "left",
                           column_names_centered = TRUE,
                           width = unit(50, "mm"),
                           col = circlize::colorRamp2(c(0,100), c("#483d8b","red1")),
                           km =4,                        # K-means clustering with 10 clusters
                           rect_gp = gpar(col = "grey99"), # Set cell border color (instead of border_color)
                           show_row_dend = FALSE,         # Hide row dendrogram to avoid km labels
                           show_column_dend = FALSE,
                           column_split = 4,                 # Ensure 3 clusters without labels
                           row_title = NULL ,
                           column_title = NULL
)
heatmap_obj_orig
pdf("results/motifbreakR/heatmap_enrichedmotif_TNBC_all_5e-4.pdf", width = 4, height =6)
draw(heatmap_obj_orig)
dev.off()

###################################################################
#                          SNP to motifs                          #
###################################################################
TNBC_overlap <- all_result_df[geneSymbol %in% TNBC_motifscan_uniq$motif.name,]
table(TNBC_overlap$geneSymbol)
# ATF2        E2F4       EOMES ERF::HOXB13       FIGLA  FOSB::JUNB FOSL1::JUNB
# FOSL2::JUNB FOSL2::JUND   MAF::NFE2       MEF2B         MGA       MYOD1     NEUROD1
# NHLH2       RUNX3        TBX1       TBX15        TBX4        TBX5       ZBTB6
# ZNF85

table(TNBC_overlap$SNP_id)
# rs10415471  rs11086068 rs112199989 rs113272353  rs11666036 rs117304397 rs117934125
# rs3745187   rs3923710   rs4255510   rs4460776   rs4754322    rs651030   rs7247493
# rs8100463   rs8108174      rs8170

SNP_orig_info<-fread('results/motifbreakR/motifbreakR_SNPinfo_TNBC.csv', header = TRUE)
TNBC_overlap$rsid <- TNBC_overlap$SNP_id
TNBC_overlap$chrpos37 <- paste0(TNBC_overlap$seqnames,":",TNBC_overlap$start)
TNBC_overlap <- left_join(TNBC_overlap, SNP_orig_info, by=c('rsid','chrpos37'))
head(TNBC_overlap)

library(highcharter)

TNBC_overlap$TNBC_Status <- ifelse(TNBC_overlap$Triple_Neg_log_or_meta > 0, "Increased TNBC", "Decreased TNBC")
TNBC_overlap$Allele_Effect <- ifelse(TNBC_overlap$alleleDiff > 0, "Enhanced", "Disputed")
TNBC_overlap$chrpos_ALT_REF <- paste(TNBC_overlap$seqnames, TNBC_overlap$start,
                                             TNBC_overlap$ALT, TNBC_overlap$REF, sep = "_")
TNBC_overlap$Gene_Allele <- paste(TNBC_overlap$geneSymbol, TNBC_overlap$Allele_Effect, sep = "_")

data <- TNBC_overlap %>% select(TNBC_Status,chrpos_ALT_REF,Gene_Allele )
#

hchart(data_to_sankey(data), "sankey", name = "TNBC") %>%
  hc_add_theme(
    hc_theme_merge(
      hc_theme_economist(),
      hc_theme(
        chart = list(
          backgroundColor = "#FFFFFF"  # Set background to white
        )
      )
    )
  )

########## manhattan plot #########################
#############################################
motifs_TNBC <- all_result_df[TNBC_unique=='Y']

BCAC_all <- fread('results/gwas/bcac_gwas_filtered.csv', header = TRUE)
BCAC_TNBC_unique <- BCAC_all%>% select(chr.iCOGs, Position.iCOGs, chrpos37,Triple_Neg_P_meta)%>%unique()
names(BCAC_TNBC_unique)<-c("chr",'pos','chrpos37',"TNBC_pval")
BCAC_TNBC_unique$TNBC_pval<-ifelse(BCAC_TNBC_unique$TNBC_pval==0,1e-300,
                                   ifelse(BCAC_TNBC_unique$TNBC_pval==1,0.99999999999999999999999,
                                          BCAC_TNBC_unique$TNBC_pval))
BCAC_TNBC_unique$chrpos37<- paste0("chr",BCAC_TNBC_unique$chrpos37)

snp_data<-fread("snp_index.txt", header = TRUE)
snp_data$chrpos37 <- paste0('chr',snp_data$CHR,':',snp_data$POS)
names(snp_data)[1] <- 'rsid'
snp_data  <- snp_data  %>% dplyr::select(chrpos37, rsid, REF, ALT)

# join rsid
BCAC_TNBC_unique <- left_join(BCAC_TNBC_unique, snp_data, by='chrpos37')
BCAC_TNBC_unique<-BCAC_TNBC_unique[!is.na(rsid)&!is.na(TNBC_pval),]

library(qqman)
png("results/motifbreakR/motifbreakR_manhattan_TNBC.png", width = 16, height = 8, units = "in", res = 150)
par(mar = c(5, 5, 1, 1))
manhattan_local(BCAC_TNBC_unique,
                chr = "chr",
                bp = "pos",
                p = "TNBC_pval",
                snp = "rsid",
                col = c("grey80", "lightblue"),
                chrlabs = NULL,
                suggestiveline =FALSE,
                genomewideline = -log10(5e-08),
                highlight = unique(motifs_TNBC$SNP_id),
                logp = TRUE,
                annotateTop = FALSE)
dev.off()

############################################################
# summerize
library(data.table)
library(dplyr)
setwd(TNBC_ROOT)

SNP_orig_info<-fread('results/motifbreakR/motifbreakR_SNPinfo_TNBC.csv', header = TRUE)
BCAC_filt <- fread('results/gwas/bcac_gwas_filtered.csv', header = TRUE)
BCAC_filt$chrpos37<-paste0("chr",BCAC_filt$chrpos37)

all_result_df<-fread( "results/motifbreakR/ALLresults_TNBCsig_annotation.csv",header = TRUE)
length(unique(all_result_df$geneSymbol))

all_result_df <- all_result_df[effect=='strong',]

motifs_TNBC <- left_join(all_result_df, BCAC_filt, by=c('chrpos37'))
motifs_TNBC$rsid <- motifs_TNBC$SNP_id
head(motifs_TNBC)
length(unique(motifs_TNBC$rsid))
# 332 SNPs
length(unique(motifs_TNBC$geneSymbol)) #motif
# 454 motifs
length(unique(motifs_TNBC$SYMBOL))

# TNBC significant
length(unique(motifs_TNBC[motifs_TNBC$TNBC_sig=='Y']$rsid))
length(unique(motifs_TNBC[motifs_TNBC$TNBC_sig=='Y']$geneSymbol))
length(unique(motifs_TNBC[motifs_TNBC$TNBC_sig=='Y']$SYMBOL))

# TNBC unique
length(unique(motifs_TNBC[TNBC_unique=='Y'&annotation_category=='Promoter']$rsid))
# 55 / 18 SNPs
length(unique(motifs_TNBC[TNBC_unique=='Y'&annotation_category=='Promoter']$geneSymbol))
# 40 / 20(promoter)
length(unique(motifs_TNBC[TNBC_unique=='Y'&annotation_category=='Promoter']$SYMBOL))
# 20 /10 promoter

table(motifs_TNBC$ALT==motifs_TNBC$EA)

# how many initial SNPs at promoters and overlap TFs
BCAC_all_Anno_df<-fread("results/tfgwas/BCAC_SNP_annotation.csv", header = TRUE)
unique_SNP_anno <- left_join(unique_SNP, BCAC_all_Anno_df, by='chrpos37') %>% unique()
length(unique(unique_SNP_anno$chrpos37))
length(unique(unique_SNP_anno[annotation_category=='Promoter']$chrpos37))

fwrite(motifs_TNBC,'results/motifbreakR/motifbreakR_TNBC_all.csv',row.names = FALSE)

# Depmap
##########################################################################
####################   25Q3 gene effect              #####################
setwd(file.path(DEPMAP_DIR, ""))
ge_df <- fread("CRISPRGeneEffect_25Q3.csv", header = TRUE)   # return a data.frame
model  <- fread("Model.csv", header = TRUE)
model_TNBC <- model %>% filter(OncotreeLineage=='Breast', ModelSubtypeFeatures%in%c('TNBC','luminal TNBC',
                                                                                    'basal_B TNBC','basal_A TNBC',
                                                                                    'basal TNBC'))
                               # PrimaryOrMetastasis=='Primary',SampleCollectionSite=='breast')

#                                                                                     'basal_B TNBC','basal_A TNBC',
#                                                                                     'basal TNBC'))
#                                                                                     'basal_B TNBC','basal_A TNBC',
#                                                                                     'basal TNBC'),
#                                PrimaryOrMetastasis=='Metastatic')
# Quick sanity check
dim(ge_df)
head(colnames(ge_df), 10)

# DepMap often uses names like "A1BG (1)" – keep the gene symbol only.
colnames(ge_df) <- sub(" \\(.*\\)", "", colnames(ge_df))

#
setwd(TNBC_ROOT)
motifs_TNBC<-fread('results/motifbreakR/motifbreakR_TNBC_all.csv', header = TRUE)
target_gene_list<-c("DDA1","MDM4")

# filter gene
ge_df_tnbc <- ge_df %>% filter(V1%in%model_TNBC$ModelID) %>% dplyr::select(V1,target_gene_list)

ge_long <- melt(
  ge_df_tnbc,
  id.vars = "V1",
  variable.name = "gene",
  value.name = "gene_effect"
)

# 4. Count dependent cell lines per gene (effect < -0.5)
dep_counts <- ge_long %>%
  group_by(gene) %>%
  summarise(n_dependent = sum(gene_effect < -0.5, na.rm = TRUE)) %>%
  ungroup()

dep_counts[dep_counts$n_dependent>0,]
# gene     n_dependent
# <fct>          <int>
#   1 DDA1              11
# 2 MDM4               1
# 3 DAXX               4
# 4 GTPBP3             1
# 5 WDR43             25
# 6 RPS18             25
# 7 TMEM185B           1
# 8 MRPL34            21
# 9 NPAT              23
# 10 USE1               6

# TNBC unique
# gene  n_dependent
# <fct>       <int>
#   1 DDA1           11
# 2 DAXX            4
# 3 MDM4            1
# 4 WDR43          25
# 5 USE1            6

# 5. Merge count back to the long data
ge_long <- merge(ge_long, dep_counts, by = "gene", all.x = TRUE)
ge_long[is.na(n_dependent), n_dependent := 0]

# 6. Create the y-axis label: "gene (n_dependent)"
ge_long[, y_label := paste0(gene, " (", n_dependent, ")")]

# Order genes by decreasing number of dependent lines (optional but nice)
ge_long[, y_label := factor(y_label,
                            levels = unique(y_label[order(-n_dependent)]))]

library(ggplot2)
ge_long[order(-ge_long$n_dependent),]%>%dplyr::filter(n_dependent>=1)
library(dplyr)
ge_long <- ge_long %>%
  group_by(gene) %>%
  mutate(median_effect = median(gene_effect, na.rm = TRUE)) %>%
  ungroup()
ge_long$y_label <- reorder(ge_long$y_label, -ge_long$median_effect)
names(ge_long)[2]<-'ModelID'
ge_long<-left_join(ge_long, model_TNBC, by='ModelID')

# Create the boxplot with reordered genes
library(ggplot2)
library(RColorBrewer)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(14)
setwd(TNBC_ROOT)
jpeg('results/motifbreakR/motifbreakR_DepMap_TNBC_unique_essential_targetgenes.jpg',width=14, height =3.5, units = 'in', res = 300)
p1<-ge_long %>%
  ggplot(aes(x=gene_effect, y=y_label)) +
  geom_boxplot(outlier.shape = NA, fill="#d3ebeb", color="navyblue")+
  geom_jitter( aes(color = CellLineName))+
  labs(x = "Gene Effect", y = "Potential target gene\n(No. of dependent cell lines)",color = "Cell line of TNBC") +# Legend title for jitter
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "red") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "lines"),           # Smaller legend keys
        legend.spacing.x = unit(0.1, "cm"))+             # Narrower gap between legend items)+
  # scale_fill_manual(values = mycolors)+
  guides(color = guide_legend(ncol = 9)) +   # 6 columns; adjust as needed
  guides(fill = "none")
p1
dev.off()

# metastatic ##################################
setwd(file.path(DEPMAP_DIR, ""))
ge_df <- fread("CRISPRGeneEffect_25Q3.csv", header = TRUE)   # return a data.frame
model  <- fread("Model.csv", header = TRUE)
model_TNBC <- model %>% filter(OncotreeLineage=='Breast', ModelSubtypeFeatures%in%c('TNBC','luminal TNBC',
                                                                                    'basal_B TNBC','basal_A TNBC',
                                                                                    'basal TNBC'),
                               PrimaryOrMetastasis=='Metastatic')
# Quick sanity check
dim(ge_df)
head(colnames(ge_df), 10)

# DepMap often uses names like "A1BG (1)" – keep the gene symbol only.
colnames(ge_df) <- sub(" \\(.*\\)", "", colnames(ge_df))

#
setwd(TNBC_ROOT)
motifs_TNBC<-fread('results/motifbreakR/motifbreakR_TNBC_all.csv', header = TRUE)
motifs_TNBC<-motifs_TNBC%>% filter(TNBC_unique=="Y")
target_gene_list<-c("DDA1","TERT","OCEL1","PIK3C2B","ABHD8","MDM4","ANKLE1","RALB","ATM")

# filter gene
ge_df_tnbc <- ge_df %>% filter(V1%in%model_TNBC$ModelID) %>% dplyr::select(V1,target_gene_list)
ge_long <- melt(
  ge_df_tnbc,
  id.vars = "V1",
  variable.name = "gene",
  value.name = "gene_effect"
)

# 4. Count dependent cell lines per gene (effect < -0.5)
dep_counts <- ge_long %>%
  group_by(gene) %>%
  summarise(n_dependent = sum(gene_effect < -0.5, na.rm = TRUE)) %>%
  ungroup()

dep_counts[dep_counts$n_dependent>0,]
# gene   n_dependent
# <fct>        <int>
#   1 DDA1             5
# 2 MDM4             1
# 3 DAXX             3
# 4 GTPBP3           1
# 5 WDR43           12
# 6 RPS18           12
# 7 MRPL34          10
# 8 NPAT            11
# 9 USE1             4

# 5. Merge count back to the long data
ge_long <- merge(ge_long, dep_counts, by = "gene", all.x = TRUE)
ge_long[is.na(n_dependent), n_dependent := 0]

# 6. Create the y-axis label: "gene (n_dependent)"
ge_long[, y_label := paste0(gene, " (", n_dependent, ")")]

# Order genes by decreasing number of dependent lines (optional but nice)
ge_long[, y_label := factor(y_label,
                            levels = unique(y_label[order(-n_dependent)]))]

library(ggplot2)
ge_long[order(-ge_long$n_dependent),]%>%dplyr::filter(n_dependent>=1)
library(dplyr)
ge_long <- ge_long %>%
  group_by(gene) %>%
  mutate(median_effect = median(gene_effect, na.rm = TRUE)) %>%
  ungroup()
ge_long$y_label <- reorder(ge_long$y_label, -ge_long$median_effect)
names(ge_long)[2]<-'ModelID'
ge_long<-left_join(ge_long, model_TNBC, by='ModelID')

mycolors <- colorRampPalette(brewer.pal(9, "Set3"))(14)
jpeg('results/motifbreakR/motifbreakR_DepMap_TNBC_essential_targetgenes_metastatic.jpg',width=12, height =4.8, units = 'in', res = 300)
p2<-ge_long %>%
  ggplot(aes(x=gene_effect, y=y_label)) +
  geom_boxplot(outlier.shape = NA, fill="grey90", color="darkblue")+
  geom_jitter( aes(color = CellLineName))+
  labs(x = "Gene Effect", y = "Potential target gene (No. of dependent cell lines)",color = "Cell line of \nmetastatic TNBC") +# Legend title for jitter
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "red") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "lines"),           # Smaller legend keys
        legend.spacing.x = unit(0.1, "cm"))+             # Narrower gap between legend items)+
  # scale_fill_manual(values = mycolors)+
  guides(color = guide_legend(ncol = 9)) +   # 6 columns; adjust as needed
  guides(fill = "none")
dev.off()

library(cowplot)
jpeg('results/motifbreakR/motifbreakR_DepMap_TNBC_essential_both_celllines.jpg',width=17, height =5, units = 'in', res = 300)
plot_grid(p1, p2, ncol = 2)
dev.off()

setwd(TNBC_ROOT)
jpeg('results/motifbreakR/motifbreakR_DepMap_TNBC_primary_breast.jpg',width=11, height =8, units = 'in', res = 140)
ge_long %>%
  ggplot(aes(x=gene_effect, y=y_label)) +
  geom_boxplot(outlier.shape = NA, fill="grey90", color="navyblue")+
  geom_jitter( aes(color = CellLineName))+
  labs(x = "Gene Effect", y = "Potential target gene (No. of dependent cell lines)",color = "Cell line of primary TNBC\n(collected in breast)") +# Legend title for jitter
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "red") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "lines"),           # Smaller legend keys
        legend.spacing.x = unit(0.1, "cm"))+             # Narrower gap between legend items)+
  # scale_fill_manual(values = mycolors)+
  guides(color = guide_legend(ncol = 9)) +   # 6 columns; adjust as needed
  guides(fill = "none")
p2
dev.off()

dep_counts[dep_counts$n_dependent>0,]
# all motifbreakR result
# gene   n_dependent
# <fct>        <int>
#   1 DDA1             6
# 2 DAXX             1
# 3 WDR43           12
# 4 RPS18           12
# 5 MRPL34          10
# 6 NPAT            11
# 7 USE1             2

setwd(TNBC_ROOT)
jpeg('results/motifbreakR/motifbreakR_DepMap_TNBC_all_motifbreakR.jpg',width=13, height =8, units = 'in', res = 140)
ge_long %>%
  ggplot(aes(x=gene_effect, y=y_label)) +
  geom_boxplot(outlier.shape = NA, fill="grey90", color="navyblue")+
  geom_jitter( alpha=0.8,aes(color = CellLineName))+
  labs(x = "Gene Effect", y = "Potential target gene (No. of dependent cell lines)",color = "Cell line of primary TNBC\n(collected in breast)") +# Legend title for jitter
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "red") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "lines"),           # Smaller legend keys
        legend.spacing.x = unit(0.1, "cm"))+             # Narrower gap between legend items)+
  # scale_fill_manual(values = mycolors)+
  guides(color = guide_legend(ncol = 9)) +   # 6 columns; adjust as needed
  guides(fill = "none")
dev.off()

# TNBC significant
# # A tibble: 5 × 2
# gene  n_dependent
# <fct>       <int>
#   1 DDA1           11
# 2 DAXX            4
# 3 MDM4            1
# 4 WDR43          25
# 5 USE1            6

setwd(TNBC_ROOT)
jpeg('results/motifbreakR/motifbreakR_DepMap_TNBC_all.jpg',width=13, height =8, units = 'in', res = 140)
ge_long %>%
  ggplot(aes(x=gene_effect, y=y_label)) +
  geom_boxplot(outlier.shape = NA, fill="grey90", color="navyblue")+
  geom_jitter( alpha=0.8,aes(color = CellLineName))+
  labs(x = "Gene Effect", y = "Potential target gene (No. of dependent cell lines)",color = "Cell line of TNBC") +# Legend title for jitter
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "red") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "lines"),           # Smaller legend keys
        legend.spacing.x = unit(0.1, "cm"))+             # Narrower gap between legend items)+
  # scale_fill_manual(values = mycolors)+
  guides(color = guide_legend(ncol = 9)) +   # 6 columns; adjust as needed
  guides(fill = "none")
dev.off()

##############################################################
# gene expression
atac_subset_qc<-readRDS("results/singlecell/atac_motif.rds")
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc)<-'RNA'

levels(atac_subset_qc)
new.cluster.ids <- c('LASP-1','Fb','LHS','LASP-2','B-Myo',
                     'Endo','T','Pc','Myl','B')

names(new.cluster.ids) <- levels(atac_subset_qc)
atac_subset_qc <- RenameIdents(atac_subset_qc, new.cluster.ids)
atac_subset_qc$celltype_RNA <- Idents(atac_subset_qc)

# new levels
new_levels <- c('LASP-1','LASP-2','LHS','B-Myo','Fb','Endo','T','B','Myl','Pc')
# Reorder the Idents based on sorted levels
Idents(atac_subset_qc) <- factor(Idents(atac_subset_qc), levels = new_levels)

jpeg('results/motifbreakR/motifbreakR_RNA_TNBC_target.jpg',width=7, height =8, units = 'in', res = 300)
#Dot plot of sig motifbreakR target genes
DotPlot(atac_subset_qc, features = c("DDA1","MDM4"),
        cols = c("grey90", "purple2")) + RotatedAxis()+
  scale_size(name = "Percent Expressed", labels = scales::percent_format(scale = 1), limits = c(0, 100))+
  theme(legend.position = "bottom",          # Move legend to bottom
        legend.box = "veritcal",             # Stack color and size legends vertically
        legend.box.just = "center",          # Center the legend box
        legend.justification = "center",
        legend.spacing.y = unit(0.3, "cm"),  # Optional: adjust spacing
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))+
  ylab("Cell type") +                       # Set y-axis label
  xlab("Genes")
dev.off()

# TCGA
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(dunn.test)
library(data.table)
library(purrr)
library(stringr)
library(Seurat)
# read RDS
seurat_TCGA<-readRDS(file.path(TCGA_DIR, "tcga_brca_RNA_seurat.rds"))
order <- c("Basal",  "LumA", "LumB", "Her2","NAT")
# Reorder the identities in the Seurat object
Idents(seurat_TCGA) <- factor(Idents(seurat_TCGA), levels = order)

count_matrix <- GetAssayData(seurat_TCGA, assay = "RNA", slot = "counts")
# Apply the filter: keep genes where at least 5 samples have counts >= 10
keep <- rowSums(count_matrix >= 10) >= 5
# Subset the Seurat object to keep only those features (genes)
seurat_TCGA <- subset(seurat_TCGA, features = rownames(count_matrix)[keep])

library(ggplot2)
library(vioplot)

jpeg('results/motifbreakR/motifbreakR_TCGA_TNBC_target.jpg',width=12, height =8, units = 'in', res = 300)
features <- c("DDA1", "MDM4")
plots <- VlnPlot(object = seurat_TCGA, features = features, combine = FALSE)
plots <- lapply(plots, function(p) p + geom_boxplot(width = 0.2, fill = "black", color = "black") +
                  stat_summary(fun = median, geom = "point", color = "white", size = 2.5) +
                  theme(legend.position = "none",
                        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)))
p <- patchwork::wrap_plots(plots)
print(p)
dev.off()

# TCGA-BRCA
data<-readRDS("results/tcga/TCGA_BRCA_data_all.rds")
metadata <- colData(data) %>% as.data.frame()

metadata$paper_BRCA_Subtype_PAM50<- ifelse(metadata$sample_type== "Solid Tissue Normal",'NAT',
                                           metadata$paper_BRCA_Subtype_PAM50)

metadata_df <- DataFrame(metadata, row.names = rownames(metadata))
colData(data) <- metadata_df

keep_samples <- metadata$sample_type %in% c("Primary Tumor","Solid Tissue Normal") &
  metadata$paper_BRCA_Subtype_PAM50 %in% c( 'Basal', 'Her2', 'LumA', 'LumB','NAT') &
  metadata$race %in%'white' & metadata$gender%in%'female'
table(keep_samples)
# keep_samples
table(metadata$paper_BRCA_Subtype_PAM50)
# Basal   Her2   LumA   LumB    NAT Normal
Tumor_data <- data[, keep_samples]
table(Tumor_data$paper_BRCA_Subtype_PAM50)
# Basal  Her2  LumA  LumB   NAT

# class: RangedSummarizedExperiment
# dim: 60660 836
# assays(6): unstranded stranded_first ... fpkm_unstrand fpkm_uq_unstrand
# rownames(60660): ENSG00000000003.15 ENSG00000000005.6 ... ENSG00000288674.1 ENSG00000288675.1
# rowData names(10): source type ... hgnc_id havana_gene
# TCGA-BH-A0H3-01A-11R-A12P-07 TCGA-E2-A1IU-01A-11R-A14D-07
# colData names(93): barcode patient ... paper_PARADIGM.Clusters paper_Pan.Gyn.Cluster

# Step 3: Extract filtered expression data and metadata
Tumor_expression_data <- assay(Tumor_data)# Raw counts, row.name is gene id
Tumor_metadata <- colData(Tumor_data) %>% as.data.frame()

# Step 4: Extract PAM50 subtypes
pam50_subtypes <- Tumor_metadata$paper_BRCA_Subtype_PAM50
table(pam50_subtypes)
# pam50_subtypes
# Basal  Her2  LumA  LumB   NAT

# DEseq2
library(DESeq2)
library(dplyr)
features <- c( "DDA1", "MDM4")
row_data <- rowData(Tumor_data)
features_ensembl <- rownames(row_data)[which(row_data$gene_name %in% features)]
ensembl_map <- data.frame(
  Gene = row_data[features_ensembl, "gene_name"],
  Ensembl = features_ensembl
)

dds <- DESeqDataSetFromMatrix(
  countData = Tumor_expression_data,
  colData = Tumor_metadata,
  design = ~ paper_BRCA_Subtype_PAM50
)
# Filter low count genes for better performance
keep <- rowSums(counts(dds) >= 1) >= 5 # remove not expressed genes, you can change the threshold
dds <- dds[keep, ]
# class: DESeqDataSet
# dim: 56633 836
# rownames(56633): ENSG00000000003.15 ENSG00000000005.6 ... ENSG00000288674.1 ENSG00000288675.1
# rowData names(0):
# TCGA-E2-A1IU-01A-11R-A14D-07
# colData names(93): barcode patient ... paper_PARADIGM.Clusters paper_Pan.Gyn.Clusters

# Step 6: Run DESeq
dds <- estimateSizeFactors(dds) #Estimate the size factors
colData(dds)$paper_BRCA_Subtype_PAM50 <- relevel(
  factor(colData(dds)$paper_BRCA_Subtype_PAM50),
  ref = "Basal"
)
dds <- DESeq(dds)

res_BasalVSLumA <- results(dds, contrast = c("paper_BRCA_Subtype_PAM50", "Basal", "LumA"))
res_BasalVSLumB <- results(dds, contrast = c("paper_BRCA_Subtype_PAM50", "Basal", "LumB"))
res_BasalVSHer2 <- results(dds, contrast = c("paper_BRCA_Subtype_PAM50", "Basal", "Her2"))
res_BasalVSNAT <- results(dds, contrast = c("paper_BRCA_Subtype_PAM50", "Basal", "NAT"))

# Convert each to data.frame and add contrast column
df_LumA <- as.data.frame(res_BasalVSLumA)%>%filter(rownames(.) %in% ensembl_map$Ensembl)
df_LumA$Contrast <- "Basal vs LumA"
df_LumA$Ensembl <- rownames(df_LumA)
df_LumA<-left_join(df_LumA,ensembl_map,by='Ensembl')

df_LumB <- as.data.frame(res_BasalVSLumB)%>%filter(rownames(.) %in% ensembl_map$Ensembl)
df_LumB$Contrast <- "Basal vs LumB"
df_LumB$Ensembl <- rownames(df_LumB)
df_LumB<-left_join(df_LumB,ensembl_map,by='Ensembl')

df_Her2 <- as.data.frame(res_BasalVSHer2)%>%filter(rownames(.) %in% ensembl_map$Ensembl)
df_Her2$Contrast <- "Basal vs Her2"
df_Her2$Ensembl <- rownames(df_Her2)
df_Her2<-left_join(df_Her2,ensembl_map,by='Ensembl')

df_NAT <- as.data.frame(res_BasalVSNAT)%>%filter(rownames(.) %in% ensembl_map$Ensembl)
df_NAT$Contrast <- "Basal vs NAT"
df_NAT$Ensembl <- rownames(df_NAT)
df_NAT<-left_join(df_NAT,ensembl_map,by='Ensembl')

# Combine all into one data.frame
combined_results <- rbind(df_LumA, df_LumB, df_Her2, df_NAT)

library(tidyr)
wide_results <- combined_results %>%
  select(Gene, Ensembl, baseMean,Contrast, log2FoldChange, lfcSE, padj) %>%
  pivot_wider(
    names_from = Contrast,
    values_from = c( log2FoldChange, lfcSE, padj),
    names_sep = "_",
    names_prefix = "_"  # Results in e.g., "log2FoldChange_Basal vs LumA"
  )
# add * to pval
padj_cols <- c("padj__Basal vs LumA", "padj__Basal vs LumB", "padj__Basal vs Her2", "padj__Basal vs NAT")

# For each padj column, create a new formatted column
for (col in padj_cols) {
  new_col <- paste0("padj_formatted__", gsub("padj__", "", col))
  wide_results <- wide_results %>%
    mutate(!!new_col := {
      padj_val <- .data[[col]]
      formatted <- ifelse(is.na(padj_val), "NA",
                          ifelse(padj_val < 0.001,
                                 sprintf("%.2e", padj_val),
                                 sprintf("%.3f", padj_val)))
      ifelse(padj_val < 0.05 & !is.na(padj_val),
             paste0(formatted, "*"),
             formatted)
    })
}
fwrite(wide_results,'results/motifbreakR/motifbreakR_TCGA_targetgenes_DEseq2.csv',row.names = FALSE)

########DESeq2 plot
library(dplyr)
combined_results <- combined_results %>%
  mutate(padj_formatted = ifelse(padj < 0.05 & !is.na(padj), "*", ""))

# Filter to the relevant genes (assuming updated_TNBC_motif_list contains the genes of interest)
plot_data <- combined_results %>%
  filter(Gene %in% features) %>%
  # Set Contrast as factor with specific order for faceting
  mutate(Contrast = factor(Contrast, levels = c("Basal vs LumA", "Basal vs LumB", "Basal vs Her2", "Basal vs NAT"))) %>%
  # Add direction for coloring
  mutate(direction = ifelse(log2FoldChange < 0, "negative", "positive"))

# Get gene order based on log2FoldChange in "Basal vs NAT" from high to low
nat_order <- plot_data %>%
  filter(Contrast == "Basal vs NAT") %>%
  arrange(log2FoldChange) %>%
  pull(Gene) %>%
  unique()

# Apply the factor order to Gene
plot_data <- plot_data %>%
  mutate(Gene = factor(Gene, levels = nat_order))

# Load ggplot2 if not already
library(ggplot2)

# Create the faceted plot with free x scales, custom colors, and 2x2 layout
p <- ggplot(plot_data, aes(x = log2FoldChange, y = Gene, fill = direction)) +
  geom_bar(stat = "identity", width = 0.8) +  # Bar length shows magnitude, direction left/right of zero
  geom_errorbarh(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), width = 0.3) +
  geom_text(aes(label = padj_formatted,
                hjust = ifelse(log2FoldChange > 0, -3, 4)),
            size = 6, color = "black") +  # Add only "*" for significance
  geom_vline(xintercept = 0)+ #linetype = "dashed", color = "gray") +  # Zero line for reference
  facet_wrap(~ Contrast, scales = "free_x", ncol = 2) +  # 2 columns for 2x2 grid, free x scales
  labs(#title = "Log2 Fold Change Across Contrasts (Basal as Reference)",
    x = "Log2 Fold Change",
    y = "Gene") +
  scale_fill_manual(values = c("negative" = "red", "positive" = "darkblue"), guide = "none") +  # Custom colors, no legend
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12, color = 'black'),
        axis.title.y = element_text(size = 12, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(size = 9, color = 'black'),
        strip.text = element_text(size = 12,  color = 'black'),
        strip.text.y = element_text(size = 12, angle = 0,  color = 'black'),
        strip.text.x = element_text(size = 12, angle = 0,  color = 'black'),
        plot.title = element_text(size = 12, face = "bold", color = 'black'),
        panel.grid = element_blank(),  # Remove all panel grid lines
        axis.line = element_line(color = "black")) +  # Add x and y axis lines
  coord_cartesian(xlim = c(min(plot_data$log2FoldChange-0.1, na.rm = TRUE),
                           max(plot_data$log2FoldChange+0.1, na.rm = TRUE)))

# Display the plot
print(p)

ggsave("results/motifbreakR/motifbreakR_deseq2_logfc_plot.png", p, width = 12, height = 8, dpi = 300)

# KW analysis
# Step 5: Normalize expression data (log2(CPM+1)) Counts per Million or Reads per Kilobase per Million
cpm <- edgeR::cpm(Tumor_expression_data, log = TRUE, prior.count = 1)
Tumor_expression_data <- cpm

gene_map <- fread(file.path(TCGA_DIR, "gencode.v36.annotation.gtf.gene.probemap"), header = TRUE)
updated_TNBC_gene_id <- gene_map[gene %in% features,id]
genes_to_test <- intersect(rownames(Tumor_expression_data), updated_TNBC_gene_id)
Tumor_expression_data.TNBC <- Tumor_expression_data[genes_to_test, ]

kw_results <- data.frame(
  id = genes_to_test,
  KW_p_value = NA
)

for (i in 1:length(genes_to_test)) {
  id <- genes_to_test[i]
  gene_expr <- Tumor_expression_data.TNBC[id, ]

  kw_test <- kruskal.test(gene_expr ~ pam50_subtypes)

  kw_results$KW_p_value[i] <- kw_test$p.value
}

# Adjust p-values for multiple testing (Benjamini-Hochberg)
kw_results$KW_p_adj <- p.adjust(kw_results$KW_p_value, method = "BH")

# Filter significant genes (adjusted p < 0.05)
sig_genes <- kw_results %>% filter(KW_p_adj < 0.05)

# Step 8: Perform Dunn's test for significant genes
dunn_results <- list()
for (id in sig_genes$id) {
  gene_expr <- Tumor_expression_data.TNBC[id, ]
  dunn_test <- dunn.test(
    x = gene_expr,
    g = pam50_subtypes,
    method = "bh"  # Benjamini-Hochberg adjustment
  )
  dunn_df <- data.frame(
    id = id,
    Comparison = dunn_test$comparisons,
    Z = dunn_test$Z,
    P_unadj = dunn_test$P,
    P_adj = dunn_test$P.adjusted
  )
  dunn_results[[id]] <- dunn_df
}

# Combine Dunn's test results
dunn_results_combined <- do.call(rbind, dunn_results)

# Step 9: Save results
kw_results <- left_join(kw_results, dplyr::select(gene_map, id, gene),by='id')
dunn_results_combined <- left_join(dunn_results_combined, dplyr::select(gene_map, id, gene),by='id')
fwrite(kw_results, "results/motifbreakR/TCGA_BRCA_PAM50_KruskalWallis_targetgene.csv", row.names = FALSE)
fwrite(dunn_results_combined, "results/motifbreakR/TCGA_BRCA_PAM50_DunnTest_targetgene.csv", row.names = FALSE)

library(dplyr)
library(tidyr)

# Reshape dunn_results_combined by pivoting Comparison into column names
dunn_results_reshaped <- dunn_results_combined %>%
  # Replace " - " with "_" in Comparison for cleaner column names
  mutate(Comparison = gsub(" - ", "_", Comparison)) %>%
  # Pivot the data to create columns for Z, P_unadj, and P_adj for each Comparison
  pivot_wider(
    id_cols = gene,
    names_from = Comparison,
    values_from = c(Z, P_unadj, P_adj),
    names_glue = "{Comparison}_{.value}"
  ) %>%
  # Rename Z to Z_score in column names
  rename_with(~gsub("_Z$", "_Z_score", .), ends_with("_Z"))

dunn_results_reshaped<- left_join(dplyr::select(kw_results, gene), dunn_results_reshaped,by='gene')
names(dunn_results_reshaped)

# Select gene and Basal-related columns
basal_columns <- c(
  "gene",
  "Basal_LumA_Z_score", "Basal_LumA_P_adj",
  "Basal_LumB_Z_score", "Basal_LumB_P_adj", "Basal_Her2_Z_score","Basal_Her2_P_adj",
  "Basal_NAT_Z_score", "Basal_NAT_P_adj"
)

# Filter the data frame and format p-values
basal_table <- dunn_results_reshaped %>%
  dplyr::select(all_of(basal_columns)) %>%
  mutate(across(contains("P_adj"), ~ifelse(.x < 0.05,
                                           paste0(format(.x, scientific = TRUE, digits = 2), "*"),
                                           format(.x, scientific = TRUE, digits = 2))),
         across(contains("Z_score"), ~round(.x, digits = 2))
  )

TNBC_table <- left_join(kw_results, basal_table, by='gene')
head(TNBC_table)
fwrite(TNBC_table, "results/motifbreakR/Dunntest_TNBC_TCGA_targetgene.csv", row.names=FALSE)
