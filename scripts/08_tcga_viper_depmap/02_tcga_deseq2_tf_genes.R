# ==============================================================================
# 02_tcga_deseq2_tf_genes.R
# TCGA-BRCA differential expression (DESeq2) of the 61 TF genes.
# Produces: Fig 4A, Table S9
# Inputs  : results/tfgwas/tfgwas_significant_TNBC_unique.csv
#           results/tcga/TCGA_BRCA_data_all.rds
#           results/tcga/TCGA_BRCA_PAM50_metadata.rds
#           results/tcga/TFgenes_TNBCmotif_DEseq2.csv
# Run     : Rscript scripts/08_tcga_viper_depmap/02_tcga_deseq2_tf_genes.R   (from the repository root)
# ==============================================================================
## 02_tcga_deseq2_tf_genes.R — Fig 4A, Table S9: TCGA-BRCA differential expression of the 61 TF genes that encode
## the 57 TNBC-specific motifs (TF-GWAS, tfgwas_significant_TNBC_unique.csv; heterodimers split).
##
## Analysis (reproduces Table S9; the published table is read from
## results/tcga/TFgenes_TNBCmotif_DEseq2.csv for verification):
##   samples : TCGA-BRCA STAR counts; Primary Tumor + Solid Tissue Normal (NAT); PAM50 Basal, Her2,
##             LumA, LumB; race "white"; gender "female" -> 836 samples
##             (Basal 116, Her2 39, LumA 446, LumB 131, NAT 104)
##   model   : DESeq2 ~ paper_BRCA_Subtype_PAM50 on ALL genes (low-count filter: >= 1 count in
##             >= 5 samples), Basal as reference; Wald test; padj = BH over all genes per contrast
##   genes   : the 57 motif names split on "::" -> 61 symbols, matched by rowData gene_name
## Sample counts before the female + white restriction are 197 / 571 / 209 / 82 / 113; the
## analysis uses the 836 samples above.
##
## The genome-wide DESeq2 fit is cached as results/tcga/TCGA_BRCA_DESeq2_BasalRef_4contrasts_allgenes.rds
## and reused by 05_expression_activity_panel.R and 08_tcga_mdm4_rps18.R. The script verifies
## that it reproduces the published table before plotting.
## Run with R 4.5.1 (see environment/README.md).
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressMessages({library(SummarizedExperiment); library(DESeq2); library(data.table); library(ggplot2)})
setwd(TNBC_ROOT)
OUT <- "results/tf_activity/"; FIG <- "results/tf_activity/figures/tcga"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
CONTRASTS <- c("Basal vs LumA", "Basal vs LumB", "Basal vs Her2", "Basal vs NAT")

## ---- 61 TF genes ----
motifs <- unique(fread("results/tfgwas/tfgwas_significant_TNBC_unique.csv")$motif.name)
TF61   <- unique(unlist(strsplit(motifs, "::")))
stopifnot(length(motifs) == 57, length(TF61) == 61)

## ---- genome-wide DESeq2 results (cached) ----
cache <- "results/tcga/TCGA_BRCA_DESeq2_BasalRef_4contrasts_allgenes.rds"
if (file.exists(cache)) {
  RES <- readRDS(cache)
} else {
  data <- readRDS("results/tcga/TCGA_BRCA_data_all.rds")
  metadata <- readRDS("results/tcga/TCGA_BRCA_PAM50_metadata.rds")      # NAT already recoded
  colData(data) <- DataFrame(metadata, row.names = rownames(metadata))
  keep_samples <- metadata$sample_type %in% c("Primary Tumor", "Solid Tissue Normal") &
    metadata$paper_BRCA_Subtype_PAM50 %in% c("Basal", "Her2", "LumA", "LumB", "NAT") &
    metadata$race %in% "white" & metadata$gender %in% "female"
  Tumor_data <- data[, keep_samples]
  print(table(Tumor_data$paper_BRCA_Subtype_PAM50))
  dds <- DESeqDataSetFromMatrix(countData = assay(Tumor_data, "unstranded"),
                                colData = as.data.frame(colData(Tumor_data)),
                                design = ~ paper_BRCA_Subtype_PAM50)
  dds <- dds[rowSums(counts(dds) >= 1) >= 5, ]
  colData(dds)$paper_BRCA_Subtype_PAM50 <- relevel(factor(colData(dds)$paper_BRCA_Subtype_PAM50), ref = "Basal")
  dds <- DESeq(dds)
  RES <- rbindlist(lapply(c("LumA", "LumB", "Her2", "NAT"), function(g) {
    r <- as.data.frame(results(dds, contrast = c("paper_BRCA_Subtype_PAM50", "Basal", g)))
    data.table(Ensembl = rownames(r), Contrast = paste("Basal vs", g), r)
  }))
  RES[, gene_name := rowData(Tumor_data)[Ensembl, "gene_name"]]
  saveRDS(RES, cache)
}

## ---- the 61 genes (matched by gene_name) ----
plot_data <- RES[gene_name %in% TF61]
setnames(plot_data, "gene_name", "Gene")
dup <- plot_data[, uniqueN(Ensembl), by = Gene][V1 > 1]
if (nrow(dup)) stop("gene symbol maps to >1 Ensembl ID: ", paste(dup$Gene, collapse = ", "))
miss <- setdiff(TF61, plot_data$Gene)
cat(sprintf("genes with DESeq2 results: %d of 61%s\n", uniqueN(plot_data$Gene),
            if (length(miss)) paste0(" | missing: ", paste(miss, collapse = ", ")) else ""))

## ---- verify against the published table (Table S9) ----
pub <- fread("results/tcga/TFgenes_TNBCmotif_DEseq2.csv")
chk <- rbindlist(lapply(CONTRASTS, function(cc) {
  a <- plot_data[Contrast == cc, .(Ensembl, lfc_new = log2FoldChange, padj_new = padj)]
  b <- pub[, .(Ensembl, lfc_pub = get(paste0("log2FoldChange__", cc)), padj_pub = get(paste0("padj__", cc)))]
  cbind(Contrast = cc, merge(a, b, by = "Ensembl"))
}))
cat(sprintf("verification vs TFgenes_TNBCmotif_DEseq2.csv: %d gene x contrast rows | max |dlog2FC| = %.2g | max |dpadj| = %.2g\n",
            nrow(chk), max(abs(chk$lfc_new - chk$lfc_pub), na.rm = TRUE),
            max(abs(chk$padj_new - chk$padj_pub), na.rm = TRUE)))

## ---- table (same layout as TFgenes_TNBCmotif_DEseq2.csv) + summary flags ----
fmt <- function(p) { f <- ifelse(is.na(p), "NA", ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.3f", p)))
                     ifelse(!is.na(p) & p < 0.05, paste0(f, "*"), f) }
W <- dcast(plot_data, Gene + Ensembl + baseMean ~ Contrast, value.var = c("log2FoldChange", "lfcSE", "padj"), sep = "__")
for (cc in CONTRASTS) W[, (paste0("padj_formatted__", cc)) := fmt(get(paste0("padj__", cc)))]
## "significantly upregulated in basal-like tumours compared with other subtypes and NAT"
## (manuscript Results) = FDR < 0.05 and log2FC > 0 in ALL FOUR contrasts
up_all <- Reduce(`&`, lapply(CONTRASTS, function(cc)
  W[[paste0("padj__", cc)]] < 0.05 & W[[paste0("log2FoldChange__", cc)]] > 0))
W[, up_in_Basal_all4 := !is.na(up_all) & up_all]
fwrite(W, file.path(OUT, "TCGA_DESeq2_61TF.csv"))
up_genes <- W[up_in_Basal_all4 == TRUE, Gene]
cat(sprintf("\nup in Basal vs all four comparators (FDR < 0.05): %d genes\n  %s\n",
            length(up_genes), paste(up_genes[order(-W[up_in_Basal_all4 == TRUE, `log2FoldChange__Basal vs NAT`])], collapse = ", ")))
MS13 <- c("E2F2","CDX1","MSANTD3","FOSL1","RUNX3","TBX4","TFAP4","ELK1","CEBPG","FIGLA","ERF","E2F4","HOXA1")
cat("  manuscript list (13) reproduced:", setequal(up_genes, MS13),
    if (!setequal(up_genes, MS13)) paste0(" | only new: ", paste(setdiff(up_genes, MS13), collapse = ","),
                                          " | only manuscript: ", paste(setdiff(MS13, up_genes), collapse = ",")) else "", "\n")

## ---- plot (Fig 4A style; free x per facet, "*" beyond the error bar, clipping off) ----
plot_data[, padj_formatted := ifelse(padj < 0.05 & !is.na(padj), "*", "")]
plot_data[, Contrast := factor(Contrast, levels = CONTRASTS)]
plot_data[, direction := ifelse(log2FoldChange < 0, "negative", "positive")]
plot_data[, star_x := ifelse(log2FoldChange > 0, log2FoldChange + lfcSE, log2FoldChange - lfcSE)]
nat_order <- plot_data[Contrast == "Basal vs NAT"][order(log2FoldChange), unique(Gene)]
plot_data[, Gene := factor(Gene, levels = nat_order)]
p <- ggplot(plot_data, aes(x = log2FoldChange, y = Gene, fill = direction)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), width = 0.3, orientation = "y") +
  geom_text(aes(x = star_x, label = padj_formatted, hjust = ifelse(log2FoldChange > 0, -0.4, 1.4)),
            vjust = 0.75, size = 4.5, color = "black") +
  geom_vline(xintercept = 0) +
  facet_wrap(~ Contrast, scales = "free_x", ncol = 4) +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  labs(x = "Log2 Fold Change", y = "TF gene",
       caption = "* FDR < 0.05 (DESeq2 Wald test, Benjamini–Hochberg-adjusted)") +
  scale_fill_manual(values = c(negative = "red", positive = "darkblue"), guide = "none") +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, color = "black"), axis.text = element_text(size = 11, color = "black"),
        strip.text = element_text(size = 12, color = "black"), panel.grid = element_blank(),
        axis.line = element_line(color = "black"), panel.spacing.x = unit(1.2, "lines"),
        plot.caption = element_text(size = 10, color = "black", hjust = 0))
f <- file.path(FIG, "TCGA_DESeq2_log2FC_61TF.png")
ggsave(f, p, width = 13, height = 1.5 + 0.2 * uniqueN(plot_data$Gene), dpi = 300, bg = "white")
cat("\nwrote", f, "and", file.path(OUT, "TCGA_DESeq2_61TF.csv"), "\n")
