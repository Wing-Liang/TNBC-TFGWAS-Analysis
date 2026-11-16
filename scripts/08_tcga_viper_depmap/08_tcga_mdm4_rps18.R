# ==============================================================================
# 08_tcga_mdm4_rps18.R
# MDM4 and RPS18 expression across PAM50 subtypes, TCGA-BRCA.
# Produces: Table S14
# Inputs  : results/tf_activity/viper_vst/steps/vst_836_blindFALSE.rds
#           results/tcga/TCGA_BRCA_PAM50_metadata.rds
#           results/tcga/TCGA_BRCA_DESeq2_BasalRef_4contrasts_allgenes.rds
# Outputs : results/target_genes/TCGA_MDM4_RPS18_expression_by_subtype.csv
#           results/target_genes/TCGA_MDM4_RPS18_DESeq2.csv
#           results/target_genes/figures/tcga/TCGA_MDM4_RPS18_boxplot_by_subtype.png
# Run     : Rscript scripts/08_tcga_viper_depmap/08_tcga_mdm4_rps18.R   (from the repository root)
# ==============================================================================
## 08_tcga_mdm4_rps18.R — MDM4 and RPS18 expression in TCGA-BRCA, female, White.
## Same 836 samples as the DESeq2 / VIPER analyses: primary tumours (PAM50 basal-like, luminal A, luminal B,
## HER2-enriched) and normal-adjacent tissue (NAT). Expression = DESeq2 VST (blind = FALSE) from the cached
## results/tf_activity/viper_vst/steps/vst_836_blindFALSE.rds; DESeq2 log2FC / FDR (genome-wide BH) from the cached
## results/tcga/TCGA_BRCA_DESeq2_BasalRef_4contrasts_allgenes.rds -- no refit.
## Outputs (results/target_genes/): data/TCGA_MDM4_RPS18_expression_by_subtype.csv, data/TCGA_MDM4_RPS18_DESeq2.csv,
##   figures/tcga/TCGA_MDM4_RPS18_boxplot_by_subtype.png (PNG only).
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
setwd(TNBC_ROOT)
ENS <- c(MDM4 = "ENSG00000198625", RPS18 = "ENSG00000231500")
V <- readRDS("results/tf_activity/viper_vst/steps/vst_836_blindFALSE.rds")$vst
meta <- readRDS("results/tcga/TCGA_BRCA_PAM50_metadata.rds")
rid <- sub("\\..*$", "", rownames(V))
stopifnot(all(ENS %in% rid), all(colnames(V) %in% rownames(meta)))
grp <- as.character(meta[colnames(V), "paper_BRCA_Subtype_PAM50"])
cat("samples by group:\n"); print(table(grp))
D <- rbindlist(lapply(names(ENS), function(g) data.table(gene = g, sample = colnames(V), group = grp,
                                                          vst = V[match(ENS[[g]], rid), ])))
LV <- c("Basal", "LumA", "LumB", "Her2", "NAT")
LAB <- c(Basal = "Basal-like", LumA = "Luminal A", LumB = "Luminal B", Her2 = "HER2-
enriched", NAT = "NAT")
NN <- table(grp)
D[, group := factor(group, levels = LV)]
S <- D[, .(n = .N, mean_vst = round(mean(vst), 3), median_vst = round(median(vst), 3)), by = .(gene, group)][order(gene, group)]
fwrite(S, "results/target_genes/TCGA_MDM4_RPS18_expression_by_subtype.csv"); print(S)
R <- readRDS("results/tcga/TCGA_BRCA_DESeq2_BasalRef_4contrasts_allgenes.rds")
R[, ensg := sub("\\..*$", "", Ensembl)]
R <- R[ensg %in% ENS][, gene := names(ENS)[match(ensg, ENS)]]
R <- R[, .(gene, Contrast, log2FC = round(log2FoldChange, 3), lfcSE = round(lfcSE, 3), FDR = signif(padj, 3))][order(gene, Contrast)]
fwrite(R, "results/target_genes/TCGA_MDM4_RPS18_DESeq2.csv"); print(R)
star <- R[, .(gene, group = sub("Basal vs ", "", Contrast),
              lab = ifelse(FDR < 0.001, "***", ifelse(FDR < 0.01, "**", ifelse(FDR < 0.05, "*", "ns"))))]
ymax <- D[, .(y = max(vst)), by = gene]; star <- merge(star, ymax, by = "gene")
star[, group := factor(group, levels = LV)]
p <- ggplot(D, aes(group, vst, fill = group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.6, linewidth = 0.3) +
  geom_text(data = star, aes(group, y + 0.15, label = lab), inherit.aes = FALSE, size = 3.2) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_x_discrete(labels = function(x) sprintf("%s\n(n = %s)", LAB[x], NN[x])) +
  scale_fill_manual(values = c(Basal = "#E64B35", LumA = "#4DBBD5", LumB = "#00A087", Her2 = "#F39B7F", NAT = "grey70"),
                    guide = "none") +
  labs(x = NULL, y = "Expression (DESeq2 VST)", title = "TCGA-BRCA, female, White",
       subtitle = "Basal-like vs each group: DESeq2 FDR * < 0.05, ** < 0.01, *** < 0.001, ns not significant") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold.italic", size = 11),
        axis.text.x = element_text(colour = "black", size = 8), plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8))
dir.create("results/target_genes/figures/tcga", showWarnings = FALSE, recursive = TRUE)
ggsave("results/target_genes/figures/tcga/TCGA_MDM4_RPS18_boxplot_by_subtype.png", p, width = 8, height = 4, dpi = 300, bg = "white")
cat("wrote results/target_genes/figures/tcga/TCGA_MDM4_RPS18_boxplot_by_subtype.png\n")
