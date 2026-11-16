# ==============================================================================
# 08_gtex_v11_genotype_violins.R
# Gene expression by genotype, GTEx v11 breast tissue.
# Produces: Fig 5D-H
# Inputs  : results/motifbreakR/motifbreakR_TNBC_all.csv
#           results/tfgwas/tfgwas_significant_TNBC_unique.csv
#           reference/promoter_snp_genes.txt
#           results/gwas/bcac_gwas_ldblock.csv
# Run     : Rscript scripts/09_motifbreakR_gtex_coloc/08_gtex_v11_genotype_violins.R   (from the repository root)
# ==============================================================================
## 08_gtex_v11_genotype_violins.R — Fig 5D-H: gene expression by genotype (GTEx-portal-style
## violins) for promoter SNPs that alter TNBC-associated motifs.
##
## Self-contained: the SNP set can be switched with MOTIF_SET.
##
## SNP SET (env var MOTIF_SET):
##   "specific"    (default) promoter SNPs altering a TNBC-SPECIFIC motif (57-motif set):
##                 18 SNPs / 20 motifs / 10 genes -- the set used for the Fig 5D-H violins
##   "significant" promoter SNPs altering ANY TNBC Bonferroni-significant motif:
##                 57 SNPs / 106 motifs / 20 genes (adds RPS18, MRPL34, POGLUT3, WDR43, ...)
##
## DATA (dbGaP-controlled, stays on the HPC cluster):
##   expression : GTEx v11 Breast_Mammary_Tissue normalized_expression.bed.gz (511 donors)
##   covariates : GTEx v11 Breast_Mammary_Tissue covariates.txt (68 covariates)
##   genotypes  : GTEx v9 WGS VCF, hg38 (all 511 v11 breast donors present)
## MODEL: expr ~ dosage + 68 covariates (what GTEx fits); the violin shows GTEx normalised
## expression by genotype and the caption carries that model's p-value.
## NOT TESTABLE: TERT (below GTEx expression threshold in breast), SNORD53 / SNORD92 (snoRNAs,
## not quantified); RALB rs4849843 is absent from the WGS VCF.
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
GTEX_DIR <- Sys.getenv("GTEX_DIR", "data/GTEx")
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
setwd(TNBC_ROOT)
GT   <- GTEX_DIR
IN   <- "results/coloc/"                  # upstream results (read only)
OUT  <- "results/target_genes/"; FIG <- "results/target_genes/figures/violin_v11"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
SET  <- Sys.getenv("MOTIF_SET", "specific"); stopifnot(SET %in% c("specific", "significant"))
BED  <- file.path(GT, "v11/GTEx_Analysis_v11_eQTL_expression_matrices/Breast_Mammary_Tissue.v11.normalized_expression.bed.gz")
CVF  <- file.path(GT, "v11/GTEx_Analysis_v11_eQTL_covariates/Breast_Mammary_Tissue.v11.covariates.txt")
VCF  <- file.path(GT, "GTEx_Analysis_2021-02-11_v9_WGS_VCF_files_WholeGenomeSeq_944Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz")
chrn <- function(x) as.integer(sub("^chr", "", as.character(x)))

## ---- style ----
## one colour per genotype CLASS (hom. reference / heterozygous / hom. alternative), so a colour
## means the same allele dose in every panel whatever the letters are (CC/CT/TT, GG/GA/AA, ...)
FILL3 <- c("#D6D8F2", "#9FA4DC", "#5E64BE")                     # purple-blue, light -> dark with dose
EDGE3 <- c("#B3B7E3", "#7E84CC", "#43489E")
BOX <- "#2F3270"
PCOL <- "black"
fmtp <- function(p) if (p < 0.001) sprintf("p value = %.2e", p) else sprintf("p value = %.3f", p)
PANEL_W <- 3.7; PANEL_H <- 5.2

## ---- SNP-gene pairs from motifbreakR (hg19), hg19 -> hg38 via the BCAC crosswalk ----
mb  <- fread("results/motifbreakR/motifbreakR_TNBC_all.csv")
u57 <- unique(fread("results/tfgwas/tfgwas_significant_TNBC_unique.csv")$motif)
p   <- mb[gsub('"', '', TNBC_sig) == "Y" & grepl("Promoter", annotation_category)]
if (SET == "specific") p <- p[providerName %in% u57]
## Intronic snoRNAs excluded: SNORD53 / SNORD92 lie inside WDR43, so their
## "promoter" SNPs are promoter only relative to the snoRNA start (intronic to WDR43), and snoRNAs
## are not reliably quantified by poly-A RNA-seq (GTEx: absent; TCGA: baseMean < 1).
SNORNA <- c("SNORD53", "SNORD92")
p <- p[!SYMBOL %in% SNORNA]
## 1q32 (MDM4 locus) restricted to the selected SNPs: the 7 MDM4 promoter SNPs;
## the PIK3C2B promoter SNPs at the same locus are not shown.
KEEP_1Q32 <- c("rs4245736", "rs4245737", "rs4252685", "rs4252687", "rs4951078", "rs4951391", "rs61817958")
p <- p[!SYMBOL %in% c("MDM4", "PIK3C2B") | SNP_id %in% KEEP_1Q32]
S <- unique(p[, .(SNP_id, gene = SYMBOL, chr = chrn(geneChr), pos19 = as.integer(Pos37),
                  distanceToTSS = as.integer(distanceToTSS))])
ens <- fread("reference/promoter_snp_genes.txt", header = FALSE, col.names = c("ensg", "gene"))
S <- merge(S, ens, by = "gene")
xw <- unique(fread("results/gwas/bcac_gwas_ldblock.csv")[, .(chr = chrn(Chr), pos19 = as.integer(Pos), pos38 = as.integer(pos38))])
S <- merge(S, xw, by = c("chr", "pos19"), all.x = TRUE)[!is.na(pos38)]
cat(sprintf("MOTIF_SET = %s | %d SNPs, %d motifs, %d genes | SNP-gene pairs with hg38: %d\n",
            SET, uniqueN(p$SNP_id), uniqueN(p$providerName), uniqueN(p$SYMBOL), nrow(S)))

## optional annotation from script 31 (QTD000211 eQTL, lead SNP, LD); NA where 31 was not run
f31 <- file.path(IN, "TNBCuniq_promoterSNPs_eQTL_and_LD.csv")
if (file.exists(f31)) S <- merge(S, fread(f31)[, .(gene, SNP_id, beta_QTD394 = beta, p_QTD394 = pvalue,
                                                   lead_rsid, r2_to_lead = r2_to_lead_eQTL)],
                                 by = c("gene", "SNP_id"), all.x = TRUE)
## published v11 reference (table only; not drawn on the figure)
v11g <- fread(file.path(IN, "GTExv11_breast_eGene_status_20genes.csv"))[, .(gene, v11_qval = qval, v11_eGene)]
v11p <- fread(file.path(IN, "GTExv11_breast_signifpairs_ourSNPs.csv"))[our_gene == TRUE,
           .(gene = promoter_gene, SNP_id, v11_pair_p = pval_nominal, v11_pair_slope = slope)]
v11p <- v11p[order(v11_pair_p)][!duplicated(paste(gene, SNP_id))]
S <- merge(merge(S, v11g, by = "gene", all.x = TRUE), v11p, by = c("gene", "SNP_id"), all.x = TRUE)

## ---- expression + covariates ----
E <- fread(cmd = paste("zcat", shQuote(BED))); setnames(E, 1:4, c("chr", "start", "end", "gene_id"))
don <- names(E)[-(1:4)]
CV <- fread(CVF); CM <- t(as.matrix(CV[, ..don])); colnames(CM) <- CV$ID
cat(sprintf("GTEx v11 breast donors: %d | covariates: %d\n", length(don), ncol(CM)))

## ---- genotypes (cached per SNP set) ----
gf <- file.path(OUT, sprintf("promoterSNPs_%s_v11.gt.tsv.gz", SET))
if (!file.exists(gf)) {
  rf <- tempfile(); writeLines(sprintf("chr%d\t%d", S$chr, S$pos38), rf)
  dl <- tempfile(); writeLines(don, dl)
  system(sprintf(paste("{ printf 'ID\\tCHROM\\tPOS\\tREF\\tALT'; printf '\\t%%s' $(cat %s); printf '\\n';",
                       "bcftools view -S %s -R %s -v snps,indels -m2 -M2 -Ou %s",
                       "| bcftools query -f '%%ID\\t%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n'; } | gzip > %s"),
                 shQuote(dl), shQuote(dl), shQuote(rf), shQuote(VCF), shQuote(gf)))
}
G <- fread(gf)
gm <- as.matrix(G[, -(1:5)])
dos <- matrix(as.integer(substr(gm, 1, 1) == "1") + as.integer(substr(gm, 3, 3) == "1"), nrow(gm))
colnames(dos) <- colnames(gm)
miss <- setdiff(don, colnames(dos)); if (length(miss)) cat("donors without genotypes (dropped):", length(miss), "\n")
don <- intersect(don, colnames(dos)); dos <- dos[, don, drop = FALSE]; CM <- CM[don, , drop = FALSE]
G[, pos38 := as.integer(POS)]

## ---- per pair: fit, plot ----
stats <- list(); plots <- list()
for (i in seq_len(nrow(S))) {
  gi  <- which(G$pos38 == S$pos38[i] & chrn(G$CHROM) == S$chr[i])
  row <- E[grepl(paste0("^", S$ensg[i]), gene_id)]
  if (!length(gi) || !nrow(row)) {
    cat(sprintf("%-8s %-12s skipped (%s)\n", S$gene[i], S$SNP_id[i],
                if (!length(gi)) "not in VCF" else "gene not in GTEx expression matrix")); next }
  gi <- gi[1]
  d <- dos[gi, ]; expr <- setNames(as.numeric(row[1, ..don]), don)
  ok <- !is.na(d) & !is.na(expr)
  ref <- G$REF[gi]; alt <- G$ALT[gi]
  lab <- c(paste0(ref, ref), paste0(ref, alt), paste0(alt, alt))
  fit <- summary(lm(expr[ok] ~ d[ok] + CM[ok, , drop = FALSE]))
  b <- coef(fit)[2, 1]; pv <- coef(fit)[2, 4]
  DT <- data.table(dosage = d[ok], expr = expr[ok])
  DT[, gt := factor(lab[dosage + 1], levels = lab)]
  n_by <- DT[, .N, by = gt][match(lab, gt)]$N; n_by[is.na(n_by)] <- 0
  DT[, gt_n := factor(sprintf("%s\n(%d)", gt, n_by[as.integer(gt)]), levels = sprintf("%s\n(%d)", lab, n_by))]
  stats[[length(stats) + 1]] <- cbind(S[i, .(gene, SNP_id, distanceToTSS)],
      data.table(ref, alt, n = sum(ok), n_RR = n_by[1], n_RA = n_by[2], n_AA = n_by[3], beta_v11 = b, p_v11 = pv),
      S[i, intersect(c("beta_QTD394", "p_QTD394", "lead_rsid", "r2_to_lead", "v11_eGene", "v11_qval",
                       "v11_pair_p", "v11_pair_slope"), names(S)), with = FALSE])
  g <- ggplot(DT, aes(gt_n, expr)) +
    geom_violin(aes(fill = gt_n, colour = gt_n), trim = FALSE, linewidth = 0.4, scale = "width", width = 0.9) +
    scale_fill_manual(values = FILL3, drop = FALSE, guide = "none") +
    scale_colour_manual(values = EDGE3, drop = FALSE, guide = "none") +
    geom_boxplot(width = 0.22, outlier.shape = NA, fill = BOX, colour = BOX, linewidth = 0.3, fatten = NULL) +
    stat_summary(fun = median, geom = "crossbar", width = 0.22, colour = "white", linewidth = 0.4) +
    labs(x = NULL, y = "Norm. expression", title = sprintf("%s ~ %s (%s>%s)", S$gene[i], S$SNP_id[i], ref, alt),
         caption = fmtp(pv)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 10, face = "bold", colour = "black", hjust = 0.5),
          plot.caption = element_text(size = 9, colour = PCOL, hjust = 0.5),
          axis.title.y = element_text(colour = "black"), axis.text = element_text(colour = "black", size = 10),
          axis.line = element_line(colour = "grey60"), axis.ticks = element_line(colour = "grey60"))
  plots[[paste(S$gene[i], S$SNP_id[i])]] <- g
  ggsave(file.path(FIG, sprintf("%s_%s.png", S$SNP_id[i], S$gene[i])), g,          # one file per SNP: <SNP>_<gene>.png
         width = PANEL_W, height = PANEL_H, dpi = 300, bg = "white")
  cat(sprintf("%-8s %-12s %s>%s n=%d (%d/%d/%d) beta=%+.3f P=%.3g\n",
              S$gene[i], S$SNP_id[i], ref, alt, sum(ok), n_by[1], n_by[2], n_by[3], b, pv))
}
ST <- rbindlist(stats, fill = TRUE); setorder(ST, p_v11)
fwrite(ST, file.path(OUT, sprintf("violin_v11_genotype_expression_%s.csv", SET)))
ord <- ST[, paste(gene, SNP_id)]
ggsave(file.path(FIG, sprintf("violin_v11_ALL_%s.png", SET)), wrap_plots(plots[ord], ncol = 4),
       width = PANEL_W * 4, height = PANEL_H * ceiling(length(ord) / 4), dpi = 200, bg = "white", limitsize = FALSE)

## ---- combined one-row figures, panels in the listed order ----
ROWS <- list(
  ## rs61817958 removed from row 1 (: not in the GTEx v8 eQTL Catalogue file used for coloc)
  MDM4_row1 = c("MDM4 rs4951391", "MDM4 rs4951078", "MDM4 rs4252685", "MDM4 rs4252687"),
  MDM4_row2 = c("MDM4 rs4245737", "MDM4 rs4245736"),
  MDM4_row3 = c("MDM4 rs4951391", "MDM4 rs4252685", "MDM4 rs4951078"))   #; rs61817958 dropped (same reason as row 1)
## per-row figure height (inches); default PANEL_H. Rows 1 and 3 shorter.
ROW_H <- c(MDM4_row1 = 4.0, MDM4_row3 = 4.0)
for (nm in names(ROWS)) {
  k <- ROWS[[nm]]; have <- k[k %in% names(plots)]
  if (length(miss <- setdiff(k, have))) cat("  ", nm, ": not available in this run:", paste(miss, collapse = ", "), "\n")
  if (!length(have)) next
  f <- file.path(FIG, sprintf("violin_v11_%s_%s.png", nm, paste(sub("^\\S+ ", "", have), collapse = "_")))
  ggsave(f, wrap_plots(plots[have], nrow = 1), width = PANEL_W * length(have),
         height = if (nm %in% names(ROW_H)) ROW_H[[nm]] else PANEL_H, dpi = 300, bg = "white")
  cat("  wrote", f, "\n")
}
cat("\n"); print(as.data.frame(ST[, .(gene, SNP_id, alleles = paste0(ref, ">", alt), n_RR, n_RA, n_AA,
                                      beta_v11 = round(beta_v11, 3), p_v11 = signif(p_v11, 3))]), row.names = FALSE)
cat("\nwrote", file.path(OUT, sprintf("violin_v11_genotype_expression_%s.csv", SET)), "and", length(plots), "figures in", FIG, "\n")
