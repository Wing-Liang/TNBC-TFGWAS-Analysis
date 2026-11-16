# ==============================================================================
# 05_coloc_harmonization_audit.R
# Audit allele harmonisation between the GWAS and eQTL data.
# Produces: Table S13 (check)
# Inputs  : results/gwas/bcac_gwas_ldblock.csv
# Run     : Rscript scripts/09_motifbreakR_gtex_coloc/05_coloc_harmonization_audit.R   (from the repository root)
# ==============================================================================
# 05_coloc_harmonization_audit.R — how many variants does the +/-300 kb coloc DROP at the
# allele-matching step, and are the dropped ones strand flips or genuine mismatches?
## ---- paths (set the environment variables; see README) ----
BCAC_DIR <- Sys.getenv("BCAC_DIR", "data/BCAC")
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({library(data.table)})
setwd(TNBC_ROOT)
chrn <- function(x) as.integer(sub("^chr","",x))
GW <- file.path(BCAC_DIR, "CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt")
gw <- fread(GW)[,.(chr=suppressWarnings(as.integer(CHR)),pos19=suppressWarnings(as.integer(position)),
                   A1=toupper(Allele1),A2=toupper(Allele2),bg=Effect,sg=StdErr)]
xw <- unique(fread("results/gwas/bcac_gwas_ldblock.csv")[,.(chr=chrn(Chr),pos19=as.integer(Pos),pos38=as.integer(pos38))])
gwh <- merge(gw,xw,by=c("chr","pos19")); gwh <- gwh[!duplicated(paste(chr,pos38))]

comp <- c(A="T",T="A",C="G",G="C")
# ABHD8 gene body hg38 chr19:17,280,000-17,300,000 -> +/-300kb window
REG <- "19:16980000-17600000"
tf <- tempfile()
system(paste("tabix", shQuote("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000211/QTD000211.all.tsv.gz"), REG, ">", shQuote(tf)))
cn <- c("mtid","chr","pos","ref","alt","variant","ma_samples","maf","pvalue","beta","se","type","ac","an","r2","mtoid","gene_id","median_tpm","rsid")
E <- fread(tf, header=FALSE, col.names=cn); E[,chr:=chrn(chr)]

for (g in c("ENSG00000131013","ENSG00000167842","ENSG00000104881")) {   # ABHD8, ANKLE1-ish, OCEL1-ish
  e <- E[mtid==g & !is.na(beta) & se>0][!duplicated(pos)]
  if (!nrow(e)) next
  m <- merge(e, gwh[chr==19], by.x="pos", by.y="pos38")
  if (!nrow(m)) next
  A1<-toupper(m$A1);A2<-toupper(m$A2);R<-toupper(m$ref);A<-toupper(m$alt)
  st <- rep("x", nrow(m)); st[A1==A&A2==R]<-"m"; st[A1==R&A2==A]<-"s"
  # of the dropped, how many are STRAND FLIPS (complement matches)?
  fl <- st=="x" & (comp[A1]==A & comp[A2]==R); fl[is.na(fl)] <- FALSE
  fs <- st=="x" & (comp[A1]==R & comp[A2]==A); fs[is.na(fs)] <- FALSE
  amb <- st=="x" & ((A1==comp[A2]) | (R==comp[A]))    # A/T or C/G strand-ambiguous
  amb[is.na(amb)] <- FALSE
  cat(sprintf("\n%s  merged=%d\n", g, nrow(m)))
  cat(sprintf("   kept  match(m)=%d  swap(s)=%d  -> used %d (%.1f%%)\n",
      sum(st=="m"), sum(st=="s"), sum(st %in% c("m","s")), 100*mean(st %in% c("m","s"))))
  cat(sprintf("   DROPPED x=%d (%.1f%%):  strand-flip=%d  strand-flip+swap=%d  strand-ambiguous=%d  other=%d\n",
      sum(st=="x"), 100*mean(st=="x"), sum(fl), sum(fs), sum(amb), sum(st=="x")-sum(fl)-sum(fs)))
}
