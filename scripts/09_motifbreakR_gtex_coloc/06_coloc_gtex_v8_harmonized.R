# ==============================================================================
# 06_coloc_gtex_v8_harmonized.R
# Colocalisation with strand-harmonised alleles (reported analysis).
# Produces: Table S13
# Inputs  : results/gwas/bcac_gwas_ldblock.csv
#           reference/promoter_snp_genes.txt
#           results/motifbreakR/motifbreakR_TNBC_all.csv
#           results/tfgwas/tfgwas_significant_TNBC_unique.csv
# Outputs : results/coloc/coloc_GTEx_bulk_strandharmonized_300kb.csv
#           results/coloc/harmonization_audit_300kb.csv
# Run     : Rscript scripts/09_motifbreakR_gtex_coloc/06_coloc_gtex_v8_harmonized.R   (from the repository root)
# ==============================================================================
## 06_coloc_gtex_v8_harmonized.R — GTEx v8 breast BULK eQTL vs TN-GWAS coloc, +/-300 kb,
## with STRAND-AWARE allele harmonisation. Identical to 04_coloc_gtex_v8.R except the allele block.
##
## 04_coloc_gtex_v8.R resolves only allele-ORDER swaps:
##     m: GWAS A1==alt & A2==ref  -> keep beta
##     s: GWAS A1==ref & A2==alt  -> flip beta
##   everything else dropped -- including STRAND FLIPS (GWAS A/G vs GTEx T/C, the same
##   variant reported on the opposite strand), which are recoverable.
## Here strand flips are resolved too. Palindromic (A/T, C/G) variants that do not match
## directly are still dropped: their strand cannot be determined from alleles alone.
## ---- paths (set the environment variables; see README) ----
BCAC_DIR <- Sys.getenv("BCAC_DIR", "data/BCAC")
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
GTEX_DIR <- Sys.getenv("GTEX_DIR", "data/GTEx")
LIFTOVER_DIR <- Sys.getenv("LIFTOVER_DIR", "data/liftover")
suppressMessages({library(data.table); library(coloc); library(rtracklayer); library(GenomicRanges)})
setwd(TNBC_ROOT)
GWAS_N<-16883L+96201L; GWAS_s<-16883L/GWAS_N; P12<-c(1e-5,5e-5); WIN<-c(300000L)
GZ<-file.path(GTEX_DIR, "eQTL_v8/QTD000211.all.tsv.gz")
chrn<-function(x) as.integer(sub("^chr","",as.character(x)))
comp<-c(A="T",T="A",C="G",G="C")

GW<-file.path(BCAC_DIR, "CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt")
gw<-fread(GW)[,.(chr=suppressWarnings(as.integer(CHR)),pos19=suppressWarnings(as.integer(position)),A1=toupper(Allele1),A2=toupper(Allele2),bg=Effect,sg=StdErr)]
xw<-unique(fread("results/gwas/bcac_gwas_ldblock.csv")[,.(chr=chrn(Chr),pos19=as.integer(Pos),pos38=as.integer(pos38))])
gwh<-merge(gw,xw,by=c("chr","pos19")); gwh<-gwh[!duplicated(paste(chr,pos38))]

ens<-fread("reference/promoter_snp_genes.txt",header=FALSE,col.names=c("ensg","gene"))
mb<-fread("results/motifbreakR/motifbreakR_TNBC_all.csv")
u57<-fread("results/tfgwas/tfgwas_significant_TNBC_unique.csv"); ma<-unique(u57$motif)
anchor<-unique(mb[providerName %in% ma & gsub('"','',TNBC_sig)=="Y" & grepl("Promoter",annotation_category)]$SYMBOL)
gpos<-unique(mb[SYMBOL %in% ens$gene,.(gene=SYMBOL,chr=chrn(geneChr),gstart=as.integer(geneStart),gend=as.integer(geneEnd))])[!duplicated(gene)]
gb<-merge(gpos,ens,by="gene")
ch<-import.chain(file.path(LIFTOVER_DIR, "hg19ToHg38.over.chain"))
g19<-GRanges(paste0("chr",gb$chr),IRanges(gb$gstart,gb$gend)); lo<-liftOver(g19,ch)
g38<-t(sapply(lo,function(x) if(length(x))c(min(start(x)),max(end(x))) else c(NA,NA))); gb$s38<-g38[,1]; gb$e38<-g38[,2]; gb<-gb[!is.na(s38)]
gb[,anchor:=ifelse(gene %in% anchor,"anchor","neighbor")]
setorder(gb,chr,s38); gb[,locus:=cumsum(c(1,as.integer(diff(s38)>2e6 | diff(chr)!=0)))]
locname<-c("1"="1q32","2"="2p23","3"="2q14","4"="5p15","5"="6p21_HLA","6"="11q22","7"="19p13.11")
cn<-c("mtid","chr","pos","ref","alt","variant","ma_samples","maf","pvalue","beta","se","type","ac","an","r2","mtoid","gene_id","median_tpm","rsid")

aud<-list()
run1<-function(gg,ensg,p12,ichr,rs,re){
  tf<-tempfile(); system(paste("tabix",shQuote(GZ),sprintf("%d:%d-%d",ichr,rs,re),">",shQuote(tf)))
  if(!file.exists(tf)||file.info(tf)$size==0) return(NULL)
  E<-fread(tf,header=FALSE,col.names=cn); E[,chr:=chrn(chr)]
  e<-E[mtid==ensg & chr==ichr & pos>=rs & pos<=re & !is.na(beta) & se>0][!duplicated(pos)]
  if(nrow(e)<20)return(NULL)
  m<-merge(e,gwh[chr==ichr],by.x="pos",by.y="pos38"); if(nrow(m)<20)return(NULL)
  A1<-toupper(m$A1);A2<-toupper(m$A2);R<-toupper(m$ref);A<-toupper(m$alt)
  st<-rep("x",nrow(m))
  st[A1==A&A2==R]<-"m"; st[A1==R&A2==A]<-"s"                      # direct
  palin<-!is.na(comp[A1]) & !is.na(A2) & comp[A1]==A2             # A/T or C/G: strand unknowable
  f1<-comp[A1]; f2<-comp[A2]
  ok<-st=="x" & !palin & !is.na(f1) & !is.na(f2)
  st[ok & f1==A & f2==R]<-"mf"; st[ok & f1==R & f2==A]<-"sf"      # strand-flipped
  if(p12==P12[1]) aud[[gg]]<<-data.table(gene=gg,merged=nrow(m),
      match=sum(st=="m"),swap=sum(st=="s"),flip=sum(st=="mf"),flip_swap=sum(st=="sf"),
      palindromic_dropped=sum(st=="x"&palin),other_dropped=sum(st=="x"&!palin),
      used_old=sum(st%in%c("m","s")),used_new=sum(st%in%c("m","s","mf","sf")))
  m[,bgh:=ifelse(st %in% c("s","sf"),-bg,bg)]
  m<-m[st%in%c("m","s","mf","sf")][is.finite(beta)&se>0&is.finite(bgh)&sg>0]; if(nrow(m)<20)return(NULL)
  de<-list(beta=m$beta,varbeta=m$se^2,type="quant",N=394L,MAF=pmin(pmax(m$maf,1e-4),0.9999),snp=as.character(m$pos))
  dg<-list(beta=m$bgh,varbeta=m$sg^2,type="cc",N=GWAS_N,s=GWAS_s,snp=as.character(m$pos))
  s<-tryCatch(suppressMessages(coloc.abf(de,dg,p1=1e-4,p2=1e-4,p12=p12))$summary,error=function(e)NULL)
  if(is.null(s))return(NULL)
  data.table(gene=gg,PP0=round(s["PP.H0.abf"],3),PP1=round(s["PP.H1.abf"],3),PP2=round(s["PP.H2.abf"],3),
             PP3=round(s["PP.H3.abf"],3),PP4=round(s["PP.H4.abf"],3),nsnps=as.integer(s["nsnps"]),p12=p12)}

out<-list()
for(i in 1:nrow(gb)) for(w in WIN) for(p in P12){
  r<-run1(gb$gene[i],gb$ensg[i],p,gb$chr[i],gb$s38[i]-w,gb$e38[i]+w)
  if(!is.null(r)){ r[,`:=`(locus=locname[as.character(gb$locus[i])],anchor=gb$anchor[i],flank_kb=w/1000L)]; out[[length(out)+1]]<-r }}
R<-rbindlist(out)
fwrite(R,"results/coloc/coloc_GTEx_bulk_strandharmonized_300kb.csv")
AU<-rbindlist(aud); fwrite(AU,"results/coloc/harmonization_audit_300kb.csv")
cat("\n=== HARMONISATION AUDIT (variants per gene) ===\n"); print(as.data.frame(AU),row.names=FALSE)
cat(sprintf("\n  total used OLD (swap only): %d | NEW (+strand flips): %d | recovered: %d (%.2f%%)\n",
    sum(AU$used_old), sum(AU$used_new), sum(AU$used_new)-sum(AU$used_old),
    100*(sum(AU$used_new)-sum(AU$used_old))/max(sum(AU$used_old),1)))
cat("\n=== POSTERIORS, strand-harmonised, p12 = 5e-5 ===\n")
print(as.data.frame(R[p12==5e-5][order(-PP4),.(gene,locus,anchor,nsnps,PP0,PP1,PP2,PP3,PP4)]),row.names=FALSE)
