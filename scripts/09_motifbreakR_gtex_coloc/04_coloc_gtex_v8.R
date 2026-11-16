# ==============================================================================
# 04_coloc_gtex_v8.R
# Colocalisation of TNBC GWAS with GTEx v8 breast cis-eQTL for promoter genes.
# Produces: Table S13
# Inputs  : results/gwas/bcac_gwas_ldblock.csv
#           reference/promoter_snp_genes.txt
#           results/motifbreakR/motifbreakR_TNBC_all.csv
#           results/tfgwas/tfgwas_significant_TNBC_unique.csv
# Outputs : results/coloc/coloc_GTEx_breast_bulk_bylocus.csv
# Run     : Rscript scripts/09_motifbreakR_gtex_coloc/04_coloc_gtex_v8.R   (from the repository root)
# ==============================================================================
## GTEx v8 breast BULK eQTL (eQTL Catalogue QTD000211, hg38) vs TN-GWAS coloc, by locus, windows 150/300/500 kb
## ---- paths (set the environment variables; see README) ----
BCAC_DIR <- Sys.getenv("BCAC_DIR", "data/BCAC")
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
LIFTOVER_DIR <- Sys.getenv("LIFTOVER_DIR", "data/liftover")
suppressMessages({library(data.table); library(coloc); library(rtracklayer); library(GenomicRanges)})
setwd(TNBC_ROOT)
GWAS_N<-16883L+96201L; GWAS_s<-16883L/GWAS_N; P12<-c(1e-5,5e-5); WIN<-c(300000L)
BASE<-"https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000211/QTD000211.all.tsv.gz"
chrn<-function(x) as.integer(sub("^chr","",as.character(x)))

## --- GWAS hg19 -> hg38 ---
GW<-file.path(BCAC_DIR, "CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt")
gw<-fread(GW)[,.(chr=suppressWarnings(as.integer(CHR)),pos19=suppressWarnings(as.integer(position)),A1=toupper(Allele1),A2=toupper(Allele2),bg=Effect,sg=StdErr)]
xw<-unique(fread("results/gwas/bcac_gwas_ldblock.csv")[,.(chr=chrn(Chr),pos19=as.integer(Pos),pos38=as.integer(pos38))])
gwh<-merge(gw,xw,by=c("chr","pos19")); gwh<-gwh[!duplicated(paste(chr,pos38))]

## --- gene bodies hg19 -> hg38, with ENSG + locus + anchor label ---
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
## locus id: same chr within 2Mb
setorder(gb,chr,s38); gb[,locus:=cumsum(c(1,as.integer(diff(s38)>2e6 | diff(chr)!=0)))]
locname<-c("1"="1q32","2"="2p23","3"="2q14","4"="5p15","5"="6p21_HLA","6"="11q22","7"="19p13.11")

## --- one tabix call, all gene windows (max window) ---
reg<-gb[,sprintf("%d:%d-%d",chr,pmax(1,s38-max(WIN)),e38+max(WIN))]
tf<-tempfile(fileext=".tsv")
cmd<-paste("tabix", shQuote(BASE), paste(reg,collapse=" "), ">", shQuote(tf))
cat("querying GTEx (",length(reg),"regions ) ...\n"); system(cmd)
cn<-c("mtid","chr","pos","ref","alt","variant","ma_samples","maf","pvalue","beta","se","type","ac","an","r2","mtoid","gene_id","median_tpm","rsid")
E<-fread(tf,header=FALSE,col.names=cn); E[,chr:=chrn(chr)]
cat("GTEx rows:",nrow(E)," unique genes:",length(unique(E$mtid)),"\n")

run1<-function(gg,ensg,p12,ichr,rs,re){
  e<-E[mtid==ensg & chr==ichr & pos>=rs & pos<=re & !is.na(beta) & se>0][!duplicated(pos)]; if(nrow(e)<20)return(NULL)
  m<-merge(e,gwh[chr==ichr],by.x="pos",by.y="pos38"); if(nrow(m)<20)return(NULL)
  A1<-toupper(m$A1);A2<-toupper(m$A2);R<-toupper(m$ref);A<-toupper(m$alt);st<-rep("x",nrow(m))
  st[A1==A&A2==R]<-"m"; st[A1==R&A2==A]<-"s"
  m[,bgh:=ifelse(st=="s",-bg,bg)]; m<-m[st%in%c("m","s")][is.finite(beta)&se>0&is.finite(bgh)&sg>0]; if(nrow(m)<20)return(NULL)
  de<-list(beta=m$beta,varbeta=m$se^2,type="quant",N=394L,MAF=pmin(pmax(m$maf,1e-4),0.9999),snp=as.character(m$pos))
  dg<-list(beta=m$bgh,varbeta=m$sg^2,type="cc",N=GWAS_N,s=GWAS_s,snp=as.character(m$pos))
  s<-tryCatch(suppressMessages(coloc.abf(de,dg,p1=1e-4,p2=1e-4,p12=p12))$summary,error=function(e)NULL);if(is.null(s))return(NULL)
  data.table(gene=gg,PP0=round(s["PP.H0.abf"],3),PP1=round(s["PP.H1.abf"],3),PP2=round(s["PP.H2.abf"],3),PP3=round(s["PP.H3.abf"],3),PP4=round(s["PP.H4.abf"],3),nsnps=as.integer(s["nsnps"]),p12=p12,win=re-rs)}

out<-list()
for(i in 1:nrow(gb)) for(w in WIN) for(p in P12){
  r<-run1(gb$gene[i],gb$ensg[i],p,gb$chr[i],gb$s38[i]-w,gb$e38[i]+w)
  if(!is.null(r)){ r[,`:=`(locus=locname[as.character(gb$locus[i])],anchor=gb$anchor[i],flank_kb=w/1000L)]; out[[length(out)+1]]<-r }
}
R<-rbindlist(out,fill=TRUE)
fwrite(R,"results/coloc/coloc_GTEx_breast_bulk_bylocus.csv")

for(pp in P12){
  cat(sprintf("\n=== GTEx breast bulk-eQTL coloc | flank=300kb, p12=%g, by locus (anchor genes *) ===\n",pp))
  pr<-R[p12==pp][order(locus,-PP4)]
  pr[,g:=ifelse(anchor=="anchor",paste0(gene,"*"),gene)]
  print(pr[,.(locus,gene=g,nsnps,PP0,PP1,PP2,PP3,PP4)],nrow=Inf)
}
cat("\nnot testable:",paste(setdiff(ens$gene,unique(R$gene)),collapse=", "),"\n")
cat("GTEX_COLOC_DONE\n")
