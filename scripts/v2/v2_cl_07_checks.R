#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table)})
options(warn=1); set.seed(20260908)
REV<-"/path/to/revision"; DM<-file.path(REV,"data/depmap24q4")
OUT<-file.path(REV,"results/v2")
msg<-function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")
sp<-function(x,y){ok<-is.finite(x)&is.finite(y);n<-sum(ok)
 ct<-suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE))
 c(rho=unname(ct$estimate),p=ct$p.value,n=n)}
W<-readRDS(file.path(DM,"cl_workspace.rds")); E<-W$E
mod<-fread(file.path(DM,"Model.csv"))
mi<-readRDS(file.path(DM,"ccle_mirna.rds")); Mm<-mi$M; rownames(Mm)<-mi$desc; Mm<-log2(Mm+1)
key<-mod[CCLEName!="",.(CCLEName,ModelID,OncotreeLineage)]; setkey(key,CCLEName)
cn<-colnames(Mm); mid<-key[cn,ModelID]
br<-which(!is.na(mid)&mid%in%W$breast_ca&key[cn,OncotreeLineage]=="Breast"); ids<-mid[br]

## dynamic range of miR-29 in the 50 breast lines
dr<-rbindlist(lapply(c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","hsa-miR-200c","hsa-miR-21"),
  function(m){v<-Mm[m,br]; data.table(miRNA=m,n=length(v),mean_log2=mean(v),sd_log2=sd(v),
    min_log2=min(v),max_log2=max(v),IQR_log2=IQR(v),frac_at_floor=mean(2^v-1 < 5))}))
fwrite(dr,file.path(OUT,"celllines_mirna_dynamic_range.csv")); print(dr)

## miR-29 vs COL1A1 restricted to lines with real COL1A1 dynamic range
hi<-ids[(2^E[ids,"COL1A1"]-1)>=10]; hidx<-br[(2^E[ids,"COL1A1"]-1)>=10]
msg("breast lines with COL1A1 >= 10 TPM and miRNA data:",length(hi))
rr<-rbindlist(lapply(c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"),function(m){
  s<-sp(Mm[m,hidx],E[hi,"COL1A1"]); data.table(miRNA=m,target="COL1A1",n=s["n"],rho=s["rho"],p=s["p"])}))
rr[,q:=p.adjust(p,"BH")]; fwrite(rr,file.path(OUT,"celllines_miR29_COL1A1_expressed_subset.csv")); print(rr)

## N11 restatement on the full 53-line CRISPR set
ce<-fread("/path/to/home/Desktop/DD/R_GPR/ESIA_BRCA/brca_eisa_pilot/depmap/CRISPRGeneEffect.csv")
setnames(ce,1,"ModelID"); gc_<-setdiff(names(ce),"ModelID"); gs<-sub(" \\(\\d+\\)$","",gc_)
brc<-intersect(mod[OncotreeLineage=="Breast",ModelID],ce$ModelID)
sel<-c("NFKB1","RELA","SP1","ETS1","COL1A1","COL3A1","MYC","TP53")
sel<-sel[sel%in%gs]
M<-as.matrix(ce[ModelID%in%brc, gc_[match(sel,gs)],with=FALSE]); colnames(M)<-sel
es<-data.table(gene=sel,n_lines=colSums(!is.na(M)),mean_chronos=colMeans(M,na.rm=TRUE),
  n_essential_lt_minus0.5=colSums(M< -0.5,na.rm=TRUE))
fwrite(es,file.path(OUT,"celllines_crispr_essentiality_53lines.csv")); print(es)
msg("done")
