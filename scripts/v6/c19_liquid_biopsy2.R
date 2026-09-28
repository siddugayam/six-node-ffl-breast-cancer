## ==========================================================================
## c19_liquid_biopsy2.R -- Part C (proper version)
## GSE73002 serum: detection + case/control discrimination for the prioritised
## miRNAs, BENCHMARKED against all 2,540 assayed miRNAs and repeated after
## within-array rank normalisation, because the raw contrast carries a large
## global array-level shift.
## ==========================================================================
suppressPackageStartupMessages({library(data.table)})
setwd("/path/to/revision")
OUT<-"results/v6"; B<-"cache/v6/cohorts/built"
TARGETS <- c("miR-21","miR-195","miR-204","miR-383","miR-124","miR-155","miR-429",
             "miR-141","miR-34a","miR-101","miR-130a","miR-29a")
RX <- c("miR-21"="^hsa-mir-21(-[0-9])?(-[35]p)?$","miR-195"="^hsa-mir-195(-[35]p)?$",
        "miR-204"="^hsa-mir-204(-[35]p)?$","miR-383"="^hsa-mir-383(-[35]p)?$",
        "miR-124"="^hsa-mir-124a?(-[123])?(-[35]p)?$","miR-155"="^hsa-mir-155(-[35]p)?$",
        "miR-429"="^hsa-mir-429(-[35]p)?$","miR-141"="^hsa-mir-141(-[35]p)?$",
        "miR-34a"="^hsa-mir-34a(-[35]p)?$","miR-101"="^hsa-mir-101(-[12])?(-[35]p)?$",
        "miR-130a"="^hsa-mir-130a(-[35]p)?$","miR-29a"="^hsa-mir-29a(-[35]p)?$")
aucf <- function(x,g){ ok<-is.finite(x); x<-x[ok]; g<-g[ok]; r<-rank(x)
  n1<-sum(g==1); n0<-sum(g==0); if(n1<3||n0<3) return(NA_real_)
  (sum(r[g==1])-n1*(n1+1)/2)/(n1*n0) }

E<-fread(file.path(B,"GSE73002_expr.csv")); M<-as.matrix(E[,-1]); rownames(M)<-E[[1]]
P<-fread(file.path(B,"GSE73002_pheno.csv")); stopifnot(identical(colnames(M),P$gsm))
grp<-P$diagnosis
modev<-apply(M,2,function(v){v<-v[is.finite(v)];t<-table(round(v,4));as.numeric(names(t)[which.max(t)])})
Mr <- round(M,4)
DET <- sweep(Mr,2,round(modev,4),">") & is.finite(M)
cat("median miRNAs above the per-array floor:", median(colSums(DET,na.rm=TRUE)), "of", nrow(M), "\n")
## within-array rank normalisation (removes the global per-array shift)
R <- apply(M,2,function(v){ r<-rank(v,na.last="keep"); r/sum(is.finite(v)) })

## dominant probe per target
tgtprobe <- sapply(TARGETS, function(tg){
  rn<-rownames(M); nm<-gsub("\\*","",tolower(gsub("miR","mir",rn)))
  hit<-which(grepl(RX[[tg]],nm) & !grepl("\\*",rownames(M)))
  if(!length(hit)) return(NA_character_)
  rownames(M)[hit[which.max(rowMeans(M[hit,,drop=FALSE],na.rm=TRUE))]] })

res<-list()
for(cmp in list(c("breast cancer","non-cancer"),
                c("breast cancer","benign breast disease"),
                c("prostate disease","non-cancer"))){
  i<-which(grp %in% cmp); g<-as.integer(grp[i]==cmp[1])
  if(length(i)<40) next
  a_raw <- apply(M[,i,drop=FALSE],1,aucf,g=g)
  a_rank<- apply(R[,i,drop=FALSE],1,aucf,g=g)
  glob <- data.table(miRNA_probe=rownames(M), auc_raw=a_raw, auc_rank=a_rank)
  fwrite(glob, file.path(OUT, sprintf("cohorts_liquidbiopsy_GSE73002_allmiRNA_%s.csv",
        gsub("[^a-z]","",gsub(" ","",cmp[1])))))
  cat(sprintf("\n--- %s vs %s (n=%d/%d) ---\n", cmp[1],cmp[2],sum(g==1),sum(g==0)))
  cat(sprintf("  ALL %d assayed miRNAs: raw AUC median %.3f (IQR %.3f-%.3f), %.1f%% with AUC>0.8\n",
      nrow(M), median(a_raw,na.rm=TRUE), quantile(a_raw,.25,na.rm=TRUE),
      quantile(a_raw,.75,na.rm=TRUE), 100*mean(a_raw>0.8,na.rm=TRUE)))
  cat(sprintf("  after within-array rank normalisation: median %.3f (IQR %.3f-%.3f), %.1f%% with AUC>0.8\n",
      median(a_rank,na.rm=TRUE), quantile(a_rank,.25,na.rm=TRUE),
      quantile(a_rank,.75,na.rm=TRUE), 100*mean(a_rank>0.8,na.rm=TRUE)))
  for(tg in TARGETS){
    pr<-tgtprobe[[tg]]; if(is.na(pr)) next
    x<-M[pr,i]; xr<-R[pr,i]; d<-DET[pr,i]
    w<-suppressWarnings(wilcox.test(xr~g))
    res[[length(res)+1]]<-data.frame(cohort="GSE73002",accession="GSE73002",matrix="serum",
      platform="3D-Gene Human miRNA V20_1.0.0 (GPL18941)", comparison=paste(cmp,collapse=" vs "),
      miRNA=tg, probe=pr, n_case=sum(g==1), n_ctrl=sum(g==0),
      det_rate_case=mean(d[g==1],na.rm=TRUE), det_rate_ctrl=mean(d[g==0],na.rm=TRUE),
      auc_raw=aucf(x,g), auc_rank=aucf(xr,g), p_rank=w$p.value,
      pctile_auc_raw_among_all=100*mean(abs(a_raw-0.5) < abs(aucf(x,g)-0.5), na.rm=TRUE),
      pctile_auc_rank_among_all=100*mean(abs(a_rank-0.5) < abs(aucf(xr,g)-0.5), na.rm=TRUE),
      global_median_auc_raw=median(a_raw,na.rm=TRUE),
      global_median_auc_rank=median(a_rank,na.rm=TRUE), stringsAsFactors=FALSE)
  }
}
RES<-rbindlist(res); RES[, q_rank:=p.adjust(p_rank,"BH"), by=comparison]
fwrite(RES, file.path(OUT,"cohorts_liquidbiopsy_GSE73002.csv"))
cat("\n=== GSE73002 targets, breast cancer vs non-cancer ===\n")
print(as.data.frame(RES[comparison=="breast cancer vs non-cancer"][order(-abs(auc_rank-0.5))][,
  .(miRNA,probe,det_case=round(det_rate_case,3),det_ctrl=round(det_rate_ctrl,3),
    AUCraw=round(auc_raw,3),AUCrank=round(auc_rank,3),pctile_rank=round(pctile_auc_rank_among_all,1),
    q=signif(q_rank,3))]))
cat("\n=== breast cancer vs benign breast disease ===\n")
print(as.data.frame(RES[comparison=="breast cancer vs benign breast disease"][,
  .(miRNA,AUCrank=round(auc_rank,3),q=signif(q_rank,3))]))
cat("\n=== prostate disease vs non-cancer (disease-specificity control) ===\n")
print(as.data.frame(RES[comparison=="prostate disease vs non-cancer"][,
  .(miRNA,AUCraw=round(auc_raw,3),AUCrank=round(auc_rank,3))]))
