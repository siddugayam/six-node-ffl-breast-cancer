## ==========================================================================
## c18_liquid_biopsy.R -- Part C
## Are the prioritised miRNAs (and hsa-miR-130a) detectable and differential in
## serum / plasma of breast cancer patients?
## ==========================================================================
suppressPackageStartupMessages({library(data.table); library(metafor)})
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
pick <- function(M,tgt){ rn<-rownames(M); nm<-tolower(gsub("miR","mir",rn))
  star<-grepl("\\*",nm); nm2<-gsub("\\*","",nm)
  hit<-which(grepl(RX[[tgt]],nm2)&!star); if(!length(hit)) hit<-which(grepl(RX[[tgt]],nm2))
  if(!length(hit)) return(NULL)
  mu<-rowMeans(M[hit,,drop=FALSE],na.rm=TRUE); mu[!is.finite(mu)]<--Inf
  if(all(!is.finite(mu))) return(NULL); i<-hit[which.max(mu)]; list(row=as.numeric(M[i,]),name=rn[i]) }
auc <- function(x,g){ r<-rank(x); n1<-sum(g==1); n0<-sum(g==0)
  if(n1<3||n0<3) return(NA); (sum(r[g==1])-n1*(n1+1)/2)/(n1*n0) }

## ---------------------------------------------------------- GSE73002 -------
E <- fread(file.path(B,"GSE73002_expr.csv")); M<-as.matrix(E[,-1]); rownames(M)<-E[[1]]
P <- fread(file.path(B,"GSE73002_pheno.csv"))
stopifnot(identical(colnames(M),P$gsm))
## 3D-Gene floors undetected miRNAs at a per-array constant: detect = above that mode
modev <- apply(M,2,function(v){ v<-v[is.finite(v)]; t<-table(round(v,4)); as.numeric(names(t)[which.max(t)]) })
nfloor <- sapply(seq_len(ncol(M)), function(j) sum(round(M[,j],4)==round(modev[j],4), na.rm=TRUE))
cat("GSE73002: median miRNAs at the per-array floor:", median(nfloor), "of", nrow(M), "\n")
DET <- sweep(M,2,modev,">")
grp <- P$diagnosis
rows<-list()
for(tg in TARGETS){
  pk<-pick(M,tg); if(is.null(pk)) next
  x<-pk$row; d<-as.numeric(DET[pk$name,])
  for(cmp in list(c("breast cancer","non-cancer"), c("breast cancer","benign breast disease"),
                  c("prostate disease","non-cancer"))){
    i <- grp %in% cmp; if(sum(i)<40) next
    g <- as.integer(grp[i]==cmp[1])
    w <- suppressWarnings(wilcox.test(x[i]~g))
    rows[[length(rows)+1]]<-data.frame(cohort="GSE73002",accession="GSE73002",
      matrix="serum", platform="3D-Gene human miRNA v21 (GPL18941)",
      comparison=paste(cmp,collapse=" vs "), miRNA=tg, probe=pk$name,
      n1=sum(g==1), n0=sum(g==0),
      det_rate_case=mean(d[i][g==1]), det_rate_ctrl=mean(d[i][g==0]),
      median_case=median(x[i][g==1]), median_ctrl=median(x[i][g==0]),
      delta_median=median(x[i][g==1])-median(x[i][g==0]),
      AUC=auc(x[i],g), p=w$p.value, stringsAsFactors=FALSE)
  }
}
G73<-rbindlist(rows); G73[, q:=p.adjust(p,"BH"), by=comparison]
fwrite(G73, file.path(OUT,"cohorts_liquidbiopsy_GSE73002.csv"))
cat("\n=== GSE73002 serum: breast cancer (n=1280) vs non-cancer (n=2686) ===\n")
print(as.data.frame(G73[comparison=="breast cancer vs non-cancer"][order(-abs(AUC-0.5))][,
  .(miRNA,probe,det_case=round(det_rate_case,3),det_ctrl=round(det_rate_ctrl,3),
    dMed=round(delta_median,3),AUC=round(AUC,3),p=signif(p,3),q=signif(q,3))]))
cat("\n=== specificity: prostate disease vs non-cancer (same platform) ===\n")
print(as.data.frame(G73[comparison=="prostate disease vs non-cancer"][,
  .(miRNA,AUC=round(AUC,3),p=signif(p,3))]))
cat("\n=== breast cancer vs benign breast disease ===\n")
print(as.data.frame(G73[comparison=="breast cancer vs benign breast disease"][,
  .(miRNA,AUC=round(AUC,3),p=signif(p,3),q=signif(q,3))]))
