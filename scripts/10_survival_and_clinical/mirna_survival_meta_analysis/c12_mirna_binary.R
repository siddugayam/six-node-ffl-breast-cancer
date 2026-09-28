## Secondary: breast-cancer miRNA cohorts whose outcome is a binary flag with no
## follow-up time -> logistic OR per SD, meta-analysed separately from the Cox HRs.
suppressPackageStartupMessages({library(metafor); library(data.table)})
setwd("/path/to/revision")
OUT <- "results/v6"; B <- "cache/v6/cohorts/built"
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
  if(all(!is.finite(mu))) return(NULL); i<-hit[which.max(mu)]
  list(row=as.numeric(M[i,]), name=rn[i]) }
ld <- function(g){ list(E=as.matrix(fread(file.path(B,paste0(g,"_expr.csv"))),rownames=1),
                        P=fread(file.path(B,paste0(g,"_pheno.csv")))) }

JOB <- list()
d<-ld("GSE57897"); keep <- d$P$tissue=="breast cancer"
JOB$GSE57897 <- list(acc="GSE57897", plat="Homo sapiens microRNA array (GPL18722)",
  E=d$E[,keep,drop=FALSE], y=as.numeric(d$P$survival[keep]), pair=NULL,
  outcome="death event flag (no follow-up time deposited)")
d<-ld("GSE103161")
JOB$GSE103161 <- list(acc="GSE103161", plat="Custom miRCURY LNA microRNA array, Exiqon 208410 (GPL23960)", E=d$E,
  y=as.numeric(d$P[["metastasis (1= yes; 0=no)"]]), pair=d$P[["pair no."]],
  outcome="distant metastasis (matched pairs, systemically untreated LN- ER+)")
d<-ld("GSE97811"); keep <- grepl("primary breast cancer", d$P$tissue) &
  d$P[["recurrence status in the first  5 years after surgery"]] %in% c("Yes","No")
JOB$GSE97811 <- list(acc="GSE97811", plat="3D-Gene Human miRNA V21_1.0.0 (GPL21263)",
  E=d$E[,keep,drop=FALSE],
  y=as.numeric(d$P[["recurrence status in the first  5 years after surgery"]][keep]=="Yes"),
  pair=NULL, outcome="5-year recurrence (women <35)")
d<-ld("GSE28321")
JOB$GSE28321 <- list(acc="GSE28321", plat="mirVANA miRNA Bioarray V2 (GPL5106)", E=d$E,
  y=as.numeric(d$P[["rfs status"]]), pair=NULL, outcome="relapse-free-survival status flag")
d<-ld("GSE40267"); cod <- d$P$cause_of_death
keep <- cod %in% c("Breast Cancer","NA","other_cause")
JOB$GSE40267 <- list(acc="GSE40267", plat="Agilent-021827 Human miRNA Microarray V3 (GPL10850)",
  E=d$E, y=as.numeric(cod=="Breast Cancer"), pair=NULL,
  outcome="breast-cancer death (no censoring time for survivors)")
d<-ld("GSE26666")
JOB$GSE26659_26666 <- list(acc="GSE26659 (= GSE26666, identical 94 GSMs)",
  plat="Agilent-019118 human miRNA 2.0 (GPL8227)", E=d$E,
  y=ifelse(d$P$relapse=="yes",1,ifelse(d$P$relapse=="no",0,NA)), pair=NULL,
  outcome="72-month relapse flag")

rows <- list()
for(cn in names(JOB)){
  j <- JOB[[cn]]; y <- j$y
  for(tg in TARGETS){
    pk <- pick(j$E, tg); if(is.null(pk)) next
    z <- as.numeric(scale(pk$row))
    dd <- data.frame(y=y, z=z); if(!is.null(j$pair)) dd$pair <- j$pair
    dd <- dd[complete.cases(dd),,drop=FALSE]
    if(nrow(dd)<25 || length(unique(dd$y))<2) next
    fit <- try(glm(y ~ z, data=dd, family=binomial()), silent=TRUE)
    if(inherits(fit,"try-error")) next
    s <- summary(fit)$coefficients
    rows[[length(rows)+1]] <- data.frame(cohort=cn, accession=j$acc, platform=j$plat,
      outcome=j$outcome, miRNA=tg, probe=pk$name, n=nrow(dd), nevent=sum(dd$y),
      OR=exp(s["z","Estimate"]), lo=exp(s["z","Estimate"]-1.96*s["z","Std. Error"]),
      hi=exp(s["z","Estimate"]+1.96*s["z","Std. Error"]),
      logOR=s["z","Estimate"], se=s["z","Std. Error"], p=s["z","Pr(>|z|)"],
      stringsAsFactors=FALSE)
  }
}
BIN <- rbindlist(rows)
BIN[, q := p.adjust(p, "BH"), by=cohort]
fwrite(BIN, file.path(OUT,"mirna_meta_binary_logistic.csv"))
out <- list()
for(tg in unique(BIN$miRNA)){
  s <- BIN[miRNA==tg]; if(nrow(s)<2) next
  m <- rma(yi=s$logOR, sei=s$se, method="DL")
  out[[length(out)+1]] <- data.frame(miRNA=tg, k=m$k, n=sum(s$n), nevent=sum(s$nevent),
    OR=exp(m$b[1]), lo=exp(m$ci.lb), hi=exp(m$ci.ub), p=m$pval, I2=m$I2, tau2=m$tau2,
    Q=m$QE, p_Q=m$QEp, stringsAsFactors=FALSE)
}
BM <- rbindlist(out); BM[, q := p.adjust(p,"BH")]
fwrite(BM, file.path(OUT,"mirna_meta_binary_pooled.csv"))
cat("\n=== binary-outcome cohorts, per-SD logistic OR, DL random effects ===\n")
print(as.data.frame(BM[order(p), .(miRNA,k,n,nevent,OR=round(OR,3),lo=round(lo,3),
      hi=round(hi,3),p=signif(p,3),q=signif(q,3),I2=round(I2,1))]))
cat("\ncohorts used:\n"); print(BIN[,.(n=max(n),nevent=max(nevent),miRNAs=.N),by=.(cohort,accession,outcome)])
