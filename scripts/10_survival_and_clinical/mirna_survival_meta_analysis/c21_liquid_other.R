## Part C (continued): three further circulating-miRNA cohorts.
suppressPackageStartupMessages({library(data.table); library(survival)})
setwd("/path/to/revision")
OUT<-"results/v6"; B<-"cache/v6/cohorts/built"
TARGETS <- c("miR-21","miR-195","miR-204","miR-383","miR-124","miR-155","miR-429",
             "miR-141","miR-34a","miR-101","miR-29a")
RX <- c("miR-21"="^hsa-mir-21(-[0-9])?(-[35]p)?$","miR-195"="^hsa-mir-195(-[35]p)?$",
        "miR-204"="^hsa-mir-204(-[35]p)?$","miR-383"="^hsa-mir-383(-[35]p)?$",
        "miR-124"="^hsa-mir-124a?(-[123])?(-[35]p)?$","miR-155"="^hsa-mir-155(-[35]p)?$",
        "miR-429"="^hsa-mir-429(-[35]p)?$","miR-141"="^hsa-mir-141(-[35]p)?$",
        "miR-34a"="^hsa-mir-34a(-[35]p)?$","miR-101"="^hsa-mir-101(-[12])?(-[35]p)?$",
        "miR-29a"="^hsa-mir-29a(-[35]p)?$")
pick<-function(M,tg){rn0<-rownames(M)
  rn<-sub("_st$","",rn0); rn<-gsub("-star","*",rn)
  nm<-gsub("\\*","",tolower(gsub("miR","mir",rn)))
  hit<-which(grepl(RX[[tg]],nm)&!grepl("\\*",rn)); if(!length(hit)) return(NULL)
  mu<-rowMeans(M[hit,,drop=FALSE],na.rm=TRUE); mu[!is.finite(mu)]<--Inf
  if(all(!is.finite(mu))) return(NULL); i<-hit[which.max(mu)]; list(x=as.numeric(M[i,]),nm=rn0[i])}
ld<-function(g) list(E={x<-fread(file.path(B,paste0(g,"_expr.csv")));m<-as.matrix(x[,-1]);rownames(m)<-x[[1]];m},
                     P=fread(file.path(B,paste0(g,"_pheno.csv"))))
rank_norm <- function(M) apply(M,2,function(v){r<-rank(v,na.last="keep"); r/sum(is.finite(v))})
rows<-list()

## ---- GSE44281 : Sister Study prospective serum, matched case/control ------
d<-ld("GSE44281"); M<-rank_norm(d$E); P<-d$P
y<-as.integer(P$status=="case"); pr<-P$pair_index
cat("GSE44281 serum: cases",sum(y),"controls",sum(1-y),"pairs",length(unique(pr)),"\n")
for(tg in TARGETS){ pk<-pick(M,tg); if(is.null(pk)) next
  dd<-data.frame(y=y,z=as.numeric(scale(pk$x)),pr=pr); dd<-dd[complete.cases(dd),]
  f<-try(clogit(y~z+strata(pr),data=dd),silent=TRUE); if(inherits(f,"try-error")) next
  s<-summary(f)$coefficients
  rows[[length(rows)+1]]<-data.frame(cohort="GSE44281",accession="GSE44281",matrix="serum",
    platform="[miRNA-2] Affymetrix Multispecies miRNA-2 array (GPL14613)",
    design="prospective nested case/control, matched pairs (Sister Study)",
    outcome="incident breast cancer", miRNA=tg, probe=pk$nm, n=nrow(dd), nevent=sum(dd$y),
    effect="OR per SD (conditional logistic)", est=exp(s["z","coef"]),
    lo=exp(s["z","coef"]-1.96*s["z","se(coef)"]), hi=exp(s["z","coef"]+1.96*s["z","se(coef)"]),
    p=s["z","Pr(>|z|)"], stringsAsFactors=FALSE) }

## ---- GSE110651 : serum, metastatic BC on eribulin, new distant metastasis --
d<-ld("GSE110651"); M<-rank_norm(d$E); P<-d$P
y<-as.integer(P[["new distant metastasis"]]=="Yes")
cat("GSE110651 serum: new distant metastasis",sum(y,na.rm=TRUE),"of",length(y),"\n")
for(tg in TARGETS){ pk<-pick(M,tg); if(is.null(pk)) next
  dd<-data.frame(y=y,z=as.numeric(scale(pk$x))); dd<-dd[complete.cases(dd),]
  f<-try(glm(y~z,data=dd,family=binomial()),silent=TRUE); if(inherits(f,"try-error")) next
  s<-summary(f)$coefficients
  rows[[length(rows)+1]]<-data.frame(cohort="GSE110651",accession="GSE110651",matrix="serum",
    platform="3D-Gene Human miRNA V21_1.0.0 (GPL21263)", design="metastatic BC on eribulin",
    outcome="new distant metastasis", miRNA=tg, probe=pk$nm, n=nrow(dd), nevent=sum(dd$y),
    effect="OR per SD (logistic)", est=exp(s["z","Estimate"]),
    lo=exp(s["z","Estimate"]-1.96*s["z","Std. Error"]),
    hi=exp(s["z","Estimate"]+1.96*s["z","Std. Error"]), p=s["z","Pr(>|z|)"], stringsAsFactors=FALSE) }

## ---- GSE118782 : plasma, cancer vs control + DFS ---------------------------
d<-ld("GSE118782"); M<-rank_norm(d$E); P<-d$P
st<-P[["sample type"]]; cat("GSE118782 plasma sample types:",paste(names(table(st)),table(st),collapse="; "),"\n")
y<-as.integer(st!="Control")
tt<-suppressWarnings(as.numeric(P$dfs_time)); ev<-suppressWarnings(as.numeric(P$dfs_status))
cat("  DFS usable:",sum(is.finite(tt)&is.finite(ev)),"events",sum(ev==1,na.rm=TRUE),"\n")
for(tg in TARGETS){ pk<-pick(M,tg); if(is.null(pk)) next
  dd<-data.frame(y=y,z=as.numeric(scale(pk$x))); dd<-dd[complete.cases(dd),]
  if(length(unique(dd$y))==2 && nrow(dd)>=20){
    f<-glm(y~z,data=dd,family=binomial()); s<-summary(f)$coefficients
    rows[[length(rows)+1]]<-data.frame(cohort="GSE118782",accession="GSE118782",matrix="plasma",
      platform="[miRNA-1] Affymetrix Multispecies miRNA-1 array (GPL8786)", design="case/control",
      outcome="breast cancer vs healthy control", miRNA=tg, probe=pk$nm, n=nrow(dd), nevent=sum(dd$y),
      effect="OR per SD (logistic)", est=exp(s["z","Estimate"]),
      lo=exp(s["z","Estimate"]-1.96*s["z","Std. Error"]),
      hi=exp(s["z","Estimate"]+1.96*s["z","Std. Error"]), p=s["z","Pr(>|z|)"], stringsAsFactors=FALSE) }
  d2<-data.frame(t=tt,e=ev,z=as.numeric(scale(pk$x))); d2<-d2[complete.cases(d2)&is.finite(d2$t)&d2$t>0,]
  if(nrow(d2)>=15 && sum(d2$e)>=5){
    f<-try(coxph(Surv(t,e)~z,data=d2),silent=TRUE); if(inherits(f,"try-error")) next
    s<-summary(f)$coefficients
    rows[[length(rows)+1]]<-data.frame(cohort="GSE118782",accession="GSE118782",matrix="plasma",
      platform="[miRNA-1] Affymetrix Multispecies miRNA-1 array (GPL8786)", design="cancer plasmas only",
      outcome="disease-free survival", miRNA=tg, probe=pk$nm, n=nrow(d2), nevent=sum(d2$e),
      effect="HR per SD (Cox)", est=exp(s["z","coef"]),
      lo=exp(s["z","coef"]-1.96*s["z","se(coef)"]), hi=exp(s["z","coef"]+1.96*s["z","se(coef)"]),
      p=s["z","Pr(>|z|)"], stringsAsFactors=FALSE) }
}
## ---- GSE68373 : plasma, metastasis vs no metastasis -----------------------
d<-ld("GSE68373"); M<-rank_norm(d$E); P<-d$P
mc <- P[["metastasis class"]]
y <- ifelse(grepl("^metastasis", mc), 1L, ifelse(mc %in% c("no metastasis","No metastasis"), 0L, NA))
cat("GSE68373 plasma metastasis class:", paste(names(table(mc)), table(mc), collapse="; "), "\n")
for(tg in TARGETS){ pk<-pick(M,tg); if(is.null(pk)) next
  dd<-data.frame(y=y,z=as.numeric(scale(pk$x))); dd<-dd[complete.cases(dd),]
  if(nrow(dd)<25 || length(unique(dd$y))<2) next
  f<-glm(y~z,data=dd,family=binomial()); s<-summary(f)$coefficients
  rows[[length(rows)+1]]<-data.frame(cohort="GSE68373",accession="GSE68373",matrix="plasma",
    platform="Agilent-031181 Human miRNA V16.0 (GPL16770)", design="plasma, metastatic vs non-metastatic",
    outcome="distant metastasis", miRNA=tg, probe=pk$nm, n=nrow(dd), nevent=sum(dd$y),
    effect="OR per SD (logistic)", est=exp(s["z","Estimate"]),
    lo=exp(s["z","Estimate"]-1.96*s["z","Std. Error"]),
    hi=exp(s["z","Estimate"]+1.96*s["z","Std. Error"]), p=s["z","Pr(>|z|)"], stringsAsFactors=FALSE) }

R<-rbindlist(rows); R[, q:=p.adjust(p,"BH"), by=.(cohort,outcome)]
fwrite(R, file.path(OUT,"cohorts_liquidbiopsy_other.csv"))
for(o in unique(R$outcome)){
  cat("\n===",o,"===\n")
  print(as.data.frame(R[outcome==o][order(p),.(cohort,miRNA,probe,n,nevent,effect,
    est=round(est,3),lo=round(lo,3),hi=round(hi,3),p=signif(p,3),q=signif(q,3))]))
}
