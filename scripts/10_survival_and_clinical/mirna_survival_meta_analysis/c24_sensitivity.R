# Sensitivity analyses of the miRNA survival meta-analysis: pooled estimates for alternative cohort sets,
# fixed- against random-effects estimates and Egger tests; writes results/v6/mirna_meta_sensitivity.csv and
# mirna_meta_FE_vs_RE_and_egger.csv.
suppressPackageStartupMessages({library(data.table); library(metafor)})
setwd("/path/to/revision"); OUT<-"results/v6"
C<-fread(file.path(OUT,"mirna_meta_percohort_cox.csv"))
P<-C[model=="univariate"&primary==1]
poolset<-function(df,label){
  o<-list()
  for(tg in unique(df$miRNA)){ s<-df[miRNA==tg]; if(nrow(s)<2) next
    m<-rma(yi=s$logHR,sei=s$se,method="DL")
    o[[length(o)+1]]<-data.frame(analysis=label,miRNA=tg,k=m$k,n=sum(s$n),nevent=sum(s$nevent),
      HR=exp(m$b[1]),lo=exp(m$ci.lb),hi=exp(m$ci.ub),p=m$pval,I2=m$I2,tau2=m$tau2,p_Q=m$QEp,
      cohorts=paste(s$cohort,collapse=";"),stringsAsFactors=FALSE) }
  d<-rbindlist(o); d[, q_BH:=p.adjust(p,"BH")]; d }
S<-rbindlist(list(
  poolset(P,"all6"),
  poolset(P[cohort!="GSE78870"],"drop_GSE78870_advanced_disease"),
  poolset(P[cohort!="TCGA_BRCA"],"drop_TCGA_discovery"),
  poolset(P[!cohort %in% c("TCGA_BRCA","GSE19783")],"drop_both_cohorts_already_in_the_paper"),
  poolset(P[cohort %in% c("GSE59829","GSE78870")],"the_two_cohorts_new_in_v6_only")))
## fixed-effect (inverse variance) comparison + Egger-style small-study check
fe<-list()
for(tg in unique(P$miRNA)){ s<-P[miRNA==tg]; if(nrow(s)<3) next
  mf<-rma(yi=s$logHR,sei=s$se,method="FE"); mr<-rma(yi=s$logHR,sei=s$se,method="DL")
  rg<-try(regtest(mr, model="lm"),silent=TRUE)
  fe[[length(fe)+1]]<-data.frame(miRNA=tg,k=mr$k,HR_RE=exp(mr$b[1]),p_RE=mr$pval,
    HR_FE=exp(mf$b[1]),p_FE=mf$pval,
    egger_p=if(inherits(rg,"try-error")) NA else rg$pval, stringsAsFactors=FALSE) }
FE<-rbindlist(fe)
fwrite(S, file.path(OUT,"mirna_meta_sensitivity.csv"))
fwrite(FE, file.path(OUT,"mirna_meta_FE_vs_RE_and_egger.csv"))
cat("=== sensitivity: pooled HR under cohort exclusions (miR-29a, miR-195, miR-204) ===\n")
print(as.data.frame(S[miRNA %in% c("miR-29a","miR-195","miR-204"),
  .(miRNA,analysis,k,n,nevent,HR=round(HR,3),lo=round(lo,3),hi=round(hi,3),
    p=signif(p,3),I2=round(I2,1))][order(miRNA,analysis)]))
cat("\n=== fixed vs random effects + Egger regression test ===\n")
print(as.data.frame(FE[order(p_RE),.(miRNA,k,HR_RE=round(HR_RE,3),p_RE=signif(p_RE,3),
  HR_FE=round(HR_FE,3),p_FE=signif(p_FE,3),egger_p=signif(egger_p,3))]))
