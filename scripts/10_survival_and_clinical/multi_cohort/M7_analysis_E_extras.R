## ==========================================================================
## M7_analysis_E_extras.R
##  E) subtype-stratified miR-130a target-activity score (basal/TNBC vs LumA)
##  + pooled tumour-vs-normal effect size for miR-130a
##  + head-to-head 3-node vs higher-order-only FFL modules in one Cox model
##  + module scores in tumour vs non-tumour breast (GTEx / TCGA normals)
## ==========================================================================
suppressPackageStartupMessages({library(data.table); library(survival); library(metafor)})
setwd("/path/to/revision")
CA <- "cache/v4/multicohort"; OUT <- "results/v4"
S  <- readRDS(file.path(CA,"genesets.rds"))
CO <- readRDS(file.path(CA,"mrna_cohorts.rds"))
MI <- readRDS(file.path(CA,"mirna_cohorts.rds"))
SC <- readRDS(file.path(CA,"module_scores.rds"))
dl <- function(y,v){ ok<-is.finite(y)&is.finite(v)&v>0; y<-y[ok]; v<-v[ok]; k<-length(y)
  if(k<2) return(NULL); w<-1/v; yF<-sum(w*y)/sum(w); Q<-sum(w*(y-yF)^2)
  C<-sum(w)-sum(w^2)/sum(w); tau2<-max(0,(Q-(k-1))/C); I2<-max(0,100*(Q-(k-1))/Q)
  ws<-1/(v+tau2); yR<-sum(ws*y)/sum(ws); se<-sqrt(1/sum(ws))
  list(k=k,est=yR,se=se,lo=yR-1.96*se,hi=yR+1.96*se,z=yR/se,p=2*pnorm(-abs(yR/se)),
       Q=Q,df=k-1,pQ=pchisq(Q,k-1,lower.tail=FALSE),I2=I2,tau2=tau2,w=100*ws/sum(ws)) }

################################################################################
## B4 -- pooled tumour-vs-normal effect for miR-130a (Hedges g)
################################################################################
cat("\n=========== B4  pooled miR-130a tumour-vs-normal ===========\n")
rows <- list()
for(nm in names(MI)){
  o <- MI[[nm]]; if(is.null(o$normals) || ncol(o$normals)==0) next
  r <- rownames(o$M)[tolower(rownames(o$M)) %in% c("hsa-mir-130a","hsa-mir-130a-3p")]
  r <- r[!grepl("\\*|-5p$", r)][1]; if(is.na(r)) next
  x <- as.numeric(o$M[r,]); y <- as.numeric(o$normals[r,])
  n1 <- sum(is.finite(x)); n2 <- sum(is.finite(y))
  s1 <- var(x, na.rm=TRUE); s2 <- var(y, na.rm=TRUE)
  sp <- sqrt(((n1-1)*s1 + (n2-1)*s2)/(n1+n2-2))
  d  <- (mean(x,na.rm=TRUE) - mean(y,na.rm=TRUE))/sp
  J  <- 1 - 3/(4*(n1+n2)-9); g <- J*d
  vg <- J^2*((n1+n2)/(n1*n2) + d^2/(2*(n1+n2-2)))
  rows[[nm]] <- data.table(cohort=nm, accession=o$accession, platform=o$platform,
    n_tumour=n1, n_normal=n2, hedges_g=g, se_g=sqrt(vg),
    lo=g-1.96*sqrt(vg), hi=g+1.96*sqrt(vg),
    p=2*pnorm(-abs(g/sqrt(vg))), wilcox_p=wilcox.test(x,y)$p.value)
}
TN <- rbindlist(rows)
m <- dl(TN$hedges_g, TN$se_g^2)
mf <- rma(yi=TN$hedges_g, vi=TN$se_g^2, method="DL")
TN[, weight_pct := m$w]
POOL <- data.table(cohort="POOLED (DerSimonian-Laird)", accession="", platform="",
  n_tumour=sum(TN$n_tumour), n_normal=sum(TN$n_normal), hedges_g=m$est, se_g=m$se,
  lo=m$lo, hi=m$hi, p=m$p, wilcox_p=NA_real_, weight_pct=100)
TNall <- rbind(TN, POOL, fill=TRUE)
TNall[cohort=="POOLED (DerSimonian-Laird)", `:=`(k=m$k, I2=m$I2, tau2=m$tau2, Q=m$Q, p_Q=m$pQ,
  metafor_max_abs_diff=max(abs(c(mf$b[1]-m$est, mf$se-m$se, mf$I2-m$I2))))]
print(TNall[, .(cohort,n_tumour,n_normal,hedges_g=round(hedges_g,3),lo=round(lo,3),
                hi=round(hi,3),p=signif(p,3),I2=round(I2,1))])
fwrite(TNall, file.path(OUT,"multicohort_B_mir130a_tumour_vs_normal_meta.csv"))
cat(sprintf("pooled Hedges g = %+.3f [%+.3f,%+.3f] p=%.3g  I2=%.1f%%  Q=%.2f (p=%.3g)  k=%d\n",
    m$est, m$lo, m$hi, m$p, m$I2, m$Q, m$pQ, m$k))

################################################################################
## D2 -- 3-node vs higher-order-only in the SAME Cox model
################################################################################
cat("\n=========== D2  3-node vs higher-order head-to-head ===========\n")
PRIMARY <- c(TCGA_BRCA="OS", METABRIC="OS", GSE96058="OS", GSE20685="OS",
             GSE58812="OS", GSE21653="DFS", GSE22219="DRFS")
h2h <- list()
for(nm in names(PRIMARY)){
  ph <- CO[[nm]]$pheno; M <- SC[[nm]]; ep <- PRIMARY[nm]
  d <- data.frame(time=suppressWarnings(as.numeric(ph[[paste0(ep,"_time")]])),
                  event=suppressWarnings(as.integer(ph[[paste0(ep,"_event")]])),
                  z3=as.numeric(scale(M$FFL_3NODE_UNION)),
                  zH=as.numeric(scale(M$FFL_HIGHER_ONLY)))
  d <- d[complete.cases(d) & d$time>0, ]
  cat(sprintf("  %-10s cor(3node, higher-only) = %+.3f\n", nm, cor(d$z3,d$zH)))
  fit <- coxph(Surv(time,event) ~ z3 + zH, data=d)
  s <- summary(fit)$coefficients
  for(v in c("z3","zH"))
    h2h[[length(h2h)+1]] <- data.table(cohort=nm, endpoint=ep, term=v,
      n=nrow(d), nevent=sum(d$event), HR=exp(s[v,1]), se_logHR=s[v,3],
      lo=exp(s[v,1]-1.96*s[v,3]), hi=exp(s[v,1]+1.96*s[v,3]), p=s[v,5],
      cor_3node_higher=cor(d$z3,d$zH))
}
H <- rbindlist(h2h)
mr <- list()
for(v in c("z3","zH")){
  s <- H[term==v]; mm <- dl(log(s$HR), s$se_logHR^2)
  mr[[v]] <- data.table(term=v, k=mm$k, n=sum(s$n), nevent=sum(s$nevent),
    HR=exp(mm$est), lo=exp(mm$lo), hi=exp(mm$hi), p=mm$p, I2=mm$I2, Q=mm$Q, p_Q=mm$pQ)
}
MR <- rbindlist(mr)
print(H[, .(cohort,term,HR=round(HR,3),lo=round(lo,3),hi=round(hi,3),p=signif(p,3))])
cat("\nmutually adjusted pooled:\n"); print(MR[, .(term,k,n,nevent,HR=round(HR,3),
  lo=round(lo,3),hi=round(hi,3),p=signif(p,3),I2=round(I2,1))])
fwrite(rbind(H, MR, fill=TRUE), file.path(OUT,"multicohort_D_3node_vs_higherorder_headtohead.csv"))

################################################################################
## E -- subtype-stratified module survival
################################################################################
cat("\n=========== E  subtype-stratified ===========\n")
## map every cohort's subtype/ER field onto the two in-vitro models
strat_of <- function(nm, ph){
  st <- as.character(ph$subtype); er <- as.character(ph$ER)
  basal <- rep(FALSE, nrow(ph)); luma <- rep(FALSE, nrow(ph))
  if(!all(is.na(st))){
    basal <- st %in% c("Basal","TNBC")
    luma  <- st %in% c("LumA")
  }
  if(!any(basal) && !all(is.na(er))) basal <- er=="Negative" & !is.na(er)
  if(!any(luma)  && !all(is.na(er))) luma  <- er=="Positive" & !is.na(er)
  list(basal=basal & !is.na(basal), luma=luma & !is.na(luma),
       src=if(!all(is.na(st))) "PAM50/subtype call" else "ER IHC proxy")
}
MODS <- c("MIR130A_ACTIVITY","MIR130A_STRONG","MIR130A_TS_ANCHOR",
          "FFL_3NODE_UNION","FFL_HIGHER_ONLY","MIR29_ECM","PROLIF")
erows <- list()
for(nm in names(PRIMARY)){
  ph <- CO[[nm]]$pheno; M <- SC[[nm]]; ep <- PRIMARY[nm]
  ss <- strat_of(nm, ph)
  cat(sprintf("  %-10s stratum source=%-20s basal/TNBC n=%4d  LumA/ER+ n=%4d\n",
      nm, ss$src, sum(ss$basal), sum(ss$luma)))
  for(lab in c("basal_TNBC","luminalA_ERpos")){
    sel <- if(lab=="basal_TNBC") ss$basal else ss$luma
    if(sum(sel) < 40) next
    for(k in MODS){
      x <- M[[k]][sel]; if(all(is.na(x))) next
      d <- data.frame(time=suppressWarnings(as.numeric(ph[[paste0(ep,"_time")]]))[sel],
                      event=suppressWarnings(as.integer(ph[[paste0(ep,"_event")]]))[sel],
                      z=as.numeric(scale(x)))
      d <- d[complete.cases(d) & d$time>0, ]
      if(nrow(d)<40 || sum(d$event)<8) next
      fit <- coxph(Surv(time,event) ~ z, data=d)
      s <- summary(fit)$coefficients["z",]
      erows[[length(erows)+1]] <- data.table(cohort=nm, stratum=lab, stratum_source=ss$src,
        module=k, endpoint=ep, n=nrow(d), nevent=sum(d$event), HR=exp(s[1]),
        se_logHR=s[3], lo=exp(s[1]-1.96*s[3]), hi=exp(s[1]+1.96*s[3]), p=s[5])
    }
  }
}
E <- rbindlist(erows)
emeta <- list()
for(lab in unique(E$stratum)) for(k in unique(E$module)){
  s <- E[stratum==lab & module==k & is.finite(se_logHR) & se_logHR>0]
  if(nrow(s)<2) next
  mm <- dl(log(s$HR), s$se_logHR^2); if(is.null(mm)) next
  emeta[[length(emeta)+1]] <- data.table(stratum=lab, module=k, k=mm$k,
    n=sum(s$n), nevent=sum(s$nevent), HR=exp(mm$est), lo=exp(mm$lo), hi=exp(mm$hi),
    p=mm$p, I2=mm$I2, Q=mm$Q, p_Q=mm$pQ, cohorts=paste(s$cohort, collapse=";"))
}
EM <- rbindlist(emeta)
EM[, q := p.adjust(p,"BH"), by=stratum]
cat("\npooled HR per SD within stratum:\n")
print(EM[, .(stratum,module,k,n,nevent,HR=round(HR,3),lo=round(lo,3),hi=round(hi,3),
             p=signif(p,3),q=signif(q,3),I2=round(I2,1))][order(stratum,p)])
fwrite(E,  file.path(OUT,"multicohort_E_subtype_cox_percohort.csv"))
fwrite(EM, file.path(OUT,"multicohort_E_subtype_meta.csv"))

## E2 -- difference in the score itself between the two strata + basal-vs-luminal
cat("\n=========== E2  score level by subtype ===========\n")
lev <- list()
for(nm in names(CO)){
  ph <- CO[[nm]]$pheno; M <- SC[[nm]]
  st <- as.character(ph$subtype)
  if(all(is.na(st)) || length(unique(na.omit(st)))<2) next
  for(k in MODS){
    x <- M[[k]]; if(all(is.na(x))) next
    z <- as.numeric(scale(x))
    b <- z[st %in% c("Basal","TNBC")]; l <- z[st %in% "LumA"]
    if(length(b)<10 || length(l)<10) next
    tt <- t.test(b,l); wt <- wilcox.test(b,l)
    lev[[length(lev)+1]] <- data.table(cohort=nm, module=k, n_basal=length(b),
      n_lumA=length(l), mean_basal=mean(b,na.rm=TRUE), mean_lumA=mean(l,na.rm=TRUE),
      diff_SD=mean(b,na.rm=TRUE)-mean(l,na.rm=TRUE), t_p=tt$p.value, wilcox_p=wt$p.value)
  }
}
LV <- rbindlist(lev)
print(LV[module %in% c("MIR130A_ACTIVITY","MIR130A_STRONG","MIR130A_TS_ANCHOR"),
         .(cohort,module,n_basal,n_lumA,diff_SD=round(diff_SD,3),wilcox_p=signif(wilcox_p,3))])
fwrite(LV, file.path(OUT,"multicohort_E_score_by_subtype.csv"))

## E3 -- measured miR-130a by subtype (TCGA + GSE19783), the direct analogue of
## the MDA-MB-231 vs MCF7 comparison
cat("\n=========== E3  measured miR-130a by subtype ===========\n")
ph <- MI$TCGA_BRCA$pheno; x <- as.numeric(MI$TCGA_BRCA$M["hsa-miR-130a-3p",])
tb <- data.table(subtype=ph$pam50, v=x)[!is.na(subtype) & subtype!=""]
agg <- tb[, .(n=.N, mean=mean(v), sd=sd(v), median=median(v)), by=subtype][order(-mean)]
print(agg)
bl <- tb[subtype %in% c("Basal","LumA")]
w <- wilcox.test(v ~ subtype, data=bl); t <- t.test(v ~ subtype, data=bl)
cat(sprintf("TCGA Basal vs LumA miR-130a-3p: diff=%+.3f log2 RPM, wilcox p=%.3g, t p=%.3g\n",
    mean(bl$v[bl$subtype=="Basal"])-mean(bl$v[bl$subtype=="LumA"]), w$p.value, t$p.value))
agg[, cohort := "TCGA_BRCA"]
fwrite(agg, file.path(OUT,"multicohort_E_mir130a_by_subtype.csv"))

################################################################################
## F -- module scores: tumour vs non-tumour breast
################################################################################
cat("\n=========== F  module score, tumour vs normal ===========\n")
frows <- list()
sc1 <- function(X, genes){ g<-intersect(genes,rownames(X)); if(length(g)<3) return(NULL)
  Z<-t(scale(t(X[g,,drop=FALSE]))); Z<-Z[is.finite(rowSums(Z)),,drop=FALSE]; colMeans(Z) }
for(k in c("MIR130A_ACTIVITY","MIR130A_STRONG","MIR130A_TS_ANCHOR","FFL_3NODE_UNION",
           "FFL_HIGHER_ONLY","MIR29_ECM","ALL_NETWORK_PROTEIN")){
  gs <- if(k=="MIR130A_ACTIVITY") S$MIR130A_ANTICORR else S[[k]]
  sgn <- if(k=="MIR130A_ACTIVITY") -1 else 1
  ## TCGA tumour vs its own matched normals (same platform)
  X <- cbind(CO$TCGA_BRCA$X, CO$TCGA_BRCA$normals)
  grp <- c(rep("tumour", ncol(CO$TCGA_BRCA$X)), rep("normal", ncol(CO$TCGA_BRCA$normals)))
  v <- sgn*sc1(X, gs)
  w <- wilcox.test(v[grp=="tumour"], v[grp=="normal"])
  frows[[length(frows)+1]] <- data.table(comparison="TCGA tumour vs TCGA normal",
    module=k, n_tumour=sum(grp=="tumour"), n_normal=sum(grp=="normal"),
    mean_tumour=mean(v[grp=="tumour"]), mean_normal=mean(v[grp=="normal"]),
    diff=mean(v[grp=="tumour"])-mean(v[grp=="normal"]), wilcox_p=w$p.value,
    note="same platform, within-study")
}
F <- rbindlist(frows)
print(F[, .(module,n_tumour,n_normal,diff=round(diff,3),wilcox_p=signif(wilcox_p,3))])
fwrite(F, file.path(OUT,"multicohort_F_module_tumour_vs_normal.csv"))
cat("\nGTEx is retained as a non-tumour baseline but is cross-study/cross-platform,\n",
    "so no tumour-vs-GTEx module contrast is reported here.\n")
