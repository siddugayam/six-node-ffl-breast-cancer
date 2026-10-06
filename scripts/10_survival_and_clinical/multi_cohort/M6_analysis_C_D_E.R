## ==========================================================================
## M6_analysis_C_D_E.R
##  C/D) Cox + random-effects meta-analysis across every mRNA cohort with an
##     outcome for the miR-29 target / ECM module and for the 3-node vs
##     higher-order FFL module scores, univariate and adjusted for
##     age/grade/size/node
##  E) subtype-stratified (basal / TNBC vs luminal A) versions
## ==========================================================================
suppressPackageStartupMessages({library(data.table); library(survival); library(metafor)})
setwd("/path/to/revision")
CA <- "cache/v4/multicohort"; OUT <- "results/v4"
S  <- readRDS(file.path(CA,"genesets.rds"))
CO <- readRDS(file.path(CA,"mrna_cohorts.rds"))
set.seed(1)
dl <- function(y,v){ ok<-is.finite(y)&is.finite(v)&v>0; y<-y[ok]; v<-v[ok]; k<-length(y)
  if(k<2) return(NULL); w<-1/v; yF<-sum(w*y)/sum(w); Q<-sum(w*(y-yF)^2)
  C<-sum(w)-sum(w^2)/sum(w); tau2<-max(0,(Q-(k-1))/C); I2<-max(0,100*(Q-(k-1))/Q)
  ws<-1/(v+tau2); yR<-sum(ws*y)/sum(ws); se<-sqrt(1/sum(ws))
  list(k=k,est=yR,se=se,lo=yR-1.96*se,hi=yR+1.96*se,z=yR/se,p=2*pnorm(-abs(yR/se)),
       Q=Q,df=k-1,pQ=pchisq(Q,k-1,lower.tail=FALSE),I2=I2,tau2=tau2,w=100*ws/sum(ws)) }

## ---------------- module scores -------------------------------------------
## mean of per-gene z-scores, computed within cohort over the tumours only.
score <- function(X, genes){
  g <- intersect(genes, rownames(X)); if(length(g)<3) return(list(v=NULL,n=length(g)))
  Z <- t(scale(t(X[g,,drop=FALSE])))
  Z <- Z[is.finite(rowSums(Z)), , drop=FALSE]
  list(v=colMeans(Z, na.rm=TRUE), n=nrow(Z))
}
SETS <- c("MIR29_TARGETS_NET","MIR29_ECM",
          "FFL_3NODE_UNION","FFL_HIGHERORDER_UNION","FFL_HIGHER_ONLY",
          "FFL_4_node","FFL_5_node","FFL_6_node","FFL_3_miR","FFL_3_TF","FFL_3_Comp",
          "ALL_NETWORK_PROTEIN","HUB_PROTEIN","PROLIF","CAF_A")
SC <- list(); COV <- list()
for(nm in names(CO)){
  o <- CO[[nm]]; X <- o$X
  s <- sapply(SETS, function(k) score(X, S[[k]])$n)
  n <- sapply(SETS, function(k){ r <- score(X,S[[k]]); if(is.null(r$v)) rep(NA_real_, ncol(X)) else r$v })
  M <- as.data.frame(n); rownames(M) <- colnames(X)
  SC[[nm]] <- M
  COV[[nm]] <- data.table(cohort=nm, set=names(s), n_genes_used=as.integer(s),
                          n_genes_in_set=sapply(SETS, function(k) length(S[[k]])))
  cat(sprintf("%-12s scored %d samples; genes used: FFL_3NODE %d/%d, FFL_HIGHER_ONLY %d/%d, MIR29_ECM %d/%d\n",
    nm, ncol(X),
    s["FFL_3NODE_UNION"], length(S$FFL_3NODE_UNION),
    s["FFL_HIGHER_ONLY"], length(S$FFL_HIGHER_ONLY), s["MIR29_ECM"], length(S$MIR29_ECM)))
}
COVd <- rbindlist(COV); fwrite(COVd, file.path(OUT,"multicohort_module_gene_coverage.csv"))
saveRDS(SC, file.path(CA,"module_scores.rds"))
ALLSETS <- SETS

################################################################################
## C/D -- Cox models for every module score in every cohort
################################################################################
cat("\n=========== C/D  module Cox models ===========\n")
EPMAP <- list(TCGA_BRCA=c("OS","PFI","DSS"), METABRIC=c("OS","RFS"), GSE96058="OS",
  GSE20685="OS", GSE58812=c("OS","DMFS"), GSE21653="DFS", GSE22219="DRFS")
covs_available <- function(ph){
  out <- character(0)
  for(v in c("age","grade","size")){
    vv <- suppressWarnings(as.numeric(ph[[v]]))
    if(sum(is.finite(vv)) > 0.6*length(vv) && length(unique(vv[is.finite(vv)]))>2) out <- c(out,v)
  }
  nd <- ph$node
  if(!all(is.na(nd))){
    ndn <- suppressWarnings(as.numeric(nd))
    if(sum(is.finite(ndn)) > 0.6*length(nd)) out <- c(out,"node")
    else if(sum(!is.na(nd)) > 0.6*length(nd) && length(unique(na.omit(nd)))>1) out <- c(out,"node")
  }
  out
}
rows <- list()
for(nm in names(EPMAP)){
  ph <- CO[[nm]]$pheno; M <- SC[[nm]]
  stopifnot(identical(rownames(M), as.character(ph$sample)))
  cv <- covs_available(ph)
  cat(sprintf("\n-- %s  covariates usable: %s\n", nm, if(length(cv)) paste(cv,collapse="+") else "none"))
  for(ep in EPMAP[[nm]]){
    tv <- suppressWarnings(as.numeric(ph[[paste0(ep,"_time")]]))
    ev <- suppressWarnings(as.integer(ph[[paste0(ep,"_event")]]))
    for(k in ALLSETS){
      x <- M[[k]]; if(all(is.na(x))) next
      for(md in c("univariate","adjusted")){
        d <- data.frame(time=tv, event=ev, z=as.numeric(scale(x)))
        use <- character(0)
        if(md=="adjusted"){
          if(!length(cv)) next
          for(v in cv){
            vv <- ph[[v]]
            d[[v]] <- if(v=="node" && !all(is.finite(suppressWarnings(as.numeric(vv))))) factor(vv) else suppressWarnings(as.numeric(vv))
          }
          use <- cv
        }
        d <- d[is.finite(d$time) & d$time>0 & is.finite(d$event) & is.finite(d$z), , drop=FALSE]
        d <- d[complete.cases(d), , drop=FALSE]
        if(nrow(d)<30 || sum(d$event)<10) next
        f <- paste0("Surv(time,event) ~ z", if(length(use)) paste0(" + ", paste(use, collapse=" + ")) else "")
        fit <- tryCatch(coxph(as.formula(f), data=d), error=function(e) NULL)
        if(is.null(fit)) next
        s <- summary(fit)$coefficients["z",]
        rows[[length(rows)+1]] <- data.table(cohort=nm, accession=CO[[nm]]$accession,
          module=k, endpoint=ep, model=md, covariates=paste(use, collapse="+"),
          n=nrow(d), nevent=sum(d$event), HR=exp(s[1]), se_logHR=s[3],
          lo=exp(s[1]-1.96*s[3]), hi=exp(s[1]+1.96*s[3]), p=s[5])
      }
    }
    cat(sprintf("   %-5s n=%d events=%d\n", ep, sum(is.finite(tv)&tv>0&is.finite(ev)), sum(ev, na.rm=TRUE)))
  }
}
CX <- rbindlist(rows)
fwrite(CX, file.path(OUT,"multicohort_CD_module_cox_all.csv"))
cat("\nmodule Cox rows:", nrow(CX), "\n")

################################################################################
## meta-analysis of every module, primary endpoint per cohort
################################################################################
cat("\n=========== C/D  meta-analysis ===========\n")
PRIMARY <- c(TCGA_BRCA="OS", METABRIC="OS", GSE96058="OS", GSE20685="OS",
             GSE58812="OS", GSE21653="DFS", GSE22219="DRFS")
meta_rows <- list()
run_meta <- function(sub, analysis, module, pool, note=""){
  if(nrow(sub)<2) return(invisible(NULL))
  y <- log(sub$HR); v <- sub$se_logHR^2
  m <- dl(y,v); if(is.null(m)) return(invisible(NULL))
  mf <- tryCatch(rma(yi=y, vi=v, method="DL"), error=function(e) NULL)
  chk <- if(is.null(mf)) NA_real_ else max(abs(c(mf$b[1]-m$est, mf$se-m$se, mf$I2-m$I2)))
  for(i in seq_len(nrow(sub)))
    meta_rows[[length(meta_rows)+1]] <<- data.table(analysis=analysis, module=module,
      pool=pool, row_type="study", cohort=sub$cohort[i], accession=sub$accession[i],
      endpoint=sub$endpoint[i], model=sub$model[i], n=sub$n[i], nevent=sub$nevent[i],
      HR=sub$HR[i], lo=sub$lo[i], hi=sub$hi[i], p=sub$p[i], weight_pct=m$w[i],
      k=NA_integer_, I2=NA_real_, tau2=NA_real_, Q=NA_real_, p_Q=NA_real_,
      metafor_max_abs_diff=NA_real_, note=note)
  meta_rows[[length(meta_rows)+1]] <<- data.table(analysis=analysis, module=module,
    pool=pool, row_type="RE_pooled", cohort="POOLED (DerSimonian-Laird)", accession="",
    endpoint="primary per cohort", model=sub$model[1], n=sum(sub$n), nevent=sum(sub$nevent),
    HR=exp(m$est), lo=exp(m$lo), hi=exp(m$hi), p=m$p, weight_pct=100, k=m$k,
    I2=m$I2, tau2=m$tau2, Q=m$Q, p_Q=m$pQ, metafor_max_abs_diff=chk, note=note)
  invisible(m)
}
for(md in c("univariate","adjusted")) for(k in ALLSETS){
  sub <- CX[module==k & model==md & endpoint==PRIMARY[cohort]]
  sub <- sub[is.finite(HR) & is.finite(se_logHR) & se_logHR>0]
  for(pool in c("all_cohorts","independent_of_TCGA")){
    s <- if(pool=="all_cohorts") sub else sub[cohort!="TCGA_BRCA"]
    run_meta(s, paste0("module_HR_primary_", md), k, pool,
             if(pool=="independent_of_TCGA") "TCGA removed (target set derived there)" else "")
  }
}
MT <- rbindlist(meta_rows)
MT[row_type=="RE_pooled", q := p.adjust(p, "BH"), by=.(analysis,pool)]
fwrite(MT, file.path(OUT,"multicohort_CD_module_meta.csv"))
P <- MT[row_type=="RE_pooled"]
cat("\n--- pooled HR per SD, univariate, all cohorts ---\n")
print(P[analysis=="module_HR_primary_univariate" & pool=="all_cohorts",
        .(module,k,n,nevent,HR=round(HR,3),lo=round(lo,3),hi=round(hi,3),
          p=signif(p,3),q=signif(q,3),I2=round(I2,1))][order(p)])
cat("\n--- pooled HR per SD, adjusted, all cohorts ---\n")
print(P[analysis=="module_HR_primary_adjusted" & pool=="all_cohorts",
        .(module,k,n,nevent,HR=round(HR,3),lo=round(lo,3),hi=round(hi,3),
          p=signif(p,3),q=signif(q,3),I2=round(I2,1))][order(p)])
cat("\n--- pooled HR, univariate, TCGA excluded ---\n")
print(P[analysis=="module_HR_primary_univariate" & pool=="independent_of_TCGA",
        .(module,k,n,nevent,HR=round(HR,3),lo=round(lo,3),hi=round(hi,3),
          p=signif(p,3),q=signif(q,3),I2=round(I2,1))][order(p)])
