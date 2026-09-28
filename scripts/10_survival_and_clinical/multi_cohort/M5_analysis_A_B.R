## ==========================================================================
## M5_analysis_A_B.R
##  A) DE-direction concordance with TCGA for the network hubs (mRNA + miRNA)
##  B) miR-130a tumour-vs-normal, Cox models in every miRNA cohort with an
##     outcome, and a DerSimonian-Laird random-effects meta-analysis of the
##     hazard ratio.
## ==========================================================================
suppressPackageStartupMessages({library(data.table); library(survival); library(metafor)})
setwd("/path/to/revision")
CA <- "cache/v4/multicohort"; OUT <- "results/v4"
S  <- readRDS(file.path(CA,"genesets.rds"))
MI <- readRDS(file.path(CA,"mirna_cohorts.rds"))
set.seed(1)

## ---------- DerSimonian-Laird random effects (own code, checked vs metafor)
dl <- function(y, v){
  ok <- is.finite(y) & is.finite(v) & v>0; y<-y[ok]; v<-v[ok]; k<-length(y)
  if(k<2) return(NULL)
  w <- 1/v; yF <- sum(w*y)/sum(w); Q <- sum(w*(y-yF)^2)
  C <- sum(w)-sum(w^2)/sum(w); tau2 <- max(0,(Q-(k-1))/C)
  I2 <- max(0,100*(Q-(k-1))/Q); ws <- 1/(v+tau2)
  yR <- sum(ws*y)/sum(ws); se <- sqrt(1/sum(ws))
  list(k=k, est=yR, se=se, lo=yR-1.96*se, hi=yR+1.96*se, z=yR/se,
       p=2*pnorm(-abs(yR/se)), Q=Q, df=k-1, pQ=pchisq(Q,k-1,lower.tail=FALSE),
       I2=I2, tau2=tau2, w=100*ws/sum(ws), fixed=yF, ok=ok)
}
resolve <- function(rn, base){
  cand <- c(base, paste0(base,"-3p"))
  hit <- rn[tolower(rn) %in% tolower(cand)]
  hit <- hit[!grepl("\\*|-5p$|-pre$", hit)]
  unique(hit)[1]
}
zs <- function(x) as.numeric(scale(x))

################################################################################
## A1 -- mRNA hub DE-direction concordance vs TCGA (re-derived, not copied)
################################################################################
cat("\n=========== A1  mRNA hub DE direction concordance ===========\n")
DEt <- fread("results/BRCA_DEX_genes.csv")
cat("TCGA DE table cols:", paste(names(DEt), collapse=","), "\n")
fc <- names(DEt)[grepl("^logFC$", names(DEt))][1]
gn <- names(DEt)[1]
cat("TCGA DE rows:", nrow(DEt), " feature column:", gn,
    " looks-like-symbol fraction:", round(mean(grepl("^[A-Za-z]", DEt[[gn]])),4), "\n")
stopifnot(mean(grepl("^[A-Za-z]", DEt[[gn]])) > 0.95)   ## Rule 4 guard
tcga_fc <- setNames(DEt[[fc]], DEt[[gn]])
qcol <- names(DEt)[grepl("adj.P.Val|^FDR$|^q", names(DEt))][1]
tcga_q  <- setNames(DEt[[qcol]], DEt[[gn]])

EG <- fread("results/external_GEO_DE.csv")
A1 <- list()
for(g in unique(EG$gse)){
  e <- EG[gse==g]
  hub <- intersect(S$HUB_PROTEIN, e$feature)
  hub <- hub[hub %in% names(tcga_fc)]
  hub <- hub[is.finite(tcga_q[hub]) & tcga_q[hub] < 0.05]
  ef <- setNames(e$logFC, e$feature)[hub]
  conc <- sum(sign(ef)==sign(tcga_fc[hub]), na.rm=TRUE)
  n <- sum(is.finite(ef))
  bt <- binom.test(conc, n, 0.5)
  ## all network protein-coding nodes as a wider reference
  nod <- intersect(S$ALL_NETWORK_PROTEIN, e$feature); nod <- nod[nod %in% names(tcga_fc)]
  nod <- nod[is.finite(tcga_q[nod]) & tcga_q[nod]<0.05]
  ef2 <- setNames(e$logFC, e$feature)[nod]
  c2 <- sum(sign(ef2)==sign(tcga_fc[nod]), na.rm=TRUE); n2 <- sum(is.finite(ef2))
  A1[[g]] <- data.table(cohort=g, layer="mRNA", feature_set=c("network_hubs","all_network_nodes"),
    n_tested=c(n,n2), n_concordant=c(conc,c2), pct=round(100*c(conc/n,c2/n2),2),
    binom_p=c(bt$p.value, binom.test(c2,n2,0.5)$p.value),
    n_normal=e$n_normal[1], n_tumour=e$n_tumour[1],
    spearman_logFC=c(cor(ef, tcga_fc[hub], method="spearman", use="complete.obs"),
                     cor(ef2, tcga_fc[nod], method="spearman", use="complete.obs")))
  cat(sprintf("  %-10s hubs %3d/%3d (%.1f%%) p=%.3g | all nodes %3d/%3d (%.1f%%) p=%.3g | rho=%.3f\n",
      g, conc, n, 100*conc/n, bt$p.value, c2, n2, 100*c2/n2,
      binom.test(c2,n2,0.5)$p.value, cor(ef2, tcga_fc[nod], method="spearman", use="complete.obs")))
}
A1 <- rbindlist(A1)

################################################################################
## A2 -- miRNA hub / miR-130a DE-direction concordance vs TCGA (NEW)
################################################################################
cat("\n=========== A2  miRNA DE direction concordance ===========\n")
DEm <- fread("results/BRCA_DEX_mirnas.csv")
cat("TCGA miRNA DE cols:", paste(names(DEm), collapse=","), " rows:", nrow(DEm), "\n")
mfc <- setNames(DEm[[grep("^logFC$", names(DEm))[1]]], DEm[[1]])
mq  <- setNames(DEm[[grep("adj.P.Val|^FDR$|^q", names(DEm))[1]]], DEm[[1]])
strip <- function(x) sub("-3p$|-5p$","", x)
A2 <- list(); B_dn <- list()
for(nm in names(MI)){
  o <- MI[[nm]]
  if(is.null(o$normals) || ncol(o$normals)==0) next
  M <- o$M; N <- o$normals
  common <- intersect(rownames(M), rownames(N))
  M <- M[common,,drop=FALSE]; N <- N[common,,drop=FALSE]
  lfc <- rowMeans(M, na.rm=TRUE) - rowMeans(N, na.rm=TRUE)
  pv  <- apply(cbind(M,N), 1, function(r){
    g <- c(rep(1,ncol(M)), rep(0,ncol(N)))
    if(sum(is.finite(r[g==1]))<3 || sum(is.finite(r[g==0]))<3) return(NA_real_)
    tryCatch(wilcox.test(r[g==1], r[g==0])$p.value, error=function(e) NA_real_)})
  qv <- p.adjust(pv,"BH")
  ## match to TCGA miRNA names, allowing the arm suffix to be absent on arrays
  key <- names(mfc); names(key) <- strip(key)
  map <- ifelse(names(lfc) %in% names(mfc), names(lfc), key[names(lfc)])
  ok <- !is.na(map) & is.finite(mfc[map]) & is.finite(lfc)
  cat(sprintf("  %-10s miRNAs measured %4d, matched to the TCGA miRNA DE table %4d\n",
      nm, length(lfc), sum(ok)))
  for(lab in c("all_matched_miRNAs","miRNA_hubs")){
    sel <- ok
    if(lab=="miRNA_hubs") sel <- ok & (map %in% S$HUB_MIRNA | strip(map) %in% strip(S$HUB_MIRNA))
    sel <- sel & is.finite(mq[map]) & mq[map] < 0.05
    n <- sum(sel); if(n<3) next
    conc <- sum(sign(lfc[sel])==sign(mfc[map[sel]]))
    bt <- binom.test(conc,n,0.5)
    A2[[length(A2)+1]] <- data.table(cohort=nm, layer="miRNA", feature_set=lab,
      n_tested=n, n_concordant=conc, pct=round(100*conc/n,2), binom_p=bt$p.value,
      n_normal=ncol(N), n_tumour=ncol(M),
      spearman_logFC=cor(lfc[sel], mfc[map[sel]], method="spearman"))
    cat(sprintf("    %-20s %3d/%3d (%.1f%%) p=%.3g rho=%.3f\n", lab, conc, n,
        100*conc/n, bt$p.value, cor(lfc[sel], mfc[map[sel]], method="spearman")))
  }
  ## ---- B: miR-130a tumour vs normal in this cohort
  r <- resolve(rownames(o$M),"hsa-miR-130a")
  if(!is.na(r)){
    x <- o$M[r,]; y <- o$normals[r,]
    wt <- wilcox.test(x,y); tt <- t.test(x,y)
    B_dn[[nm]] <- data.table(cohort=nm, accession=o$accession, platform=o$platform,
      probe=r, n_tumour=length(x), n_normal=length(y),
      mean_tumour=mean(x,na.rm=TRUE), mean_normal=mean(y,na.rm=TRUE),
      log2FC=mean(x,na.rm=TRUE)-mean(y,na.rm=TRUE),
      t_p=tt$p.value, wilcox_p=wt$p.value,
      cohens_d=(mean(x,na.rm=TRUE)-mean(y,na.rm=TRUE))/sd(c(x,y),na.rm=TRUE))
    cat(sprintf("    miR-130a T-vs-N: log2FC %+.3f  wilcox p=%.3g  (%d vs %d)\n",
        B_dn[[nm]]$log2FC, wt$p.value, length(x), length(y)))
  }
}
A2 <- rbindlist(A2)
A <- rbind(A1, A2, fill=TRUE)
fwrite(A, file.path(OUT,"multicohort_A_DE_direction_concordance.csv"))
B1 <- rbindlist(B_dn)
fwrite(B1, file.path(OUT,"multicohort_B_mir130a_tumour_vs_normal.csv"))
cat("\nmiR-130a tumour-vs-normal table:\n"); print(B1[, .(cohort,n_tumour,n_normal,log2FC,wilcox_p)])

################################################################################
## B2 -- Cox models for miR-130a in every miRNA cohort with an outcome
################################################################################
cat("\n=========== B2  miR-130a Cox models ===========\n")
COVS <- list(TCGA_BRCA=c("age"), GSE19783=character(0), GSE22216=c("age","grade","size","node"),
             GSE37405=c("age","size","node"))
cox_rows <- list()
for(nm in names(MI)){
  o <- MI[[nm]]; if(!length(o$endpoints)) next
  r <- resolve(rownames(o$M),"hsa-miR-130a"); if(is.na(r)) next
  ph <- o$pheno; x <- as.numeric(o$M[r,])
  ## standardise within cohort (and within platform for GSE37405)
  if(nm=="GSE37405"){
    z <- rep(NA_real_, length(x))
    for(pl in unique(ph$platform)){ i <- ph$platform==pl; z[i] <- zs(x[i]) }
  } else z <- zs(x)
  for(ep in o$endpoints){
    if(ep=="relapse72m_binary"){
      d <- data.frame(y=ph$relapse_bin, z=z)
      d <- d[is.finite(d$y) & is.finite(d$z),]
      fit <- glm(y ~ z, data=d, family=binomial)
      s <- summary(fit)$coefficients["z",]
      cox_rows[[length(cox_rows)+1]] <- data.table(cohort=nm, accession=o$accession,
        endpoint=ep, model="univariate", n=nrow(d), nevent=sum(d$y),
        HR=exp(s[1]), se_logHR=s[2], lo=exp(s[1]-1.96*s[2]), hi=exp(s[1]+1.96*s[2]),
        p=s[4], scale="odds ratio per SD (logistic, no time variable)")
      cat(sprintf("  %-10s %-18s OR=%.3f [%.3f,%.3f] p=%.3g (n=%d)\n", nm, ep,
          exp(s[1]), exp(s[1]-1.96*s[2]), exp(s[1]+1.96*s[2]), s[4], nrow(d)))
      next
    }
    tv <- ph[[paste0(ep,"_time")]]; ev <- ph[[paste0(ep,"_event")]]
    if(is.null(tv)) next
    base <- data.frame(time=as.numeric(tv), event=as.integer(ev), z=z)
    if(nm=="GSE37405") base$strat <- ph$platform
    for(md in c("univariate","adjusted")){
      d <- base; cv <- COVS[[nm]]
      if(md=="adjusted" && length(cv)){
        for(v in cv){
          vv <- suppressWarnings(as.numeric(ph[[v]]))
          if(sum(is.finite(vv)) > 0.6*length(vv)) d[[v]] <- vv
        }
      }
      cv <- intersect(cv, names(d))
      if(md=="adjusted" && !length(cv)) next
      d <- d[is.finite(d$time) & d$time>0 & is.finite(d$event) & is.finite(d$z),]
      d <- d[complete.cases(d),]
      if(nrow(d)<20 || sum(d$event)<5) next
      f <- paste0("Surv(time,event) ~ z", if(md=="adjusted") paste0(" + ", paste(cv,collapse=" + ")) else "",
                  if(nm=="GSE37405") " + strata(strat)" else "")
      fit <- tryCatch(coxph(as.formula(f), data=d), error=function(e) NULL)
      if(is.null(fit)) next
      s <- summary(fit)$coefficients["z",]
      cox_rows[[length(cox_rows)+1]] <- data.table(cohort=nm, accession=o$accession,
        endpoint=ep, model=md, n=nrow(d), nevent=sum(d$event),
        HR=exp(s[1]), se_logHR=s[3], lo=exp(s[1]-1.96*s[3]), hi=exp(s[1]+1.96*s[3]),
        p=s[5], scale=paste0("HR per SD of miR-130a", if(md=="adjusted") paste0("; adj ", paste(cv,collapse="+")) else ""))
      cat(sprintf("  %-10s %-6s %-11s HR=%.3f [%.3f,%.3f] p=%.3g (n=%d, events=%d)\n",
          nm, ep, md, exp(s[1]), exp(s[1]-1.96*s[3]), exp(s[1]+1.96*s[3]), s[5], nrow(d), sum(d$event)))
    }
  }
}
CX <- rbindlist(cox_rows)
fwrite(CX, file.path(OUT,"multicohort_B_mir130a_cox_all.csv"))

################################################################################
## B3 -- random-effects meta-analysis of the miR-130a hazard ratio
################################################################################
cat("\n=========== B3  miR-130a HR meta-analysis ===========\n")
PRIMARY <- c(TCGA_BRCA="OS", GSE19783="OS", GSE22216="DRFS", GSE37405="RFS")
meta_rows <- list()
run_meta <- function(sub, label, note=""){
  if(nrow(sub)<2){ cat("  <2 studies for", label, "\n"); return(NULL) }
  y <- log(sub$HR); v <- sub$se_logHR^2
  m <- dl(y,v); if(is.null(m)) return(NULL)
  mf <- tryCatch(rma(yi=y, vi=v, method="DL"), error=function(e) NULL)
  chk <- if(is.null(mf)) NA_real_ else max(abs(c(mf$b[1]-m$est, mf$se-m$se, mf$I2-m$I2)))
  for(i in seq_len(nrow(sub)))
    meta_rows[[length(meta_rows)+1]] <<- data.table(analysis=label, row_type="study",
      cohort=sub$cohort[i], accession=sub$accession[i], endpoint=sub$endpoint[i],
      model=sub$model[i], n=sub$n[i], nevent=sub$nevent[i],
      HR=sub$HR[i], lo=sub$lo[i], hi=sub$hi[i], p=sub$p[i],
      weight_pct=m$w[i], k=NA_integer_, I2=NA_real_, tau2=NA_real_, Q=NA_real_,
      p_Q=NA_real_, metafor_max_abs_diff=NA_real_, note=note)
  meta_rows[[length(meta_rows)+1]] <<- data.table(analysis=label, row_type="RE_pooled",
    cohort="POOLED (DerSimonian-Laird)", accession="", endpoint="primary per cohort",
    model=sub$model[1], n=sum(sub$n), nevent=sum(sub$nevent),
    HR=exp(m$est), lo=exp(m$lo), hi=exp(m$hi), p=m$p, weight_pct=100,
    k=m$k, I2=m$I2, tau2=m$tau2, Q=m$Q, p_Q=m$pQ, metafor_max_abs_diff=chk, note=note)
  cat(sprintf("  %-34s k=%d  pooled HR=%.3f [%.3f,%.3f] p=%.4g  I2=%.1f%%  Q=%.2f (p=%.3g)  tau2=%.4f  metafor-diff=%.2e\n",
      label, m$k, exp(m$est), exp(m$lo), exp(m$hi), m$p, m$I2, m$Q, m$pQ, m$tau2, chk))
  invisible(m)
}
prim <- CX[model=="univariate" & endpoint==PRIMARY[cohort]]
cat("primary-endpoint univariate rows:\n"); print(prim[, .(cohort,endpoint,n,nevent,HR,lo,hi,p)])
run_meta(prim, "mir130a_HR_primary_univariate", "one endpoint per cohort")
run_meta(prim[cohort!="TCGA_BRCA"], "mir130a_HR_primary_univariate_noTCGA", "TCGA (discovery) removed")
adj <- CX[model=="adjusted" & endpoint==PRIMARY[cohort]]
run_meta(adj, "mir130a_HR_primary_adjusted", "clinically adjusted where covariates exist")
## sensitivity: TCGA PFI instead of OS (PFI is the endpoint the manuscript quotes)
alt <- rbind(CX[cohort=="TCGA_BRCA" & endpoint=="PFI" & model=="univariate"],
             prim[cohort!="TCGA_BRCA"])
run_meta(alt, "mir130a_HR_primary_univariate_TCGA_PFI", "TCGA contributes PFI, not OS")
MT <- rbindlist(meta_rows)
fwrite(MT, file.path(OUT,"multicohort_B_mir130a_meta.csv"))
cat("\nWROTE", file.path(OUT,"multicohort_B_mir130a_meta.csv"), nrow(MT), "rows\n")
