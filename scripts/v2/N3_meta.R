## ==========================================================================
## N3_meta.R -- random-effects (DerSimonian-Laird) meta-analysis across the
## independent cohorts of the ETS1->COL1A1 / NFKB1->COL1A1 correlations
## (unadjusted and CAF-adjusted) and of COL1A1~COL3A1.
## Fisher-z scale; own DL implementation cross-checked against metafor::rma.
## Output: results/v2/newcohorts_meta.csv (forest-plot ready)
## ==========================================================================
suppressPackageStartupMessages({library(metafor)})
setwd("/path/to/revision")
OUT <- "results/v2"
E <- read.csv(file.path(OUT,"newcohorts_per_edge.csv"), stringsAsFactors=FALSE)
E <- E[E$class %in% c("gene_gene","TF_target"), ]
cat("per-edge rows available for meta-analysis:", nrow(E), "\n")
cat("cohorts:", paste(unique(E$cohort), collapse=", "), "\n")

NORMAL_COHORTS <- c("GTEx_breast")   ## non-tumour; excluded from primary pool

fz  <- function(r) 0.5*log((1+r)/(1-r))
ifz <- function(z) (exp(2*z)-1)/(exp(2*z)+1)

## ---- DerSimonian-Laird random-effects, implemented here -------------------
dl <- function(y, v){
  ok <- is.finite(y) & is.finite(v) & v > 0
  y <- y[ok]; v <- v[ok]; k <- length(y)
  if(k < 2) return(NULL)
  w  <- 1/v
  yF <- sum(w*y)/sum(w)
  Q  <- sum(w*(y-yF)^2)
  C  <- sum(w) - sum(w^2)/sum(w)
  tau2 <- max(0, (Q - (k-1))/C)
  I2 <- max(0, 100*(Q-(k-1))/Q)
  ws <- 1/(v+tau2)
  yR <- sum(ws*y)/sum(ws)
  se <- sqrt(1/sum(ws))
  list(k=k, est=yR, se=se, lo=yR-1.96*se, hi=yR+1.96*se,
       z=yR/se, p=2*pnorm(-abs(yR/se)),
       Q=Q, df=k-1, pQ=pchisq(Q, k-1, lower.tail=FALSE), I2=I2, tau2=tau2,
       w=100*ws/sum(ws), fixed=yF)
}

rows <- list()
add <- function(...) rows[[length(rows)+1]] <<- data.frame(..., stringsAsFactors=FALSE)

EDGES <- c("COL1A1~COL3A1","ETS1->COL1A1","NFKB1->COL1A1","RELA->COL1A1","SP1->COL1A1")
ADJ   <- c(raw="rho_raw", cafA="rho_cafA", cafB="rho_cafB")
SE    <- c(raw="se_raw",  cafA="se_cafA",  cafB="se_cafB")

for(ed in EDGES) for(aj in names(ADJ)){
  sub <- E[E$edge == ed, ]
  sub <- sub[is.finite(sub[[ADJ[aj]]]) & is.finite(sub[[SE[aj]]]), ]
  if(!nrow(sub)) next
  for(pool in c("tumour_only","all_cohorts","independent_tumour_only")){
    s <- switch(pool,
      tumour_only             = sub[!(sub$cohort %in% NORMAL_COHORTS), ],
      all_cohorts             = sub,
      ## discovery cohort (TCGA) AND the non-tumour baseline removed: this is
      ## the pool that answers "does it replicate in cohorts we did not use to
      ## derive the finding?"
      independent_tumour_only = sub[!(sub$cohort %in% c(NORMAL_COHORTS,"TCGA_BRCA")), ])
    if(nrow(s) < 2) next
    y <- fz(s[[ADJ[aj]]]); v <- s[[SE[aj]]]^2
    m <- dl(y, v); if(is.null(m)) next
    ## cross-check with metafor
    mf <- tryCatch(rma(yi=y, vi=v, method="DL"), error=function(e) NULL)
    chk <- if(is.null(mf)) NA else max(abs(c(mf$b[1]-m$est, mf$se-m$se, mf$I2-m$I2)))
    ## forest rows (one per cohort)
    for(i in seq_len(nrow(s)))
      add(edge=ed, adjustment=aj, pool=pool, row_type="study",
          cohort=s$cohort[i], accession=s$accession[i], n=s$n[i],
          rho=s[[ADJ[aj]]][i], fisher_z=y[i], se_z=s[[SE[aj]]][i],
          ci_lo_rho=ifz(y[i]-1.96*s[[SE[aj]]][i]),
          ci_hi_rho=ifz(y[i]+1.96*s[[SE[aj]]][i]),
          weight_pct=m$w[i], k=m$k, I2=NA, tau2=NA, Q=NA, p_Q=NA,
          pooled_z=NA, p_pooled=NA, metafor_max_abs_diff=NA)
    add(edge=ed, adjustment=aj, pool=pool, row_type="RE_pooled",
        cohort="POOLED (DerSimonian-Laird)", accession="", n=sum(s$n),
        rho=ifz(m$est), fisher_z=m$est, se_z=m$se,
        ci_lo_rho=ifz(m$lo), ci_hi_rho=ifz(m$hi), weight_pct=100,
        k=m$k, I2=m$I2, tau2=m$tau2, Q=m$Q, p_Q=m$pQ,
        pooled_z=m$z, p_pooled=m$p, metafor_max_abs_diff=chk)
    if(pool %in% c("tumour_only","independent_tumour_only"))
      cat(sprintf("%-14s %-5s %-24s k=%d  pooled rho=%+.3f [%+.3f,%+.3f]  z=%+.2f p=%.3g  I2=%.1f%%  Q=%.1f (p=%.3g)  metafor-diff=%.2e\n",
                  ed, aj, pool, m$k, ifz(m$est), ifz(m$lo), ifz(m$hi), m$z, m$p, m$I2, m$Q, m$pQ, chk))
  }
}
## ==========================================================================
## Additional: random-effects meta-analysis of the Cox log-hazard-ratios
## (a) the two module scores, (b) each of the 50 protein-coding hubs.
## One endpoint per cohort (the primary one) so that studies stay independent.
## ==========================================================================
EC <- read.csv(file.path(OUT,"newcohorts_per_edge.csv"), stringsAsFactors=FALSE)
COX <- EC[EC$class %in% c("cox_hub","cox_module"), ]
PRIMARY <- c(TCGA_BRCA="OS", METABRIC="OS", GSE96058="OS", GSE20685="OS",
             GSE58812="OS", GSE21653="DFS", GSE22219="DRFS")
COX$ep <- COX$target
COX <- COX[!is.na(PRIMARY[COX$cohort]) & COX$ep == PRIMARY[COX$cohort], ]
cat("\nCox rows in the primary-endpoint meta-analysis:", nrow(COX),
    " cohorts:", paste(unique(COX$cohort), collapse=","), "\n")

cox_rows <- list()
for(feat in unique(COX$source)){
  s <- COX[COX$source==feat & is.finite(COX$rho_raw) & is.finite(COX$se_raw) & COX$se_raw>0, ]
  if(nrow(s) < 2) next
  y <- log(s$rho_raw); v <- s$se_raw^2
  m <- dl(y, v); if(is.null(m)) next
  for(i in seq_len(nrow(s)))
    cox_rows[[length(cox_rows)+1]] <- data.frame(
      edge=paste0("COXMETA_",feat), adjustment="univariate", pool="tumour_only",
      row_type="study", cohort=s$cohort[i], accession=s$accession[i], n=s$n[i],
      rho=s$rho_raw[i], fisher_z=y[i], se_z=s$se_raw[i],
      ci_lo_rho=exp(y[i]-1.96*s$se_raw[i]), ci_hi_rho=exp(y[i]+1.96*s$se_raw[i]),
      weight_pct=m$w[i], k=m$k, I2=NA, tau2=NA, Q=NA, p_Q=NA,
      pooled_z=NA, p_pooled=NA, metafor_max_abs_diff=NA, stringsAsFactors=FALSE)
  cox_rows[[length(cox_rows)+1]] <- data.frame(
    edge=paste0("COXMETA_",feat), adjustment="univariate", pool="tumour_only",
    row_type="RE_pooled", cohort="POOLED (DerSimonian-Laird)", accession="",
    n=sum(s$n), rho=exp(m$est), fisher_z=m$est, se_z=m$se,
    ci_lo_rho=exp(m$lo), ci_hi_rho=exp(m$hi), weight_pct=100,
    k=m$k, I2=m$I2, tau2=m$tau2, Q=m$Q, p_Q=m$pQ,
    pooled_z=m$z, p_pooled=m$p, metafor_max_abs_diff=NA, stringsAsFactors=FALSE)
}
CX <- do.call(rbind, cox_rows)
pool_rows <- CX$row_type=="RE_pooled"
CX$q_pooled <- NA
hubfeat <- !grepl("_module$", CX$edge)
idx <- which(pool_rows & hubfeat)
CX$q_pooled[idx] <- p.adjust(CX$p_pooled[idx], "BH")
cat("hub features meta-analysed:", length(idx),
    " pooled FDR<0.05:", sum(CX$q_pooled[idx] < 0.05, na.rm=TRUE), "\n")
pr <- CX[pool_rows, c("edge","k","rho","ci_lo_rho","ci_hi_rho","I2","p_pooled","q_pooled")]
pr <- pr[order(pr$p_pooled), ]
cat("\ntop 12 pooled Cox associations (HR per SD, primary endpoint per cohort):\n")
print(head(pr, 12), row.names=FALSE, digits=3)
cat("\nmodule rows:\n"); print(pr[grepl("_module$", pr$edge), ], row.names=FALSE, digits=3)

M <- do.call(rbind, rows)
M$q_pooled <- NA
M <- rbind(M, CX)
write.csv(M, file.path(OUT,"newcohorts_meta.csv"), row.names=FALSE)
cat("\nWROTE", file.path(OUT,"newcohorts_meta.csv"), nrow(M), "rows\n")
