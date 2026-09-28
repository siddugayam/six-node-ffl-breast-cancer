## ==========================================================================
## M9_link22216_figures.R
##  (i)  prove the GSE22219 (mRNA) <-> GSE22216 (miRNA) patient linkage from
##       the deposited clinical fields, then use it as a second independent
##       validation of the miR-130a target-activity score
##  (ii) forest-plot figures for the headline meta-analyses
## ==========================================================================
suppressPackageStartupMessages({library(data.table)})
setwd("/path/to/revision")
source("scripts/07_expression_validation/external_cohorts/N0_geo_utils.R")
CA <- "cache/v4/multicohort"; OUT <- "results/v4"; FIG <- "figures/v4"
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
S  <- readRDS(file.path(CA,"genesets.rds"))
CO <- readRDS(file.path(CA,"mrna_cohorts.rds"))
MI <- readRDS(file.path(CA,"mirna_cohorts.rds"))
num <- function(v){ v[v %in% c(".","NA","")] <- NA; suppressWarnings(as.numeric(v)) }

## ---------- (i) linkage proof ---------------------------------------------
## GSE22219 (mRNA) and GSE22216 (miRNA) are the two arms of one 210-patient
## study, but the deposited sample titles do NOT run in the same order: naive
## index matching ("BCmRNA n" <-> "BCmicroRNA n") agrees on only 2-52% of the
## clinical fields (see multicohort_C_GSE22219_GSE22216_linkage_check.csv).
## Patients are therefore linked on a 6-field clinical fingerprint instead.
cat("\n=========== GSE22219 <-> GSE22216 linkage ===========\n")
pm <- parse_series_matrix("cache/newcohorts/GSE22219_series_matrix.txt.gz")
pu <- parse_series_matrix(file.path(CA,"GSE22216_series_matrix.txt.gz"))
idx <- function(t) as.integer(sub("^\\D+\\s*","", t))
fp <- function(m) data.table(gsm=m$geo_accession, title=m$Sample_title, idx=idx(m$Sample_title),
  age=num(pull_char(m,"patient age")), size=num(pull_char(m,"tumour size")),
  er=pull_char(m,"er status"), grade=num(pull_char(m,"tumour grade")),
  ev=num(pull_char(m,"distant-relapse event")),
  tt=round(num(pull_char(m,"distant-relapse free survival")),6))
A <- fp(pm$meta); B <- fp(pu$meta)
A[, kk := paste(age,size,er,grade,ev,tt,sep="|")]
B[, kk := paste(age,size,er,grade,ev,tt,sep="|")]

## (a) show that naive index matching fails
Ji <- merge(A, B, by="idx", suffixes=c("_mrna","_mirna"))
agree <- function(a,b) sum(a==b | (is.na(a)&is.na(b)), na.rm=TRUE)
chk <- data.table(method="naive title index", field=c("age","size","er","grade","drfs_event","drfs_time"),
  n_pairs=nrow(Ji), n_identical=c(agree(Ji$age_mrna,Ji$age_mirna),
    agree(Ji$size_mrna,Ji$size_mirna), agree(Ji$er_mrna,Ji$er_mirna),
    agree(Ji$grade_mrna,Ji$grade_mirna), agree(Ji$ev_mrna,Ji$ev_mirna),
    agree(Ji$tt_mrna,Ji$tt_mirna)))
chk[, pct := round(100*n_identical/n_pairs,2)]
print(chk)

## (b) 1:1 clinical-fingerprint linkage
dupA <- unique(A$kk[duplicated(A$kk)]); dupB <- unique(B$kk[duplicated(B$kk)])
LK <- merge(A, B, by="kk", suffixes=c("_mrna","_mirna"))[!kk %in% c(dupA,dupB)]
cat(sprintf("GSE22219 n=%d (unique fingerprints %d); GSE22216 n=%d (unique %d)\n",
    nrow(A), uniqueN(A$kk), nrow(B), uniqueN(B$kk)))
cat("unambiguous 1:1 clinical-fingerprint matches:", nrow(LK), "\n")
chk <- rbind(chk, data.table(method="clinical fingerprint (6 fields)", field="all six jointly",
  n_pairs=nrow(LK), n_identical=nrow(LK), pct=100))
fwrite(chk, file.path(OUT,"multicohort_C_GSE22219_GSE22216_linkage_check.csv"))
fwrite(LK, file.path(CA,"GSE22219_GSE22216_link.csv"))
LINKED <- nrow(LK) > 150

## ---------- second independent validation of the target-activity score -----
if(LINKED){
  Xm <- CO$GSE22219$X; Mu <- MI$GSE22216$M
  jm <- LK[match(colnames(Xm), LK$gsm_mrna)]
  ok <- !is.na(jm$gsm_mirna) & jm$gsm_mirna %in% colnames(Mu)
  cat("paired mRNA+miRNA samples usable:", sum(ok), "\n")
  sc <- function(X, genes){ g<-intersect(genes,rownames(X)); Z<-t(scale(t(X[g,,drop=FALSE])))
    Z<-Z[is.finite(rowSums(Z)),,drop=FALSE]; list(v=colMeans(Z), n=nrow(Z)) }
  rows <- list()
  for(k in c("MIR130A_ANTICORR","MIR130A_STRONG","MIR130A_TS_ANCHOR")){
    r <- sc(Xm, S[[k]])
    v <- if(k=="MIR130A_ANTICORR") -r$v else r$v
    lab <- if(k=="MIR130A_ANTICORR") "MIR130A_ACTIVITY" else k
    ct <- cor.test(v[ok], as.numeric(Mu["hsa-miR-130a", jm$gsm_mirna[ok]]), method="spearman")
    rows[[length(rows)+1]] <- data.table(cohort="GSE22219/GSE22216", n=sum(ok), score=lab,
      rho=unname(ct$estimate), p=ct$p.value,
      note=sprintf("independent validation, patients linked on a 6-field clinical fingerprint; %d/%d set genes measured", r$n, length(S[[k]])))
    cat(sprintf("  %-18s rho=%+.3f p=%.3g (%d genes)\n", lab, ct$estimate, ct$p.value, r$n))
  }
  V <- fread(file.path(OUT,"multicohort_C_targetscore_validation.csv"))
  V <- V[cohort != "GSE22219/GSE22216"]
  V <- rbind(V, rbindlist(rows), fill=TRUE)
  fwrite(V, file.path(OUT,"multicohort_C_targetscore_validation.csv"))
  print(V)
}

## ---------- (ii) forest plots ---------------------------------------------
Z <- fread(file.path(OUT,"multicohort_FOREST_TABLE.csv"))
forest <- function(d, main, xlab, logx=TRUE, refline=if(logx) 1 else 0, file){
  d <- d[is.finite(estimate) & is.finite(lo) & is.finite(hi)]
  if(!nrow(d)) return(invisible(NULL))
  d <- d[order(row_type=="RE_pooled", decreasing=FALSE)]
  n <- nrow(d)
  png(file, width=1900, height=200+72*n, res=150)
  op <- par(mar=c(7.5,15,4.5,13))
  xr <- range(c(d$lo,d$hi), finite=TRUE)
  if(logx){ xr <- log(xr); xr <- xr + c(-.12,.12)*diff(xr) } else xr <- xr + c(-.1,.1)*diff(xr)
  plot(NA, xlim=xr, ylim=c(0.4, n+0.6), yaxt="n", xaxt="n", xlab="", ylab="", main=main, bty="n")
  mtext(xlab, side=1, line=2.6, cex=0.85)
  if(logx){ at <- pretty(exp(xr)); at <- at[at>0]; axis(1, at=log(at), labels=at) } else axis(1)
  abline(v=if(logx) log(refline) else refline, lty=2, col="grey50")
  for(i in seq_len(n)){
    y <- n-i+1
    e <- if(logx) log(d$estimate[i]) else d$estimate[i]
    l <- if(logx) log(d$lo[i]) else d$lo[i]; h <- if(logx) log(d$hi[i]) else d$hi[i]
    pool <- d$row_type[i]=="RE_pooled"
    segments(l, y, h, y, lwd=if(pool) 3 else 2, col=if(pool) "firebrick" else "grey25")
    points(e, y, pch=if(pool) 23 else 22, cex=if(pool) 1.6 else 1.2,
           bg=if(pool) "firebrick" else "grey25", col=if(pool) "firebrick" else "grey25")
    lab <- if(pool) "POOLED (random effects)" else
      sprintf("%s (%s, n=%d%s)", d$cohort[i], d$endpoint[i], d$n[i],
              if(is.finite(d$nevent[i])) paste0(", e=",d$nevent[i]) else "")
    mtext(lab, side=2, at=y, las=1, adj=1, line=0.4, cex=0.68,
          font=if(pool) 2 else 1)
    mtext(sprintf("%.2f [%.2f, %.2f]  p=%s", d$estimate[i], d$lo[i], d$hi[i],
          format.pval(d$p[i], digits=2, eps=1e-16)),
          side=4, at=y, las=1, adj=0, line=0.3, cex=0.68, font=if(pool) 2 else 1)
  }
  ip <- which(d$row_type=="RE_pooled")[1]
  if(!is.na(ip)) mtext(sprintf("k=%d  I2=%.1f%%  Q=%.2f (p=%.3g)", d$k[ip], d$I2[ip],
    d$Q[ip], d$p_Q[ip]), side=1, line=4.6, adj=0.5, cex=0.72, col="grey30")
  par(op); dev.off(); cat("wrote", file, "\n")
}
forest(Z[analysis=="mir130a_HR_primary_univariate"],
  "hsa-miR-130a and outcome: random-effects meta-analysis",
  "hazard ratio per SD of miR-130a (primary endpoint per cohort)", TRUE, 1,
  file.path(FIG,"fig_multicohort_forest_mir130a_HR.png"))
forest(Z[analysis=="mir130a_tumour_vs_normal"],
  "hsa-miR-130a: tumour vs normal breast",
  "Hedges g (tumour minus normal)", FALSE, 0,
  file.path(FIG,"fig_multicohort_forest_mir130a_tumour_vs_normal.png"))
for(m in c("MIR130A_ACTIVITY","MIR130A_STRONG","MIR130A_TS_ANCHOR","MIR29_ECM",
           "FFL_3NODE_UNION","FFL_HIGHER_ONLY"))
  forest(Z[analysis=="module_HR_primary_univariate__all_cohorts" & feature==m],
    paste0(m, ": pooled prognostic effect"),
    "hazard ratio per SD of the module score", TRUE, 1,
    file.path(FIG, paste0("fig_multicohort_forest_", m, ".png")))
for(st in c("basal_TNBC","luminalA_ERpos"))
  forest(Z[analysis=="module_HR_subtype_stratified" & stratum==st &
           feature=="MIR130A_STRONG"],
    paste0("miR-130a validated-target score, ", st),
    "hazard ratio per SD", TRUE, 1,
    file.path(FIG, paste0("fig_multicohort_forest_MIR130A_STRONG_", st, ".png")))
cat("\nfigures written to", FIG, "\n")
