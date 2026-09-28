#!/usr/bin/env Rscript
## 49_neoadj_meanz_sensitivity.R
## Sensitivity analysis: repeat the pCR association using the *mean-of-z* score that Farmer
## et al. actually used ("the average of the 50 genes"), instead of ssGSEA, and check that
## the small five-gene exemplar module behaves the same way under both scorings.
suppressPackageStartupMessages({library(data.table); library(matrixStats)})
BASE <- "/path/to/revision"
source(file.path(BASE, "scripts/v5/46_helpers.R"))
LOG <- file.path(BASE, "logs/v5/49_neoadj_meanz_sensitivity.log")
logf <- function(...) { m <- paste0(...); cat(m,"\n"); cat(m,"\n",file=LOG,append=TRUE) }
cat("", file = LOG)
sets  <- readRDS(file.path(BASE, "cache/v5/farmer_signature_sets.rds"))
store <- readRDS(file.path(BASE, "cache/v5/farmer_neoadjuvant_store.rds"))
PLAT <- c(GSE25066="GPL96", GSE20194="GPL96", GSE41998="GPL571", GSE22093="GPL96",
          GSE23988="GPL96", GSE32646="GPL570", GSE42822="GPL96")
KEY <- c("FARMER_STROMAL","COLLAGEN_PAIR","MS_exemplar_module","FFL_3node",
         "FFL_higher_order_only","FFL_6node_all","WEST_DTF_FIBROMATOSIS","WINSLOW_STROMAL_SIG1",
         "CAF_scRNA_50","NABA_CORE_MATRISOME","HALLMARK_EMT","CHANG_WOUND_UP_VANTVEER")
auc_fun <- function(s, y) { r <- rank(s); n1 <- sum(y==1); n0 <- sum(y==0)
  (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }
rows <- list(); k <- 0L; comp <- list(); kc <- 0L
for (g in names(PLAT)) {
  s <- read_series(g); X <- s$X
  if (max(X, na.rm=TRUE) > 60) { X[X < 1] <- 1; X <- log2(X) }
  X <- X[rowSums(is.na(X)) == 0, , drop=FALSE]
  Y <- collapse_to_symbol(X, annot_map(PLAT[[g]])); Y <- Y[rowSds(Y) > 0, , drop=FALSE]
  Z <- (Y - rowMeans(Y)) / rowSds(Y)
  y <- store[[g]]$y; er <- store[[g]]$er
  for (nm in KEY) {
    gg <- intersect(sets[[nm]], rownames(Z)); if (length(gg) < 2) next
    x <- as.numeric(scale(colMeans(Z[gg,,drop=FALSE])))
    ok <- !is.na(y) & is.finite(x)
    if (sum(ok) < 40 || sum(y[ok]==1) < 5) next
    f <- glm(y[ok] ~ x[ok], family=binomial); cs <- summary(f)$coefficients
    k <- k + 1L
    rows[[k]] <- data.table(cohort=g, score=nm, scoring="meanZ", n_genes=length(gg),
      n=sum(ok), n_pCR=sum(y[ok]==1), logOR=cs[2,1], se=cs[2,2], OR=exp(cs[2,1]),
      OR_lo=exp(cs[2,1]-1.96*cs[2,2]), OR_hi=exp(cs[2,1]+1.96*cs[2,2]), p=cs[2,4],
      AUC=auc_fun(x[ok], y[ok]))
    kc <- kc + 1L
    comp[[kc]] <- data.table(cohort=g, score=nm,
      rho_meanZ_vs_ssGSEA = cor(x, as.numeric(store[[g]]$SC[nm, ]), method="spearman"))
  }
  logf("scored ", g, " (mean-z), genes=", nrow(Z))
}
R <- rbindlist(rows)
fwrite(R, file.path(BASE, "results/v5/farmer_neoadjuvant_pcr_meanz.csv"))
CMP <- rbindlist(comp)
fwrite(CMP, file.path(BASE, "results/v5/farmer_scoring_method_agreement.csv"))
meta_one <- function(d) { w <- 1/d$se^2; mu <- sum(w*d$logOR)/sum(w); se <- sqrt(1/sum(w))
  Q <- sum(w*(d$logOR-mu)^2); df <- nrow(d)-1
  tau2 <- max(0,(Q-df)/(sum(w)-sum(w^2)/sum(w))); wr <- 1/(d$se^2+tau2)
  mur <- sum(wr*d$logOR)/sum(wr); ser <- sqrt(1/sum(wr))
  data.table(k=nrow(d), n=sum(d$n), OR_fixed=exp(mu), lo=exp(mu-1.96*se), hi=exp(mu+1.96*se),
    p_fixed=2*pnorm(-abs(mu/se)), OR_random=exp(mur), p_random=2*pnorm(-abs(mur/ser)),
    I2=max(0,100*(Q-df)/Q)) }
M <- R[, meta_one(.SD), by=score]; setorder(M, OR_fixed)
fwrite(M, file.path(BASE, "results/v5/farmer_neoadjuvant_meta_meanz.csv"))
logf("\n=== mean-z scoring: meta-analysis of pCR odds per +1 SD ===")
for (i in seq_len(nrow(M))) { r <- M[i]
  logf(sprintf("  %-26s k=%d n=%4d | fixed OR=%.3f (%.3f-%.3f) p=%.3g | random OR=%.3f p=%.3g | I2=%.0f%%",
    r$score, r$k, r$n, r$OR_fixed, r$lo, r$hi, r$p_fixed, r$OR_random, r$p_random, r$I2)) }
logf("\n=== agreement between mean-z and ssGSEA scorings (Spearman rho) ===")
a <- CMP[, .(median_rho=median(rho_meanZ_vs_ssGSEA), min_rho=min(rho_meanZ_vs_ssGSEA)), by=score]
for (i in seq_len(nrow(a))) logf(sprintf("  %-26s median rho=%.3f  min=%.3f", a$score[i], a$median_rho[i], a$min_rho[i]))
logf("DONE 49")
