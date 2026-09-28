#!/usr/bin/env Rscript
## 48_neoadj_stratified.R
## Farmer's original result was obtained in ER-NEGATIVE tumours receiving FEC.  This script
##  (a) repeats the pCR analysis within ER-negative and ER-positive strata,
##  (b) runs the head-to-head logistic models (Farmer score vs FFL module score) in the
##      largest single cohort, GSE25066, and
##  (c) reports a GSE25066-only result, since the seven series overlap in patients.
suppressPackageStartupMessages({library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v5/48_neoadj_stratified.log")
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG)
store <- readRDS(file.path(BASE, "cache/v5/farmer_neoadjuvant_store.rds"))
KEY <- c("FARMER_STROMAL","COLLAGEN_PAIR","MS_exemplar_module","FFL_3node",
         "FFL_higher_order_only","FFL_6node_all","WEST_DTF_FIBROMATOSIS","WINSLOW_STROMAL_SIG1",
         "CAF_scRNA_50","ESTIMATE_STROMAL_noCOL","NABA_CORE_MATRISOME","HALLMARK_EMT",
         "HALLMARK_TGF_BETA","FINAK_SDPP","TRIULZI_ECM","CHANG_WOUND_UP_VANTVEER")
auc_fun <- function(score, y) { r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

rows <- list(); k <- 0L
for (g in names(store)) {
  s <- store[[g]]
  for (st in c("ER_negative", "ER_positive")) {
    sel <- if (st == "ER_negative") s$er == "neg" else s$er == "pos"
    sel[is.na(sel)] <- FALSE
    for (nm in intersect(KEY, rownames(s$SC))) {
      x <- as.numeric(scale(s$SC[nm, ]))
      ok <- sel & !is.na(s$y) & is.finite(x)
      if (sum(ok) < 30 || sum(s$y[ok] == 1) < 5 || sum(s$y[ok] == 0) < 5) next
      f <- glm(s$y[ok] ~ x[ok], family = binomial); cs <- summary(f)$coefficients
      k <- k + 1L
      rows[[k]] <- data.table(cohort = g, stratum = st, score = nm, n = sum(ok),
        n_pCR = sum(s$y[ok] == 1), logOR = cs[2,1], se = cs[2,2], OR = exp(cs[2,1]),
        OR_lo = exp(cs[2,1] - 1.96*cs[2,2]), OR_hi = exp(cs[2,1] + 1.96*cs[2,2]),
        p = cs[2,4], AUC = auc_fun(x[ok], s$y[ok]))
    }
  }
}
R <- rbindlist(rows)
fwrite(R, file.path(BASE, "results/v5/farmer_neoadjuvant_pcr_ERstratified.csv"))
meta_one <- function(d) { w <- 1/d$se^2; mu <- sum(w*d$logOR)/sum(w); se <- sqrt(1/sum(w))
  Q <- sum(w*(d$logOR-mu)^2); df <- nrow(d)-1
  tau2 <- max(0,(Q-df)/(sum(w)-sum(w^2)/sum(w))); wr <- 1/(d$se^2+tau2)
  mur <- sum(wr*d$logOR)/sum(wr); ser <- sqrt(1/sum(wr))
  data.table(k=nrow(d), n=sum(d$n), n_pCR=sum(d$n_pCR), OR_fixed=exp(mu),
    lo=exp(mu-1.96*se), hi=exp(mu+1.96*se), p_fixed=2*pnorm(-abs(mu/se)),
    OR_random=exp(mur), p_random=2*pnorm(-abs(mur/ser)), I2=max(0,100*(Q-df)/Q)) }
M <- R[, meta_one(.SD), by = .(score, stratum)]
setorder(M, stratum, OR_fixed)
fwrite(M, file.path(BASE, "results/v5/farmer_neoadjuvant_meta_ERstratified.csv"))
logf("=== ER-stratified meta-analysis, OR of pCR per +1 SD (OR<1 = chemoresistance) ===")
for (st in unique(M$stratum)) { logf("\n-- ", st, " --")
  for (i in which(M$stratum == st)) { r <- M[i]
    logf(sprintf("  %-26s k=%d n=%3d pCR=%3d | fixed OR=%.3f (%.3f-%.3f) p=%.3g | random OR=%.3f p=%.3g | I2=%.0f%%",
      r$score, r$k, r$n, r$n_pCR, r$OR_fixed, r$lo, r$hi, r$p_fixed, r$OR_random, r$p_random, r$I2)) } }

## ------------------------------------- (b) head-to-head in GSE25066 --------------------
s <- store$GSE25066
y <- s$y; er <- s$er
h2h <- list(); h <- 0L
for (mod in c("FFL_3node","FFL_higher_order_only","FFL_6node_all","MS_exemplar_module","COLLAGEN_PAIR")) {
  f <- as.numeric(scale(s$SC["FARMER_STROMAL", ])); m <- as.numeric(scale(s$SC[mod, ]))
  ok <- !is.na(y)
  d <- data.frame(y = y[ok], f = f[ok], m = m[ok], er = er[ok])
  d <- d[complete.cases(d), ]
  g0 <- glm(y ~ f, d, family = binomial); g1 <- glm(y ~ m, d, family = binomial)
  g2 <- glm(y ~ f + m, d, family = binomial); g3 <- glm(y ~ f + m + er, d, family = binomial)
  h <- h + 1L
  h2h[[h]] <- data.table(cohort = "GSE25066", module = mod, n = nrow(d), n_pCR = sum(d$y == 1),
    OR_farmer_alone = exp(coef(g0)[2]), p_farmer_alone = summary(g0)$coefficients[2,4],
    OR_module_alone = exp(coef(g1)[2]), p_module_alone = summary(g1)$coefficients[2,4],
    OR_farmer_adj   = exp(coef(g2)["f"]), p_farmer_adj = summary(g2)$coefficients["f",4],
    OR_module_adj   = exp(coef(g2)["m"]), p_module_adj = summary(g2)$coefficients["m",4],
    LRT_p_module_adds = anova(g0, g2, test = "Chisq")[2,"Pr(>Chi)"],
    LRT_p_farmer_adds = anova(g1, g2, test = "Chisq")[2,"Pr(>Chi)"],
    OR_farmer_ERadj = exp(coef(g3)["f"]), p_farmer_ERadj = summary(g3)$coefficients["f",4],
    OR_module_ERadj = exp(coef(g3)["m"]), p_module_ERadj = summary(g3)$coefficients["m",4])
}
H <- rbindlist(h2h)
fwrite(H, file.path(BASE, "results/v5/farmer_neoadjuvant_headtohead_GSE25066.csv"))
logf("\n=== GSE25066 head-to-head: Farmer stromal score vs FFL module scores for pCR ===")
for (i in seq_len(nrow(H))) { r <- H[i]
  logf(sprintf("  %-22s Farmer alone OR=%.3f p=%.3g | module alone OR=%.3f p=%.3g | together: Farmer OR=%.3f p=%.3g, module OR=%.3f p=%.3g | LRT module adds p=%.3g",
    r$module, r$OR_farmer_alone, r$p_farmer_alone, r$OR_module_alone, r$p_module_alone,
    r$OR_farmer_adj, r$p_farmer_adj, r$OR_module_adj, r$p_module_adj, r$LRT_p_module_adds)) }
logf("DONE 48")
