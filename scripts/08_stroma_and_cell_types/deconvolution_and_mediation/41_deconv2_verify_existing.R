#!/usr/bin/env Rscript
## 41_deconv2_verify_existing.R
## INDEPENDENT re-derivation of the deposited deconvolution results, with code that
## deliberately does NOT use the mediation package and does not read any of the
## deposited intermediate objects except the per-sample cell-type estimates.
##  (i)  re-derive the pairwise Spearman matrix between fibroblast/stromal estimates
##  (ii) re-derive ACME / ADE / proportion mediated / Sobel z from first principles
##       (for a linear-linear mediation model ACME = a*b and ADE = c' exactly), with
##       my own nonparametric bootstrap, and compare with results/v3/deconv_mediation_by_method.csv
## Outputs: results/v3/deconv_verify_mediation_agreement.csv
##          results/v3/deconv_verify_fibcorr_agreement.csv
suppressPackageStartupMessages({library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/41_deconv2_verify.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }

gexp <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
mexp <- readRDS(file.path(BASE, "data/brca_mirna_expr_canonical.rds"))
stopifnot(!any(duplicated(rownames(gexp))))
logf("gene matrix ", nrow(gexp), " x ", ncol(gexp), "; duplicated rownames = ",
     sum(duplicated(rownames(gexp))), " (RULE 4 guard)")
logf("fraction of gene rownames that are pure digits: ",
     round(mean(grepl("^[0-9]+$", rownames(gexp))), 5))

FIB <- readRDS(file.path(BASE, "cache/v3/deconv/fibroblast_estimates_all.rds"))
logf("fibroblast estimate matrix: ", nrow(FIB), " samples x ", ncol(FIB), " estimates")

## ---------------------------------------------------------------- (i) correlations ----
dep <- as.matrix(read.csv(file.path(BASE, "results/v3/deconv_fibroblast_correlation_spearman.csv"),
                          row.names = 1, check.names = FALSE))
common <- intersect(colnames(FIB), colnames(dep))
mine <- cor(FIB[, common], method = "spearman")
delta <- max(abs(mine - dep[common, common]))
logf("(i) recomputed Spearman matrix over ", length(common), " estimates; max |mine - deposited| = ",
     signif(delta, 3))
fwrite(data.frame(n_estimates = length(common), max_abs_diff = delta,
                  min_offdiag_rho = min(mine[upper.tri(mine)]),
                  median_offdiag_rho = median(mine[upper.tri(mine)]),
                  max_offdiag_rho = max(mine[upper.tri(mine)])),
       file.path(BASE, "results/v3/deconv_verify_fibcorr_agreement.csv"))

## ---------------------------------------------------------------- (ii) mediation ----
ns <- function(v) qnorm((rank(v) - 0.5) / length(v))

med_closed <- function(x, m, y, B = 1000L, seed = 11L) {
  d <- data.frame(x = ns(x), m = ns(m), y = ns(y))
  f <- function(dd) {
    mm <- lm(m ~ x, data = dd); om <- lm(y ~ x + m, data = dd)
    c(a = unname(coef(mm)["x"]), b = unname(coef(om)["m"]), cp = unname(coef(om)["x"]))
  }
  p <- f(d)
  ACME <- p["a"] * p["b"]; ADE <- p["cp"]; TOT <- ACME + ADE
  ## analytic Sobel
  mm <- lm(m ~ x, data = d); om <- lm(y ~ x + m, data = d)
  sa <- summary(mm)$coef["x", 2]; sb <- summary(om)$coef["m", 2]
  z  <- p["a"] * p["b"] / sqrt(p["b"]^2 * sa^2 + p["a"]^2 * sb^2)
  ## nonparametric bootstrap (percentile)
  set.seed(seed); n <- nrow(d)
  bs <- vapply(seq_len(B), function(i) {
    ii <- sample.int(n, n, replace = TRUE); q <- f(d[ii, ])
    c(q["a"] * q["b"], q["cp"], q["a"] * q["b"] / (q["a"] * q["b"] + q["cp"]))
  }, numeric(3))
  ## two-sided bootstrap p (proportion of resamples on the other side of 0, x2)
  bp <- function(v, obs) { pv <- 2 * min(mean(v <= 0), mean(v >= 0)); min(1, max(pv, 1 / B)) }
  c(ACME = unname(ACME), ADE = unname(ADE), total = unname(TOT),
    prop_mediated = unname(ACME / TOT),
    acme_lo = unname(quantile(bs[1, ], 0.025)), acme_hi = unname(quantile(bs[1, ], 0.975)),
    acme_p = bp(bs[1, ], ACME),
    ade_lo = unname(quantile(bs[2, ], 0.025)), ade_hi = unname(quantile(bs[2, ], 0.975)),
    ade_p = bp(bs[2, ], ADE),
    sobel_z = unname(z), sobel_p = unname(2 * pnorm(-abs(z))),
    a_path = unname(p["a"]), b_path = unname(p["b"]), cprime = unname(p["cp"]))
}

DEP <- fread(file.path(BASE, "results/v3/deconv_mediation_by_method.csv"))
logf("deposited mediation rows = ", nrow(DEP), "; axes = ", uniqueN(paste(DEP$x, DEP$y)),
     "; mediators = ", uniqueN(DEP$mediator))

sG <- rownames(FIB)
sB <- sort(intersect(sG, colnames(mexp)))
logf("gene-assay n = ", length(sG), "; gene+miRNA n = ", length(sB))

res <- vector("list", nrow(DEP))
for (i in seq_len(nrow(DEP))) {
  r <- DEP[i]
  ismir <- grepl("^hsa-", r$x)
  ss <- if (ismir) sB else sG
  X <- if (ismir) mexp[r$x, ss] else gexp[r$x, ss]
  Y <- gexp[r$y, ss]; M <- FIB[ss, r$mediator]
  v <- med_closed(X, M, Y, B = 1000L)
  res[[i]] <- data.table(x = r$x, y = r$y, mediator = r$mediator, n = length(ss),
                         t(v), dep_ACME = r$ACME, dep_ADE = r$ADE,
                         dep_prop = r$prop_mediated, dep_sobel_z = r$sobel_z,
                         dep_sobel_p = r$sobel_p, dep_acme_p = r$acme_p,
                         dep_n = r$n)
  if (i %% 25 == 0) logf("  ... ", i, "/", nrow(DEP))
}
R <- rbindlist(res)
R[, `:=`(d_ACME = ACME - dep_ACME, d_ADE = ADE - dep_ADE,
         d_prop = prop_mediated - dep_prop, d_sobelz = sobel_z - dep_sobel_z,
         n_match = n == dep_n)]
fwrite(R, file.path(BASE, "results/v3/deconv_verify_mediation_agreement.csv"))
logf("\n=== AGREEMENT WITH DEPOSITED MEDIATION (n = ", nrow(R), " models) ===")
logf("max |ACME_mine - ACME_deposited| = ", signif(max(abs(R$d_ACME)), 3))
logf("max |ADE_mine  - ADE_deposited|  = ", signif(max(abs(R$d_ADE)), 3))
logf("max |propmed difference|         = ", signif(max(abs(R$d_prop)), 3))
logf("max |Sobel z difference|         = ", signif(max(abs(R$d_sobelz)), 3))
logf("sample sizes identical in all models: ", all(R$n_match))
logf("sign(ACME) agrees in ", sum(sign(R$ACME) == sign(R$dep_ACME)), "/", nrow(R))
logf("Sobel p < 0.05 in mine: ", sum(R$sobel_p < 0.05), " ; deposited: ", sum(R$dep_sobel_p < 0.05))
logf("inconsistent (sign ACME != sign ADE) mine: ", sum(sign(R$ACME) != sign(R$ADE)),
     " ; deposited: ", sum(DEP$inconsistent))
logf("DONE 41")
