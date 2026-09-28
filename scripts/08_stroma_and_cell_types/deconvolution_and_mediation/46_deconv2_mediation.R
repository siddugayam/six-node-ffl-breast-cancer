#!/usr/bin/env Rscript
## 46_deconv2_mediation.R -- THE DECISIVE RE-TEST, extended.
## Repeat the causal-mediation analysis of the TF -> collagen and miR-29 -| collagen
## edges once per fibroblast/stromal estimate, now including
##   * five new algorithm families run on the same Wu breast single-cell signature
##     (qprog, RLS, OLS, DWLS, dtangle) and their collagen-free counterparts,
##   * two NON-TRANSCRIPTOMIC mediators (ABSOLUTE DNA purity, pathologist slide score),
## on top of the 23 estimates used in the first pass.
## Model: normal-score transform of exposure / mediator / outcome; m ~ x ; y ~ x + m.
## ACME = a*b, ADE = c' (exact for the linear-linear model), percentile bootstrap CIs
## from 2,000 nonparametric resamples, plus the analytic Sobel test.
## Verified in script 41 to reproduce mediation::mediate to 5.6e-16 on the first pass.
suppressPackageStartupMessages({library(data.table); library(parallel)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/46_deconv2_mediation.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
B <- 2000L; CORES <- 8L

gexp <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
mexp <- readRDS(file.path(BASE, "data/brca_mirna_expr_canonical.rds"))
stopifnot(!any(duplicated(rownames(gexp))), !any(duplicated(rownames(mexp))))
logf("RULE-4 guard: duplicated gene rownames = 0, duplicated miRNA rownames = 0; ",
     "gene symbols non-numeric fraction = ", round(mean(!grepl("^[0-9]+$", rownames(gexp))), 5))

FIB <- readRDS(file.path(BASE, "cache/v3/deconv2/fibroblast_estimates_extended.rds"))
logf("mediators available: ", ncol(FIB))
SS <- rownames(FIB)
sB <- sort(intersect(SS, colnames(mexp)))
logf("gene-assay n = ", length(SS), "; gene+miRNA n = ", length(sB))

ns <- function(v) qnorm((rank(v) - 0.5) / length(v))

run1 <- function(xv, mv, yv, seed = 11L) {
  ok <- is.finite(xv) & is.finite(mv) & is.finite(yv)
  xv <- xv[ok]; mv <- mv[ok]; yv <- yv[ok]; n <- length(xv)
  if (n < 50) return(NULL)
  d <- data.frame(x = ns(xv), m = ns(mv), y = ns(yv))
  X1 <- cbind(1, d$x); X2 <- cbind(1, d$x, d$m)
  f <- function(dd) {
    a <- qr.solve(cbind(1, dd$x), dd$m)[2]
    cf <- qr.solve(cbind(1, dd$x, dd$m), dd$y)
    c(a = a, cp = cf[2], b = cf[3])
  }
  p <- f(d)
  ACME <- p["a"] * p["b"]; ADE <- p["cp"]; TOT <- ACME + ADE
  mm <- lm(m ~ x, data = d); om <- lm(y ~ x + m, data = d)
  sa <- summary(mm)$coef["x", 2]; sb <- summary(om)$coef["m", 2]
  z  <- p["a"] * p["b"] / sqrt(p["b"]^2 * sa^2 + p["a"]^2 * sb^2)
  set.seed(seed)
  bs <- matrix(NA_real_, 3, B)
  for (i in seq_len(B)) {
    ii <- sample.int(n, n, replace = TRUE); q <- f(d[ii, ])
    ac <- q["a"] * q["b"]; bs[, i] <- c(ac, q["cp"], ac / (ac + q["cp"]))
  }
  bp <- function(v) min(1, max(2 * min(mean(v <= 0), mean(v >= 0)), 1 / B))
  data.table(n = n, rho_raw = cor(xv, yv, method = "spearman"),
    rho_x_med = cor(xv, mv, method = "spearman"),
    rho_y_med = cor(yv, mv, method = "spearman"),
    rho_partial = cor(residuals(lm(ns(xv) ~ ns(mv))), residuals(lm(ns(yv) ~ ns(mv))), method = "spearman"),
    total = unname(TOT), ACME = unname(ACME),
    acme_lo = quantile(bs[1, ], .025), acme_hi = quantile(bs[1, ], .975), acme_p = bp(bs[1, ]),
    ADE = unname(ADE),
    ade_lo = quantile(bs[2, ], .025), ade_hi = quantile(bs[2, ], .975), ade_p = bp(bs[2, ]),
    prop_mediated = unname(ACME / TOT),
    prop_lo = quantile(bs[3, ], .025, na.rm = TRUE), prop_hi = quantile(bs[3, ], .975, na.rm = TRUE),
    a_path = unname(p["a"]), b_path = unname(p["b"]), cprime = unname(p["cp"]),
    sobel_z = unname(z), sobel_p = unname(2 * pnorm(-abs(z))),
    inconsistent = sign(ACME) != sign(ADE))
}

gene_jobs  <- list(c("ETS1","COL1A1"), c("NFKB1","COL1A1"), c("ETS1","COL3A1"),
                   c("NFKB1","COL3A1"), c("SP1","COL1A1"), c("RELA","COL1A1"))
mirna_jobs <- list(c("hsa-miR-29a","COL1A1"), c("hsa-miR-29a","COL3A1"),
                   c("hsa-miR-29b","COL1A1"), c("hsa-miR-29b","COL3A1"),
                   c("hsa-miR-29c","COL1A1"), c("hsa-miR-29c","COL3A1"))
JOBS <- c(lapply(gene_jobs,  function(j) list(x = j[1], y = j[2], assay = "gene",  ss = SS)),
          lapply(mirna_jobs, function(j) list(x = j[1], y = j[2], assay = "miRNA", ss = sB)))
for (j in JOBS) stopifnot(j$y %in% rownames(gexp))
GRID <- expand.grid(job = seq_along(JOBS), med = colnames(FIB), stringsAsFactors = FALSE)
logf("mediation models to fit: ", nrow(GRID), " (", length(JOBS), " axes x ", ncol(FIB), " mediators)")

t0 <- Sys.time()
out <- mclapply(seq_len(nrow(GRID)), function(i) {
  j <- JOBS[[GRID$job[i]]]; med <- GRID$med[i]
  ss <- j$ss
  xv <- if (j$assay == "miRNA") mexp[j$x, ss] else gexp[j$x, ss]
  r <- run1(as.numeric(xv), as.numeric(FIB[ss, med]), as.numeric(gexp[j$y, ss]))
  if (is.null(r)) return(NULL)
  cbind(data.table(x = j$x, y = j$y, mediator = med), r)
}, mc.cores = CORES)
nbad <- sum(sapply(out, is.null)); logf("models returning NULL: ", nbad)
R <- rbindlist(out[!sapply(out, is.null)])
logf("elapsed ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
fwrite(R, file.path(BASE, "results/v3/deconv_mediation_extended.csv"))
logf("WROTE results/v3/deconv_mediation_extended.csv  rows = ", nrow(R))

## classification of each (axis, mediator) result
R[, verdict := fifelse(sobel_p >= 0.05, "no mediation",
                fifelse(inconsistent, "INCONSISTENT (ACME and ADE opposite sign)",
                 fifelse(abs(ADE) < 0.05 | ade_p >= 0.05, "complete mediation",
                         "partial, same-signed mediation")))]
fwrite(R[, .N, by = .(x, y, verdict)][order(x, y, -N)],
       file.path(BASE, "results/v3/deconv_mediation_extended_verdicts.csv"))

NONRNA <- c("ABSOLUTE_nontumour", "PATH_stroma_frozen", "PATH_stroma_allslides",
            "METH_EpiDISH_Fib", "METH_EpiDISH_Epi", "METH_EpiDISH_IC")
S <- R[, .(n_methods = .N, n_med_RNA = sum(!mediator %in% NONRNA),
           ACME_min = min(ACME), ACME_med = median(ACME), ACME_max = max(ACME),
           ADE_min = min(ADE), ADE_med = median(ADE), ADE_max = max(ADE),
           propmed_min = min(prop_mediated), propmed_med = median(prop_mediated), propmed_max = max(prop_mediated),
           n_sobel_sig = sum(sobel_p < 0.05), n_acme_sig = sum(acme_p < 0.05),
           n_inconsistent = sum(inconsistent),
           n_ADE_null = sum(ade_p >= 0.05),
           n_ADE_sig_same_sign = sum(ade_p < 0.05 & !inconsistent)), by = .(x, y)]
fwrite(S, file.path(BASE, "results/v3/deconv_mediation_extended_summary.csv"))
logf("\n================ METHOD-INDEPENDENCE SUMMARY ================")
for (i in seq_len(nrow(S))) { r <- S[i]
  logf(sprintf("%-12s -> %-6s | ACME %+.3f [%+.3f,%+.3f] | ADE %+.3f [%+.3f,%+.3f] | propmed %+.2f..%+.2f | Sobel sig %d/%d | INCONSISTENT %d/%d",
    r$x, r$y, r$ACME_med, r$ACME_min, r$ACME_max, r$ADE_med, r$ADE_min, r$ADE_max,
    r$propmed_min, r$propmed_max, r$n_sobel_sig, r$n_methods, r$n_inconsistent, r$n_methods)) }

logf("\n================ PER-MEDIATOR DETAIL ================")
for (ax in unique(paste(R$x, R$y, sep = " -> "))) {
  logf("\n--- ", ax, " ---")
  sub <- R[paste(x, y, sep = " -> ") == ax][order(ACME)]
  logf(sprintf("%-34s %5s %7s %8s %8s %8s %9s %10s  %s", "mediator", "n", "rho", "rho|M",
               "ACME", "ADE", "prop_med", "sobel_p", "verdict"))
  for (i in seq_len(nrow(sub))) { r <- sub[i]
    logf(sprintf("%-34s %5d %+7.3f %+8.3f %+8.3f %+8.3f %+9.3f %10.2g  %s",
      r$mediator, r$n, r$rho_raw, r$rho_partial, r$ACME, r$ADE, r$prop_mediated, r$sobel_p, r$verdict)) }
}

## wide tables for the manuscript
for (v in c("ACME", "ADE", "prop_mediated", "sobel_p")) {
  W <- dcast(R, x + y ~ mediator, value.var = v)
  fwrite(W, file.path(BASE, paste0("results/v3/deconv_mediation_extended_", v, "_wide.csv")))
}
logf("DONE 46")
