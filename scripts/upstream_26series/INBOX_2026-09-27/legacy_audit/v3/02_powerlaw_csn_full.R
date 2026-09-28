#!/usr/bin/env Rscript
# (A) Is the degree distribution scale-free? Clauset, Shalizi & Newman (2009)
# SIAM Review 51:661-703, implemented in poweRlaw (Gillespie 2015 J Stat Softw 64:1).
#
# Procedure, exactly as CSN prescribe:
#   1. MLE of alpha for every candidate x_min; choose x_min minimising the KS distance.
#   2. Goodness of fit: bootstrap p-value from 1,000 synthetic data sets drawn from the
#      fitted power law below/above x_min. p <= 0.1 => power law is RULED OUT.
#   3. Likelihood-ratio (Vuong) tests of the power law against discrete log-normal,
#      exponential and Poisson alternatives, all fitted above the SAME x_min.
suppressPackageStartupMessages({library(poweRlaw)})
set.seed(20260908)
REV <- "/path/to/revision"
RES <- "/path/to/revision/INBOX_2026-09-27/legacy_audit/v3/out_full_hash0"
FIG <- "/path/to/revision/INBOX_2026-09-27/legacy_audit/v3/out_full_hash0/fig"
NBOOT <- as.integer(Sys.getenv("NBOOT", "1000"))
NCORE <- as.integer(Sys.getenv("NCORE", "4"))

dg <- read.csv(file.path(RES, "systems_degree_table.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(dg) == 587)
cat("degree table rows:", nrow(dg), "\n")

fit_one <- function(x, label) {
  x <- x[x > 0]
  cat("\n=====", label, " n =", length(x), " max =", max(x), "\n")
  m_pl <- displ$new(x)
  est <- estimate_xmin(m_pl, xmax = max(x))
  m_pl$setXmin(est)
  alpha <- m_pl$pars; xmin <- m_pl$xmin
  ntail <- sum(x >= xmin)
  cat(sprintf("  power law: alpha = %.4f  xmin = %d  ntail = %d  KS = %.4f\n",
              alpha, xmin, ntail, est$gof))
  bs <- bootstrap_p(m_pl, no_of_sims = NBOOT, threads = NCORE, seed = 1)
  p_gof <- bs$p
  cat(sprintf("  CSN bootstrap goodness-of-fit p = %.4f  (%d sims)\n", p_gof, NBOOT))

  cmpr <- function(alt_class, altname) {
    m_alt <- alt_class$new(x); m_alt$setXmin(m_pl$getXmin())
    m_alt$setPars(estimate_pars(m_alt))
    cp <- compare_distributions(m_pl, m_alt)
    cat(sprintf("  vs %-12s  LR = %+8.3f  p_two = %.4g  (%s)\n", altname,
                cp$test_statistic, cp$p_two_sided,
                if (cp$p_two_sided > 0.1) "indistinguishable"
                else if (cp$test_statistic > 0) "power law favoured" else "alternative favoured"))
    data.frame(alternative = altname, LR = cp$test_statistic,
               p_one_sided = cp$p_one_sided, p_two_sided = cp$p_two_sided)
  }
  cmps <- rbind(cmpr(dislnorm, "lognormal"), cmpr(disexp, "exponential"),
                cmpr(dispois, "poisson"))
  cmps$distribution <- label
  list(summary = data.frame(distribution = label, n = length(x), n_tail = ntail,
                            alpha = alpha, xmin = xmin, KS = est$gof,
                            bootstrap_p = p_gof, nboot = NBOOT,
                            scale_free_supported = p_gof > 0.1),
       cmp = cmps, model = m_pl, x = x)
}

res <- list()
res$total <- fit_one(dg$degree_total, "degree_total")
res$out   <- fit_one(dg$degree_out,   "degree_out")
res$inn   <- fit_one(dg$degree_in,    "degree_in")

summ <- do.call(rbind, lapply(res, `[[`, "summary"))
cmp  <- do.call(rbind, lapply(res, `[[`, "cmp"))
write.csv(summ, file.path(RES, "systems_powerlaw_csn_fits.csv"), row.names = FALSE)
write.csv(cmp,  file.path(RES, "systems_powerlaw_csn_lrt.csv"),  row.names = FALSE)
cat("\n--- SUMMARY ---\n"); print(summ); cat("\n"); print(cmp)

pdf(file.path(FIG, "Fig_v3_A_powerlaw_ccdf.pdf"), width = 11, height = 4)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1))
for (nm in names(res)) {
  m <- res[[nm]]$model
  plot(m, pch = 16, cex = .6, xlab = "degree k", ylab = "P(K >= k)",
       main = sprintf("%s\nalpha=%.2f xmin=%d p=%.3f", res[[nm]]$summary$distribution,
                      res[[nm]]$summary$alpha, res[[nm]]$summary$xmin,
                      res[[nm]]$summary$bootstrap_p))
  lines(m, col = "red", lwd = 2)
  m_ln <- dislnorm$new(res[[nm]]$x); m_ln$setXmin(m$getXmin()); m_ln$setPars(estimate_pars(m_ln))
  lines(m_ln, col = "blue", lwd = 2, lty = 2)
  legend("bottomleft", c("power law", "lognormal"), col = c("red", "blue"),
         lty = c(1, 2), lwd = 2, bty = "n", cex = .8)
}
dev.off()
png(file.path(FIG, "Fig_v3_A_powerlaw_ccdf.png"), width = 1650, height = 600, res = 150)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1))
for (nm in names(res)) {
  m <- res[[nm]]$model
  plot(m, pch = 16, cex = .6, xlab = "degree k", ylab = "P(K >= k)",
       main = sprintf("%s\nalpha=%.2f xmin=%d p=%.3f", res[[nm]]$summary$distribution,
                      res[[nm]]$summary$alpha, res[[nm]]$summary$xmin,
                      res[[nm]]$summary$bootstrap_p))
  lines(m, col = "red", lwd = 2)
  m_ln <- dislnorm$new(res[[nm]]$x); m_ln$setXmin(m$getXmin()); m_ln$setPars(estimate_pars(m_ln))
  lines(m_ln, col = "blue", lwd = 2, lty = 2)
}
dev.off()
cat("figures written\n")
