# Cross-implementation check of the Clauset-Shalizi-Newman fit: is the R/poweRlaw result
# (alpha=3.46, xmin=33) and the python/powerlaw result (alpha=2.94, xmin=24) the same fit
# evaluated at a different x_min, or a real disagreement?
suppressPackageStartupMessages(library(poweRlaw))
set.seed(20260908)
RES <- "/path/to/revision/results/v3"
dg <- read.csv(file.path(RES, "systems_degree_table.csv"))
x <- dg$degree_total; x <- x[x > 0]
cat("n =", length(x), " max =", max(x), "\n")
out <- list()
for (xm in c(24, 33)) {
  m <- displ$new(x); m$setXmin(xm); m$setPars(estimate_pars(m))
  g <- get_KS_statistic(m)
  bs <- bootstrap_p(m, no_of_sims = 500, threads = 4, seed = 2)
  ln <- dislnorm$new(x); ln$setXmin(xm); ln$setPars(estimate_pars(ln))
  cp <- compare_distributions(m, ln)
  cat(sprintf("xmin=%d  alpha=%.4f  ntail=%d  KS=%.4f  bootstrap_p=%.3f  LR_vs_lognormal=%+.3f p=%.4f\n",
              xm, m$pars, sum(x >= xm), g, bs$p, cp$test_statistic, cp$p_two_sided))
  out[[as.character(xm)]] <- data.frame(xmin = xm, alpha = m$pars, n_tail = sum(x >= xm),
    KS = g, bootstrap_p = bs$p, LR_vs_lognormal = cp$test_statistic, p_LR = cp$p_two_sided)
}
# and the free-xmin fit WITHOUT the xmax restriction used in the original run
m <- displ$new(x); est <- estimate_xmin(m); m$setXmin(est)
cat(sprintf("free xmin (no xmax argument): alpha=%.4f xmin=%d KS=%.4f\n", m$pars, m$xmin, est$gof))
res <- do.call(rbind, out)
res$free_xmin_alpha <- m$pars; res$free_xmin <- m$xmin
write.csv(res, file.path(RES, "systems_powerlaw_crosscheck.csv"), row.names = FALSE)
print(res)
