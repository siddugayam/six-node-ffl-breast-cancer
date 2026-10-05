#!/usr/bin/env Rscript
# Random-effects pooling (metafor REML) and meta-regression; a copy of the helper of analyses/perturbation_tests.
# usage: Rscript rma_pool.R <in.tsv> <out.tsv>
# in.tsv columns: group, yi, sei [, mod]   (one row per dataset).  Without 'mod': pooled estimate per group.
# With 'mod' (0/1): per group, rma(yi, sei, mods = ~ mod) -> the coefficient of mod (z-test) and the intercept.
suppressMessages(library(metafor))
a <- commandArgs(TRUE); d <- read.delim(a[1], stringsAsFactors = FALSE)
out <- list()
for (g in unique(d$group)) {
  x <- d[d$group == g & is.finite(d$yi) & is.finite(d$sei) & d$sei > 0, ]
  if (nrow(x) == 0) { out[[g]] <- data.frame(group = g, k = 0, term = "pooled", est = NA, se = NA, ci_lb = NA, ci_ub = NA, p = NA, tau2 = NA, I2 = NA); next }
  if (!"mod" %in% names(x)) {
    f <- rma(yi = x$yi, sei = x$sei, method = "REML")
    out[[g]] <- data.frame(group = g, k = nrow(x), term = "pooled", est = f$b[1], se = f$se, ci_lb = f$ci.lb, ci_ub = f$ci.ub, p = f$pval, tau2 = f$tau2, I2 = f$I2)
  } else {
    if (length(unique(x$mod)) < 2) { out[[g]] <- data.frame(group = g, k = nrow(x), term = "mod", est = NA, se = NA, ci_lb = NA, ci_ub = NA, p = NA, tau2 = NA, I2 = NA); next }
    f <- rma(yi = x$yi, sei = x$sei, mods = ~ x$mod, method = "REML")
    out[[g]] <- data.frame(group = g, k = nrow(x), term = c("intercept", "mod"), est = f$b[, 1], se = f$se, ci_lb = f$ci.lb, ci_ub = f$ci.ub, p = f$pval, tau2 = f$tau2, I2 = f$I2)
  }
}
write.table(do.call(rbind, out), a[2], sep = "\t", quote = FALSE, row.names = FALSE)
