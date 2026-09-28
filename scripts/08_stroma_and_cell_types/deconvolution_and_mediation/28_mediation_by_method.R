#!/usr/bin/env Rscript
## 28_mediation_by_method.R -- THE DECISIVE RE-TEST.
## Repeat the causal-mediation analysis of the TF -> collagen and miR-29 -| collagen
## edges once per deconvolution method, using THAT METHOD'S fibroblast/stromal
## estimate as the single mediator. Identical model, transform and bootstrap to
## scripts/10_survival_and_clinical/reactive_stroma_and_neoadjuvant/44_farmer_mediation.R, so the numbers are directly comparable to the
## deposited results/v5/farmer_mediation.csv.
##   normal-score transform of exposure, mediator and outcome
##   mediator model  m ~ x        (OLS)
##   outcome  model  y ~ x + m    (OLS)
##   mediation::mediate(boot = TRUE, sims = 1000, percentile CI) + analytic Sobel test
suppressPackageStartupMessages({library(mediation); library(parallel); library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/28_mediation.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
SIMS <- 1000L; CORES <- 8L
logf("mediation bootstrap sims = ", SIMS, " (boot=TRUE, percentile CI); cores = ", CORES)

gexp <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
mexp <- readRDS(file.path(BASE, "data/brca_mirna_expr_canonical.rds"))
FIB  <- readRDS(file.path(BASE, "cache/v3/deconv/fibroblast_estimates_all.rds"))
sG <- rownames(FIB)                       # 1,097 primary tumours with a gene assay
sB <- sort(intersect(sG, colnames(mexp))) # subset that also has a miRNA assay
logf("gene-assay tumours n = ", length(sG), " ; both-assay tumours n = ", length(sB))
logf("mediators (fibroblast/stromal estimates): ", ncol(FIB))

ns <- function(v) { n <- length(v); qnorm((rank(v) - 0.5) / n) }

run_med <- function(xname, yname, medname, ss, xassay) {
  X <- if (xassay == "miRNA") mexp[xname, ss] else gexp[xname, ss]
  Y <- gexp[yname, ss]
  M <- FIB[ss, medname]
  d <- data.frame(x = ns(X), y = ns(Y), m = ns(M))
  mm <- lm(m ~ x, data = d); om <- lm(y ~ x + m, data = d)
  set.seed(7)
  fit <- mediate(mm, om, treat = "x", mediator = "m", boot = TRUE, sims = SIMS, boot.ci.type = "perc")
  a <- coef(mm)["x"]; b <- coef(om)["m"]
  sa <- summary(mm)$coef["x", 2]; sb <- summary(om)$coef["m", 2]
  z  <- a * b / sqrt(b^2 * sa^2 + a^2 * sb^2)
  data.frame(x = xname, y = yname, mediator = medname, n = length(ss),
    rho_raw = cor(X, Y, method = "spearman"),
    rho_x_mediator = cor(X, M, method = "spearman"),
    rho_y_mediator = cor(Y, M, method = "spearman"),
    rho_partial_given_mediator = cor(residuals(lm(ns(X) ~ ns(M))), residuals(lm(ns(Y) ~ ns(M))), method = "spearman"),
    total = fit$tau.coef, total_p = fit$tau.p,
    ACME = fit$d0, acme_lo = fit$d0.ci[1], acme_hi = fit$d0.ci[2], acme_p = fit$d0.p,
    ADE = fit$z0, ade_lo = fit$z0.ci[1], ade_hi = fit$z0.ci[2], ade_p = fit$z0.p,
    prop_mediated = fit$n0, prop_lo = fit$n0.ci[1], prop_hi = fit$n0.ci[2], prop_p = fit$n0.p,
    a_path = a, b_path = b, cprime = coef(om)["x"], sobel_z = z, sobel_p = 2 * pnorm(-abs(z)),
    inconsistent = sign(fit$d0) != sign(fit$z0),
    stringsAsFactors = FALSE)
}

gene_jobs  <- list(c("ETS1","COL1A1"), c("NFKB1","COL1A1"), c("ETS1","COL3A1"),
                   c("NFKB1","COL3A1"), c("SP1","COL1A1"), c("RELA","COL1A1"))
mirna_jobs <- list(c("hsa-miR-29a","COL1A1"), c("hsa-miR-29a","COL3A1"),
                   c("hsa-miR-29b","COL1A1"), c("hsa-miR-29b","COL3A1"),
                   c("hsa-miR-29c","COL1A1"), c("hsa-miR-29c","COL3A1"))
JOBS <- c(lapply(gene_jobs, function(j) list(x = j[1], y = j[2], assay = "gene",  ss = sG)),
          lapply(mirna_jobs, function(j) list(x = j[1], y = j[2], assay = "miRNA", ss = sB)))
GRID <- expand.grid(job = seq_along(JOBS), med = colnames(FIB), stringsAsFactors = FALSE)
logf("total mediation models to fit: ", nrow(GRID))

t0 <- Sys.time()
out <- mclapply(seq_len(nrow(GRID)), function(i) {
  j <- JOBS[[GRID$job[i]]]
  tryCatch(run_med(j$x, j$y, GRID$med[i], j$ss, j$assay),
           error = function(e) { message("FAIL ", j$x, "->", j$y, " / ", GRID$med[i], ": ", conditionMessage(e)); NULL })
}, mc.cores = CORES)
bad <- sum(sapply(out, is.null)); logf("failed models: ", bad)
R <- rbindlist(out[!sapply(out, is.null)])
logf("elapsed ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
fwrite(R, file.path(BASE, "results/v3/deconv_mediation_by_method.csv"))
logf("WROTE results/v3/deconv_mediation_by_method.csv rows = ", nrow(R))

## ------------------------------------------------------------------ readable log ----
for (ax in unique(paste(R$x, R$y, sep = " -> "))) {
  logf("\n=== ", ax, " ===")
  sub <- R[paste(x, y, sep = " -> ") == ax][order(-ACME)]
  logf(sprintf("%-32s %7s %7s %8s %8s %9s %9s %10s %s", "mediator", "rho", "rho|M", "total",
               "ACME", "ADE", "prop_med", "sobel_p", "inconsistent"))
  for (i in seq_len(nrow(sub))) { r <- sub[i]
    logf(sprintf("%-32s %+7.3f %+7.3f %+8.3f %+8.3f %+9.3f %+9.3f %10.3g %s",
      r$mediator, r$rho_raw, r$rho_partial_given_mediator, r$total, r$ACME, r$ADE,
      r$prop_mediated, r$sobel_p, ifelse(r$inconsistent, "YES", "."))) }
}

## ---------------------------------------------------- method-independence summary ----
S <- R[, .(n_methods = .N,
           ACME_min = min(ACME), ACME_max = max(ACME), ACME_median = median(ACME),
           ADE_min = min(ADE), ADE_max = max(ADE), ADE_median = median(ADE),
           n_ACME_sig = sum(acme_p < 0.05), n_sobel_sig = sum(sobel_p < 0.05),
           n_inconsistent = sum(inconsistent),
           n_ACME_null_and_ADE_sig = sum(acme_p >= 0.05 & ade_p < 0.05),
           propmed_min = min(prop_mediated), propmed_max = max(prop_mediated)),
         by = .(x, y)]
fwrite(S, file.path(BASE, "results/v3/deconv_mediation_summary.csv"))
logf("\n=== METHOD-INDEPENDENCE SUMMARY (across all fibroblast estimates) ===")
for (i in seq_len(nrow(S))) { r <- S[i]
  logf(sprintf("%-12s -> %-7s  ACME %+.3f [%+.3f,%+.3f]  ADE %+.3f [%+.3f,%+.3f]  ACME sig %d/%d  Sobel sig %d/%d  INCONSISTENT %d/%d",
    r$x, r$y, r$ACME_median, r$ACME_min, r$ACME_max, r$ADE_median, r$ADE_min, r$ADE_max,
    r$n_ACME_sig, r$n_methods, r$n_sobel_sig, r$n_methods, r$n_inconsistent, r$n_methods)) }
W <- dcast(R, x + y ~ mediator, value.var = "prop_mediated")
fwrite(W, file.path(BASE, "results/v3/deconv_mediation_propmediated_wide.csv"))
W2 <- dcast(R, x + y ~ mediator, value.var = "ACME"); fwrite(W2, file.path(BASE, "results/v3/deconv_mediation_ACME_wide.csv"))
W3 <- dcast(R, x + y ~ mediator, value.var = "ADE");  fwrite(W3, file.path(BASE, "results/v3/deconv_mediation_ADE_wide.csv"))
logf("DONE 28")
