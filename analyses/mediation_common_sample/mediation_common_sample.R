## mediation_common_sample.R -- the stromal mediation of the twelve paths of 46 (ETS1, NFKB1 and hsa-miR-29a/b/c, each
## with COL1A1 and COL3A1; SP1 and RELA with COL1A1; Table 3 shows seven of them), repeated for each of the 42 distinct
## stromal estimates on one common sample of tumours
## (both assays and every estimate available), so that each path's total effect is the same for every estimate
## (see SETTINGS.md). The model function run1() is copied unchanged from
## scripts/08_stroma_and_cell_types/deconvolution_and_mediation/46_deconv_extended_mediation.R; the script first checks
## that it reproduces that script's published rows on each estimate's own tumours.
## Root of the repository: environment variable FFL_ROOT. No sample-level values are written.
## usage: FFL_ROOT=<repository checkout> Rscript mediation_common_sample.R
suppressPackageStartupMessages({library(data.table); library(parallel)})
BASE <- Sys.getenv("FFL_ROOT"); stopifnot(nzchar(BASE))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
B <- 2000L; CORES <- 4L

gexp <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
mexp <- readRDS(file.path(BASE, "data/brca_mirna_expr_canonical.rds"))
## the matrix that 46 reads; its CSV export in results/ is rounded to 15 significant digits, which adds ties in the
## MCPcounter columns and changes their ranks
FIB <- readRDS(file.path(BASE, "cache/systems_and_dynamics/deconv_extended/fibroblast_estimates_extended.rds"))
SS <- rownames(FIB); sB <- sort(intersect(SS, colnames(mexp)))

## ---- run1(), verbatim from 46 ------------------------------------------------------------------------------------
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

## ---- estimates: 42 distinct ---------------------------------------------------------------------------------------
DROP <- c("MCPcounter_Fibroblasts_recomputed", "NNLS_Wu641_full_CAFs", "QPROG_Wu641_full_CAFs")
for (d in DROP) { same <- setdiff(colnames(FIB)[apply(FIB, 2, function(v) isTRUE(all.equal(v, FIB[, d], tolerance = 0)))], d)
  cat(d, "equals", paste(same, collapse = ";"), "\n"); stopifnot(length(same) >= 1) }
EST <- setdiff(colnames(FIB), DROP); stopifnot(length(EST) == 42)

## the twelve paths of 46, in its order
ARMS <- c(lapply(list(c("ETS1","COL1A1"), c("NFKB1","COL1A1"), c("ETS1","COL3A1"), c("NFKB1","COL3A1"), c("SP1","COL1A1"),
                      c("RELA","COL1A1")), function(j) list(x = j[1], y = j[2], assay = "gene")),
          lapply(list(c("hsa-miR-29a","COL1A1"), c("hsa-miR-29a","COL3A1"), c("hsa-miR-29b","COL1A1"), c("hsa-miR-29b","COL3A1"),
                      c("hsa-miR-29c","COL1A1"), c("hsa-miR-29c","COL3A1")), function(j) list(x = j[1], y = j[2], assay = "miRNA")))
TABLE3 <- data.table(x = c("ETS1","NFKB1","SP1","RELA","hsa-miR-29a","hsa-miR-29a","hsa-miR-29b"),
                     y = c("COL1A1","COL1A1","COL1A1","COL1A1","COL3A1","COL1A1","COL3A1"), withexemplar_inconsistent = c(23, 32, 38, 22, 6, 6, 3))
xval <- function(a, ss) as.numeric(if (a$assay == "miRNA") mexp[a$x, ss] else gexp[a$x, ss])

## ---- control: published rows on each estimate's own tumours (46 used all gene-assay tumours for the TF paths) ----
CTRL <- rbindlist(mclapply(seq_len(length(ARMS) * length(EST)), function(i) {
  a <- ARMS[[(i - 1) %/% length(EST) + 1]]; med <- EST[(i - 1) %% length(EST) + 1]
  ss <- if (a$assay == "miRNA") sB else SS
  r <- run1(xval(a, ss), as.numeric(FIB[ss, med]), as.numeric(gexp[a$y, ss]))
  cbind(data.table(x = a$x, y = a$y, mediator = med), r[, .(n, total, ACME, acme_lo, acme_hi, ADE, ade_lo, ade_hi)])
}, mc.cores = CORES))
PUB <- fread(file.path(BASE, "results/systems_and_dynamics/deconv_mediation_extended.csv"))
J <- merge(CTRL, PUB[, .(x, y, mediator, n_p = n, total_p = total, ACME_p = ACME, acme_lo_p = acme_lo, acme_hi_p = acme_hi,
                         ADE_p = ADE, ade_lo_p = ade_lo, ade_hi_p = ade_hi)], by = c("x", "y", "mediator"))
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12))
ctrl <- data.table(models_compared = nrow(J), models_expected = nrow(CTRL), n_equal = all(J$n == J$n_p),
                   max_rel_total = rel(J$total, J$total_p), max_rel_ACME = rel(J$ACME, J$ACME_p), max_rel_ADE = rel(J$ADE, J$ADE_p),
                   max_rel_CI = max(rel(J$acme_lo, J$acme_lo_p), rel(J$acme_hi, J$acme_hi_p), rel(J$ade_lo, J$ade_lo_p), rel(J$ade_hi, J$ade_hi_p)))
ctrl[, reproduced := models_compared == models_expected & n_equal & max_rel_total < 1e-9 & max_rel_ACME < 1e-9 & max_rel_ADE < 1e-9 & max_rel_CI < 1e-9]
fwrite(ctrl, file.path(HERE, "control_reproduction.csv")); print(ctrl); stopifnot(ctrl$reproduced)

## ---- common sample ------------------------------------------------------------------------------------------------
okx <- sapply(ARMS, function(a) is.finite(xval(a, sB)))
common <- sB[apply(is.finite(FIB[sB, EST, drop = FALSE]), 1, all) & apply(okx, 1, all) & is.finite(gexp["COL1A1", sB]) &
              is.finite(gexp["COL3A1", sB])]
cat("common sample: tumours with both assays and all 42 estimates:", length(common), "\n")
R <- rbindlist(mclapply(seq_len(length(ARMS) * length(EST)), function(i) {
  a <- ARMS[[(i - 1) %/% length(EST) + 1]]; med <- EST[(i - 1) %% length(EST) + 1]
  r <- run1(xval(a, common), as.numeric(FIB[common, med]), as.numeric(gexp[a$y, common]))
  cbind(data.table(x = a$x, y = a$y, mediator = med), r)
}, mc.cores = CORES))
stopifnot(all(R$n == length(common)))
tot <- R[, .(spread = max(total) - min(total)), by = .(x, y)]; print(tot); stopifnot(all(tot$spread < 1e-12))
R[, ADE_CI_excludes_zero := ade_lo > 0 | ade_hi < 0]
fwrite(R, file.path(HERE, "mediation_common_sample_per_estimate.csv"))
## published counts on all available tumours, for the same 42 estimates
P42 <- PUB[mediator %in% EST]
P42[, ADE_CI_excludes_zero := ade_lo > 0 | ade_hi < 0]
S <- merge(
  R[, .(n_common = n[1], total_common = total[1], ACME_median = median(ACME), ACME_min = min(ACME), ACME_max = max(ACME),
        ADE_median = median(ADE), ADE_min = min(ADE), ADE_max = max(ADE), propmed_median = median(prop_mediated),
        propmed_min = min(prop_mediated), propmed_max = max(prop_mediated), n_estimates = .N, n_inconsistent = sum(inconsistent),
        n_inconsistent_ADE_CI_excl0 = sum(inconsistent & ADE_CI_excludes_zero), n_sobel_p05 = sum(sobel_p < 0.05)), by = .(x, y)],
  P42[, .(published_n_range = paste(min(n), max(n), sep = "-"), published_total_range = paste(signif(min(total), 4), signif(max(total), 4), sep = " to "),
          published_ACME_median = median(ACME), published_ADE_median = median(ADE), published_propmed_median = median(prop_mediated),
          published_n_inconsistent = sum(inconsistent), published_n_inconsistent_ADE_CI_excl0 = sum(inconsistent & ADE_CI_excludes_zero)), by = .(x, y)],
  by = c("x", "y"))
S <- S[match(paste(sapply(ARMS, `[[`, "x"), sapply(ARMS, `[[`, "y")), paste(S$x, S$y))]
S <- merge(S, TABLE3, by = c("x", "y"), all.x = TRUE, sort = FALSE)
S[, in_withexemplar := !is.na(withexemplar_inconsistent)]
chk3 <- S[in_withexemplar == TRUE, all(published_n_inconsistent == withexemplar_inconsistent)]
cat("published inconsistent counts equal Table 3 for its seven rows:", chk3, "\n"); stopifnot(chk3)
fwrite(S, file.path(HERE, "mediation_common_sample_summary.csv")); print(S)
cat("DONE\n")
