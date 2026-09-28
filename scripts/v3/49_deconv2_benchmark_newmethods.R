#!/usr/bin/env Rscript
## 49_deconv2_benchmark_newmethods.R
## An INDEPENDENT accuracy benchmark that includes the new algorithm families.
## Pseudo-bulk mixtures of known composition are built from the Wu et al. 2021 atlas
## with a PATIENT-level hold-out: the signature matrix and the InstaPrism reference are
## built from one half of the patients, the mixtures are simulated from the other half.
## Caveat (same as the first pass): simulated mixtures carry no bulk technical noise,
## no ambient RNA and no cell-type mRNA-content differences, so absolute accuracy is an
## upper bound and only the ranking between methods should be used.
suppressPackageStartupMessages({library(data.table); library(nnls); library(MASS);
                                library(limSolve); library(quadprog); library(matrixStats)})
.libPaths(c(path.expand("~/Rlib_deconv2"), .libPaths()))
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/49_deconv2_benchmark.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
D <- file.path(BASE, "cache/v3/deconv")
set.seed(20260910)

pc <- as.matrix(read.delim(file.path(D, "wu2021_sum_umi_by_patient_celltype.tsv"),
                           row.names = 1, check.names = FALSE))
bulkgenes <- rownames(readRDS(file.path(D, "bulk_log2.rds")))
pat <- sub("\\|.*$", "", colnames(pc)); ct <- sub("^.*\\|", "", colnames(pc))
cts <- sort(unique(ct)); pats <- sort(unique(pat))
logf("pseudobulk profiles: ", ncol(pc), " from ", length(pats), " patients, ", length(cts), " cell types")
trainP <- sort(sample(pats, floor(length(pats) / 2))); testP <- setdiff(pats, trainP)
logf("train patients ", length(trainP), " / test patients ", length(testP))
cpm <- sweep(pc, 2, colSums(pc), "/") * 1e6

NSIM <- 300
P <- matrix(rgamma(NSIM * length(cts), 0.7, 1), NSIM, length(cts)); P <- P / rowSums(P)
colnames(P) <- cts
MIX <- matrix(0, nrow(cpm), NSIM, dimnames = list(rownames(cpm), paste0("sim", seq_len(NSIM))))
for (i in seq_len(NSIM)) for (k in seq_along(cts)) {
  cand <- which(pat %in% testP & ct == cts[k])
  if (!length(cand)) { P[i, k] <- 0; next }
  MIX[, i] <- MIX[, i] + P[i, k] * cpm[, sample(cand, 1)]
}
P <- P / rowSums(P)
MIX <- sweep(MIX, 2, colSums(MIX), "/") * 1e6
MIX <- MIX[intersect(rownames(MIX), bulkgenes), , drop = FALSE]
LMIX <- log2(MIX + 1)
logf("simulated ", NSIM, " mixtures on ", nrow(MIX), " genes; true CAF fraction median ",
     round(median(P[, "CAFs"]), 3))

## ------ train-half reference and signature matrix (marker genes by fold-change) ------
trmask <- pat %in% trainP
REF <- sapply(cts, function(k) rowMeans(cpm[, trmask & ct == k, drop = FALSE]))
REF <- REF[intersect(rownames(REF), rownames(MIX)), , drop = FALSE]
logf("train reference: ", nrow(REF), " genes x ", ncol(REF), " types")
pick <- unique(unlist(lapply(seq_len(ncol(REF)), function(k) {
  v <- REF[, k]; o <- rowMaxs(REF[, -k, drop = FALSE])
  fc <- log2((v + 1) / (o + 1))
  names(sort(fc[v > 20], decreasing = TRUE))[1:80]
})))
pick <- pick[!is.na(pick)]
SIG <- REF[pick, , drop = FALSE]
logf("train signature: ", nrow(SIG), " marker genes; COL* among them: ",
     paste(grep("^COL", rownames(SIG), value = TRUE), collapse = ", "))

S <- SIG; Y <- MIX[rownames(SIG), ]; K <- ncol(S)
norm1 <- function(w) { w[!is.finite(w) | w < 0] <- 0; if (sum(w) <= 0) rep(NA_real_, length(w)) else w / sum(w) }

EST <- list()
EST$NNLS  <- t(apply(Y, 2, function(b) norm1(nnls::nnls(S, b)$x)))
EST$OLS   <- t(apply(Y, 2, function(b) norm1(qr.solve(S, b))))
EST$RLS   <- t(apply(Y, 2, function(b) norm1(tryCatch(coef(MASS::rlm(S, b, maxit = 100)),
                                                      error = function(e) rep(NA_real_, K)))))
sc <- max(S)
EST$QPROG <- t(apply(Y / sc, 2, function(b) norm1(tryCatch(
  limSolve::lsei(A = S / sc, B = b, E = matrix(1, 1, K), F = 1, G = diag(K), H = rep(0, K), type = 2)$X,
  error = function(e) rep(NA_real_, K)))))

## DWLS (same implementation as script 43)
solveOLSint <- function(S, b) { D <- t(S) %*% S; d <- t(S) %*% b; s <- norm(D, "2")
  quadprog::solve.QP(D / s, d / s, cbind(diag(ncol(S))), rep(0, ncol(S)))$solution }
solveDWLSj <- function(S, b, sol, j) { mult <- 2^(j - 1)
  w <- as.vector((1 / (S %*% sol))^2); w[!is.finite(w)] <- max(w[is.finite(w)])
  w <- w / min(w); w[w > mult] <- mult
  SW <- S * sqrt(w); bW <- b * sqrt(w); D <- t(SW) %*% SW; d <- t(SW) %*% bW; s <- norm(D, "2")
  quadprog::solve.QP(D / s, d / s, cbind(diag(ncol(S))), rep(0, ncol(S)))$solution }
EST$DWLS <- t(apply(Y, 2, function(b) {
  sol <- tryCatch(solveOLSint(S, b), error = function(e) rep(1 / K, K)); sol[sol <= 0] <- 1e-8
  it <- 0; ch <- 1
  while (ch > 0.01 && it < 100) { new <- tryCatch(solveDWLSj(S, b, sol, 4L), error = function(e) sol)
    new[new <= 0] <- 1e-10; avg <- rowMeans(cbind(new, sol, sol, sol))
    ch <- sqrt(sum((avg - sol)^2)); sol <- avg; it <- it + 1 }
  norm1(sol) }))

## dtangle against the train-half reference
{
  library(dtangle)
  lref <- log2(REF + 1); gg <- rownames(lref)
  Yd <- t(LMIX[gg, , drop = FALSE]); Rd <- t(lref)
  ps <- lapply(seq_len(nrow(Rd)), function(i) nrow(Yd) + i)
  mk <- dtangle::find_markers(Y = rbind(Yd, Rd), pure_samples = ps, data_type = "rna-seq",
                              marker_method = "ratio")
  nmk <- unlist(lapply(mk$L, function(v) max(5L, floor(length(v) * 0.05))))
  dd <- dtangle::dtangle(Y = rbind(Yd, Rd), pure_samples = ps, n_markers = nmk, markers = mk$L,
                         data_type = "rna-seq")
  M <- dd$estimates[seq_len(nrow(Yd)), , drop = FALSE]; colnames(M) <- colnames(REF)
  EST$DTANGLE <- M
}
## MCP-counter and EPIC and InstaPrism, for continuity with the first-pass benchmark
EST$MCPcounter <- tryCatch({ library(MCPcounter)
  t(MCPcounter::MCPcounter.estimate(LMIX, featuresType = "HUGO_symbols")) }, error = function(e) NULL)
EST$EPIC <- tryCatch({ library(EPIC); EPIC::EPIC(bulk = MIX, reference = "TRef",
  withOtherCells = TRUE)$cellFractions }, error = function(e) NULL)
EST$InstaPrism <- tryCatch({ library(InstaPrism)
  gi <- intersect(rownames(pc), rownames(MIX))
  keep <- trmask
  rp <- InstaPrism::refPrepare(sc_Expr = pc[gi, keep], cell.type.labels = ct[keep],
                               cell.state.labels = colnames(pc)[keep])
  ip <- InstaPrism::InstaPrism(bulk_Expr = MIX[gi, ], refPhi_cs = rp, n.iter = 100,
                               n.core = 6, verbose = FALSE)
  t(ip@Post.ini.ct@theta) }, error = function(e) { logf("InstaPrism failed: ", conditionMessage(e)); NULL })

trueCAF <- P[, "CAFs"]
res <- rbindlist(lapply(names(EST), function(nm) {
  M <- EST[[nm]]; if (is.null(M)) return(NULL)
  cn <- colnames(M)
  k <- which(grepl("CAF|[Ff]ibroblast", cn))[1]
  if (is.na(k)) return(NULL)
  v <- M[, k]
  isfrac <- !is.null(cn) && all(abs(rowSums(M, na.rm = TRUE) - 1) < 0.05, na.rm = TRUE)
  data.table(method = nm, feature = cn[k],
             spearman_vs_true_CAF = cor(v, trueCAF, method = "spearman", use = "complete.obs"),
             pearson_vs_true_CAF  = cor(v, trueCAF, use = "complete.obs"),
             rmse_if_fraction = if (isfrac) sqrt(mean((v - trueCAF)^2, na.rm = TRUE)) else NA_real_,
             bias_if_fraction = if (isfrac) mean(v - trueCAF, na.rm = TRUE) else NA_real_,
             spearman_vs_true_CancerEpi = cor(v, P[, "Cancer Epithelial"], method = "spearman", use = "complete.obs"),
             spearman_vs_true_PVL = cor(v, P[, "PVL"], method = "spearman", use = "complete.obs"),
             n_failed = sum(!is.finite(v)))
}))[order(-spearman_vs_true_CAF)]
fwrite(res, file.path(BASE, "results/v3/deconv_benchmark_newmethods.csv"))
logf("\n=== ACCURACY OF THE CAF ESTIMATE ON HELD-OUT SIMULATED MIXTURES (n = ", NSIM, ") ===")
for (i in seq_len(nrow(res))) { r <- res[i]
  logf(sprintf("%-12s rho %+0.3f  r %+0.3f  RMSE %s  bias %s  (rho vs true CancerEpi %+0.3f)",
    r$method, r$spearman_vs_true_CAF, r$pearson_vs_true_CAF,
    ifelse(is.na(r$rmse_if_fraction), "  -  ", sprintf("%.4f", r$rmse_if_fraction)),
    ifelse(is.na(r$bias_if_fraction), "  -  ", sprintf("%+.4f", r$bias_if_fraction)),
    r$spearman_vs_true_CancerEpi)) }
saveRDS(list(P = P, EST = EST), file.path(BASE, "cache/v3/deconv2/sim_benchmark2.rds"))
logf("DONE 49")
