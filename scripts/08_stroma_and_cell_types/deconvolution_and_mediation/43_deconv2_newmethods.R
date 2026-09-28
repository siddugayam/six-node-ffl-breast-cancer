#!/usr/bin/env Rscript
## 43_deconv2_newmethods.R -- ADD independent deconvolution methods that were not in the
## first pass, so that the mediation re-test spans several algorithm FAMILIES and not
## only several signature sets.
##
## Families now represented:
##   marker/enrichment  : MCP-counter, xCell, ssGSEA, ConsensusTME, ESTIMATE   (already run)
##   nu-SVR             : CIBERSORT LM22, CIBERSORTx-style on the Wu atlas      (already run)
##   constrained LS     : EPIC, quanTIseq, NNLS                                 (already run)
##   Bayesian           : InstaPrism / BayesPrism                               (already run)
##   NEW: dtangle (linear log-scale marker model, Hunt et al. 2019)
##   NEW: DWLS  (dampened weighted least squares, Tsoucas et al. 2019)          -- own impl.
##   NEW: qprog (equality-constrained quadratic programming, sum-to-1)          -- limSolve
##   NEW: RLS   (robust iteratively-reweighted least squares, MASS::rlm)
##   NEW: TIMER and ABIS via the immunedeconv wrapper
##   NEW: immunedeconv-wrapper reruns of EPIC / MCP-counter / quanTIseq as a
##        cross-implementation control on the first pass
## Every method that can see fibroblasts is run against the SAME Wu et al. 2021 breast
## single-cell signature matrix, so the algorithm is the only thing that differs.
suppressPackageStartupMessages({library(data.table); library(parallel); library(MASS)})
.libPaths(c(path.expand("~/Rlib_deconv2"), .libPaths()))
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/43_deconv2_newmethods.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
OUT  <- file.path(BASE, "cache/v3/deconv2/out"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CORES <- 8L
ONLY <- Sys.getenv("ONLY", "")

LOG2 <- readRDS(file.path(BASE, "cache/v3/deconv/bulk_log2.rds"))
TPM  <- readRDS(file.path(BASE, "cache/v3/deconv/bulk_tpm.rds"))
SS   <- colnames(LOG2)
logf("bulk: ", nrow(LOG2), " genes x ", ncol(LOG2), " primary tumours")

SIGWU <- as.matrix(read.delim(file.path(BASE, "cache/v3/deconv/wu2021_signature_matrix.txt"),
                              row.names = 1, check.names = FALSE))
LM22  <- as.matrix(read.delim(file.path(BASE, "data/deconv/LM22.txt"),
                              row.names = 1, check.names = FALSE))
REFWU <- as.matrix(read.delim(file.path(BASE, "cache/v3/deconv/wu2021_meancpm_by_celltype.tsv"),
                              row.names = 1, check.names = FALSE))
logf("Wu signature matrix ", nrow(SIGWU), " x ", ncol(SIGWU), "; Wu full reference ",
     nrow(REFWU), " x ", ncol(REFWU), "; LM22 ", nrow(LM22), " x ", ncol(LM22))

status <- list()
save_res <- function(name, M, note) {
  M <- as.matrix(M)
  if (!all(rownames(M) %in% SS)) M <- t(M)
  M <- M[SS, , drop = FALSE]
  saveRDS(M, file.path(OUT, paste0(name, ".rds")))
  logf("SAVED ", name, ": ", ncol(M), " features [", paste(colnames(M), collapse = ", "), "]")
  status[[name]] <<- data.frame(method = name, ok = TRUE, n_features = ncol(M),
                                n_samples = nrow(M), note = note, stringsAsFactors = FALSE)
}
attempt <- function(name, f, note = "") {
  if (nzchar(ONLY) && !grepl(ONLY, name)) return(invisible())
  logf("\n---- ", name, " ----"); t0 <- Sys.time()
  r <- tryCatch(f(), error = function(e) { logf("FAILED: ", conditionMessage(e)); NULL })
  if (is.null(r)) { status[[name]] <<- data.frame(method = name, ok = FALSE, n_features = NA,
      n_samples = NA, note = paste("failed:", note), stringsAsFactors = FALSE); return(invisible()) }
  save_res(name, r, note)
  logf("  elapsed ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
}

## ============================================================ constrained LS family ====
## shared helper: align a signature matrix to the bulk
align <- function(SIG, BULK) {
  gi <- intersect(rownames(SIG), rownames(BULK))
  gi <- gi[apply(SIG[gi, , drop = FALSE], 1, function(v) all(is.finite(v))) &
           rowSums(BULK[gi, , drop = FALSE]) > 0]
  list(S = SIG[gi, , drop = FALSE], Y = BULK[gi, , drop = FALSE], n = length(gi))
}

## --- qprog: min ||Sw - y||^2 s.t. w >= 0, sum(w) = 1  (limSolve::lsei) ---------------
attempt("QPROG_Wu", function() {
  library(limSolve)
  a <- align(SIGWU, TPM); logf("  genes used = ", a$n)
  k <- ncol(a$S)
  ## NOTE: run serially and on a rescaled signature. Unscaled (CPM-magnitude) input
  ## makes t(S)S badly conditioned and solve.QP reports "constraints are inconsistent"
  ## for 428/1097 samples; dividing S and Y by max(S) fixes it (0 failures). Forked
  ## mclapply additionally corrupted solve.QP through the threaded BLAS, so this one
  ## method is single-threaded on purpose.
  sc <- max(a$S); S2 <- a$S / sc; Y2 <- a$Y / sc
  M <- t(simplify2array(lapply(seq_len(ncol(Y2)), function(j) {
    s <- tryCatch(limSolve::lsei(A = S2, B = Y2[, j], E = matrix(1, 1, k), F = 1,
                                 G = diag(k), H = rep(0, k), type = 2)$X,
                  error = function(e) rep(NA_real_, k))
    s / sum(s)
  })))
  logf("  samples with no QP solution: ", sum(!complete.cases(M)))
  dimnames(M) <- list(colnames(a$Y), colnames(a$S)); M
}, "limSolve::lsei, non-negative + sum-to-one, Wu breast signature, TPM/max(S) rescaled")

## --- RLS: MASS::rlm (Huber) with negatives clipped, renormalised ----------------------
attempt("RLS_Wu", function() {
  a <- align(SIGWU, TPM); logf("  genes used = ", a$n)
  M <- t(simplify2array(mclapply(seq_len(ncol(a$Y)), function(j) {
    co <- tryCatch(coef(MASS::rlm(a$S, a$Y[, j], maxit = 100)),
                   error = function(e) rep(NA_real_, ncol(a$S)))
    co[!is.finite(co) | co < 0] <- 0
    if (sum(co) <= 0) rep(NA_real_, ncol(a$S)) else co / sum(co)
  }, mc.cores = CORES)))
  dimnames(M) <- list(colnames(a$Y), colnames(a$S)); M
}, "MASS::rlm Huber robust regression, negatives clipped, renormalised, Wu signature")

## --- OLS: plain least squares, clipped, renormalised (the crudest baseline) ----------
attempt("OLS_Wu", function() {
  a <- align(SIGWU, TPM); logf("  genes used = ", a$n)
  B <- qr.solve(a$S, a$Y)
  B[B < 0] <- 0
  M <- t(sweep(B, 2, colSums(B), "/"))
  dimnames(M) <- list(colnames(a$Y), colnames(a$S)); M
}, "unconstrained OLS, negatives clipped, renormalised, Wu signature")

## --- DWLS: dampened weighted least squares (Tsoucas et al. 2019) --------------------
## Own implementation of the published algorithm. Deviation recorded: the dampening
## constant search uses 25 random half-gene subsets per candidate j (paper uses 100)
## to make 1,097 samples tractable; the constant is searched per sample as published.
attempt("DWLS_Wu", function() {
  library(quadprog)
  a <- align(SIGWU, TPM); logf("  genes used = ", a$n)
  S <- a$S; K <- ncol(S)
  solveOLSint <- function(S, b) {
    D <- t(S) %*% S; d <- t(S) %*% b
    sc <- norm(D, "2")
    quadprog::solve.QP(D / sc, d / sc, cbind(diag(K)), rep(0, K))$solution
  }
  solveDWLSj <- function(S, b, sol, j) {
    mult <- 2^(j - 1)
    ws <- as.vector((1 / (S %*% sol))^2); ws[!is.finite(ws)] <- max(ws[is.finite(ws)])
    wsS <- ws / min(ws); wsS[wsS > mult] <- mult
    W <- wsS
    SW <- S * sqrt(W); bW <- b * sqrt(W)
    D <- t(SW) %*% SW; d <- t(SW) %*% bW
    sc <- norm(D, "2")
    quadprog::solve.QP(D / sc, d / sc, cbind(diag(K)), rep(0, K))$solution
  }
  findJ <- function(S, b, sol, nrep = 25L) {
    ws <- as.vector((1 / (S %*% sol))^2); ws[!is.finite(ws)] <- NA
    wsS <- ws / min(ws, na.rm = TRUE)
    jmax <- max(2L, ceiling(log2(max(wsS, na.rm = TRUE))))
    jmax <- min(jmax, 20L)
    sds <- numeric(jmax)
    ng <- nrow(S)
    for (j in seq_len(jmax)) {
      mult <- 2^(j - 1); wd <- wsS; wd[is.na(wd) | wd > mult] <- mult
      sols <- matrix(NA_real_, K, nrep)
      for (i in seq_len(nrep)) {
        set.seed(i)
        sub <- sample.int(ng, floor(ng * 0.5))
        fit <- tryCatch(lm.wfit(S[sub, , drop = FALSE], b[sub], w = wd[sub]), error = function(e) NULL)
        if (!is.null(fit)) { cf <- fit$coefficients
          if (sum(cf) != 0) sols[, i] <- cf * sum(sol) / sum(cf) }
      }
      sds[j] <- mean(apply(sols, 1, sd, na.rm = TRUE)^2, na.rm = TRUE)
    }
    which.min(sds)
  }
  M <- t(simplify2array(mclapply(seq_len(ncol(a$Y)), function(jx) {
    b <- a$Y[, jx]
    sol <- tryCatch(solveOLSint(S, b), error = function(e) rep(1 / K, K))
    if (any(!is.finite(sol)) || sum(sol) <= 0) sol <- rep(1 / K, K)
    sol[sol <= 0] <- 1e-8
    jj <- tryCatch(findJ(S, b, sol), error = function(e) 1L)
    it <- 0L; ch <- 1
    while (ch > 0.01 && it < 100L) {
      new <- tryCatch(solveDWLSj(S, b, sol, jj), error = function(e) sol)
      new[new <= 0] <- 1e-10
      avg <- rowMeans(cbind(new, sol, sol, sol))
      ch <- sqrt(sum((avg - sol)^2)); sol <- avg; it <- it + 1L
    }
    sol / sum(sol)
  }, mc.cores = CORES)))
  dimnames(M) <- list(colnames(a$Y), colnames(SIGWU)); M
}, "own implementation of Tsoucas 2019 DWLS (25-rep dampening search), Wu signature")

## ============================================================ dtangle family ==========
attempt("DTANGLE_Wu", function() {
  library(dtangle)
  ref <- log2(REFWU + 1)
  gi <- intersect(rownames(ref), rownames(LOG2))
  Y <- t(LOG2[gi, , drop = FALSE]); R <- t(ref[gi, , drop = FALSE])
  logf("  genes shared with Wu reference = ", length(gi))
  mk <- dtangle::find_markers(Y = rbind(Y, R), pure_samples = lapply(seq_len(nrow(R)),
        function(i) nrow(Y) + i), data_type = "rna-seq", marker_method = "ratio")
  nmk <- lapply(mk$L, function(v) max(5L, floor(length(v) * 0.05)))
  d <- dtangle::dtangle(Y = rbind(Y, R), pure_samples = lapply(seq_len(nrow(R)),
        function(i) nrow(Y) + i), n_markers = unlist(nmk), markers = mk$L,
        data_type = "rna-seq")
  P <- d$estimates[seq_len(nrow(Y)), , drop = FALSE]
  colnames(P) <- rownames(ref) [seq_len(ncol(P))]
  colnames(P) <- colnames(REFWU)
  rownames(P) <- rownames(Y); P
}, "dtangle (Hunt 2019), log2 scale, Wu breast atlas reference, ratio markers top 5%")

attempt("DTANGLE_LM22", function() {
  library(dtangle)
  ref <- log2(LM22 + 1)
  gi <- intersect(rownames(ref), rownames(LOG2))
  Y <- t(LOG2[gi, , drop = FALSE]); R <- t(ref[gi, , drop = FALSE])
  logf("  genes shared with LM22 = ", length(gi))
  ps <- lapply(seq_len(nrow(R)), function(i) nrow(Y) + i)
  mk <- dtangle::find_markers(Y = rbind(Y, R), pure_samples = ps,
        data_type = "rna-seq", marker_method = "ratio")
  nmk <- unlist(lapply(mk$L, function(v) max(5L, floor(length(v) * 0.05))))
  d <- dtangle::dtangle(Y = rbind(Y, R), pure_samples = ps, n_markers = nmk,
        markers = mk$L, data_type = "rna-seq")
  P <- d$estimates[seq_len(nrow(Y)), , drop = FALSE]
  colnames(P) <- colnames(LM22); rownames(P) <- rownames(Y); P
}, "dtangle, log2 scale, LM22 immune reference")

## ============================================================ immunedeconv wrapper ====
attempt("IDECONV_TIMER", function() {
  library(immunedeconv)
  r <- immunedeconv::deconvolute(TPM, "timer", indications = rep("brca", ncol(TPM)))
  M <- as.matrix(r[, -1]); rownames(M) <- r$cell_type; t(M)
}, "immunedeconv wrapper, TIMER, indication brca")

attempt("IDECONV_ABIS", function() {
  library(immunedeconv)
  r <- immunedeconv::deconvolute(TPM, "abis")
  M <- as.matrix(r[, -1]); rownames(M) <- r$cell_type; t(M)
}, "immunedeconv wrapper, ABIS")

attempt("IDECONV_EPIC", function() {
  library(immunedeconv)
  r <- immunedeconv::deconvolute(TPM, "epic", tumor = TRUE)
  M <- as.matrix(r[, -1]); rownames(M) <- r$cell_type; t(M)
}, "immunedeconv wrapper, EPIC (cross-implementation control on the direct EPIC call)")

attempt("IDECONV_MCP", function() {
  library(immunedeconv)
  r <- immunedeconv::deconvolute(TPM, "mcp_counter")
  M <- as.matrix(r[, -1]); rownames(M) <- r$cell_type; t(M)
}, "immunedeconv wrapper, MCP-counter (cross-implementation control)")

attempt("IDECONV_QUANTISEQ", function() {
  library(immunedeconv)
  r <- immunedeconv::deconvolute(TPM, "quantiseq", tumor = TRUE)
  M <- as.matrix(r[, -1]); rownames(M) <- r$cell_type; t(M)
}, "immunedeconv wrapper, quanTIseq (cross-implementation control)")

st <- rbindlist(status)
fwrite(st, file.path(BASE, paste0("results/v3/deconv_newmethod_status", if (nzchar(ONLY)) paste0("_", ONLY) else "", ".csv")))
logf("\n=== NEW METHOD STATUS ===")
for (i in seq_len(nrow(st))) logf(sprintf("%-18s ok=%-5s features=%s", st$method[i], st$ok[i], st$n_features[i]))
logf("DONE 43")
