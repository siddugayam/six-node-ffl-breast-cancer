#!/usr/bin/env Rscript
## 23_cibersort_impl.R -- faithful open re-implementation of the CIBERSORT
## nu-support-vector-regression deconvolution of Newman et al. 2015 (Nat Methods
## 12:453-457).  Sourced by later scripts; running it directly executes the
## validation suite and writes results/v3/deconv_cibersort_validation.csv.
##
## Why a re-implementation: the original CIBERSORT.R and CIBERSORTx web tool are
## behind a Stanford registration/licence wall that requires the account holder to
## click a confirmation link. The algorithm itself is fully published; this file
## reproduces it step for step (global standardisation of the signature matrix,
## per-sample standardisation of the mixture, nu in {0.25,0.5,0.75}, linear kernel,
## un-scaled, coefficients clipped at zero and renormalised to sum 1, model chosen
## by minimum RMSE) and is validated below against mixtures with known composition.
suppressPackageStartupMessages({library(e1071); library(parallel); library(preprocessCore)})

cibersort_core <- function(X, y, nus = c(0.25, 0.5, 0.75)) {
  mods <- lapply(nus, function(nu)
    e1071::svm(X, y, type = "nu-regression", kernel = "linear", nu = nu, scale = FALSE))
  rm <- numeric(length(mods)); cr <- numeric(length(mods))
  for (t in seq_along(mods)) {
    w <- t(mods[[t]]$coefs) %*% mods[[t]]$SV
    w[w < 0] <- 0
    if (sum(w) == 0) { rm[t] <- Inf; cr[t] <- NA; next }
    w <- w / sum(w)
    k <- as.vector(X %*% as.vector(w))
    rm[t] <- sqrt(mean((k - y)^2)); cr[t] <- cor(k, y)
  }
  mn <- which.min(rm)
  q <- t(mods[[mn]]$coefs) %*% mods[[mn]]$SV
  q[q < 0] <- 0
  list(w = as.vector(q / sum(q)), rmse = rm[mn], r = cr[mn], nu = nus[mn])
}

## Y : gene x sample mixture matrix, LINEAR space, rownames = symbols
## X : gene x celltype signature matrix, LINEAR space
run_cibersort <- function(X, Y, QN = FALSE, perm = 0, cores = 8, seed = 1) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  X <- X[order(rownames(X)), , drop = FALSE]
  Y <- Y[order(rownames(Y)), , drop = FALSE]
  ## base-R quantile normalisation (mean of the sorted columns, applied by rank).
  ## preprocessCore::normalize.quantiles is avoided because its OpenMP build conflicts
  ## with the threaded BLAS in this R installation (pthread_create returns 22).
  if (QN) { cn <- colnames(Y); rn <- rownames(Y)
            ref <- rowMeans(apply(Y, 2, sort))
            Y <- apply(Y, 2, function(v) {
                   o <- order(v); z <- numeric(length(v)); z[o] <- ref
                   ## ties get the mean of the reference values spanning the tied block,
                   ## which is what preprocessCore does
                   ave(z, v)
                 })
            colnames(Y) <- cn; rownames(Y) <- rn }
  gi <- intersect(rownames(X), rownames(Y))
  Xs <- X[gi, , drop = FALSE]; Ys <- Y[gi, , drop = FALSE]
  Xs <- (Xs - mean(Xs)) / sd(as.vector(Xs))
  fit1 <- function(j) {
    y <- Ys[, j]; sdy <- sd(y)
    if (!is.finite(sdy) || sdy == 0) return(c(rep(NA_real_, ncol(Xs)), NA, NA, NA))
    y <- (y - mean(y)) / sdy
    r <- cibersort_core(Xs, y)
    c(r$w, r$r, r$rmse, r$nu)
  }
  set.seed(seed)
  out <- parallel::mclapply(seq_len(ncol(Ys)), fit1, mc.cores = cores)
  M <- do.call(rbind, out)
  colnames(M) <- c(colnames(Xs), "correlation", "RMSE", "nu_used")
  rownames(M) <- colnames(Ys)
  pv <- rep(NA_real_, ncol(Ys))
  if (perm > 0) {
    ## doPerm exactly as published: draw length(signature genes) values at random from
    ## the pool of all mixture values, standardise, deconvolve, keep the fit correlation
    pool <- as.vector(Ys)
    set.seed(seed + 99)
    nulld <- unlist(parallel::mclapply(seq_len(perm), function(i) {
      yr <- pool[sample(length(pool), nrow(Xs))]
      s <- sd(yr); if (!is.finite(s) || s == 0) return(NA_real_)
      cibersort_core(Xs, (yr - mean(yr)) / s)$r
    }, mc.cores = cores))
    nulld <- sort(nulld[is.finite(nulld)])
    pv <- sapply(M[, "correlation"], function(r) (1 + sum(nulld >= r)) / (1 + length(nulld)))
    attr(pv, "null_median") <- median(nulld); attr(pv, "null_max") <- max(nulld)
  }
  list(fractions = M[, seq_len(ncol(Xs)), drop = FALSE],
       stats = data.frame(sample = rownames(M), correlation = M[, "correlation"],
                          RMSE = M[, "RMSE"], nu_used = M[, "nu_used"], p_value = pv,
                          row.names = NULL, stringsAsFactors = FALSE),
       n_genes_used = length(gi))
}

## ------------------------------------------------------------------- validation ------
if (sys.nframe() == 0) {
  BASE <- "/path/to/revision"
  LOG  <- file.path(BASE, "logs/v3/23_cibersort_validation.log"); cat("", file = LOG)
  logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
  LM22 <- as.matrix(read.delim(file.path(BASE, "data/deconv/LM22.txt"), row.names = 1, check.names = FALSE))
  logf("LM22 loaded: ", nrow(LM22), " genes x ", ncol(LM22), " cell types")

  ## (1) pure profiles: each LM22 column deconvolved against LM22 should return itself
  r1 <- run_cibersort(LM22, LM22, QN = FALSE, cores = 8)
  self <- diag(r1$fractions[colnames(LM22), colnames(LM22)])
  logf("PURE-PROFILE recovery: median self-fraction ", round(median(self), 3),
       "  min ", round(min(self), 3), "  n with self-fraction the largest: ",
       sum(apply(r1$fractions, 1, which.max) == seq_len(ncol(LM22))), "/", ncol(LM22))

  ## (2) synthetic mixtures with known proportions (Dirichlet), noise-free and noisy
  set.seed(42); K <- ncol(LM22); NSIM <- 60
  P <- matrix(rgamma(NSIM * K, 0.4, 1), NSIM, K); P <- P / rowSums(P)
  colnames(P) <- colnames(LM22)
  MIX <- LM22 %*% t(P)
  colnames(MIX) <- paste0("sim", seq_len(NSIM))
  r2 <- run_cibersort(LM22, MIX, QN = FALSE, cores = 8)
  E2 <- r2$fractions[colnames(MIX), colnames(LM22)]
  MIXN <- MIX * 2^matrix(rnorm(length(MIX), 0, 0.30), nrow(MIX))       # ~lognormal noise
  r3 <- run_cibersort(LM22, MIXN, QN = FALSE, cores = 8)
  E3 <- r3$fractions[colnames(MIXN), colnames(LM22)]
  res <- data.frame(
    test = c("pure_profiles_self_fraction_median", "noisefree_pearson_r_all_cells",
             "noisefree_rmse", "noisy30pct_pearson_r_all_cells", "noisy30pct_rmse",
             "noisefree_spearman_percell_median", "noisy30pct_spearman_percell_median"),
    value = c(median(self), cor(as.vector(P), as.vector(E2)),
              sqrt(mean((as.vector(P) - as.vector(E2))^2)),
              cor(as.vector(P), as.vector(E3)),
              sqrt(mean((as.vector(P) - as.vector(E3))^2)),
              median(sapply(seq_len(K), function(k) cor(P[, k], E2[, k], method = "spearman"))),
              median(sapply(seq_len(K), function(k) cor(P[, k], E3[, k], method = "spearman")))))
  for (i in seq_len(nrow(res))) logf(sprintf("%-38s %.4f", res$test[i], res$value[i]))
  write.csv(res, file.path(BASE, "results/v3/deconv_cibersort_validation.csv"), row.names = FALSE)
  logf("WROTE results/v3/deconv_cibersort_validation.csv")
}
