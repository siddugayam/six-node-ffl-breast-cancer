#!/usr/bin/env Rscript
## 50_deconv2_matched_collagen_control.R
## The initial collagen-free Wu signature was built by re-running the whole
## marker-selection sweep with COL* excluded; that matrix has 841 genes and was chosen
## at G = 100, whereas the collagen-containing matrix has 641 genes and was chosen at
## G = 75 (results/v3/deconv_wu_signature_condition_numbers.csv). The two therefore
## differ in size as well as in collagen content.
## This script builds strictly matched controls instead:
## the identical 641-gene signature with
##   (a) its 8 COL* genes deleted        -> 633 genes
##   (b) only COL1A1 deleted             -> 640 genes  (COL3A1 is not in this matrix)
## and re-run the same estimators. Any difference can then only come from the collagens.
suppressPackageStartupMessages({library(data.table); library(nnls); library(limSolve)})
.libPaths(c(path.expand("~/Rlib_deconv2"), .libPaths()))
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/50_deconv2_matched_control.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
OUT <- file.path(BASE, "cache/v3/deconv2/out")

SIG <- as.matrix(read.delim(file.path(BASE, "cache/v3/deconv/wu2021_signature_matrix.txt"),
                            row.names = 1, check.names = FALSE))
TPM <- readRDS(file.path(BASE, "cache/v3/deconv/bulk_tpm.rds"))
SS  <- colnames(TPM)
colg <- grep("^COL[0-9]", rownames(SIG), value = TRUE)
logf("base signature: ", nrow(SIG), " genes; COL* present: ", paste(colg, collapse = ", "))
logf("COL1A1 in signature: ", "COL1A1" %in% rownames(SIG),
     " ; COL3A1 in signature: ", "COL3A1" %in% rownames(SIG))

variants <- list(
  full        = SIG,
  minusAllCOL = SIG[!rownames(SIG) %in% colg, , drop = FALSE],
  minusCOL1A1 = SIG[rownames(SIG) != "COL1A1", , drop = FALSE])
for (nm in names(variants)) logf("variant ", nm, ": ", nrow(variants[[nm]]), " genes")

run_nnls <- function(S, Y) {
  gi <- intersect(rownames(S), rownames(Y)); S <- S[gi, ]; Y <- Y[gi, ]
  t(simplify2array(lapply(seq_len(ncol(Y)), function(j) { w <- nnls::nnls(S, Y[, j])$x; w / sum(w) })))
}
run_qp <- function(S, Y) {
  gi <- intersect(rownames(S), rownames(Y)); S <- S[gi, ]; Y <- Y[gi, ]; K <- ncol(S)
  sc <- max(S)
  t(simplify2array(lapply(seq_len(ncol(Y)), function(j)
    tryCatch({ s <- limSolve::lsei(A = S / sc, B = Y[, j] / sc, E = matrix(1, 1, K), F = 1,
                                   G = diag(K), H = rep(0, K), type = 2)$X; s / sum(s) },
             error = function(e) rep(NA_real_, K)))))
}
res <- list()
for (nm in names(variants)) {
  S <- variants[[nm]]
  for (alg in c("NNLS", "QPROG")) {
    M <- if (alg == "NNLS") run_nnls(S, TPM) else run_qp(S, TPM)
    dimnames(M) <- list(SS, colnames(S))
    tag <- paste0(alg, "_Wu641_", nm)
    saveRDS(M, file.path(OUT, paste0(tag, ".rds")))
    res[[tag]] <- M[, "CAFs"]
    logf(sprintf("%-28s CAF mean %.4f sd %.4f  failures %d", tag, mean(M[, "CAFs"], na.rm = TRUE),
                 sd(M[, "CAFs"], na.rm = TRUE), sum(!complete.cases(M))))
  }
}
CM <- do.call(cbind, res)
C <- cor(CM, method = "spearman", use = "pairwise.complete.obs")
logf("\n=== Spearman between the matched signature variants (CAF fraction) ===")
print(round(C, 4)); cat(capture.output(print(round(C, 4))), sep = "\n", file = LOG, append = TRUE)
fwrite(data.table(estimate = rownames(C), as.data.frame(C)),
       file.path(BASE, "results/v3/deconv_matched_collagen_control_corr.csv"))

## how much does the collagen content of the signature actually move the CAF estimate?
g <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
col1 <- g["COL1A1", SS]; col3 <- g["COL3A1", SS]
tab <- rbindlist(lapply(names(res), function(nm) data.table(estimate = nm,
  rho_COL1A1 = cor(res[[nm]], col1, method = "spearman", use = "complete.obs"),
  rho_COL3A1 = cor(res[[nm]], col3, method = "spearman", use = "complete.obs"))))
fwrite(tab, file.path(BASE, "results/v3/deconv_matched_collagen_control_vs_collagen.csv"))
logf("\n=== correlation of each variant's CAF estimate with the collagen outcomes ===")
for (i in seq_len(nrow(tab))) logf(sprintf("%-28s rho(COL1A1) %+0.3f   rho(COL3A1) %+0.3f",
  tab$estimate[i], tab$rho_COL1A1[i], tab$rho_COL3A1[i]))
logf("DONE 50")
