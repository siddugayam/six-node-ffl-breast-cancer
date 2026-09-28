#!/usr/bin/env Rscript
## 25_run_deconvolution.R -- run every deconvolution method that installed, on the
## 1,097 TCGA-BRCA primary tumours. Each method is wrapped in tryCatch so a failure
## is recorded rather than silently dropped.
suppressPackageStartupMessages({library(matrixStats); library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/25_run_deconvolution.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
D <- file.path(BASE, "cache/v3/deconv"); dir.create(file.path(D, "out"), showWarnings = FALSE)

LOG2 <- readRDS(file.path(D, "bulk_log2.rds"))
TPM  <- readRDS(file.path(D, "bulk_tpm.rds"))
SS   <- colnames(LOG2)
logf("bulk: ", nrow(LOG2), " genes x ", ncol(LOG2), " primary tumours")
status <- list()
save_res <- function(name, M, note = "") {
  M <- as.matrix(M)
  stopifnot(all(rownames(M) == SS) || all(colnames(M) == SS))
  if (!all(rownames(M) == SS)) M <- t(M)
  M <- M[SS, , drop = FALSE]
  saveRDS(M, file.path(D, "out", paste0(name, ".rds")))
  logf("SAVED ", name, ": ", ncol(M), " features  [", paste(head(colnames(M), 30), collapse = ", "), "]")
  status[[name]] <<- data.frame(method = name, ok = TRUE, n_features = ncol(M),
                                n_samples = nrow(M), note = note, stringsAsFactors = FALSE)
}
attempt <- function(name, expr, note = "") {
  logf("\n---- ", name, " ----")
  t0 <- Sys.time()
  r <- tryCatch(eval(expr), error = function(e) { logf("FAILED: ", conditionMessage(e)); NULL })
  if (is.null(r)) { status[[name]] <<- data.frame(method = name, ok = FALSE, n_features = NA,
                       n_samples = NA, note = "failed", stringsAsFactors = FALSE); return(invisible()) }
  save_res(name, r, note)
  logf("  elapsed ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
}

## ---------------------------------------------------------------- 1. ESTIMATE ----
attempt("ESTIMATE", quote({
  library(estimate)
  fin <- file.path(D, "bulk_log2.txt")
  f1 <- file.path(D, "estimate_filtered.gct"); f2 <- file.path(D, "estimate_score.gct")
  estimate::filterCommonGenes(input.f = fin, output.f = f1, id = "GeneSymbol")
  estimate::estimateScore(input.ds = f1, output.ds = f2, platform = "illumina")
  e <- read.delim(f2, skip = 2, row.names = 1, check.names = FALSE)[, -1, drop = FALSE]
  e <- t(as.matrix(e))
  rownames(e) <- gsub("\\.", "-", rownames(e))
  e
}), note = "estimate pkg, illumina platform, log2 input")

## ---------------------------------------------------------------- 2. xCell -------
attempt("xCell", quote({
  library(xCell)
  x <- xCell::xCellAnalysis(LOG2, rnaseq = TRUE, parallel.sz = 8)
  t(x)
}), note = "xCell 64 types + 3 aggregate scores, rnaseq=TRUE")

## ---------------------------------------------------------------- 3. MCP-counter -
attempt("MCPcounter", quote({
  library(MCPcounter)
  m <- MCPcounter::MCPcounter.estimate(LOG2, featuresType = "HUGO_symbols")
  t(m)
}), note = "log2 input, HUGO symbols")

## ---------------------------------------------------------------- 4. EPIC --------
attempt("EPIC", quote({
  library(EPIC)
  e <- EPIC::EPIC(bulk = TPM, reference = "TRef", withOtherCells = TRUE)
  e$cellFractions
}), note = "TRef tumour reference (CAFs + endothelial + immune), TPM-scaled input")

attempt("EPIC_BRef", quote({
  library(EPIC)
  e <- EPIC::EPIC(bulk = TPM, reference = "BRef", withOtherCells = TRUE)
  e$cellFractions
}), note = "BRef blood reference, immune only")

## ---------------------------------------------------------------- 5. quanTIseq ---
attempt("quanTIseq", quote({
  library(quantiseqr)
  q <- quantiseqr::run_quantiseq(expression_data = TPM, signature_matrix = "TIL10",
                                 is_arraydata = FALSE, is_tumordata = TRUE, scale_mRNA = TRUE)
  rn <- q$Sample; q$Sample <- NULL
  M <- as.matrix(q); rownames(M) <- rn; M
}), note = "TIL10, is_tumordata=TRUE, mRNA scaling on")

## ---------------------------------------------------------------- 6. ConsensusTME -
attempt("ConsensusTME", quote({
  library(ConsensusTME)
  c1 <- ConsensusTME::consensusTMEAnalysis(as.matrix(LOG2), cancer = "BRCA", statMethod = "ssgsea")
  t(c1)
}), note = "BRCA gene sets, ssGSEA")

## ---------------------------------------------------------------- 7. CIBERSORT ---
attempt("CIBERSORT_LM22", quote({
  source(file.path(BASE, "scripts/v3/23_cibersort_impl.R"))
  LM22 <- as.matrix(read.delim(file.path(BASE, "data/deconv/LM22.txt"), row.names = 1, check.names = FALSE))
  r <- run_cibersort(LM22, TPM, QN = FALSE, perm = 200, cores = 8)
  saveRDS(r$stats, file.path(D, "out", "CIBERSORT_LM22_stats.rds"))
  logf("  CIBERSORT genes used = ", r$n_genes_used, "; median mixture r = ",
       round(median(r$stats$correlation, na.rm = TRUE), 3),
       "; samples with p<0.05 = ", sum(r$stats$p_value < 0.05, na.rm = TRUE), "/", nrow(r$stats))
  r$fractions
}), note = "own validated nu-SVR re-implementation, LM22, QN=FALSE (RNA-seq)")

## ---------------------------------------------------- 8. Wu-atlas custom signature -
for (tag in c("", "_nocollagen")) {
  attempt(paste0("CBSX_Wu", tag), quote({
    source(file.path(BASE, "scripts/v3/23_cibersort_impl.R"))
    SIG <- as.matrix(read.delim(file.path(D, paste0("wu2021_signature_matrix", tag, ".txt")),
                                row.names = 1, check.names = FALSE))
    r <- run_cibersort(SIG, TPM, QN = FALSE, perm = 200, cores = 8)
    saveRDS(r$stats, file.path(D, "out", paste0("CBSX_Wu", tag, "_stats.rds")))
    logf("  genes used = ", r$n_genes_used, "; median mixture r = ",
         round(median(r$stats$correlation, na.rm = TRUE), 3))
    r$fractions
  }), note = paste0("nu-SVR against Wu et al. 2021 breast atlas signature", tag))
}

attempt("NNLS_Wu", quote({
  library(nnls)
  SIG <- as.matrix(read.delim(file.path(D, "wu2021_signature_matrix.txt"), row.names = 1, check.names = FALSE))
  gi <- intersect(rownames(SIG), rownames(TPM))
  X <- SIG[gi, ]; Y <- TPM[gi, ]
  M <- t(sapply(seq_len(ncol(Y)), function(j) { w <- nnls::nnls(X, Y[, j])$x; w / sum(w) }))
  dimnames(M) <- list(colnames(Y), colnames(X))
  M
}), note = "non-negative least squares against the same Wu signature matrix")

st <- rbindlist(status)
fwrite(st, file.path(BASE, "results/v3/deconv_method_status.csv"))
logf("\n=== METHOD STATUS ===")
for (i in seq_len(nrow(st))) logf(sprintf("%-22s ok=%-5s features=%s", st$method[i], st$ok[i], st$n_features[i]))
logf("DONE 25")
