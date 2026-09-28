#!/usr/bin/env Rscript
## 26_instaprism_bayesprism.R
## BayesPrism-family deconvolution (InstaPrism: the exact BayesPrism model solved
## with a deterministic fixed-point iteration instead of Gibbs sampling) using the
## Wu et al. 2021 breast-tumour atlas as reference.
## Cell states = 187 patient x celltype pseudobulks with >=20 cells; cell types =
## the 9 major types. A collagen-free variant is also produced.
suppressPackageStartupMessages({library(InstaPrism); library(Matrix)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/26_instaprism.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
D <- file.path(BASE, "cache/v3/deconv")

pc <- as.matrix(read.delim(file.path(D, "wu2021_sum_umi_by_patient_celltype.tsv"), row.names = 1, check.names = FALSE))
TPM <- readRDS(file.path(D, "bulk_tpm.rds"))
logf("reference: ", nrow(pc), " genes x ", ncol(pc), " patient-x-celltype pseudobulks")
cs <- colnames(pc); ct <- sub("^.*\\|", "", cs)
logf("cell types: ", paste(sort(unique(ct)), collapse = ", "))
logf("cell states per type: ", paste(sprintf("%s=%d", names(table(ct)), table(ct)), collapse = ", "))

run_ip <- function(ref, bulk, tag) {
  rp <- InstaPrism::refPrepare(sc_Expr = ref, cell.type.labels = ct, cell.state.labels = cs)
  t0 <- Sys.time()
  ip <- InstaPrism::InstaPrism(bulk_Expr = bulk, refPhi_cs = rp, n.iter = 100, n.core = 8, verbose = FALSE)
  th <- t(ip@Post.ini.ct@theta)                         # sample x celltype
  logf(tag, ": ", nrow(th), " samples x ", ncol(th), " cell types; elapsed ",
       round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s")
  logf("  median CAF fraction = ", round(median(th[, "CAFs"]), 4),
       "  median Cancer Epithelial = ", round(median(th[, "Cancer Epithelial"]), 4))
  ## also keep the updated (Gibbs-refined) theta if InstaPrism produced one
  saveRDS(th, file.path(D, "out", paste0(tag, ".rds")))
  th
}

gi <- intersect(rownames(pc), rownames(TPM))
logf("genes shared reference/bulk: ", length(gi))
A <- run_ip(pc[gi, ], TPM[gi, ], "InstaPrism_Wu")

colg <- grep("^COL[0-9]", gi, value = TRUE)
gi2 <- setdiff(gi, colg)
logf("collagen genes removed for the clean run: ", length(colg))
B <- run_ip(pc[gi2, ], TPM[gi2, ], "InstaPrism_Wu_nocollagen")

logf("Spearman CAF(full) vs CAF(collagen-free) = ",
     round(cor(A[, "CAFs"], B[, "CAFs"], method = "spearman"), 4))
logf("DONE 26")
