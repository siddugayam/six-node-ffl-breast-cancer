#!/usr/bin/env Rscript
## 33_cibersort_lm22_variants.R -- CIBERSORT/LM22 goodness-of-fit under both
## quantile-normalisation settings, with the published doPerm null.
suppressPackageStartupMessages({library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/33_cibersort_variants.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
source(file.path(BASE, "scripts/v3/23_cibersort_impl.R"))
D <- file.path(BASE, "cache/v3/deconv")
TPM  <- readRDS(file.path(D, "bulk_tpm.rds"))
LM22 <- as.matrix(read.delim(file.path(BASE, "data/deconv/LM22.txt"), row.names = 1, check.names = FALSE))
out <- list()
for (qn in c(FALSE, TRUE)) {
  r <- run_cibersort(LM22, TPM, QN = qn, perm = 500, cores = 4)
  tag <- paste0("LM22_QN", qn)
  logf(tag, ": genes used ", r$n_genes_used, " ; median fit r = ", round(median(r$stats$correlation), 3),
       " ; IQR ", round(quantile(r$stats$correlation, .25), 3), "-", round(quantile(r$stats$correlation, .75), 3),
       " ; samples p<0.05 = ", sum(r$stats$p_value < 0.05), "/", nrow(r$stats))
  saveRDS(r$fractions, file.path(D, "out", paste0("CIBERSORT_", tag, ".rds")))
  out[[tag]] <- data.table(variant = tag, n_genes = r$n_genes_used,
                           median_r = median(r$stats$correlation),
                           q25_r = quantile(r$stats$correlation, .25), q75_r = quantile(r$stats$correlation, .75),
                           n_p_lt_05 = sum(r$stats$p_value < 0.05), n = nrow(r$stats))
}
A <- readRDS(file.path(D, "out", "CIBERSORT_LM22_QNFALSE.rds"))
B <- readRDS(file.path(D, "out", "CIBERSORT_LM22_QNTRUE.rds"))
cc <- sapply(colnames(A), function(j) cor(A[, j], B[, j], method = "spearman"))
logf("\nSpearman QN=FALSE vs QN=TRUE per cell type:")
for (j in names(sort(cc))) logf(sprintf("  %-36s %+.3f", j, cc[j]))
fwrite(rbindlist(out), file.path(BASE, "results/v3/deconv_cibersort_lm22_variants.csv"))
fwrite(data.table(cell_type = names(cc), spearman_QNfalse_vs_QNtrue = as.numeric(cc)),
       file.path(BASE, "results/v3/deconv_cibersort_lm22_QN_agreement.csv"))
logf("DONE 33")
