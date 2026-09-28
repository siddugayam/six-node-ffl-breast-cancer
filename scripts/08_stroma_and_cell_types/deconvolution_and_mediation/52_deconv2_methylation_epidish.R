#!/usr/bin/env Rscript
## 52_deconv2_methylation_epidish.R
## A FIBROBLAST ESTIMATE FROM A DIFFERENT ASSAY.
## Every RNA-based mediator shares its measurement channel with COL1A1 mRNA, so a
## EpiDISH (Teschendorff et al. 2017) estimates epithelial / fibroblast / immune-cell
## fractions from Illumina 450k DNA METHYLATION using the centEpiFibIC.m reference.
## The TCGA-BRCA 450k matrix already cached for this project (Xena
## cache/multiomics/HumanMethylation450.gz) covers most of the same tumours, so this
## gives a fibroblast fraction measured on DNA, not RNA, in the same samples.
## Output: cache/v3/deconv2/epidish_fractions.rds
##         results/v3/deconv_methylation_epidish.csv
suppressPackageStartupMessages({library(data.table); library(EpiDISH)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/52_deconv2_epidish.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }

data("centEpiFibIC.m", package = "EpiDISH")
ref <- centEpiFibIC.m
logf("EpiDISH reference centEpiFibIC.m: ", nrow(ref), " probes x ", ncol(ref),
     " cell types -> ", paste(colnames(ref), collapse = ", "))

MFILE <- file.path(BASE, "cache/multiomics/HumanMethylation450.gz")
PRB <- file.path(BASE, "cache/v3/deconv2/epidish_probes.txt")
writeLines(rownames(ref), PRB)
logf("extracting ", nrow(ref), " reference probes from ", basename(MFILE))
cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==FNR{p[$1];next} FNR==1 || ($1 in p)' %s -",
               shQuote(MFILE), shQuote(PRB))
t0 <- Sys.time()
M <- fread(cmd = cmd, sep = "\t", header = TRUE)
logf("extraction took ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s; ",
     nrow(M), " probe rows x ", ncol(M) - 1, " samples")
setnames(M, 1, "probe")
bm <- as.matrix(M[, -1]); rownames(bm) <- M$probe
mode(bm) <- "numeric"
logf("beta range ", paste(round(range(bm, na.rm = TRUE), 3), collapse = " .. "),
     "; missing values ", sum(is.na(bm)))
keep <- rowMeans(is.na(bm)) < 0.2
bm <- bm[keep, , drop = FALSE]
logf("probes with <20% missing kept: ", nrow(bm))
for (i in which(rowSums(is.na(bm)) > 0)) bm[i, is.na(bm[i, ])] <- mean(bm[i, ], na.rm = TRUE)
logf("probes shared with the EpiDISH reference: ", length(intersect(rownames(bm), rownames(ref))))

t0 <- Sys.time()
res <- EpiDISH::epidish(beta.m = bm, ref.m = ref, method = "RPC")
logf("epidish RPC elapsed ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
F <- res$estF
logf("fractions: ", nrow(F), " samples x ", ncol(F), " -> ", paste(colnames(F), collapse = ", "))
logf("median fractions: ", paste(sprintf("%s=%.3f", colnames(F), apply(F, 2, median)), collapse = ", "))

SS <- readLines(file.path(BASE, "cache/v3/deconv/samples_primary_tumour.txt"))
rn <- rownames(F)
ov <- intersect(SS, rn)
logf("overlap with our 1,097 primary tumours: ", length(ov))
out <- matrix(NA_real_, length(SS), ncol(F), dimnames = list(SS, colnames(F)))
out[ov, ] <- F[ov, ]
saveRDS(out, file.path(BASE, "cache/v3/deconv2/epidish_fractions.rds"))
fwrite(data.table(sample = SS, as.data.frame(out)),
       file.path(BASE, "results/v3/deconv_methylation_epidish.csv"))

## how does the DNA-based fibroblast fraction compare with the RNA-based ones?
FIB <- readRDS(file.path(BASE, "cache/v3/deconv2/fibroblast_estimates_extended.rds"))
fibcol <- grep("Fib", colnames(F), value = TRUE)[1]
v <- out[, fibcol]
tab <- rbindlist(lapply(colnames(FIB), function(m) data.table(rna_estimate = m,
  rho = cor(v, FIB[, m], method = "spearman", use = "complete.obs"),
  n = sum(!is.na(v) & !is.na(FIB[, m])))))[order(-rho)]
fwrite(tab, file.path(BASE, "results/v3/deconv_methylation_vs_rna_agreement.csv"))
logf("\n=== DNA-methylation fibroblast fraction (", fibcol, ") vs each RNA estimate ===")
for (i in seq_len(nrow(tab))) logf(sprintf("  rho %+0.3f  (n = %4d)  %s", tab$rho[i], tab$n[i], tab$rna_estimate[i]))
logf("median rho across RNA estimates = ", round(median(tab$rho), 3))

g <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
for (k in colnames(F)) {
  a <- out[, k]
  logf(sprintf("EpiDISH %-6s vs COL1A1 rho %+0.3f ; vs COL3A1 rho %+0.3f ; vs ABSOLUTE purity rho %s",
    k, cor(a, g["COL1A1", SS], method = "spearman", use = "complete.obs"),
    cor(a, g["COL3A1", SS], method = "spearman", use = "complete.obs"),
    { p <- readRDS(file.path(BASE, "cache/v3/deconv2/extra_covariates.rds"))$purity
      sprintf("%+0.3f", cor(a, p[SS], method = "spearman", use = "complete.obs")) }))
}
logf("DONE 52")
