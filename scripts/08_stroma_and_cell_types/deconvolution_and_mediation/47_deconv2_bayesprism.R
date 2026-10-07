#!/usr/bin/env Rscript
## 47_deconv2_bayesprism.R
## Run BayesPrism (Danko-Lab, v2.2.3) against the Wu et al. 2021 breast
## single-cell atlas, and compare it with the InstaPrism estimates.
## BayesPrism's Gibbs sampler on 1,097 bulk samples x ~16k genes is expensive, so it is
## run on a random subset of tumours (default 200, seed fixed) and used as a validation
## of InstaPrism, which is an analytic reimplementation of the same model and was run on
## all 1,097. Whatever the subset shows about InstaPrism is reported as found.
## Reference construction: one "cell state" per (patient x cell type) group from the Wu
## atlas -- the design BayesPrism recommends -- summed UMI counts, cell types as in the
## atlas. Output: cache/v3/deconv2/out/BAYESPRISM_Wu_subset.rds
##               results/v3/deconv_bayesprism_vs_instaprism.csv
suppressPackageStartupMessages({library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/47_deconv2_bayesprism.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
NSUB <- as.integer(Sys.getenv("NSUB", "200"))
CORES <- 8L

library(BayesPrism)
logf("BayesPrism version ", as.character(packageVersion("BayesPrism")))

REF <- as.matrix(read.delim(file.path(BASE, "cache/v3/deconv/wu2021_sum_umi_by_patient_celltype.tsv"),
                            row.names = 1, check.names = FALSE))
logf("reference (genes x patient|celltype): ", nrow(REF), " x ", ncol(REF))
state <- colnames(REF)
ctype <- sub("^.*\\|", "", state)
logf("cell types: ", paste(sort(unique(ctype)), collapse = ", "))
logf("cell states per type: ", paste(sprintf("%s=%d", names(table(ctype)), as.integer(table(ctype))), collapse = ", "))

TPM  <- readRDS(file.path(BASE, "cache/v3/deconv/bulk_tpm.rds"))
set.seed(20260910)
sub <- sort(sample(colnames(TPM), min(NSUB, ncol(TPM))))
logf("bulk subset: ", length(sub), " tumours")

gi <- intersect(rownames(REF), rownames(TPM))
logf("genes shared reference/bulk: ", length(gi))
bk <- t(TPM[gi, sub, drop = FALSE])           # samples x genes
sc <- t(REF[gi, , drop = FALSE])              # cell states x genes
mode(sc) <- "numeric"

## BayesPrism's own preprocessing: drop outlier/uninformative genes
sc.filt <- BayesPrism::cleanup.genes(input = sc, input.type = "count.matrix",
             species = "hs", gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
             exp.cells = 1)
logf("reference after cleanup.genes: ", ncol(sc.filt), " genes")
gi2 <- intersect(colnames(sc.filt), colnames(bk))
bk2 <- bk[, gi2, drop = FALSE]; sc2 <- sc.filt[, gi2, drop = FALSE]
logf("genes entering the prism: ", length(gi2))

t0 <- Sys.time()
pr <- BayesPrism::new.prism(reference = sc2, mixture = bk2, input.type = "count.matrix",
        cell.type.labels = ctype, cell.state.labels = state, key = "Cancer Epithelial",
        outlier.cut = 0.01, outlier.fraction = 0.1)
res <- BayesPrism::run.prism(prism = pr, n.cores = CORES)
theta <- BayesPrism::get.fraction(bp = res, which.theta = "final", state.or.type = "type")
logf("run.prism elapsed ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
logf("theta: ", nrow(theta), " x ", ncol(theta), " -> ", paste(colnames(theta), collapse = ", "))
saveRDS(theta, file.path(BASE, "cache/v3/deconv2/out/BAYESPRISM_Wu_subset.rds"))

IP <- readRDS(file.path(BASE, "cache/v3/deconv/out/InstaPrism_Wu.rds"))
logf("InstaPrism matrix ", nrow(IP), " x ", ncol(IP), " -> ", paste(colnames(IP), collapse = ", "))
common <- intersect(colnames(theta), colnames(IP))
ss <- intersect(rownames(theta), rownames(IP))
logf("comparing on ", length(ss), " tumours and ", length(common), " cell types")
agr <- rbindlist(lapply(common, function(k) data.table(cell_type = k,
        spearman = cor(theta[ss, k], IP[ss, k], method = "spearman"),
        pearson  = cor(theta[ss, k], IP[ss, k]),
        mean_bayesprism = mean(theta[ss, k]), mean_instaprism = mean(IP[ss, k]))))
fwrite(agr, file.path(BASE, "results/v3/deconv_bayesprism_vs_instaprism.csv"))
logf("\n=== BayesPrism vs InstaPrism (same reference, same tumours) ===")
for (i in seq_len(nrow(agr))) logf(sprintf("%-20s rho = %+0.3f  r = %+0.3f   mean %.4f vs %.4f",
  agr$cell_type[i], agr$spearman[i], agr$pearson[i], agr$mean_bayesprism[i], agr$mean_instaprism[i]))
logf("median Spearman across cell types = ", round(median(agr$spearman), 3))
logf("DONE 47")
