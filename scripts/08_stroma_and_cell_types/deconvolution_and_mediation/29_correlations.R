#!/usr/bin/env Rscript
## 29_correlations.R
## (B) pairwise Spearman correlation matrix between every fibroblast/stromal estimate
## (D) correlation of every cell-type fraction from every method with COL1A1, COL3A1
##     and the miR-29 family (and, for context, ETS1/NFKB1)
suppressPackageStartupMessages({library(data.table); library(matrixStats)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/29_correlations.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }

gexp <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
mexp <- readRDS(file.path(BASE, "data/brca_mirna_expr_canonical.rds"))
FIB  <- readRDS(file.path(BASE, "cache/v3/deconv/fibroblast_estimates_all.rds"))
long <- fread(file.path(BASE, "results/v3/deconv_all_celltype_estimates.csv"))
SS   <- rownames(FIB); SB <- sort(intersect(SS, colnames(mexp)))
logf("tumours: gene assay ", length(SS), " ; both assays ", length(SB))

## ------------------------------------------------------------------- (B) ------------
C  <- cor(FIB, method = "spearman")
CP <- cor(FIB, method = "pearson")
fwrite(data.table(estimate = rownames(C), as.data.table(C)),
       file.path(BASE, "results/v3/deconv_fibroblast_correlation_spearman.csv"))
fwrite(data.table(estimate = rownames(CP), as.data.table(CP)),
       file.path(BASE, "results/v3/deconv_fibroblast_correlation_pearson.csv"))
ut <- which(upper.tri(C), arr.ind = TRUE)
PW <- data.table(a = rownames(C)[ut[, 1]], b = colnames(C)[ut[, 2]],
                 spearman = C[ut], pearson = CP[ut])
PW[, `:=`(p_spearman = sapply(seq_len(.N), function(i)
  cor.test(FIB[, a[i]], FIB[, b[i]], method = "spearman", exact = FALSE)$p.value))]
setorder(PW, -spearman)
fwrite(PW, file.path(BASE, "results/v3/deconv_fibroblast_correlation_pairs.csv"))
logf("\n=== (B) PAIRWISE SPEARMAN BETWEEN FIBROBLAST/STROMAL ESTIMATES (", ncol(FIB), " estimates, ",
     nrow(PW), " pairs) ===")
logf("median |rho| = ", round(median(abs(PW$spearman)), 3), "  min = ", round(min(PW$spearman), 3),
     "  max = ", round(max(PW$spearman), 3))
logf("pairs with rho < 0.5: ", sum(PW$spearman < 0.5), " ; rho < 0.3: ", sum(PW$spearman < 0.3),
     " ; rho < 0 : ", sum(PW$spearman < 0))
logf("\nTOP 12 pairs:");    for (i in 1:12) logf(sprintf("  %+.3f  %s  ~  %s", PW$spearman[i], PW$a[i], PW$b[i]))
logf("\nBOTTOM 12 pairs:"); for (i in (nrow(PW)-11):nrow(PW)) logf(sprintf("  %+.3f  %s  ~  %s", PW$spearman[i], PW$a[i], PW$b[i]))
## mean correlation of each estimate with all the others = how typical it is
mm <- sapply(colnames(C), function(j) mean(C[j, setdiff(colnames(C), j)]))
logf("\nmean rho with all other estimates (a low value flags an outlier method):")
for (j in names(sort(mm))) logf(sprintf("  %-32s %+.3f", j, mm[j]))
fwrite(data.table(estimate = names(mm), mean_rho_with_others = as.numeric(mm))[order(mean_rho_with_others)],
       file.path(BASE, "results/v3/deconv_fibroblast_typicality.csv"))
## hierarchical clustering of the estimates
hc <- hclust(as.dist(1 - C), method = "average")
ct <- cutree(hc, k = 3)
logf("\naverage-linkage clustering on 1-rho, k=3:")
for (k in 1:3) logf("  cluster ", k, ": ", paste(names(ct)[ct == k], collapse = ", "))
fwrite(data.table(estimate = names(ct), cluster = as.integer(ct)),
       file.path(BASE, "results/v3/deconv_fibroblast_clusters.csv"))

## ------------------------------------------------------------------- (D) ------------
TARG_G <- c("COL1A1","COL3A1","COL1A2","FN1","ACTA2","FAP","PDGFRB","ETS1","NFKB1","SP1","RELA","EPCAM","PTPRC","ESR1")
TARG_M <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c")
mir29  <- colMeans(mexp[c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"), SB, drop = FALSE])
CFV <- readRDS(file.path(BASE, "cache/v3/deconv/fibroblast_estimates_collagenfree.rds"))
long <- rbind(long, data.table(method = "collagen_free",
  feature = rep(colnames(CFV), each = nrow(CFV)), sample = rep(rownames(CFV), ncol(CFV)),
  value = as.vector(CFV)))
wide <- dcast(long, sample ~ method + feature, value.var = "value", sep = "::")
setkey(wide, sample); wide <- wide[SS]
FEATS <- setdiff(names(wide), "sample")
logf("\n=== (D) correlation of every cell-type estimate with markers; features = ", length(FEATS), " ===")
res <- rbindlist(lapply(FEATS, function(f) {
  v <- wide[[f]]
  a <- rbindlist(lapply(TARG_G, function(g) {
    ct <- suppressWarnings(cor.test(v, gexp[g, SS], method = "spearman", exact = FALSE))
    data.table(feature = f, target = g, target_type = "gene", n = length(SS),
               rho = unname(ct$estimate), p = ct$p.value) }))
  vb <- wide[[f]][match(SB, SS)]
  b <- rbindlist(lapply(TARG_M, function(g) {
    ct <- suppressWarnings(cor.test(vb, mexp[g, SB], method = "spearman", exact = FALSE))
    data.table(feature = f, target = g, target_type = "miRNA", n = length(SB),
               rho = unname(ct$estimate), p = ct$p.value) }))
  ct <- suppressWarnings(cor.test(vb, mir29, method = "spearman", exact = FALSE))
  c3 <- data.table(feature = f, target = "miR-29 family (mean)", target_type = "miRNA",
                   n = length(SB), rho = unname(ct$estimate), p = ct$p.value)
  rbind(a, b, c3)
}))
res[, `:=`(method = sub("::.*$", "", feature), cell_feature = sub("^.*::", "", feature))]
res[, fdr := p.adjust(p, "BH"), by = target]
setcolorder(res, c("method","cell_feature","target","target_type","n","rho","p","fdr"))
res[, feature := NULL]
setorder(res, target, -rho)
fwrite(res, file.path(BASE, "results/v3/deconv_celltype_target_correlations.csv"))
logf("WROTE deconv_celltype_target_correlations.csv rows = ", nrow(res))
for (tg in c("COL1A1","COL3A1","miR-29 family (mean)","hsa-miR-29a")) {
  s <- res[target == tg][order(-rho)]
  logf("\n--- ", tg, " : top 12 positively correlated cell-type estimates ---")
  for (i in 1:12) logf(sprintf("  %+.3f (FDR %8.2g)  %s :: %s", s$rho[i], s$fdr[i], s$method[i], s$cell_feature[i]))
  logf("--- ", tg, " : top 8 negatively correlated ---")
  for (i in (nrow(s)-7):nrow(s)) logf(sprintf("  %+.3f (FDR %8.2g)  %s :: %s", s$rho[i], s$fdr[i], s$method[i], s$cell_feature[i]))
}
logf("DONE 29")
