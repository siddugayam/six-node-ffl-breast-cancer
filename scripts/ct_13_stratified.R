#!/usr/bin/env Rscript
# Stratified check: partial correlation on a covariate as collinear as the CAF score
# (rho 0.88 with COL1A1) can over-adjust. Repeat the key tests WITHIN tertiles of CAF
# content, where composition is approximately constant, as a non-parametric alternative.
set.seed(42)
RES  <- "/path/to/revision/results/multiomics/"
DATA <- "/path/to/revision/data/"
CACHE<- "/path/to/revision/cache/celltype/"
expr <- readRDS(paste0(DATA,"brca_gene_expr.rds")); mir <- readRDS(paste0(DATA,"brca_mirna_expr_canonical.rds"))
sc <- readRDS(paste0(CACHE,"tcga_stromal_scores.rds")); rownames(sc) <- sc$sample
ph <- readRDS(paste0(DATA,"brca_pheno.rds")); tum <- ph$sample[ph$sample_type=="Primary Tumor"]
sg <- intersect(intersect(colnames(expr), tum), sc$sample); sb <- intersect(sg, colnames(mir))
getv <- function(n, s) if (n %in% rownames(expr)) as.numeric(expr[n,s]) else as.numeric(mir[n,s])
pairs <- data.frame(
 source=c("COL1A1","NFKB1","RELA","SP1","ETS1","NFKB1","RELA","SP1","ETS1",
          "hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","hsa-miR-29a","hsa-miR-29b","hsa-miR-29c",
          "hsa-let-7b","hsa-let-7e","hsa-miR-101"),
 target=c("COL3A1","COL1A1","COL1A1","COL1A1","COL1A1","COL3A1","COL3A1","COL3A1","COL3A1",
          "COL1A1","COL1A1","COL1A1","COL3A1","COL3A1","COL3A1","COL3A1","COL3A1","EZH2"),
 stringsAsFactors=FALSE)
out <- list()
for (i in seq_len(nrow(pairs))) {
  a <- pairs$source[i]; b <- pairs$target[i]
  s <- if (grepl("^hsa-", a)) sb else sg
  x <- getv(a,s); y <- getv(b,s); z <- sc[s,"score_CAF_scRNA"]
  ter <- cut(z, quantile(z, c(0,1/3,2/3,1)), include.lowest=TRUE, labels=c("low","mid","high"))
  rr <- sapply(levels(ter), function(L){ k <- ter==L
    ct <- suppressWarnings(cor.test(x[k], y[k], method="spearman", exact=FALSE)); c(ct$estimate, ct$p.value, sum(k)) })
  # Fisher combination of the three within-stratum p-values
  chi <- -2*sum(log(pmax(rr[2,], 1e-300))); pf <- pchisq(chi, df=6, lower.tail=FALSE)
  all_s <- suppressWarnings(cor.test(x, y, method="spearman", exact=FALSE))
  out[[i]] <- data.frame(source=a, target=b, n_total=length(s), rho_all=unname(all_s$estimate),
    p_all=all_s$p.value,
    rho_CAFlow=rr[1,1], p_CAFlow=rr[2,1], n_low=rr[3,1],
    rho_CAFmid=rr[1,2], p_CAFmid=rr[2,2], n_mid=rr[3,2],
    rho_CAFhigh=rr[1,3], p_CAFhigh=rr[2,3], n_high=rr[3,3],
    mean_within_stratum_rho=mean(rr[1,]), fisher_combined_p=pf,
    all_three_strata_same_sign_as_whole=all(sign(rr[1,])==sign(all_s$estimate)),
    stringsAsFactors=FALSE)
}
res <- do.call(rbind, out)
write.csv(res, paste0(RES,"stromal_stratified_correlations.csv"), row.names=FALSE)
cat("rows:", nrow(res), "\n")
print(res[,c("source","target","rho_all","rho_CAFlow","rho_CAFmid","rho_CAFhigh",
             "mean_within_stratum_rho","fisher_combined_p","all_three_strata_same_sign_as_whole")],
      row.names=FALSE, digits=3)
