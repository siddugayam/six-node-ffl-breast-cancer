#!/usr/bin/env Rscript
# Does the collagen signal survive adjustment for stromal / CAF content?
suppressMessages({library(ppcor)})
set.seed(42)
RES  <- "/path/to/revision/results/multiomics/"
DATA <- "/path/to/revision/data/"
CACHE<- "/path/to/revision/cache/celltype/"

expr <- readRDS(paste0(DATA,"brca_gene_expr.rds"))
mir  <- readRDS(paste0(DATA,"brca_mirna_expr_canonical.rds"))
sc   <- readRDS(paste0(CACHE,"tcga_stromal_scores.rds"))
rownames(sc) <- sc$sample
covs <- c("ssgsea_stromal","ssgsea_stromal_noCOL","score_CAF_marker","score_CAF_scRNA")

s_gene <- intersect(colnames(expr), sc$sample)
s_both <- intersect(s_gene, colnames(mir))
cat("tumours with expression + stromal score:", length(s_gene), "\n")
cat("tumours with expression + miRNA + stromal score:", length(s_both), "\n")

## partial Spearman for one pair, given one covariate vector
pspear <- function(x, y, z) {
  ok <- complete.cases(x,y,z); if (sum(ok) < 20) return(c(NA,NA,NA))
  r <- suppressWarnings(pcor.test(x[ok], y[ok], z[ok], method="spearman"))
  c(r$estimate, r$p.value, sum(ok))
}
pspear_multi <- function(x, y, Z) {
  ok <- complete.cases(x,y,Z); if (sum(ok) < 20) return(c(NA,NA,NA))
  r <- suppressWarnings(pcor.test(x[ok], y[ok], Z[ok,,drop=FALSE], method="spearman"))
  c(r$estimate, r$p.value, sum(ok))
}
getv <- function(node, samples) {
  if (node %in% rownames(expr)) return(as.numeric(expr[node, samples]))
  if (node %in% rownames(mir))  return(as.numeric(mir[node, samples]))
  rep(NA_real_, length(samples))
}

## ---------------- pairs to test ----------------
named <- read.csv("/path/to/revision/results/named_axis_correlations.csv",
                  stringsAsFactors=FALSE)[,c("source","target")]
edges <- read.delim(paste0(DATA,"canonical_edges.tsv"), stringsAsFactors=FALSE)
coll  <- edges[edges$target %in% c("COL1A1","COL3A1"), c("source","target")]
extra <- data.frame(source=c("NFKB1","RELA","SP1","ETS1","COL1A1","EZH2","MYC","TP53"),
                    target=c("COL3A1","COL3A1","COL3A1","COL3A1","COL3A1","COL1A1","COL1A1","COL1A1"),
                    stringsAsFactors=FALSE)
pairs <- unique(rbind(named, coll, extra))
pairs <- pairs[pairs$source != pairs$target, ]
cat("pairs to test:", nrow(pairs), "\n")

out <- list()
for (i in seq_len(nrow(pairs))) {
  a <- pairs$source[i]; b <- pairs$target[i]
  smp <- if (grepl("^hsa-", a) || grepl("^hsa-", b)) s_both else s_gene
  x <- getv(a, smp); y <- getv(b, smp)
  if (all(is.na(x)) || all(is.na(y))) { cat("SKIP (not measured):", a, "->", b, "\n"); next }
  ok <- complete.cases(x,y)
  base <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
  row <- data.frame(source=a, target=b, n=sum(ok),
                    rho_unadjusted=unname(base$estimate), p_unadjusted=base$p.value,
                    stringsAsFactors=FALSE)
  for (cv in covs) {
    r <- pspear(x, y, sc[smp, cv])
    row[[paste0("rho_adj_",cv)]] <- r[1]; row[[paste0("p_adj_",cv)]] <- r[2]
  }
  Z <- sc[smp, c("score_CAF_scRNA","score_immune","score_epithelial","score_endothelial")]
  r <- pspear_multi(x, y, Z)
  row$rho_adj_CAF_immune_epi_endo <- r[1]; row$p_adj_CAF_immune_epi_endo <- r[2]
  row$pct_rho_retained_CAFscRNA <- 100*r[1]/row$rho_unadjusted
  out[[length(out)+1]] <- row
}
res <- do.call(rbind, out)
res$pct_rho_retained_CAFscRNA <- round(100*res$rho_adj_score_CAF_scRNA/res$rho_unadjusted, 1)
write.csv(res, paste0(RES,"stromal_adjusted_correlations.csv"), row.names=FALSE)
cat("\nstromal_adjusted_correlations.csv rows:", nrow(res), "\n\n")
show <- res[, c("source","target","n","rho_unadjusted","p_unadjusted",
                "rho_adj_score_CAF_scRNA","p_adj_score_CAF_scRNA",
                "rho_adj_ssgsea_stromal_noCOL","p_adj_ssgsea_stromal_noCOL",
                "rho_adj_CAF_immune_epi_endo","p_adj_CAF_immune_epi_endo")]
print(format(show, digits=3), row.names=FALSE)

## ---------------- every node vs stromal content ----------------
nodes <- read.delim(paste0(DATA,"canonical_nodes.tsv"), stringsAsFactors=FALSE)
nv <- list()
for (i in seq_len(nrow(nodes))) {
  nm <- nodes$name[i]
  smp <- if (nodes$type[i]=="miRNA") s_both else s_gene
  v <- getv(nm, smp); if (all(is.na(v))) next
  z <- sc[smp,"score_CAF_scRNA"]; z2 <- sc[smp,"score_epithelial"]
  ct1 <- suppressWarnings(cor.test(v, z, method="spearman", exact=FALSE))
  ct2 <- suppressWarnings(cor.test(v, z2, method="spearman", exact=FALSE))
  nv[[length(nv)+1]] <- data.frame(node=nm, type=nodes$type[i], n=sum(complete.cases(v,z)),
      rho_vs_CAF_score=unname(ct1$estimate), p_vs_CAF=ct1$p.value,
      rho_vs_epithelial_score=unname(ct2$estimate), p_vs_epithelial=ct2$p.value,
      stringsAsFactors=FALSE)
}
nvv <- do.call(rbind, nv)
nvv$fdr_vs_CAF <- p.adjust(nvv$p_vs_CAF, "BH")
nvv <- nvv[order(-nvv$rho_vs_CAF_score),]
write.csv(nvv, paste0(RES,"node_vs_stromal_correlation.csv"), row.names=FALSE)
cat("\nnode_vs_stromal_correlation.csv rows:", nrow(nvv), "\n")
cat("nodes tested:", table(nvv$type), "\n")
cat("\nTop 15 CAF-tracking nodes:\n"); print(head(nvv[,c("node","type","rho_vs_CAF_score","fdr_vs_CAF")],15), row.names=FALSE)
cat("\nFeatured nodes:\n")
print(nvv[nvv$node %in% c("COL1A1","COL3A1","NFKB1","RELA","SP1","ETS1","VEGFA","EZH2","TP53","MYC",
   "hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","hsa-miR-101","hsa-let-7b","hsa-let-7e",
   "hsa-miR-145","hsa-miR-200c","hsa-miR-124"), c("node","type","n","rho_vs_CAF_score","p_vs_CAF","fdr_vs_CAF","rho_vs_epithelial_score")], row.names=FALSE)
saveRDS(list(s_gene=s_gene, s_both=s_both), paste0(CACHE,"ct06_samples.rds"))
