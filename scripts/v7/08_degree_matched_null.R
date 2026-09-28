#!/usr/bin/env Rscript
# v7 / 08 : the ECM-concentration test of script 07 repeated with a DEGREE-MATCHED
#           null. High-degree nodes are preferentially swallowed by large modules,
#           so a uniform null over-states concentration; this is the conservative
#           version. Also compiles the headline agreement numbers.
suppressMessages({library(igraph); library(data.table)})
setwd("/path/to/revision")
set.seed(20260912)

g <- readRDS("results/v7/cm2_graph_undirected.rds"); nm <- V(g)$name; N <- length(nm)
P <- readRDS("results/v7/cm2_all_partitions.rds"); s <- fread("results/v7/cm2_run_summary.csv")
nodes <- fread("data/canonical_nodes.tsv"); ntype <- setNames(nodes$type, nodes$name)
prot <- nm[ntype[nm] %in% c("Gene", "TF")]
deg <- degree(g)
EMT <- intersect(unique(fread("results/v7/cm2_genesets_msigdb.csv")[
  gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", gene_symbol]), nm)
ECM <- intersect(c("COL1A1","COL3A1","COL7A1","COL18A1","FN1","SPARC","POSTN","LOX",
                   "MMP2","MMP14","THBS1","THBS2","TGFB2","TGFBI","SERPINE1","ACTA2",
                   "TAGLN","VIM","CDH11","FBN1","CTGF","PDGFRB","ITGA5","ITGB1"), nm)

## degree strata over protein-coding nodes
brk <- unique(quantile(deg[prot], probs = seq(0, 1, 0.1)))
bin <- cut(deg[prot], breaks = brk, include.lowest = TRUE, labels = FALSE)
names(bin) <- prot
pool <- split(prot, bin)
draw_matched <- function(set) {
  want <- table(bin[set])
  unlist(lapply(names(want), function(b) sample(pool[[b]], want[[b]], replace = FALSE)))
}
conc_matched <- function(m, set, B = 5000) {
  obs <- max(table(m[set]))
  null <- replicate(B, max(table(m[draw_matched(set)])))
  c(obs = obs, exp = mean(null), p = (sum(null >= obs) + 1) / (B + 1))
}

runs <- s[degenerate == FALSE & k_clusters <= 30, run]
res <- rbindlist(lapply(runs, function(k) {
  m <- P[[k]]
  e <- conc_matched(m, EMT); c_ <- conc_matched(m, ECM)
  data.table(run = k, modularity = s[run == k, modularity],
             EMT_obs = e["obs"], EMT_exp = round(e["exp"], 1), EMT_p_degmatched = e["p"],
             ECM_obs = c_["obs"], ECM_exp = round(c_["exp"], 1), ECM_p_degmatched = c_["p"])
}))
res[, `:=`(EMT_q = p.adjust(EMT_p_degmatched, "BH"), ECM_q = p.adjust(ECM_p_degmatched, "BH"))]
setorder(res, -modularity)
fwrite(res, "results/v7/cm2_concentration_degreematched.csv")
cat("== degree-matched concentration null ==\n"); print(res, nrows = 30)

## headline agreement numbers
ag <- fread("results/v7/cm2_prioritisation_agreement.csv")
cat(sprintf("\ninformative partitions: %d\n", nrow(ag)))
cat(sprintf("partitions with ANY module enriched for the 30 at BH q<0.05: %d (%.0f%%)\n",
            sum(ag$hyper_q < 0.05), 100 * mean(ag$hyper_q < 0.05)))
cat(sprintf("partitions where the 30 co-cluster more than random (perm p<0.05): %d (%.0f%%)\n",
            sum(ag$concentration_p < 0.05), 100 * mean(ag$concentration_p < 0.05)))
cat("best single overlap across all informative partitions:\n")
print(ag[order(hyper_p)][1:5, .(run, best_cluster_size, prior30_in_best, hyper_p, hyper_q)])

## same for the modularity-based subset only
mb <- ag[run %in% runs]
cat(sprintf("\nmodularity-based partitions only (n=%d): q<0.05 in %d; concentration p<0.05 in %d\n",
            nrow(mb), sum(mb$hyper_q < 0.05), sum(mb$concentration_p < 0.05)))
print(mb[order(hyper_p)][, .(run, best_cluster_size, prior30_in_best,
                             hyper_p = signif(hyper_p,3), hyper_q = signif(hyper_q,3),
                             concentration_p)])
