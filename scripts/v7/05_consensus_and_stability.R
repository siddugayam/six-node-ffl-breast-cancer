#!/usr/bin/env Rscript
# v7 / 05 : degeneracy filtering, stochastic-seed stability, consensus
#           co-clustering, and the two headline tests:
#             does the stroma module exist, and does it match the prioritised 30?
suppressMessages({library(igraph); library(data.table); library(mclust)})
setwd("/path/to/revision")

g  <- readRDS("results/v7/cm2_graph_undirected.rds"); nm <- V(g)$name; N <- length(nm)
P  <- readRDS("results/v7/cm2_all_partitions.rds")
s  <- fread("results/v7/cm2_run_summary.csv")
nodes <- fread("data/canonical_nodes.tsv"); ntype <- setNames(nodes$type, nodes$name)
prior30 <- fread("results/v5/tables/Table4_prioritised_30.csv")$name
EMT <- unique(fread("results/v7/cm2_genesets_msigdb.csv")[
  gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", gene_symbol])
EMT_net <- intersect(EMT, nm)
prot <- nm[ntype[nm] %in% c("Gene", "TF")]

## ---- degeneracy rule, fixed in advance ----------------------------------
## A partition is uninformative if it is one cluster, all singletons, or has a
## giant component holding >60% of the network. Everything else is "informative".
s[, degenerate := k_clusters < 3 | largest > 0.60 * N | k_clusters > 0.90 * N]
fwrite(s, "results/v7/cm2_run_summary.csv")
info <- s[degenerate == FALSE, run]
cat(sprintf("informative partitions: %d of %d\n", length(info), nrow(s)))
cat(paste(info, collapse = ", "), "\n\n")

## ---- one representative run per distinct algorithm x implementation ------
rep_runs <- c(
  "MCL (igraph impl., I=2.5)"        = "MCL_I2.5",
  "MCL (clusterMaker2, I=2.5)"       = "cm2_mcl_I25",
  "Louvain (igraph, r=1)"            = "Louvain_r1.00",
  "Louvain (clusterMaker2)"          = "cm2_louvain_cm2",
  "Walktrap (igraph, steps=4)"       = "Walktrap_s4",
  "Infomap (igraph, 10 trials)"      = "Infomap_t10",
  "Infomap (clusterMaker2)"          = "cm2_infomap_cm2_t10",
  "Label propagation (igraph)"       = "LabelProp_best",
  "Label propagation (clusterMaker2)"= "cm2_labelprop_cm2",
  "Fast-greedy (igraph)"             = "FastGreedy",
  "Fast-greedy (clusterMaker2)"      = "cm2_fastgreedy_cm2",
  "GLay (clusterMaker2)"             = "cm2_glay_cm2",
  "Leading eigenvector (igraph)"     = "LeadingEigen",
  "Leiden (igraph, r=1)"             = "Leiden_mod_r1.0",
  "Girvan-Newman (igraph)"           = "EdgeBetweenness")
rep_runs <- rep_runs[rep_runs %in% names(P)]

trio <- c("COL1A1", "COL3A1", "hsa-miR-29a")
tab <- rbindlist(lapply(seq_along(rep_runs), function(i) {
  k <- rep_runs[i]; m <- P[[k]]
  cl <- m["COL1A1"]; mem <- names(m)[m == cl]
  deg <- s[run == k, degenerate]
  pm <- intersect(mem, prot)
  data.table(
    method = names(rep_runs)[i], run = k,
    k_clusters = s[run == k, k_clusters], modularity = round(s[run == k, modularity], 4),
    largest = s[run == k, largest], degenerate = deg,
    COL1A1_module_size = length(mem),
    COL3A1_with_COL1A1 = unname(m["COL3A1"] == cl),
    miR29a_with_COL1A1 = unname(m["hsa-miR-29a"] == cl),
    all_three_together = unname(m["COL3A1"] == cl & m["hsa-miR-29a"] == cl),
    EMT_hits = length(intersect(mem, EMT_net)),
    EMT_p = signif(phyper(length(intersect(mem, EMT_net)) - 1, length(EMT_net),
                          length(prot) - length(EMT_net), length(pm),
                          lower.tail = FALSE), 3),
    prior30_hits = length(intersect(mem, prior30)),
    prior30_p = signif(phyper(length(intersect(mem, prior30)) - 1, length(prior30),
                              N - length(prior30), length(mem), lower.tail = FALSE), 3))
}))
fwrite(tab, "results/v7/cm2_TABLE_method_battery.csv")
print(tab[, .(method, k_clusters, modularity, degenerate, COL1A1_module_size,
              COL3A1_with_COL1A1, miR29a_with_COL1A1, EMT_hits, EMT_p,
              prior30_hits, prior30_p)])

## ---- consensus co-clustering over INFORMATIVE partitions only -----------
infor <- intersect(unname(rep_runs), info)
cat("\ninformative representative runs used for consensus:", length(infor), "\n")
pairs <- data.table(a = c("COL1A1","COL1A1","COL3A1"), b = c("COL3A1","hsa-miR-29a","hsa-miR-29a"))
pairs[, n_together := mapply(function(x, y) sum(sapply(infor, function(k) P[[k]][x] == P[[k]][y])), a, b)]
pairs[, n_runs := length(infor)]
pairs[, runs := mapply(function(x, y) paste(infor[sapply(infor, function(k) P[[k]][x] == P[[k]][y])],
                                            collapse = ";"), a, b)]
fwrite(pairs, "results/v7/cm2_trio_cocluster_informative.csv")
cat("\n== trio co-clustering, informative partitions only ==\n"); print(pairs[, .(a, b, n_together, n_runs, runs)])

## null expectation for a pair co-clustering by chance in each partition
null_rate <- sapply(infor, function(k) { m <- P[[k]]; tb <- table(m); sum(tb * (tb - 1)) / (N * (N - 1)) })
cat("\nchance co-clustering rate per informative partition:\n"); print(round(null_rate, 4))
cat("mean chance rate:", round(mean(null_rate), 4), "\n")

## ---- stochastic stability of the one algorithm family that does partition -
seeds <- readRDS("results/v7/cm2_seed_replicates.rds")
lou <- seeds$louvain
ari_lou <- outer(seq_along(lou), seq_along(lou), Vectorize(function(i, j)
  adjustedRandIndex(lou[[i]], lou[[j]])))
cat(sprintf("\nLouvain, 50 random seeds: mean pairwise ARI %.3f (range %.3f-%.3f), k range %d-%d\n",
            mean(ari_lou[upper.tri(ari_lou)]), min(ari_lou[upper.tri(ari_lou)]),
            max(ari_lou[upper.tri(ari_lou)]),
            min(sapply(lou, function(m) length(unique(m)))),
            max(sapply(lou, function(m) length(unique(m))))))
iC1 <- match("COL1A1", nm); iC3 <- match("COL3A1", nm); i29 <- match("hsa-miR-29a", nm)
lou_tog <- data.table(
  pair = c("COL1A1-COL3A1", "COL1A1-miR-29a", "COL3A1-miR-29a", "all three"),
  frac_of_50_seeds = c(mean(sapply(lou, function(m) m[iC1] == m[iC3])),
                       mean(sapply(lou, function(m) m[iC1] == m[i29])),
                       mean(sapply(lou, function(m) m[iC3] == m[i29])),
                       mean(sapply(lou, function(m) m[iC1] == m[iC3] & m[iC1] == m[i29]))))
fwrite(lou_tog, "results/v7/cm2_louvain_seed_stability.csv")
cat("\n== Louvain across 50 seeds ==\n"); print(lou_tog)

## ---- global consensus matrix over informative partitions, then cluster it -
C <- matrix(0, N, N, dimnames = list(nm, nm))
for (k in info) { m <- P[[k]]; C <- C + outer(m, m, "==") }
C <- C / length(info)
saveRDS(C, "results/v7/cm2_consensus_matrix.rds")
gc_ <- graph_from_adjacency_matrix(C >= 0.5, mode = "undirected", diag = FALSE)
cons <- components(gc_)$membership
cat(sprintf("\nconsensus at >=50%% of %d informative partitions: %d components, largest %d, singletons %d\n",
            length(info), max(cons), max(table(cons)), sum(table(cons) == 1)))
cons_dt <- data.table(name = nm, consensus_cluster = cons, type = ntype[nm])
fwrite(cons_dt, "results/v7/cm2_consensus_clusters.csv")
cat("consensus cluster of COL1A1 / COL3A1 / miR-29a:",
    cons["COL1A1"], cons["COL3A1"], cons["hsa-miR-29a"], "\n")
cat("COL1A1 consensus-cluster members:",
    paste(nm[cons == cons["COL1A1"]], collapse = " "), "\n")
cat("consensus co-clustering fraction COL1A1-COL3A1:", round(C["COL1A1","COL3A1"],3),
    " COL1A1-miR29a:", round(C["COL1A1","hsa-miR-29a"],3),
    " COL3A1-miR29a:", round(C["COL3A1","hsa-miR-29a"],3), "\n")
