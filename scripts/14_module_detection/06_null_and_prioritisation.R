#!/usr/bin/env Rscript
# v7 / 06 : (a) characterise the consensus modules;
#           (b) degree-preserving rewiring null for the COL1A1/COL3A1/miR-29a triad;
#           (c) agreement between every informative partition and the 30
#               prioritised nodes (hypergeometric, and a concentration permutation test).
suppressMessages({library(igraph); library(data.table)})
setwd("/path/to/revision")
set.seed(20260912)

g <- readRDS("results/v7/cm2_graph_undirected.rds"); nm <- V(g)$name; N <- length(nm)
P <- readRDS("results/v7/cm2_all_partitions.rds")
s <- fread("results/v7/cm2_run_summary.csv"); info <- s[degenerate == FALSE, run]
nodes <- fread("data/canonical_nodes.tsv"); ntype <- setNames(nodes$type, nodes$name)
prior30 <- fread("results/v5/tables/Table4_prioritised_30.csv")$name
EMT_net <- intersect(unique(fread("results/v7/cm2_genesets_msigdb.csv")[
  gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", gene_symbol]), nm)
prot <- nm[ntype[nm] %in% c("Gene", "TF")]
cons <- fread("results/v7/cm2_consensus_clusters.csv")
cm <- setNames(cons$consensus_cluster, cons$name)

## ---- (a) consensus modules ----------------------------------------------
tb <- table(cm); keep <- as.integer(names(tb)[tb >= 5])
consmod <- rbindlist(lapply(keep, function(cl) {
  mem <- names(cm)[cm == cl]; pm <- intersect(mem, prot)
  sg <- induced_subgraph(g, mem)
  data.table(cluster = cl, size = length(mem),
             internal_edges = ecount(sg), density = edge_density(sg),
             n_gene = sum(ntype[mem] == "Gene"), n_TF = sum(ntype[mem] == "TF"),
             n_miRNA = sum(ntype[mem] == "miRNA"),
             EMT_hits = length(intersect(mem, EMT_net)),
             EMT_p = phyper(length(intersect(mem, EMT_net)) - 1, length(EMT_net),
                            length(prot) - length(EMT_net), length(pm), lower.tail = FALSE),
             prior30_hits = length(intersect(mem, prior30)),
             prior30_p = phyper(length(intersect(mem, prior30)) - 1, length(prior30),
                                N - length(prior30), length(mem), lower.tail = FALSE),
             members = paste(mem, collapse = " "))
}))
consmod[, EMT_q := p.adjust(EMT_p, "BH")][, prior30_q := p.adjust(prior30_p, "BH")]
setorder(consmod, -size)
fwrite(consmod, "results/v7/cm2_consensus_modules.csv")
cat("== consensus modules (>=5 nodes) ==\n")
print(consmod[, .(cluster, size, density = round(density, 3), n_gene, n_TF, n_miRNA,
                  EMT_hits, EMT_p = signif(EMT_p, 3), EMT_q = signif(EMT_q, 3),
                  prior30_hits, prior30_p = signif(prior30_p, 3))])
cat("\n")
for (i in seq_len(min(6, nrow(consmod))))
  cat(sprintf("module %d (n=%d): %s\n\n", consmod$cluster[i], consmod$size[i],
              consmod$members[i]))

## ---- (b) degree-preserving rewiring null --------------------------------
## Does the triad co-cluster more often than in a network with identical degree
## sequence but randomised wiring? Louvain and fast-greedy are used because they
## are the only members of the battery that produce a non-degenerate partition
## cheaply enough to permute.
trio <- c("COL1A1", "COL3A1", "hsa-miR-29a")
obs <- sapply(c("Louvain_r1.00", "FastGreedy", "Leiden_mod_r1.0", "cm2_glay_cm2"),
              function(k) { m <- P[[k]]; as.integer(m[trio[1]] == m[trio[2]] & m[trio[1]] == m[trio[3]]) })
B <- 500
null_lou <- integer(B); null_fg <- integer(B); null_q_lou <- numeric(B)
for (b in seq_len(B)) {
  gr <- rewire(g, keeping_degseq(niter = 10 * ecount(g)))
  ml <- cluster_louvain(gr); mm <- membership(ml)
  null_lou[b] <- as.integer(mm[trio[1]] == mm[trio[2]] & mm[trio[1]] == mm[trio[3]])
  null_q_lou[b] <- modularity(ml)
  mf <- membership(cluster_fast_greedy(gr))
  null_fg[b] <- as.integer(mf[trio[1]] == mf[trio[2]] & mf[trio[1]] == mf[trio[3]])
}
nulldt <- data.table(
  test = c("Louvain: all three in one module", "Fast-greedy: all three in one module",
           "Louvain modularity"),
  observed = c(obs["Louvain_r1.00"], obs["FastGreedy"], s[run == "Louvain_r1.00", modularity]),
  null_mean = c(mean(null_lou), mean(null_fg), mean(null_q_lou)),
  null_sd = c(sd(null_lou), sd(null_fg), sd(null_q_lou)),
  p_empirical = c((sum(null_lou >= obs["Louvain_r1.00"]) + 1) / (B + 1),
                  (sum(null_fg  >= obs["FastGreedy"])     + 1) / (B + 1),
                  (sum(null_q_lou >= s[run == "Louvain_r1.00", modularity]) + 1) / (B + 1)),
  B = B)
fwrite(nulldt, "results/v7/cm2_rewiring_null.csv")
cat("\n== degree-preserving rewiring null (B=500) ==\n"); print(nulldt)

## ---- (c) agreement with the 30 prioritised nodes -------------------------
agree <- rbindlist(lapply(info, function(k) {
  m <- P[[k]]
  per <- rbindlist(lapply(sort(unique(m)), function(cl) {
    mem <- names(m)[m == cl]; h <- length(intersect(mem, prior30))
    data.table(cluster = cl, size = length(mem), hits = h,
               p = phyper(h - 1, length(prior30), N - length(prior30),
                          length(mem), lower.tail = FALSE))
  }))
  per[, q := p.adjust(p, "BH")]
  best <- per[order(p)][1]
  # concentration: observed co-clustering pairs among the 30 vs. permutation
  obs_pairs <- sum(outer(m[prior30], m[prior30], "==")) - length(prior30)
  perm <- replicate(2000, { idx <- sample(nm, length(prior30))
    sum(outer(m[idx], m[idx], "==")) - length(prior30) })
  data.table(run = k, k_clusters = length(unique(m)),
             best_cluster = best$cluster, best_cluster_size = best$size,
             prior30_in_best = best$hits, hyper_p = best$p, hyper_q = best$q,
             n_clusters_with_any = per[hits > 0, .N],
             cocluster_pairs_obs = obs_pairs / 2,
             cocluster_pairs_null = mean(perm) / 2,
             concentration_p = (sum(perm >= obs_pairs) + 1) / 2001)
}))
setorder(agree, hyper_p)
fwrite(agree, "results/v7/cm2_prioritisation_agreement.csv")
cat("\n== agreement with the 30 prioritised nodes (informative partitions) ==\n")
print(agree[, .(run, k_clusters, best_cluster_size, prior30_in_best,
                hyper_p = signif(hyper_p, 3), hyper_q = signif(hyper_q, 3),
                n_clusters_with_any, cocluster_pairs_obs,
                cocluster_pairs_null = round(cocluster_pairs_null, 1),
                concentration_p)], nrows = 50)
