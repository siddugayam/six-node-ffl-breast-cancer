## Community-detection claims of Supplementary Note S6 / Table S9, with and without the 30 legacy
## GeneMANIA miRNA-miRNA edges.  Partitions come from scripts/14_module_detection/04_igraph_battery.R run unchanged in
## sandbox_<mode>/ (setup_sandbox.sh); this script repeats the logic of 05_analysis.R (run summary),
## 06_consensus_and_stability.R (degeneracy rule), 08_stroma_module_test.R (COL1A1-module membership),
## 07_null_and_prioritisation.R (b: 500-rewiring null; c: prioritised-30 overlap) on those partitions.
## The four clusterMaker2 (Cytoscape) runs can only be included for the full network, from the stored
## partition file; Cytoscape is not available here to repeat them without the legacy edges.
## usage: Rscript cd_audit.R <full|nolegacy>
suppressMessages({library(igraph); library(data.table)})
mode <- commandArgs(TRUE)[1]; stopifnot(mode %in% c("full", "nolegacy"))
REV  <- "/path/to/revision"
HERE <- file.path(REV, "analyses/analysed_network_reruns/v7"); S <- file.path(HERE, paste0("sandbox_", mode))
g <- readRDS(file.path(S, "results/v7/cm2_graph_undirected.rds")); nm <- V(g)$name; N <- length(nm)
nodes <- fread(file.path(REV, "data/canonical_nodes.tsv")); ntype <- setNames(nodes$type, nodes$name)
prior30 <- fread(file.path(REV, "results/v5/tables/Table4_prioritised_30.csv"))$name
out <- list(mode = mode, nodes = N, undirected_edges = ecount(g))

ig <- fread(file.path(S, "results/v7/cm2_igraph_partitions.csv")); setkey(ig, name); ig <- ig[nm]
P <- as.list(ig[, -"name"])
if (mode == "full") {
  ig0 <- fread(file.path(REV, "results/v7/cm2_igraph_partitions.csv")); setkey(ig0, name); ig0 <- ig0[nm]
  same <- vapply(names(P), function(k) identical(as.integer(P[[k]]), as.integer(ig0[[k]])), TRUE)
  out$igraph_partitions_identical_to_stored <- sprintf("%d of %d", sum(same), length(same))
  if (!all(same)) out$igraph_runs_not_reproduced <- paste(names(same)[!same], collapse = ";")
  cy <- fread(file.path(REV, "results/v7/cm2_cytoscape_partitions.tsv"), na.strings = c("", "NA", "None"))
  setnames(cy, "name", "node"); setkey(cy, node); cy <- cy[nm]
  for (cc in setdiff(names(cy), c("SUID", "node", "selected", "shared name"))) {
    v <- cy[[cc]]; na <- is.na(v); v <- as.character(v); v[na] <- paste0("singleton_", seq_len(sum(na)))
    P[[paste0("cm2_", cc)]] <- as.integer(factor(v))
  }
}
for (k in names(P)) names(P[[k]]) <- nm

summ <- rbindlist(lapply(names(P), function(k) {
  m <- P[[k]]; tb <- sort(table(m), decreasing = TRUE)
  data.table(run = k, k_clusters = length(tb), modularity = modularity(g, m), largest = as.integer(tb[1]))
}))
summ[, degenerate := k_clusters < 3 | largest > 0.60 * N | k_clusters > 0.90 * N]
runs <- summ[degenerate == FALSE & k_clusters <= 30, run]
trio <- c("COL1A1", "COL3A1", "hsa-miR-29a"); TFarm <- c("ETS1", "NFKB1", "RELA", "SP1")
tab <- rbindlist(lapply(runs, function(k) {
  m <- P[[k]]; cl <- m["COL1A1"]
  data.table(run = k, cytoscape = startsWith(k, "cm2_"), k_clusters = summ[run == k, k_clusters],
             modularity = round(summ[run == k, modularity], 4), COL1A1_module_size = sum(m == cl),
             COL3A1_in = unname(m["COL3A1"] == cl), miR29a_in = unname(m["hsa-miR-29a"] == cl),
             miR29b_in = unname(m["hsa-miR-29b"] == cl), miR29c_in = unname(m["hsa-miR-29c"] == cl),
             trio_together = all(m[trio] == cl), TFarm_all_in = all(m[TFarm] == cl),
             TFarm_any_in = any(m[TFarm] == cl))
}))
setorder(tab, -modularity)
fwrite(tab, file.path(HERE, sprintf("cd_trio_by_run_%s.csv", mode)))
ig_tab <- tab[cytoscape == FALSE]
out$usable_runs <- nrow(tab); out$usable_runs_igraph <- nrow(ig_tab)
out$trio_together_all_runs <- sprintf("%d of %d", sum(tab$trio_together), nrow(tab))
out$trio_together_igraph_runs <- sprintf("%d of %d", sum(ig_tab$trio_together), nrow(ig_tab))
out$miR29a_in_COL1A1_module_igraph <- sprintf("%d of %d", sum(ig_tab$miR29a_in), nrow(ig_tab))
out$TFarm_all_with_COL1A1_igraph <- sprintf("%d of %d", sum(ig_tab$TFarm_all_in), nrow(ig_tab))

## prioritised-30 overlap, best module per run (06 part c, hypergeometric)
agree <- rbindlist(lapply(runs, function(k) {
  m <- P[[k]]
  per <- rbindlist(lapply(sort(unique(m)), function(cl) {
    mem <- names(m)[m == cl]; h <- length(intersect(mem, prior30))
    data.table(cluster = cl, size = length(mem), hits = h, has_COL1A1 = "COL1A1" %in% mem,
               has_VEGFA = "VEGFA" %in% mem, has_miR200 = any(grepl("^hsa-miR-200", mem)),
               p = phyper(h - 1, length(prior30), N - length(prior30), length(mem), lower.tail = FALSE))
  }))
  b <- per[order(p)][1]
  data.table(run = k, best_size = b$size, prior30_in_best = b$hits, hyper_p = signif(b$p, 3),
             best_has_COL1A1 = b$has_COL1A1, best_has_VEGFA = b$has_VEGFA, best_has_miR200 = b$has_miR200)
}))
fwrite(agree, file.path(HERE, sprintf("cd_prior30_best_module_%s.csv", mode)))

## degree-preserving rewiring null for the trio (06 part b), same seed and design
set.seed(20260912)
B <- 500; null_lou <- integer(B); null_fg <- integer(B)
for (b in seq_len(B)) {
  gr <- rewire(g, keeping_degseq(niter = 10 * ecount(g)))
  mm <- membership(cluster_louvain(gr)); null_lou[b] <- as.integer(mm[trio[1]] == mm[trio[2]] & mm[trio[1]] == mm[trio[3]])
  mf <- membership(cluster_fast_greedy(gr)); null_fg[b] <- as.integer(mf[trio[1]] == mf[trio[2]] & mf[trio[1]] == mf[trio[3]])
}
obs <- sapply(c("Louvain_r1.00", "FastGreedy"), function(k) { m <- P[[k]]; as.integer(all(m[trio] == m[trio[1]])) })
out$rewiring_null <- list(
  louvain = list(observed = obs[["Louvain_r1.00"]], null_mean = mean(null_lou), p = (sum(null_lou >= obs[["Louvain_r1.00"]]) + 1) / (B + 1)),
  fastgreedy = list(observed = obs[["FastGreedy"]], null_mean = mean(null_fg), p = (sum(null_fg >= obs[["FastGreedy"]]) + 1) / (B + 1)))
jsonlite::write_json(out, file.path(HERE, sprintf("cd_audit_%s.json", mode)), auto_unbox = TRUE, pretty = TRUE)
print(tab[, .(run, k_clusters, modularity, COL1A1_module_size, COL3A1_in, miR29a_in, trio_together, TFarm_all_in)], nrows = 40)
print(agree, nrows = 40)
str(out)
