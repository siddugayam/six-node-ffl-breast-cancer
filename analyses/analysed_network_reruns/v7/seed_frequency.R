## How often do COL1A1, COL3A1 and miR-29a share a module across the stochastic seeds that
## 04_igraph_battery.R stores (50 Louvain, 20 Infomap, 100 label-propagation), with and without the
## 30 legacy edges.  Also the size of the shared module, since a 145-node module is not a specific one.
suppressMessages({library(igraph); library(data.table)})
HERE <- "/path/to/revision/analyses/analysed_network_reruns/v7"
trio <- c("COL1A1", "COL3A1", "hsa-miR-29a")
res <- rbindlist(lapply(c("full", "nolegacy"), function(mode) {
  S <- file.path(HERE, paste0("sandbox_", mode))
  g <- readRDS(file.path(S, "results/v7/cm2_graph_undirected.rds")); nm <- V(g)$name
  R <- readRDS(file.path(S, "results/v7/cm2_seed_replicates.rds"))
  rbindlist(lapply(names(R), function(alg) {
    x <- t(sapply(R[[alg]], function(m) { names(m) <- nm; cl <- m["COL1A1"]
      c(together = all(m[trio] == cl), size = sum(m == cl), k = length(unique(m))) }))
    data.table(network = mode, algorithm = alg, seeds = nrow(x), trio_together = sum(x[, "together"]),
               median_COL1A1_module_size_when_together = if (any(x[, "together"] == 1)) median(x[x[, "together"] == 1, "size"]) else NA_real_,
               median_k = median(x[, "k"]))
  }))
}))
fwrite(res, file.path(HERE, "cd_seed_frequency.csv")); print(res)
