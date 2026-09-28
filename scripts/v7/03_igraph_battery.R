#!/usr/bin/env Rscript
# v7 / 03 : the same battery run independently in igraph 2.x, so that every
#           clusterMaker2 result has a second implementation behind it, and so
#           that walktrap (absent from clusterMaker2 2.3.4) is covered.
suppressMessages({library(igraph); library(data.table)})
setwd("/path/to/revision")
source("scripts/v7/01_cm2_mcl.R")
set.seed(20260912)

g  <- readRDS("results/v7/cm2_graph_undirected.rds")
nm <- V(g)$name
A  <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
memb_list <- list(); meta <- list()

add <- function(label, algorithm, params, m, extra = list()) {
  # igraph returns $membership unnamed (vertex order); MCL returns it named.
  if (is.null(names(m))) { stopifnot(length(m) == length(nm)); names(m) <- nm }
  m <- as.integer(factor(m[nm]))          # canonical relabelling, node order fixed
  stopifnot(!anyNA(m))
  names(m) <- nm
  memb_list[[label]] <<- m
  meta[[label]] <<- data.table(run = label, algorithm = algorithm, params = params,
                               k = length(unique(m)),
                               modularity = modularity(g, m),
                               largest = max(table(m)),
                               singletons = sum(table(m) == 1),
                               note = if (length(extra)) paste(names(extra), unlist(extra),
                                                               sep = "=", collapse = ";") else "")
}

## ---- MCL (van Dongen), own implementation -------------------------------
for (infl in c(1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0)) {
  r <- mcl_vandongen(A, inflation = infl, prune = 1e-6)
  add(sprintf("MCL_I%.1f", infl), "MCL", sprintf("inflation=%.1f;prune=1e-6", infl),
      r$membership, list(iter = r$iterations, converged = r$converged))
  cat(sprintf("MCL I=%.1f  k=%d  iters=%d converged=%s\n", infl,
              length(unique(r$membership)), r$iterations, r$converged))
}
## pruning sensitivity at the clusterMaker2 default inflation
for (p in c(1e-4, 1e-8, 1e-15)) {
  r <- mcl_vandongen(A, inflation = 2.5, prune = p)
  add(sprintf("MCL_I2.5_p%g", p), "MCL", sprintf("inflation=2.5;prune=%g", p),
      r$membership, list(iter = r$iterations))
}

## ---- Louvain / multilevel, resolution sweep ------------------------------
for (res in c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0)) {
  m <- cluster_louvain(g, resolution = res)$membership
  add(sprintf("Louvain_r%.2f", res), "Louvain", sprintf("resolution=%.2f", res), m)
}
## Louvain seed stability at the default resolution
lou_seeds <- lapply(1:50, function(s) { set.seed(s); cluster_louvain(g)$membership })

## ---- walktrap, step sweep ------------------------------------------------
for (st in c(2, 3, 4, 5, 6, 8)) {
  wt <- cluster_walktrap(g, steps = st)
  add(sprintf("Walktrap_s%d", st), "Walktrap", sprintf("steps=%d", st), membership(wt))
}

## ---- infomap -------------------------------------------------------------
for (tr in c(1, 10, 100)) {
  m <- cluster_infomap(g, nb.trials = tr)$membership
  add(sprintf("Infomap_t%d", tr), "Infomap", sprintf("nb.trials=%d", tr), m)
}
inf_seeds <- lapply(1:20, function(s) { set.seed(s); cluster_infomap(g, nb.trials = 10)$membership })

## ---- label propagation ---------------------------------------------------
lp_seeds <- lapply(1:100, function(s) { set.seed(s); cluster_label_prop(g)$membership })
lp_k <- sapply(lp_seeds, function(m) length(unique(m)))
lp_mod <- sapply(lp_seeds, function(m) modularity(g, m))
add("LabelProp_best", "LabelProp", "100 restarts, best modularity",
    lp_seeds[[which.max(lp_mod)]], list(k_range = paste(range(lp_k), collapse = "-")))
add("LabelProp_median", "LabelProp", "100 restarts, median-modularity run",
    lp_seeds[[order(lp_mod)[ceiling(length(lp_mod)/2)]]])

## ---- fast greedy ---------------------------------------------------------
fg <- cluster_fast_greedy(g)
add("FastGreedy", "FastGreedy", "modularity-optimal cut", membership(fg))
for (k in c(3, 5, 8, 10, 20)) add(sprintf("FastGreedy_k%d", k), "FastGreedy",
                                  sprintf("cut_at k=%d", k), cut_at(fg, no = k))

## ---- extras present in clusterMaker2 -------------------------------------
add("LeadingEigen", "LeadingEigen", "default", membership(cluster_leading_eigen(g)))
for (res in c(0.5, 1.0, 2.0))
  add(sprintf("Leiden_mod_r%.1f", res), "Leiden",
      sprintf("objective=modularity;resolution=%.1f", res),
      cluster_leiden(g, objective_function = "modularity", resolution = res,
                     n_iterations = 10)$membership)
add("EdgeBetweenness", "GirvanNewman", "modularity-optimal cut",
    membership(cluster_edge_betweenness(g)))

## ---- save ----------------------------------------------------------------
M <- as.data.table(do.call(cbind, memb_list))
M[, name := nm]
setcolorder(M, "name")
fwrite(M, "results/v7/cm2_igraph_partitions.csv")
fwrite(rbindlist(meta), "results/v7/cm2_igraph_runs.csv")

saveRDS(list(louvain = lou_seeds, infomap = inf_seeds, labelprop = lp_seeds),
        "results/v7/cm2_seed_replicates.rds")
cat("\nlabel propagation over 100 restarts: k in [", paste(range(lp_k), collapse = ","),
    "], modularity in [", sprintf("%.4f,%.4f", min(lp_mod), max(lp_mod)), "]\n")
cat("wrote results/v7/cm2_igraph_partitions.csv and cm2_igraph_runs.csv\n")
