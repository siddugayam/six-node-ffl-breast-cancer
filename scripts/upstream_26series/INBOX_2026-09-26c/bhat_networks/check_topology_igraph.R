## Cross-check of topology_mcc.py with igraph (independent implementation): MCC from igraph's maximal cliques,
## betweenness normalised within the connected component, closeness 1/mean distance, local clustering.
suppressPackageStartupMessages({library(igraph); library(data.table)})
HERE <- "/path/to/revision/INBOX_2026-09-26c/bhat_networks"
for (tag in c("nolegacy", "nolegacy_nostring")) {
  sif <- fread(file.path(HERE, "sif", paste0(tag, "_composite_3to6_merged.sif")), header = FALSE, col.names = c("s", "lab", "t"))
  g <- simplify(graph_from_data_frame(sif[, .(s, t)], directed = FALSE))
  mc <- max_cliques(g); mcc <- setNames(numeric(vcount(g)), V(g)$name)
  for (C in mc) { v <- names(C); mcc[v] <- mcc[v] + factorial(length(C) - 1) }
  comp <- components(g)$membership; bet <- setNames(numeric(vcount(g)), V(g)$name)
  for (k in unique(comp)) { sg <- induced_subgraph(g, which(comp == k))
    if (vcount(sg) > 2) bet[V(sg)$name] <- betweenness(sg, directed = FALSE, normalized = TRUE) }
  D <- distances(g); D[is.infinite(D)] <- NA; diag(D) <- NA
  clo <- 1 / rowMeans(D, na.rm = TRUE)
  clu <- transitivity(g, type = "local", isolates = "zero"); names(clu) <- V(g)$name
  py <- fread(file.path(HERE, paste0("topology_", tag, "_composite_3to6_merged.csv")))
  i <- match(py$node, V(g)$name)
  cat(sprintf("%s: %d nodes | max |diff| MCC %.3g, betweenness %.3g, closeness %.3g, clustering %.3g, neighbours %d\n",
              tag, nrow(py), max(abs(py$MCC - mcc[i])), max(abs(py$betweenness - bet[i])),
              max(abs(py$closeness - clo[i])), max(abs(py$clustering - clu[i])),
              max(abs(py$neighbours - degree(g)[i]))))
}
