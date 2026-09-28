#!/usr/bin/env Rscript
# v7 / 00 : build the undirected simple graph used by every community-detection
#           algorithm in the clusterMaker2 battery, and export it for Cytoscape.
#
# Design decision (stated in the report): every algorithm in the battery
# (MCL, Louvain, walktrap, infomap, label propagation, fast-greedy) is defined
# on an UNDIRECTED graph. The canonical network is directed and signed, so we
# collapse it: direction dropped, reciprocal pairs merged, edges unweighted
# (all |sign| = 1). No self-loops, no multi-edges. This is exactly what one
# would get by loading canonical_edges.tsv into Cytoscape and clustering it.
suppressMessages({library(igraph); library(data.table)})
setwd("/path/to/revision")

e <- fread("data/canonical_edges.tsv")
n <- fread("data/canonical_nodes.tsv")

gd <- graph_from_data_frame(e[, .(source, target)], directed = TRUE,
                            vertices = n[, .(name, type)])
gu <- simplify(as_undirected(gd, mode = "collapse"),
               remove.multiple = TRUE, remove.loops = TRUE)
E(gu)$weight <- 1

stopifnot(components(gu)$no == 1)
cat(sprintf("undirected simple graph: %d nodes, %d edges, 1 component, density %.4f\n",
            vcount(gu), ecount(gu), edge_density(gu)))

saveRDS(gu, "results/v7/cm2_graph_undirected.rds")

# SIF for Cytoscape / clusterMaker2
el <- as.data.table(as_edgelist(gu))
setnames(el, c("s", "t"))
el[, i := "interacts"]
fwrite(el[, .(s, i, t)], "results/v7/cm2_network.sif", sep = "\t", col.names = FALSE)
# node attribute table (type) for Cytoscape
fwrite(n[, .(name, type)], "results/v7/cm2_nodes.tsv", sep = "\t")
cat("wrote results/v7/cm2_network.sif\n")
