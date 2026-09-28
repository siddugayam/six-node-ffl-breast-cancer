#!/usr/bin/env Rscript
# (A cont.) Independent cross-check, in R/igraph, of the topology numbers that the Python
# analysis (01_architecture.py, 14_regulator_subnetwork.py) reports with networkx.
# Different library, different code base: if the two disagree, section A is wrong.
suppressPackageStartupMessages(library(igraph))

REV <- "/path/to/revision"
RES <- file.path(REV, "results", "v3")

E <- read.delim(file.path(REV, "data", "canonical_edges.tsv"), stringsAsFactors = FALSE)
N <- read.delim(file.path(REV, "data", "canonical_nodes.tsv"), stringsAsFactors = FALSE)
stopifnot(!any(duplicated(N$name)))
g <- graph_from_data_frame(E[, c("source", "target")], directed = TRUE, vertices = N)
cat(sprintf("nodes %d  edges %d  simple=%s\n", vcount(g), ecount(g), is_simple(g)))

scc <- components(g, mode = "strong")
wcc <- components(g, mode = "weak")
cat(sprintf("SCC: %d components, largest %d, singletons %d\n",
            scc$no, max(scc$csize), sum(scc$csize == 1)))
cat(sprintf("WCC: %d components, largest %d\n", wcc$no, max(wcc$csize)))
cat(sprintf("reciprocity %.6f  transitivity(undirected) %.6f\n",
            reciprocity(g), transitivity(as_undirected(g, mode = "collapse"), type = "global")))

# bow-tie relative to the largest SCC
big <- which(scc$membership == which.max(scc$csize))
rep_v <- big[1]
reach_fwd <- subcomponent(g, rep_v, mode = "out")
reach_bwd <- subcomponent(g, rep_v, mode = "in")
IN <- setdiff(reach_bwd, big)
OUT <- setdiff(reach_fwd, big)
other <- setdiff(V(g), union(union(IN, OUT), big))
cat(sprintf("bow-tie  IN %d  CORE %d  OUT %d  other %d\n",
            length(IN), length(big), length(OUT), length(other)))
tab <- table(V(g)$type[big])
cat("core composition:", paste(names(tab), tab, collapse = "  "), "\n")

# degree summary
deg <- degree(g, mode = "all"); din <- degree(g, mode = "in"); dout <- degree(g, mode = "out")
cat(sprintf("degree: mean %.4f max %d ; in max %d ; out max %d ; out==0: %d ; in==0: %d\n",
            mean(deg), max(deg), max(din), max(dout), sum(dout == 0), sum(din == 0)))
cat(sprintf("out-degree 0 by type: %s\n",
            paste(names(table(V(g)$type[dout == 0])), table(V(g)$type[dout == 0]),
                  collapse = "  ")))

# regulator-only subgraph
reg <- V(g)[V(g)$type %in% c("TF", "miRNA")]
gr <- induced_subgraph(g, reg)
sccr <- components(gr, mode = "strong")
cat(sprintf("\nregulators only: %d nodes %d edges; SCC %d, largest %d (%.1f%%), reciprocity %.4f\n",
            vcount(gr), ecount(gr), sccr$no, max(sccr$csize),
            100 * max(sccr$csize) / vcount(gr), reciprocity(gr)))
# fraction of edges lying on no cycle == flow hierarchy (Luo & Magee 2011)
sm <- sccr$membership
fh <- mean(sm[ends(gr, E(gr))[, 1]] != sm[ends(gr, E(gr))[, 2]])
sm2 <- scc$membership
fh2 <- mean(sm2[ends(g, E(g))[, 1]] != sm2[ends(g, E(g))[, 2]])
cat(sprintf("flow hierarchy: full %.6f   regulators only %.6f\n", fh2, fh))

out <- data.frame(
  metric = c("nodes", "edges", "n_SCC", "largest_SCC", "n_WCC", "reciprocity",
             "transitivity_undirected", "bowtie_IN", "bowtie_CORE", "bowtie_OUT",
             "n_out_degree_zero", "n_in_degree_zero", "flow_hierarchy_full",
             "regulators_N", "regulators_M", "regulators_largest_SCC",
             "regulators_reciprocity", "flow_hierarchy_regulators"),
  igraph_value = c(vcount(g), ecount(g), scc$no, max(scc$csize), wcc$no, reciprocity(g),
                   transitivity(as_undirected(g, mode = "collapse"), type = "global"),
                   length(IN), length(big), length(OUT), sum(dout == 0), sum(din == 0), fh2,
                   vcount(gr), ecount(gr), max(sccr$csize), reciprocity(gr), fh),
  stringsAsFactors = FALSE)
write.csv(out, file.path(RES, "systems_igraph_crosscheck.csv"), row.names = FALSE)
print(out)
cat("\nwrote systems_igraph_crosscheck.csv\n")
