## MCODE on the six-node composite FFL network, as in Bhat et al. (2024) section 2.8: ClusterMaker MCODE
## defaults with the haircut method, k-core 2 and max depth 100 (degree cutoff 2, node score cutoff 0.2, no
## fluff), run through scripts/14_module_detection/mcode.R, the Cytoscape-faithful port.  Then edge betweenness inside the
## top-ranked cluster.  The network is undirected; reciprocal M<->T and T<->T pairs collapse to one edge.
## Node order, which breaks ties between equal node scores, is first appearance in the SIF, the order in
## which Cytoscape creates nodes when it imports the file.  100 random node orders test that tie-break.
REV  <- "/path/to/revision"
HERE <- file.path(REV, "analyses/six_node_pattern_networks")
source(file.path(REV, "scripts/14_module_detection/mcode.R"))
suppressPackageStartupMessages(library(igraph))
set.seed(20260927)

for (tag in c("nolegacy", "nolegacy_nostring")) {
  sif  <- read.delim(file.path(HERE, "sif", paste0(tag, "_6node_Composite_FFL.sif")), header = FALSE,
                     col.names = c("s", "lab", "t"), stringsAsFactors = FALSE)
  attr <- read.delim(file.path(HERE, "sif", paste0(tag, "_node_attributes.tsv")), stringsAsFactors = FALSE)
  ord  <- unique(as.vector(t(as.matrix(sif[, c("s", "t")]))))
  g <- graph_from_data_frame(sif[, c("s", "t")], directed = FALSE, vertices = data.frame(name = ord))
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  V(g)$type <- attr$type[match(V(g)$name, attr$node)]
  cat(sprintf("\n== %s six-node composite network: %d nodes, %d SIF edges, %d undirected edges\n",
              tag, vcount(g), nrow(sif), ecount(g)))

  info <- mcode_score_graph(g)
  cl <- mcode_find_clusters(g, info, degree_cutoff = 2L, k_core = 2L, max_depth = 100L,
                            node_score_cutoff = 0.2, haircut = TRUE, fluff = FALSE)
  nm <- V(g)$name; ty <- V(g)$type
  memb <- function(v, t) paste(sort(nm[v][ty[v] == t]), collapse = ";")
  tab <- do.call(rbind, lapply(cl, function(x) data.frame(
    rank = x$rank, score = round(x$score, 4), nodes = x$size, edges = x$edges, seed = nm[x$seed],
    n_TF = sum(ty[x$nodes] == "TF"), n_miRNA = sum(ty[x$nodes] == "miRNA"), n_Gene = sum(ty[x$nodes] == "Gene"),
    TF = memb(x$nodes, "TF"), miRNA = memb(x$nodes, "miRNA"), Gene = memb(x$nodes, "Gene"))))
  write.csv(tab, file.path(HERE, paste0("mcode_", tag, "_six_node_composite_clusters.csv")), row.names = FALSE)
  write.csv(data.frame(node = nm, type = ty, node_score = info$score, core_level = info$core_level,
                       core_density = info$core_density, closed_neighbourhood_size = info$num_nb)[order(-info$score, nm), ],
            file.path(HERE, paste0("mcode_", tag, "_six_node_composite_node_scores.csv")), row.names = FALSE)
  cat(sprintf("   %d clusters; top 5:\n", length(cl)))
  print(head(tab[, c("rank", "score", "nodes", "edges", "seed", "n_TF", "n_miRNA", "n_Gene")], 5), row.names = FALSE)
  if (length(cl) == 0) next
  top <- cl[[1]]
  cat("   top cluster TF:", tab$TF[1], "\n   miRNA:", tab$miRNA[1], "\n   Gene:", tab$Gene[1], "\n")

  ## top-cluster subnetwork: SIF edges with both ends in the cluster; edge betweenness on its undirected form
  keep <- nm[top$nodes]
  sub_sif <- sif[sif$s %in% keep & sif$t %in% keep, ]
  write.table(sub_sif, file.path(HERE, "sif", paste0(tag, "_mcode_top_cluster.sif")), sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  sg <- induced_subgraph(g, top$nodes)
  eb <- edge_betweenness(sg, directed = FALSE)
  el <- as_edgelist(sg)
  lab <- apply(el, 1, function(p) paste(sort(unique(sub_sif$lab[(sub_sif$s == p[1] & sub_sif$t == p[2]) |
                                                                (sub_sif$s == p[2] & sub_sif$t == p[1])])), collapse = "+"))
  ebt <- data.frame(node_a = el[, 1], node_b = el[, 2], type_a = ty[match(el[, 1], nm)], type_b = ty[match(el[, 2], nm)],
                    sif_labels = lab, edge_betweenness_pairs_once = eb, edge_betweenness_ordered_pairs = 2 * eb)
  ebt <- ebt[order(-ebt$edge_betweenness_pairs_once, ebt$node_a, ebt$node_b), ]
  write.csv(ebt, file.path(HERE, paste0("mcode_", tag, "_top_cluster_edge_betweenness.csv")), row.names = FALSE)
  cat("   top-cluster edge betweenness, highest 8 (pairs counted once):\n")
  print(head(ebt[, c("node_a", "node_b", "sif_labels", "edge_betweenness_pairs_once")], 8), row.names = FALSE)

  ## tie-break sensitivity: same parameters, 100 random node orders
  same <- 0L; sizes <- integer(0)
  for (i in 1:100) {
    ci <- mcode_find_clusters(g, info, seed_order = sample(vcount(g)))
    same <- same + as.integer(setequal(ci[[1]]$nodes, top$nodes)); sizes <- c(sizes, ci[[1]]$size)
  }
  cat(sprintf("   top cluster identical in %d of 100 random node orders (top-cluster size range %d-%d)\n",
              same, min(sizes), max(sizes)))
  writeLines(sprintf("%s\t%d\t%d\t%d", tag, same, min(sizes), max(sizes)),
             file.path(HERE, paste0("mcode_", tag, "_tiebreak_100_orders.tsv")))
}
