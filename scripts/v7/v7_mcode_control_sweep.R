## v7 -- the negative controls, swept across the full node-score-cutoff range,
## so that the negative result cannot be a single-parameter accident.
source("scripts/v7/mcode.R")
suppressPackageStartupMessages(library(igraph))
set.seed(20260912)
OUT <- "results/v7"; p <- function(x) file.path(OUT, x)
nodes <- read.delim("data/canonical_nodes.tsv", stringsAsFactors = FALSE)
edges <- read.delim("data/canonical_edges.tsv", stringsAsFactors = FALSE)
L  <- read.delim("data/layer_miRNA_miRNA.tsv", stringsAsFactors = FALSE)
ev <- read.delim("data/edge_evidence_tier.tsv", stringsAsFactors = FALSE)
mm <- read.csv("results/ffl_module_membership.csv", stringsAsFactors = FALSE)
ms_ex <- mm$node[mm$in_MS_exemplar_module]
mir29 <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"); collagen <- c("COL1A1","COL3A1")
k <- function(a,b) paste(pmin(a,b), pmax(a,b))
hyper <- function(hit,K,n,N) if (hit==0) 1 else phyper(hit-1,K,N-K,n,lower.tail=FALSE)
build <- function(el) { gd <- graph_from_data_frame(el[,c("source","target")], directed=TRUE,
    vertices=nodes[,c("name","type")])
  g <- simplify(as_undirected(gd, mode="collapse"), remove.multiple=TRUE, remove.loops=TRUE)
  V(g)$type <- nodes$type[match(V(g)$name, nodes$name)]; g }

is_mm <- edges$edge_type == "miRNA_miRNA"
genomic <- k(L$miRNA_1, L$miRNA_2)[L$threshold_10kb]
supported_mm <- is_mm & (k(edges$source, edges$target) %in% genomic)
strongk <- k(ev$source, ev$target)[ev$tier == "strong"]
keepstrong <- (k(edges$source,edges$target) %in% strongk) | (edges$edge_type != "miRNA_target")
gen_el <- L[L$threshold_10kb & L$miRNA_1 %in% nodes$name & L$miRNA_2 %in% nodes$name, c("miRNA_1","miRNA_2")]
names(gen_el) <- c("source","target")

nets <- list(
  as_published = edges,
  drop_all_miRNA_miRNA = edges[!is_mm, ],
  keep_only_coclustered = edges[!is_mm | supported_mm, ],
  consistent_genomic_layer = rbind(edges[!is_mm, c("source","target")], gen_el),
  strong_miRNA_target_only = edges[keepstrong, ],
  strong_and_no_miRNA_miRNA = edges[keepstrong & !is_mm, ])

grid <- seq(0.10, 0.50, by = 0.05)
res <- do.call(rbind, lapply(names(nets), function(v) {
  g <- build(nets[[v]]); nm <- V(g)$name; N <- vcount(g)
  info <- mcode_score_graph(g)
  do.call(rbind, lapply(grid, function(q) {
    cl <- mcode_find_clusters(g, info, node_score_cutoff = q)
    if (!length(cl)) return(NULL)
    mem <- lapply(cl, function(c) nm[c$nodes])
    axis <- which(sapply(mem, function(v) all(collagen %in% v) && any(mir29 %in% v)))
    ovl <- sapply(mem, function(v) sum(ms_ex %in% v))
    bi <- which.max(ovl)
    data.frame(network = v, nsc = q, n_clusters = length(cl),
      top_size = length(mem[[1]]), top_score = round(cl[[1]]$score, 3),
      axis_found = length(axis) > 0,
      axis_rank = if (length(axis)) axis[1] else NA,
      axis_size = if (length(axis)) length(mem[[axis[1]]]) else NA,
      best_exemplar_overlap = ovl[bi], best_exemplar_cluster_size = length(mem[[bi]]),
      best_exemplar_p = signif(hyper(ovl[bi], length(ms_ex), length(mem[[bi]]), N), 3),
      best_exemplar_members = paste(intersect(ms_ex, mem[[bi]]), collapse=";"),
      stringsAsFactors = FALSE) })) }))
write.csv(res, p("mcode_18_control_cutoff_sweep.csv"), row.names = FALSE)
cat("== CONTROL NETWORKS x NODE SCORE CUTOFF 0.10-0.50 ==\n")
print(res[, c("network","nsc","n_clusters","top_size","top_score","axis_found",
              "axis_rank","best_exemplar_overlap","best_exemplar_p")], row.names = FALSE)
cat("\n-- axis recovered in how many of the 9 cutoffs? --\n")
print(tapply(res$axis_found, res$network, sum))
cat("\n-- best exemplar members recovered per network (over all cutoffs) --\n")
for (v in names(nets)) { d <- res[res$network==v,]
  cat(sprintf("%-28s max overlap %d/10 : %s\n", v, max(d$best_exemplar_overlap),
      d$best_exemplar_members[which.max(d$best_exemplar_overlap)])) }
cat("\nDONE\n")
