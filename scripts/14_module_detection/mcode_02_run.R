## ---------------------------------------------------------------------------
## MCODE (Bader & Hogue 2003) as an independent convergent-validity test
## of the higher-order reactive-stroma module.
## Implementation: scripts/14_module_detection/mcode.R, a faithful port of the Cytoscape MCODE
## app's Java source (validated against that app's own unit tests).
## Outputs: results/v7/mcode_*
## ---------------------------------------------------------------------------
source("scripts/14_module_detection/mcode.R")
suppressPackageStartupMessages({library(igraph)})
set.seed(20260912)
OUT <- "results/v7"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
p <- function(x) file.path(OUT, x)

## ---- 1. network -----------------------------------------------------------
nodes <- read.delim("data/canonical_nodes.tsv", stringsAsFactors = FALSE)
edges <- read.delim("data/canonical_edges.tsv", stringsAsFactors = FALSE)
stopifnot(all(edges$source %in% nodes$name), all(edges$target %in% nodes$name))

## undirected projection: drop direction, collapse reciprocal/multi edges, no loops
gd <- graph_from_data_frame(edges[, c("source","target")], directed = TRUE,
                            vertices = nodes[, c("name","type")])
g  <- simplify(as_undirected(gd, mode = "collapse"),
               remove.multiple = TRUE, remove.loops = TRUE)
V(g)$type <- nodes$type[match(V(g)$name, nodes$name)]
cat(sprintf("undirected projection: %d nodes, %d edges, density %.4f, components %d\n",
            vcount(g), ecount(g), edge_density(g), components(g)$no))
nm <- V(g)$name; ty <- V(g)$type

net_summary <- data.frame(
  n_nodes = vcount(g), n_edges_directed = nrow(edges), n_edges_undirected = ecount(g),
  density = edge_density(g), n_components = components(g)$no,
  largest_component = max(components(g)$csize),
  mean_degree = mean(degree(g)), median_degree = median(degree(g)),
  max_degree = max(degree(g)), n_degree_lt2 = sum(degree(g) < 2),
  global_transitivity = transitivity(g, type = "global"))
write.csv(net_summary, p("mcode_01_network_summary.csv"), row.names = FALSE)

## ---- 2. reference sets ----------------------------------------------------
mm  <- read.csv("results/ffl_module_membership.csv", stringsAsFactors = FALSE)
ms_exemplar <- mm$node[mm$in_MS_exemplar_module]          # the paper's module
higher_only <- mm$node[mm$higher_order_only]
pri30 <- read.csv("results/v5/tables/Table4_prioritised_30.csv", stringsAsFactors = FALSE)$name
sigs <- readRDS("cache/v5/farmer_signature_sets.rds")
farmer50   <- intersect(sigs$FARMER_STROMAL, nm)
caf50      <- intersect(sigs$CAF_scRNA_50, nm)
matrisome  <- intersect(sigs$NABA_CORE_MATRISOME, nm)
stromal_union <- intersect(unique(unlist(sigs[c("FARMER_STROMAL","WEST_DTF_FIBROMATOSIS",
   "WINSLOW_STROMAL_SIG1","FINAK_SDPP","FARMER2005_STROMAL_CLUSTER4","CAF_scRNA_50",
   "NABA_CORE_MATRISOME","ESTIMATE_STROMAL")])), nm)
mir29 <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c")
collagen <- c("COL1A1","COL3A1")
## collagen axis: the two collagens, their TFs, miR-29
collagen_axis <- intersect(c(collagen, "ETS1","NFKB1","RELA","SP1", mir29,
                             "hsa-miR-101","EZH2"), nm)
refsets <- list(MS_exemplar_module = ms_exemplar, collagen_axis = collagen_axis,
                FARMER_STROMAL = farmer50, CAF_scRNA_50 = caf50,
                NABA_CORE_MATRISOME = matrisome, stromal_union = stromal_union,
                prioritised_30 = pri30, higher_order_only = higher_only)
write.csv(data.frame(set = names(refsets),
                     n_in_network = sapply(refsets, length),
                     members = sapply(refsets, paste, collapse = ";")),
          p("mcode_02_reference_sets.csv"), row.names = FALSE)

## ---- 3. helpers -----------------------------------------------------------
hyper <- function(hit, set_size, draw, N)          # P(X >= hit)
  if (hit == 0) 1 else phyper(hit - 1, set_size, N - set_size, draw, lower.tail = FALSE)

describe <- function(cl, g, tag, params) {
  if (length(cl) == 0)
    return(data.frame(run = tag, rank = NA, seed = NA, size = 0)[0, ])
  do.call(rbind, lapply(cl, function(c) {
    v <- nm[c$nodes]
    data.frame(
      run = tag, params = params, rank = c$rank, seed = nm[c$seed],
      size = c$size, edges = c$edges,
      density = if (c$size > 1) c$edges / (c$size * (c$size - 1) / 2) else 0,
      score = c$score,
      n_TF = sum(ty[c$nodes] == "TF"), n_miRNA = sum(ty[c$nodes] == "miRNA"),
      n_Gene = sum(ty[c$nodes] == "Gene"),
      frac_of_network = c$size / vcount(g),
      COL1A1 = "COL1A1" %in% v, COL3A1 = "COL3A1" %in% v,
      n_miR29 = sum(mir29 %in% v),
      collagens_and_mir29_together = all(collagen %in% v) && any(mir29 %in% v),
      n_MS_exemplar = sum(ms_exemplar %in% v),
      p_MS_exemplar = hyper(sum(ms_exemplar %in% v), length(ms_exemplar), c$size, vcount(g)),
      n_prioritised30 = sum(pri30 %in% v),
      p_prioritised30 = hyper(sum(pri30 %in% v), length(pri30), c$size, vcount(g)),
      n_farmer50 = sum(farmer50 %in% v),
      p_farmer50 = hyper(sum(farmer50 %in% v), length(farmer50), c$size, vcount(g)),
      n_stromal_union = sum(stromal_union %in% v),
      p_stromal_union = hyper(sum(stromal_union %in% v), length(stromal_union), c$size, vcount(g)),
      members = paste(v, collapse = ";"), stringsAsFactors = FALSE)
  }))
}

## ---- 4. MCODE at Cytoscape defaults --------------------------------------
info <- mcode_score_graph(g, degree_cutoff = 2L, include_loops = FALSE)
ns <- data.frame(node = nm, type = ty, degree = degree(g),
                 mcode_node_score = info$score, core_level = info$core_level,
                 core_density = info$core_density, nbhd_density = info$density,
                 in_prioritised30 = nm %in% pri30, in_MS_exemplar = nm %in% ms_exemplar)
ns <- ns[order(-ns$mcode_node_score), ]
write.csv(ns, p("mcode_03_node_scores.csv"), row.names = FALSE)
cat("node score: max", max(info$score), " n>0:", sum(info$score > 0),
    " n distinct:", length(unique(info$score)), "\n")

cl_def <- mcode_find_clusters(g, info, degree_cutoff = 2L, k_core = 2L,
                              max_depth = 100L, node_score_cutoff = 0.2,
                              haircut = TRUE, fluff = FALSE)
res_def <- describe(cl_def, g, "default",
                    "deg=2;nsc=0.2;kcore=2;depth=100;haircut=T;fluff=F")
write.csv(res_def, p("mcode_04_clusters_default.csv"), row.names = FALSE)
cat("\n== DEFAULTS ==\n"); print(res_def[, c("rank","seed","size","score","n_TF","n_miRNA",
  "n_Gene","COL1A1","COL3A1","n_miR29","n_MS_exemplar","n_prioritised30","n_farmer50")])

## ---- 5. sensitivity: node score cutoff 0.1 - 0.5 -------------------------
grid <- seq(0.10, 0.50, by = 0.05)
sens <- do.call(rbind, lapply(grid, function(q) {
  cl <- mcode_find_clusters(g, info, node_score_cutoff = q)
  describe(cl, g, sprintf("nsc=%.2f", q), sprintf("deg=2;nsc=%.2f;kcore=2;depth=100;haircut=T;fluff=F", q))
}))
write.csv(sens, p("mcode_05_sensitivity_nodescorecutoff.csv"), row.names = FALSE)
cat("\n== NODE SCORE CUTOFF SWEEP (top cluster per setting) ==\n")
print(do.call(rbind, lapply(split(sens, sens$run), function(d) {
  d1 <- d[d$rank == 1, ]
  data.frame(run = d1$run, n_clusters = nrow(d), top_size = d1$size,
             top_score = round(d1$score, 3), top_frac_net = round(d1$frac_of_network, 3),
             COL1A1 = d1$COL1A1, COL3A1 = d1$COL3A1, miR29 = d1$n_miR29,
             any_cluster_has_axis = any(d$collagens_and_mir29_together))})))

## ---- 6. sensitivity: other parameters ------------------------------------
oth <- list()
for (dc in c(2L, 3L, 4L, 5L)) {
  i2 <- mcode_score_graph(g, degree_cutoff = dc)
  oth[[sprintf("degcut=%d", dc)]] <- describe(
    mcode_find_clusters(g, i2, degree_cutoff = dc), g, sprintf("degcut=%d", dc),
    sprintf("deg=%d;nsc=0.2;kcore=2;depth=100", dc))
}
for (kc in 2:6) oth[[sprintf("kcore=%d", kc)]] <- describe(
  mcode_find_clusters(g, info, k_core = kc), g, sprintf("kcore=%d", kc),
  sprintf("deg=2;nsc=0.2;kcore=%d;depth=100", kc))
for (md in c(1L, 2L, 3L, 5L, 10L, 100L)) oth[[sprintf("depth=%d", md)]] <- describe(
  mcode_find_clusters(g, info, max_depth = md), g, sprintf("depth=%d", md),
  sprintf("deg=2;nsc=0.2;kcore=2;depth=%d", md))
oth[["haircut=F"]] <- describe(mcode_find_clusters(g, info, haircut = FALSE), g,
                               "haircut=F", "deg=2;nsc=0.2;kcore=2;depth=100;haircut=F")
oth[["fluff=T"]]   <- describe(mcode_find_clusters(g, info, fluff = TRUE), g,
                               "fluff=T", "deg=2;nsc=0.2;kcore=2;depth=100;haircut=T;fluff=T;fdt=0.1")
other <- do.call(rbind, oth)
write.csv(other, p("mcode_06_sensitivity_other_params.csv"), row.names = FALSE)
cat("\n== OTHER PARAMETERS (top cluster per setting) ==\n")
print(do.call(rbind, lapply(split(other, factor(other$run, levels = unique(other$run))),
  function(d) { d1 <- d[d$rank == 1, ]
    data.frame(run = d1$run, n_clusters = nrow(d), top_size = d1$size,
               top_score = round(d1$score,3), top_frac = round(d1$frac_of_network,3),
               COL1A1 = d1$COL1A1, COL3A1 = d1$COL3A1, miR29 = d1$n_miR29,
               axis_any = any(d$collagens_and_mir29_together)) })))

## ---- 7. tie-break sensitivity (seed order) -------------------------------
B <- 100
perm <- do.call(rbind, lapply(seq_len(B), function(b) {
  cl <- mcode_find_clusters(g, info, seed_order = sample.int(vcount(g)))
  d <- describe(cl, g, sprintf("perm%03d", b), "defaults, permuted tie-break")
  d1 <- d[d$rank == 1, ]
  data.frame(b = b, n_clusters = nrow(d), top_size = d1$size, top_score = d1$score,
             top_has_COL1A1 = d1$COL1A1, top_has_COL3A1 = d1$COL3A1,
             top_n_miR29 = d1$n_miR29,
             any_axis = any(d$collagens_and_mir29_together),
             axis_cluster_size = if (any(d$collagens_and_mir29_together))
               d$size[which(d$collagens_and_mir29_together)[1]] else NA_integer_)
}))
write.csv(perm, p("mcode_07_seedorder_permutations.csv"), row.names = FALSE)
cat("\n== TIE-BREAK PERMUTATIONS (n=100) ==\n")
print(summary(perm[, c("n_clusters","top_size","top_score")]))
cat("fraction of permutations where some cluster holds both collagens + a miR-29:",
    mean(perm$any_axis), "\n")

## ---- 8. strong-evidence subnetwork ---------------------------------------
ev <- read.delim("data/edge_evidence_tier.tsv", stringsAsFactors = FALSE)
strong <- ev[ev$tier == "strong", c("source","target")]
keep_key <- paste(pmin(strong$source, strong$target), pmax(strong$source, strong$target))
nonmir <- edges[edges$edge_type != "miRNA_target", c("source","target")]
keep_key <- unique(c(keep_key, paste(pmin(nonmir$source, nonmir$target),
                                     pmax(nonmir$source, nonmir$target))))
el <- as_edgelist(g)
gs <- subgraph.edges(g, which(paste(pmin(el[,1], el[,2]), pmax(el[,1], el[,2])) %in% keep_key),
                     delete.vertices = FALSE)
cat(sprintf("\nstrong-evidence projection: %d nodes, %d edges\n", vcount(gs), ecount(gs)))
info_s <- mcode_score_graph(gs)
nm_all <- nm
res_strong <- local({
  nmx <- V(gs)$name; tyx <- V(gs)$type
  cl <- mcode_find_clusters(gs, info_s)
  nm <<- nmx; ty <<- tyx
  d <- describe(cl, gs, "strong_evidence_only", "deg=2;nsc=0.2;kcore=2;depth=100;haircut=T")
  nm <<- nm_all; ty <<- nodes$type[match(nm_all, nodes$name)]
  d })
write.csv(res_strong, p("mcode_08_clusters_strong_evidence.csv"), row.names = FALSE)
cat("\n== STRONG-EVIDENCE SUBNETWORK ==\n")
print(head(res_strong[, c("rank","seed","size","score","n_TF","n_miRNA","n_Gene",
  "COL1A1","COL3A1","n_miR29","collagens_and_mir29_together","n_prioritised30")], 15))

## ---- 9. independent implementation: ProNet::mcode -------------------------
if (requireNamespace("ProNet", quietly = TRUE)) {
  pn <- ProNet::mcode(g, vwp = 0.2, haircut = TRUE, fluff = FALSE, loops = FALSE)
  pnd <- do.call(rbind, lapply(seq_along(pn$COMPLEX), function(i) {
    v <- nm[pn$COMPLEX[[i]]]
    data.frame(rank = i, size = length(v), score = pn$score[i],
               n_TF = sum(ty[pn$COMPLEX[[i]]] == "TF"),
               n_miRNA = sum(ty[pn$COMPLEX[[i]]] == "miRNA"),
               n_Gene = sum(ty[pn$COMPLEX[[i]]] == "Gene"),
               COL1A1 = "COL1A1" %in% v, COL3A1 = "COL3A1" %in% v,
               n_miR29 = sum(mir29 %in% v),
               collagens_and_mir29_together = all(collagen %in% v) && any(mir29 %in% v),
               n_MS_exemplar = sum(ms_exemplar %in% v),
               n_prioritised30 = sum(pri30 %in% v),
               p_prioritised30 = hyper(sum(pri30 %in% v), length(pri30), length(v), vcount(g)),
               n_farmer50 = sum(farmer50 %in% v),
               members = paste(v, collapse = ";"), stringsAsFactors = FALSE) }))
  write.csv(pnd, p("mcode_09_pronet_crosscheck.csv"), row.names = FALSE)
  cat("\n== ProNet::mcode 1.0.0 (paper-pseudocode variant), vwp=0.2 ==\n")
  print(pnd[, c("rank","size","score","n_TF","n_miRNA","n_Gene","COL1A1","COL3A1",
                "n_miR29","collagens_and_mir29_together","n_prioritised30")])
}

save(g, info, cl_def, res_def, sens, other, perm, refsets,
     file = p("mcode_workspace.RData"))
cat("\nDONE\n")
