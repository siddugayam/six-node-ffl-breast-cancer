#!/usr/bin/env Rscript
## ============================================================================
## 09_ffl_hubs.R
## Define "high-centrality hubs in the FFLs" on the REPAIRED canonical network
##
## Two independent hub definitions are produced:
##   (A) CENTRALITY hubs   -- top decile of total degree AND/OR betweenness
##   (B) FFL-PARTICIPATION hubs -- exhaustive enumeration of every 3-node
##       feed-forward loop (a->b, b->c, a->c) in the directed canonical
##       network; a node is a hub if it is in the top decile of the number of
##       distinct FFLs it participates in.
## Both are written out so the overlap test can be reported both ways.
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

ROOT <- "/path/to/revision"
LOG  <- file.path(ROOT, "logs", "exir_agent.log")
say <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | 09 | ", paste0(..., collapse = ""))
  cat(msg, "\n", file = LOG, append = TRUE); cat(msg, "\n"); flush.console()
}
say("=== 09_ffl_hubs.R START ===")

edges <- fread(file.path(ROOT, "data", "canonical_edges.tsv"))
nodes <- fread(file.path(ROOT, "data", "canonical_nodes.tsv"))
say("edges ", nrow(edges), "  nodes ", nrow(nodes))

## ---- directed graph, de-duplicated -----------------------------------------
el <- unique(edges[, .(source, target)])
say("unique directed (source,target) pairs: ", nrow(el))
g <- graph_from_data_frame(as.data.frame(el),
                           directed = TRUE,
                           vertices = as.data.frame(nodes[, .(name, type)]))
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
say("igraph: V=", vcount(g), " E=", ecount(g))

## ---- (A) centrality ---------------------------------------------------------
deg_all <- degree(g, mode = "all")
deg_in  <- degree(g, mode = "in")
deg_out <- degree(g, mode = "out")
btw     <- betweenness(g, directed = TRUE, normalized = FALSE)

## ---- (B) exhaustive 3-node FFL enumeration ---------------------------------
succ <- lapply(V(g), function(v) as.integer(neighbors(g, v, mode = "out")))
names(succ) <- V(g)$name
succ_set <- lapply(succ, function(x) x)

n_ffl <- 0L
part  <- integer(vcount(g)); names(part) <- V(g)$name
role_src <- integer(vcount(g)); names(role_src) <- V(g)$name
role_mid <- integer(vcount(g)); names(role_mid) <- V(g)$name
role_snk <- integer(vcount(g)); names(role_snk) <- V(g)$name

for (a in seq_len(vcount(g))) {
  sa <- succ_set[[a]]
  if (length(sa) < 2L) next
  sa_lookup <- logical(vcount(g)); sa_lookup[sa] <- TRUE
  for (b in sa) {
    sb <- succ_set[[b]]
    if (!length(sb)) next
    cs <- sb[sa_lookup[sb]]
    cs <- cs[cs != a & cs != b]
    if (!length(cs)) next
    n_ffl <- n_ffl + length(cs)
    role_src[a] <- role_src[a] + length(cs)
    role_mid[b] <- role_mid[b] + length(cs)
    role_snk[cs] <- role_snk[cs] + 1L
  }
}
part <- role_src + role_mid + role_snk
say("exhaustive 3-node FFLs (a->b, b->c, a->c) in the canonical directed network: ", n_ffl)
say("nodes participating in >=1 FFL: ", sum(part > 0))

## ---- hub tables -------------------------------------------------------------
hub <- data.table(node       = V(g)$name,
                  node_type  = V(g)$type,
                  degree_all = as.integer(deg_all),
                  degree_in  = as.integer(deg_in),
                  degree_out = as.integer(deg_out),
                  betweenness = as.numeric(btw),
                  ffl_participation = as.integer(part),
                  ffl_as_source = as.integer(role_src),
                  ffl_as_mid    = as.integer(role_mid),
                  ffl_as_sink   = as.integer(role_snk))

q90_deg <- quantile(hub$degree_all,        0.90)
q90_btw <- quantile(hub$betweenness,       0.90)
q90_ffl <- quantile(hub$ffl_participation, 0.90)
say("90th percentile thresholds: degree ", q90_deg, "  betweenness ",
    round(q90_btw, 1), "  FFL-participation ", q90_ffl)

hub[, hub_degree      := degree_all        >= q90_deg]
hub[, hub_betweenness := betweenness       >= q90_btw]
hub[, hub_centrality  := hub_degree | hub_betweenness]
hub[, hub_ffl         := ffl_participation >= q90_ffl]
setorder(hub, -ffl_participation, -degree_all)

fwrite(hub, file.path(ROOT, "results", "ffl_network_hubs.csv"))
say("hub_degree=", sum(hub$hub_degree), " hub_betweenness=", sum(hub$hub_betweenness),
    " hub_centrality(union)=", sum(hub$hub_centrality), " hub_ffl=", sum(hub$hub_ffl))
say("wrote results/ffl_network_hubs.csv  rows=", nrow(hub))
say("=== 09 DONE ===")
