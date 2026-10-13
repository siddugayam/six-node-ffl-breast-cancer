## ---------------------------------------------------------------------------
## mcode.R -- a faithful R/igraph port of the Cytoscape MCODE app
##            (Bader & Hogue, BMC Bioinformatics 2003;4:2)
##
## Ported line-by-line from the reference Java implementation
##   BaderLab/MCODE, src/main/java/ca/utoronto/tdccbr/mcode/internal/model/
##     MCODEAlgorithm.java   (Gary Bader; LGPL)
##   commit fetched 2026-09-12 from https://github.com/BaderLab/MCODE
##
## Why a port and not a package call: there is no MCODE implementation on CRAN
## or Bioconductor.  The only R implementation, ProNet::mcode (CRAN archive),
## transcribes the *paper pseudocode*, in which the node-score threshold is
## re-referenced to the current node at every recursion step, and it omits the
## degree cutoff, the k-core filter and the max-depth parameter.  The Cytoscape
## app -- the de-facto standard the question asks about -- keeps the threshold
## fixed at the seed's score and implements all four parameters.  Both are run;
## this file is the Cytoscape semantics.
##
## Defaults, from MCODEParameters.setDefaultParams():
##   includeLoops = FALSE, degreeCutoff = 2, kCore = 2, maxDepthFromStart = 100,
##   nodeScoreCutoff = 0.2, fluff = FALSE, haircut = TRUE,
##   fluffNodeDensityCutoff = 0.1
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(igraph))

## calcDensity(): MCODEAlgorithm.java L837-857.  Multi-edges merged, loops
## optionally counted; possibleEdgeNum uses Java integer division.
mcode_density <- function(g, include_loops = FALSE) {
  n <- vcount(g)
  if (n == 0) return(0)
  gg <- simplify(g, remove.multiple = TRUE, remove.loops = !include_loops)
  actual <- ecount(gg)
  possible <- if (include_loops) (n * (n + 1)) %/% 2 else (n * (n - 1)) %/% 2
  if (possible == 0) 0 else actual / possible
}

## getKCore(): L889-936.  Iteratively drop vertices of degree < k until
## convergence; returns NULL when nothing survives.  Equivalent to the set of
## vertices with coreness >= k, computed inside the supplied subgraph.
mcode_kcore <- function(g, k, include_loops = FALSE) {
  if (vcount(g) == 0) return(NULL)
  gg <- simplify(g, remove.multiple = TRUE, remove.loops = !include_loops)
  keep <- which(coreness(gg, mode = "all") >= k)
  if (length(keep) == 0) return(NULL)
  induced_subgraph(g, keep)
}
mcode_kcore_vids <- function(g, k, include_loops = FALSE) {
  if (vcount(g) == 0) return(integer(0))
  gg <- simplify(g, remove.multiple = TRUE, remove.loops = !include_loops)
  which(coreness(gg, mode = "all") >= k)
}

## getHighestKCore(): L946-971.  Returns the largest k with a non-empty k-core
## and that core.  k = 0 when the graph has no edges.
mcode_highest_kcore <- function(g, include_loops = FALSE) {
  if (vcount(g) == 0) return(list(k = 0L, core = NULL))
  gg <- simplify(g, remove.multiple = TRUE, remove.loops = !include_loops)
  cn <- coreness(gg, mode = "all")
  k <- max(cn)
  if (k == 0) return(list(k = 0L, core = NULL))
  list(k = as.integer(k), core = induced_subgraph(g, which(cn >= k)))
}

## calcNodeInfo(): L575-652  +  scoreNode(): L539-546
## Note the two exact behaviours that matter downstream:
##  (i) a vertex with < 2 neighbours gets numNodeNeighbors = 0 and an EMPTY
##      neighbour array, so it can never expand a cluster and never seeds one;
## (ii) nodeNeighbors is the CLOSED neighbourhood: self first, then neighbours
##      in ascending vertex order -- this fixes the DFS visiting order.
mcode_score_graph <- function(g, degree_cutoff = 2L, include_loops = FALSE,
                              verbose = FALSE) {
  n <- vcount(g)
  adj <- as_adj_list(g, mode = "all")
  density <- numeric(n); core_level <- integer(n); core_density <- numeric(n)
  num_nb <- integer(n); nb_list <- vector("list", n)

  for (i in seq_len(n)) {
    nb <- sort(unique(setdiff(as.integer(adj[[i]]), i)))
    if (length(nb) < 2L) {                       # L589-600
      if (length(nb) == 1L) {
        core_level[i] <- 1L; core_density[i] <- 1.0; density[i] <- 1.0
      }
      num_nb[i] <- 0L; nb_list[[i]] <- integer(0)
      next
    }
    nbhd <- c(i, nb)                             # L609-619 (self first)
    sub  <- induced_subgraph(g, nbhd)
    density[i] <- mcode_density(sub, include_loops)
    num_nb[i]  <- length(nbhd)
    hk <- mcode_highest_kcore(sub, include_loops)
    core_level[i] <- hk$k
    core_density[i] <- if (is.null(hk$core)) 0 else mcode_density(hk$core, include_loops)
    nb_list[[i]] <- nbhd
  }
  score <- ifelse(num_nb > degree_cutoff, core_density * core_level, 0)  # L540
  list(score = score, density = density, core_level = core_level,
       core_density = core_density, num_nb = num_nb, nb = nb_list)
}

## getClusterCoreInternal(): L696-742.  Depth-first; the threshold is
## startNodeScore*(1 - nodeScoreCutoff) and startNodeScore is the SEED's score,
## held fixed through the whole recursion.  A node is marked seen before the
## depth test, exactly as in the Java.
mcode_cluster_core <- function(seed, info, seen, node_score_cutoff, max_depth) {
  cluster <- integer(0)
  thr <- info$score[seed] - info$score[seed] * node_score_cutoff
  rec <- function(v, depth) {
    if (isTRUE(seen[[as.character(v)]])) return(invisible(NULL))
    assign(as.character(v), TRUE, envir = seen)
    if (depth > max_depth) return(invisible(NULL))
    nbs <- info$nb[[v]]
    if (length(nbs) == 0L) return(invisible(NULL))
    for (u in nbs) {
      if (!isTRUE(seen[[as.character(u)]]) && info$score[u] >= thr) {
        if (!(u %in% cluster)) cluster <<- c(cluster, u)
        rec(u, depth + 1L)
      }
    }
    invisible(NULL)
  }
  rec(seed, 1L)
  list(cluster = cluster, seen = seen)
}

## fluffClusterBoundary(): L753-787.  Uses each candidate's own neighbourhood
## density; fluffed nodes are NOT added to the global seen map.
mcode_fluff <- function(cluster, seen, info, fluff_density_cutoff) {
  add <- integer(0); internal <- new.env(hash = TRUE, parent = emptyenv())
  for (v in cluster) {
    for (u in info$nb[[v]]) {
      if (!isTRUE(seen[[as.character(u)]]) &&
          !isTRUE(internal[[as.character(u)]]) &&
          info$density[u] > fluff_density_cutoff) {
        add <- c(add, u); assign(as.character(u), TRUE, envir = internal)
      }
    }
  }
  c(cluster, add)
}

## findClusters(): L333-461, incl. filterCluster (L795-803), haircutCluster
## (L812-826), scoreCluster (L555-564) and rank() (L978-995).
## seed_order: the order in which equally scored nodes are tried.  Java uses the
## network's node order; pass a permutation to test tie-break sensitivity.
mcode_find_clusters <- function(g, info = NULL,
                                degree_cutoff = 2L, k_core = 2L,
                                max_depth = 100L, node_score_cutoff = 0.2,
                                haircut = TRUE, fluff = FALSE,
                                fluff_density_cutoff = 0.1,
                                include_loops = FALSE,
                                seed_order = NULL) {
  if (is.null(info))
    info <- mcode_score_graph(g, degree_cutoff, include_loops)
  n <- vcount(g)
  base_order <- if (is.null(seed_order)) seq_len(n) else seed_order
  ord <- base_order[order(-info$score[base_order], seq_along(base_order))]

  seen <- new.env(hash = TRUE, parent = emptyenv())
  clusters <- list()
  for (s in ord) {
    if (isTRUE(seen[[as.character(s)]])) next
    res <- mcode_cluster_core(s, info, seen, node_score_cutoff, max_depth)
    seen <- res$seen
    cc <- res$cluster
    if (length(cc) == 0L) next                        # L390
    if (!(s %in% cc)) cc <- c(cc, s)                  # L392-393
    sub <- induced_subgraph(g, cc)
    if (length(mcode_kcore_vids(sub, k_core, include_loops)) == 0L) next   # filterCluster
    if (haircut) {                                    # L399-400
      keep <- mcode_kcore_vids(sub, 2L, include_loops)
      if (length(keep) > 0L) cc <- cc[keep]
    }
    if (fluff) cc <- mcode_fluff(cc, seen, info, fluff_density_cutoff)
    sub <- induced_subgraph(g, cc)
    clusters[[length(clusters) + 1L]] <- list(
      seed = s, nodes = cc,
      score = mcode_density(sub, include_loops) * vcount(sub),
      size = vcount(sub), edges = ecount(simplify(sub)))
  }
  if (length(clusters) == 0L) return(clusters)
  clusters <- clusters[order(-vapply(clusters, function(x) x$score, 0),
                             seq_along(clusters))]
  for (i in seq_along(clusters)) clusters[[i]]$rank <- i
  clusters
}
