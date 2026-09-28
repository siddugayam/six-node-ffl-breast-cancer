## jam_common.R -- shared machinery for the jActiveModules re-implementation
## Ideker, Ozier, Schwikowski & Siegel (2002) Bioinformatics 18:S233-S240.
suppressPackageStartupMessages({
  library(data.table); library(igraph); library(Rcpp); library(matrixStats)
})
BASE <- "/path/to/revision"
OUT  <- file.path(BASE, "results/v7")
SCR  <- file.path(BASE, "scripts/v7")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

jam_log <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

## ---------------------------------------------------------------- network --
jam_graph <- function() {
  nodes <- fread(file.path(BASE, "data/canonical_nodes.tsv"))
  ed    <- fread(file.path(BASE, "data/canonical_edges.tsv"))
  ed    <- ed[source != target]
  g <- graph_from_data_frame(unique(ed[, .(source, target)]), directed = FALSE,
                             vertices = nodes[, .(name, type)])
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  g
}

## ------------------------------------------------------------- node scores --
## z_i = Phi^-1(1 - p_i), computed in the lower tail for numerical accuracy.
p_to_z <- function(p, cap = 30) {
  z <- qnorm(p, lower.tail = FALSE)
  z[!is.finite(z) & !is.na(z) & z > 0] <-  cap
  z[!is.finite(z) & !is.na(z) & z < 0] <- -cap
  pmin(pmax(z, -cap), cap)
}

## score vector for the network's nodes
## pcol   : "adj.P.Val" (default, as specified) or "P.Value"
## mirna  : TRUE  -> use the matched TCGA-BRCA miRNA differential expression
##          FALSE -> miRNA nodes get z = 0 (gene table only, as literally given)
## signed : TRUE  -> z is multiplied by sign(logFC)
jam_scores <- function(g, pcol = "adj.P.Val", mirna = TRUE, signed = FALSE) {
  nm  <- V(g)$name; ty <- V(g)$type
  gde <- fread(file.path(BASE, "results/v2/v2_DE_genes.csv"))
  setnames(gde, "feature", "id")
  mde <- fread(file.path(BASE, "results/BRCA_DEX_mirnas.csv"))
  setnames(mde, "feature", "id")
  de  <- rbind(gde[, .(id, logFC, P.Value, adj.P.Val)], mde[, .(id, logFC, P.Value, adj.P.Val)])
  de  <- de[!duplicated(id)]
  i   <- match(nm, de$id)
  p   <- de[[pcol]][i]; fc <- de$logFC[i]
  if (!mirna) { p[ty == "miRNA"] <- NA; fc[ty == "miRNA"] <- NA }
  z <- p_to_z(p)
  if (signed) z <- z * sign(fc)
  z[is.na(z)] <- 0                      # no expression evidence -> no contribution
  names(z) <- nm
  list(z = z, p = p, logFC = fc, has_data = !is.na(p))
}

## ------------------------------------------------------------- calibration --
## Ideker's correction: mu_k and sigma_k of the aggregate z over randomly
## sampled node sets of size k.  `mode`:
##   "sets"      Monte-Carlo random node sets, sampled WITHOUT replacement (the
##               published calibration; carries the finite-population factor)
##   "indep"     random nodes sampled WITH replacement (no finite-population
##               factor; sigma_k = sd(z) for every k)
##   "connected" randomly grown CONNECTED subnetworks of size k (stricter)
jam_calibrate <- function(z, g = NULL, mode = c("sets", "indep", "connected"),
                          B = 10000, seed = 1) {
  mode <- match.arg(mode); n <- length(z); set.seed(seed)
  mu <- numeric(n + 1); sd_ <- numeric(n + 1)
  if (mode %in% c("sets", "indep")) {
    ## One random permutation of the node scores supplies, through its running
    ## sum, one uniformly sampled set of every size k = 1..N at once.
    rep_ <- (mode == "indep")
    M <- vapply(seq_len(B), function(b)
      cumsum(if (rep_) sample(z, n, replace = TRUE) else sample(z)), numeric(n))
    M <- M / sqrt(seq_len(n))
    mu[-1]  <- rowMeans(M)
    sd_[-1] <- matrixStats::rowSds(M)
  } else {
    al <- jam_adj(g)
    M  <- grow_null_cpp(al$flat, al$start, al$len, z[al$order], B, n)
    mu[-1]  <- colMeans(M)
    sd_[-1] <- matrixStats::colSds(M)
  }
  sd_[sd_ < 1e-9] <- 1e-9
  mu[1] <- 0; sd_[1] <- 1e-9           # index 1 == size 0, never used
  list(mu = mu, sd = sd_, mode = mode, B = B)
}

## 0-based CSR adjacency in the vertex order of g
jam_adj <- function(g) {
  al <- as_adj_list(g, mode = "all")
  len <- vapply(al, length, 1L)
  list(flat = as.integer(unlist(al)) - 1L,
       start = as.integer(c(0L, cumsum(len)[-length(len)])),
       len = as.integer(len), order = V(g)$name)
}

## score one node set with the calibrated score
jam_score_set <- function(nodes, z, cal) {
  k <- length(nodes); if (k == 0) return(NA_real_)
  (sum(z[nodes]) / sqrt(k) - cal$mu[k + 1]) / cal$sd[k + 1]
}

## ------------------------------------------------------------------ search --
## One annealing run.  Returns the highest-scoring connected component found.
jam_run <- function(g, z, cal, seed, iters = 2e5, T0 = 1, T1 = 1e-4,
                    allowed = NULL, p_init = 0.5, trace_every = 1000) {
  al <- jam_adj(g); n <- length(z)
  stopifnot(identical(al$order, names(z)))
  if (is.null(allowed)) allowed <- rep(TRUE, n)
  set.seed(seed)
  init <- runif(n) < p_init
  r <- sa_search_cpp(al$flat, al$start, al$len, as.numeric(z), cal$mu, cal$sd,
                     init, allowed, as.integer(iters), T0, T1, as.integer(trace_every))
  bs <- r$best_state
  sub <- induced_subgraph(g, which(bs))
  cc  <- components(sub)
  mem <- split(V(sub)$name, cc$membership)
  sc  <- vapply(mem, jam_score_set, 0.0, z = z, cal = cal)
  ord <- order(sc, decreasing = TRUE)
  fsub <- induced_subgraph(g, which(r$final_state))
  fcc  <- components(fsub)
  fmem <- split(V(fsub)$name, fcc$membership)
  fsc  <- vapply(fmem, jam_score_set, 0.0, z = z, cal = cal)
  ford <- order(fsc, decreasing = TRUE)
  list(seed = seed, best_score = r$best_score, final_score = r$final_score,
       modules = mem[ord], scores = sc[ord],
       final_modules = fmem[ford], final_scores = fsc[ford],
       n_active_best = sum(bs), n_active_final = sum(r$final_state),
       accept_rate = r$n_accept / iters, trace = r$trace, iters = iters)
}

## Sequential module extraction ("peeling"): find the top component, remove its
## nodes, search again.  Needed because after annealing the final active set is
## one giant component plus singletons, so components 2..5 of a single run are
## not modules in any useful sense.
jam_run_peel <- function(g, z, cal, seed, n_mod = 5, ...) {
  n <- length(z); allowed <- rep(TRUE, n); out <- list()
  for (m in seq_len(n_mod)) {
    r <- jam_run(g, z, cal, seed = seed * 1000 + m, allowed = allowed, ...)
    if (length(r$modules) == 0) break
    mod <- r$modules[[1]]
    out[[m]] <- list(rank = m, nodes = mod, score = unname(r$scores[1]),
                     size = length(mod), seed = seed,
                     accept_rate = r$accept_rate, n_active_best = r$n_active_best)
    allowed[match(mod, names(z))] <- FALSE
  }
  out
}
