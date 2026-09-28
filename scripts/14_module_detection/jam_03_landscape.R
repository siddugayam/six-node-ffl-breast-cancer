#!/usr/bin/env Rscript
## jam_03_landscape.R -- the shape of the Ideker objective on this network.
## For every size k, the highest-scoring set of size k ignoring connectivity is
## the top-k nodes by z.  Because the network is a single dense component these
## sets are almost always connected, so they are (near-)global optima of the
## search and let us ask (a) where the objective actually peaks under each null
## and (b) how close simulated annealing got.
source("/path/to/revision/scripts/14_module_detection/jam_common.R")
sourceCpp(file.path(SCR, "jam_engine.cpp"))
R <- readRDS(file.path(BASE, "cache/v7/jam/jam_runs.rds"))
g <- R$g; z <- R$scores$primary$z; nodes <- V(g)$name
ord <- names(sort(z, decreasing = TRUE))

land <- rbindlist(lapply(names(R$cals), function(cn) {
  cal <- R$cals[[cn]]
  zz <- if (cn == "genes") R$scores$genes$z else if (cn == "rawp") R$scores$rawp$z else
        if (cn == "signed") R$scores$signed$z else z
  o <- names(sort(zz, decreasing = TRUE))
  cs <- cumsum(zz[o])
  sc <- (cs / sqrt(seq_along(o)) - cal$mu[-1]) / cal$sd[-1]
  data.table(null_model = cn, k = seq_along(o), score_topk = sc)
}))
fwrite(land, file.path(OUT, "jam_score_landscape_topk.csv"))
best <- land[, .(k_opt = k[which.max(score_topk)], score_opt = max(score_topk)), by = null_model]
## connectivity of the optimal top-k set
best[, n_components := vapply(k_opt, function(k) as.integer(components(induced_subgraph(g, ord[1:k]))$no), 0L)]
print(best)
fwrite(best, file.path(OUT, "jam_score_landscape_optima.csv"))

## diagnostic: does the collapse to a single node under the independent-sampling
## null persist when the annealing is run 10x longer?
cal <- R$cals$indep
diag <- rbindlist(lapply(c(2e5, 2e6), function(it) rbindlist(lapply(1:3, function(s) {
  r <- jam_run(g, z, cal, seed = s, iters = it)
  data.table(null_model = "indep", iters = it, seed = s,
             size = length(r$modules[[1]]), score = unname(r$scores[1]))
}))))
print(diag)
fwrite(diag, file.path(OUT, "jam_indep_null_iteration_diagnostic.csv"))
jam_log("landscape done")
