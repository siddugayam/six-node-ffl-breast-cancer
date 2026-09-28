#!/usr/bin/env Rscript
## ============================================================================
## 14_overlap_sweep.R
## Robustness of the ExIR-driver / FFL-hub overlap to the (arbitrary) choice of
## how many top-ranked ExIR drivers count as "drivers". Sweeps the driver-set
## size from 25 to 4311 and recomputes the hypergeometric p in the universe of
## canonical-network nodes that ExIR could rank.
## ============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- "/path/to/revision"
LOG  <- file.path(ROOT, "logs", "exir_agent.log")
say <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | 14 | ", paste0(..., collapse = ""))
  cat(msg, "\n", file = LOG, append = TRUE); cat(msg, "\n"); flush.console()
}
say("=== 14 START ===")

cls  <- fread(file.path(ROOT, "results", "exir_classification.csv"))
hub  <- fread(file.path(ROOT, "results", "ffl_network_hubs.csv"))
meta <- readRDS(file.path(ROOT, "data", "exir_inputs_meta.rds"))
nodes<- fread(file.path(ROOT, "data", "canonical_nodes.tsv"))

drv <- cls[class == "Driver"][order(rank)]
U_all <- unique(c(cls$feature))                       # every feature ExIR ranked
U_net <- intersect(nodes$name, U_all)
hub_u <- hub[hub_centrality == TRUE | hub_ffl == TRUE, node]
say("universe(all ExIR-ranked features)=", length(U_all),
    "  universe(network nodes)=", length(U_net),
    "  hubs(union)=", length(hub_u), "  hubs in network universe=", length(intersect(hub_u, U_net)))

sweep <- rbindlist(lapply(c(25, 50, 70, 100, 200, 300, 500, 1000, 2000, 4311), function(k) {
  A <- head(drv$feature, k)
  rbindlist(lapply(list(list("network_nodes", U_net), list("all_ExIR_features", U_all)),
    function(u) {
      U <- u[[2]]
      a <- intersect(A, U); b <- intersect(hub_u, U); ov <- intersect(a, b)
      e <- length(a) * length(b) / length(U)
      data.table(top_k_drivers = k, universe_name = u[[1]], universe = length(U),
                 n_drivers_in_U = length(a), n_hubs_in_U = length(b), overlap = length(ov),
                 expected = e, fold = length(ov)/max(e,1e-12),
                 p_hyper = phyper(length(ov)-1, length(b), length(U)-length(b), length(a),
                                  lower.tail = FALSE),
                 overlap_features = paste(sort(ov), collapse = ";"))
    }))
}))
fwrite(sweep, file.path(ROOT, "results", "exir_overlap_sweep.csv"))
for (i in seq_len(nrow(sweep))) {
  say(sprintf("top%-5d %-18s |U|=%5d A=%4d B=%3d ov=%3d exp=%6.2f fold=%.2f p=%.3g",
      sweep$top_k_drivers[i], sweep$universe_name[i], sweep$universe[i],
      sweep$n_drivers_in_U[i], sweep$n_hubs_in_U[i], sweep$overlap[i],
      sweep$expected[i], sweep$fold[i], sweep$p_hyper[i]))
}
say("wrote results/exir_overlap_sweep.csv rows=", nrow(sweep))
say("min p across the whole sweep: ", format(min(sweep$p_hyper), digits = 3))
say("=== 14 DONE ===")
