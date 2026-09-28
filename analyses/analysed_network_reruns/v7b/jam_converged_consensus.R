## jActiveModules consensus restricted to the primary seeds whose rank-1 module converged (>= 100
## nodes), because without the legacy edges one of the ten seeds stalls at a 2-node module and the
## all-seed consensus rule (jam_02_characterise.R) then collapses to a single node.  Everything else is
## as in jam_02 (hypergeometric overlap with the prioritised 30, all 587 nodes as universe) and jam_05
## (z-rank-stratified null, set.seed(1), B = 10,000).  Run for both networks; with all ten seeds
## converged (full network) it must return the stored consensus of 205 nodes.
suppressMessages({library(data.table)})
HERE <- "/path/to/revision/analyses/analysed_network_reruns/v7b"
res <- rbindlist(lapply(c("full", "nolegacy"), function(mode) {
  R <- file.path(HERE, paste0("sandbox_", mode), "results/v7")
  mods <- fread(file.path(R, "jam_modules_all_runs.csv")); ns <- fread(file.path(R, "jam_node_scores.csv"))
  p30 <- fread(file.path(HERE, paste0("sandbox_", mode), "results/v5/tables/Table4_prioritised_30.csv"))$name
  top <- mods[variant == "primary" & rank == 1]
  conv <- top[size >= 100]
  L <- lapply(conv$members, function(x) strsplit(x, ";", fixed = TRUE)[[1]])
  cons <- Reduce(intersect, L); N <- nrow(ns)
  h <- length(intersect(cons, p30))
  p_h <- phyper(h - 1, length(p30), N - length(p30), length(cons), lower.tail = FALSE)
  set.seed(1); B <- 10000
  zr <- rank(ns$z); names(zr) <- ns$node
  str <- cut(zr, breaks = quantile(zr, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  tab <- table(str[match(p30, ns$node)])
  null <- replicate(B, { pick <- unlist(lapply(names(tab), function(s) sample(ns$node[str == as.integer(s)], tab[[s]])))
                         length(intersect(pick, cons)) })
  mir29 <- sapply(c("hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c"), function(x) sum(vapply(L, function(l) x %in% l, TRUE)))
  data.table(network = mode, primary_seeds = nrow(top), converged_seeds = nrow(conv),
             stalled_seed_sizes = paste(top[size < 100, size], collapse = ";"),
             consensus_size = length(cons), prior30_in_consensus = h, hyper_p = signif(p_h, 3),
             zmatched_null_mean = round(mean(null), 2), p_zmatched = round((sum(null >= h) + 1) / (B + 1), 4),
             COL1A1_in = "COL1A1" %in% cons, COL3A1_in = "COL3A1" %in% cons,
             miR29abc_converged_seeds_selecting = paste(mir29, collapse = "/"))
}))
fwrite(res, file.path(HERE, "jam_converged_consensus.csv")); print(res)
