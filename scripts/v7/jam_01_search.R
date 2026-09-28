#!/usr/bin/env Rscript
## jam_01_search.R -- jActiveModules (Ideker et al. 2002) active-subnetwork
## search on the canonical miRNA-TF network, primary run + sensitivity grid.
source("/path/to/revision/scripts/v7/jam_common.R")
sourceCpp(file.path(SCR, "jam_engine.cpp"))
set.seed(20260908)

g  <- jam_graph()
jam_log("network: ", vcount(g), " nodes, ", ecount(g), " undirected edges, ",
        components(g)$no, " component(s)")

## ------------------------------------------------------------ node scores --
sc_primary <- jam_scores(g, "adj.P.Val", mirna = TRUE,  signed = FALSE)
sc_genes   <- jam_scores(g, "adj.P.Val", mirna = FALSE, signed = FALSE)
sc_rawp    <- jam_scores(g, "P.Value",   mirna = TRUE,  signed = FALSE)
sc_signed  <- jam_scores(g, "adj.P.Val", mirna = TRUE,  signed = TRUE)

fwrite(data.table(node = V(g)$name, type = V(g)$type, degree = degree(g),
                  adj_P = sc_primary$p, logFC = sc_primary$logFC,
                  z = sc_primary$z, z_signed = sc_signed$z,
                  z_rawp = sc_rawp$z, has_DE = sc_primary$has_data),
       file.path(OUT, "jam_node_scores.csv"))

## ----------------------------------------------------------- calibrations --
cal_sets  <- jam_calibrate(sc_primary$z, mode = "sets",      B = 10000)
cal_indep <- jam_calibrate(sc_primary$z, mode = "indep",     B = 10000)
cal_conn  <- jam_calibrate(sc_primary$z, g, mode = "connected", B = 10000)
cal_genes <- jam_calibrate(sc_genes$z,  mode = "sets", B = 10000)
cal_rawp  <- jam_calibrate(sc_rawp$z,   mode = "sets", B = 10000)
cal_sign  <- jam_calibrate(sc_signed$z, mode = "sets", B = 10000)

## transcriptome background: random gene sets drawn from all measured features
gde <- fread(file.path(BASE, "results/v2/v2_DE_genes.csv"))
mde <- fread(file.path(BASE, "results/BRCA_DEX_mirnas.csv"))
zall <- p_to_z(c(gde$adj.P.Val, mde$adj.P.Val))
set.seed(11); n <- vcount(g)
M <- vapply(seq_len(10000), function(b) cumsum(sample(zall, n)), numeric(n)) / sqrt(seq_len(n))
cal_tx <- list(mu = c(0, rowMeans(M)), sd = pmax(c(1e-9, matrixStats::rowSds(M)), 1e-9),
               mode = "transcriptome", B = 10000)
jam_log("z over all ", length(zall), " measured features: mean ", round(mean(zall), 2),
        " vs network nodes ", round(mean(sc_primary$z), 2))

cal_tab <- rbindlist(lapply(list(sets = cal_sets, indep = cal_indep, connected = cal_conn,
                                 transcriptome = cal_tx),
  function(cl) data.table(k = 1:n, mu_k = cl$mu[-1], sd_k = cl$sd[-1])), idcol = "null_model")
fwrite(cal_tab, file.path(OUT, "jam_calibration.csv"))

## ------------------------------------------------------------------- runs --
run_variant <- function(tag, z, cal, seeds, n_mod, iters, T0 = 1, T1 = 1e-4) {
  jam_log("variant ", tag, ": ", length(seeds), " seeds x ", n_mod, " modules, ",
          format(iters, scientific = FALSE), " iterations, null=", cal$mode)
  res <- list()
  for (s in seeds) {
    t0 <- Sys.time()
    pk <- jam_run_peel(g, z, cal, seed = s, n_mod = n_mod, iters = iters, T0 = T0, T1 = T1)
    for (m in pk) res[[length(res) + 1]] <- data.table(
      variant = tag, seed = s, rank = m$rank, size = m$size, score = m$score,
      accept_rate = m$accept_rate, n_active_best = m$n_active_best,
      members = paste(sort(m$nodes), collapse = ";"))
    jam_log("   seed ", s, ": top size ", pk[[1]]$size, " score ",
            round(pk[[1]]$score, 2), "  (", round(difftime(Sys.time(), t0, units = "secs"), 1), "s)")
  }
  rbindlist(res)
}

SEEDS_P <- 1:10; SEEDS_S <- 1:5
all <- list()
all[["primary"]] <- run_variant("primary", sc_primary$z, cal_sets, SEEDS_P, 5, 2e5)
all[["S1_genes_only"]]   <- run_variant("S1_genes_only",   sc_genes$z,   cal_genes, SEEDS_S, 3, 2e5)
all[["S2_raw_p"]]        <- run_variant("S2_raw_p",        sc_rawp$z,    cal_rawp,  SEEDS_S, 3, 2e5)
all[["S3_indep_null"]]   <- run_variant("S3_indep_null",   sc_primary$z, cal_indep, SEEDS_S, 3, 2e5)
all[["S4_connected_null"]] <- run_variant("S4_connected_null", sc_primary$z, cal_conn, SEEDS_S, 3, 2e5)
all[["S5a_iters_2e4"]]   <- run_variant("S5a_iters_2e4",   sc_primary$z, cal_sets,  SEEDS_S, 3, 2e4)
all[["S5b_iters_1e6"]]   <- run_variant("S5b_iters_1e6",   sc_primary$z, cal_sets,  SEEDS_S, 3, 1e6)
all[["S6_signed_z"]]     <- run_variant("S6_signed_z",     sc_signed$z,  cal_sign,  SEEDS_S, 3, 2e5)
all[["S7_transcriptome_null"]] <- run_variant("S7_transcriptome_null", sc_primary$z, cal_tx, SEEDS_S, 3, 2e5)

mods <- rbindlist(all)
fwrite(mods, file.path(OUT, "jam_modules_all_runs.csv"))
saveRDS(list(g = g, scores = list(primary = sc_primary, genes = sc_genes,
                                  rawp = sc_rawp, signed = sc_signed),
             cals = list(sets = cal_sets, indep = cal_indep, connected = cal_conn,
                         genes = cal_genes, rawp = cal_rawp, signed = cal_sign,
                         transcriptome = cal_tx),
             modules = mods),
        file.path(BASE, "cache/v7/jam/jam_runs.rds"))
jam_log("DONE: ", nrow(mods), " modules over ", uniqueN(mods$variant), " variants")
