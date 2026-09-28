#!/usr/bin/env Rscript
## jam_05_circularity_check.R -- the agreement between the active module and the
## paper's prioritised nodes is not independent evidence: the prioritisation
## score contains a differential-expression domain built from the same
## TCGA-BRCA contrast that supplies the jActiveModules node scores.
source("/path/to/revision/scripts/14_module_detection/jam_common.R")
ns <- fread(file.path(OUT, "jam_node_scores.csv"))
s4 <- fread(file.path(BASE, "results/v5/tables/TableS4_node_prioritisation_full.csv"))
m  <- merge(ns, s4[, .(node = name, priority, T_topology, E_evidence, D_expression,
                       R_replication, M_multiomic, C_clinical, F_functional)], by = "node")
f  <- fread(file.path(OUT, "jam_node_selection_frequency.csv"))
m  <- merge(m, f[, .(node, freq_top_module)], by = "node")
cc <- rbindlist(lapply(c("priority", "T_topology", "E_evidence", "D_expression",
                         "R_replication", "M_multiomic", "C_clinical", "F_functional"),
  function(v) {
    x <- suppressWarnings(as.numeric(m[[v]]))
    ok <- is.finite(x) & is.finite(m$z)
    data.table(domain = v, n = sum(ok),
               rho_with_z = if (sum(ok) > 3) cor(x[ok], m$z[ok], method = "spearman") else NA_real_,
               rho_with_membership = if (sum(ok) > 3)
                 cor(x[ok], m$freq_top_module[ok], method = "spearman") else NA_real_)
  }))
fwrite(cc, file.path(OUT, "jam_prioritisation_circularity.csv"))
print(cc)
## partial check: overlap of the module with the prioritised 30 after matching on z
p30 <- fread(file.path(BASE, "results/v5/tables/Table4_prioritised_30.csv"))$name
cons <- f[n_seeds_in_top_module == 10]$node
set.seed(1); B <- 10000
zr <- rank(ns$z); names(zr) <- ns$node
obs <- length(intersect(p30, cons))
## z-matched random 30-node sets: sample nodes with the same z-rank strata
str <- cut(zr, breaks = quantile(zr, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
tab <- table(str[match(p30, ns$node)])
null <- replicate(B, {
  pick <- unlist(lapply(names(tab), function(s)
    sample(ns$node[str == as.integer(s)], tab[[s]])))
  length(intersect(pick, cons))
})
res <- data.table(observed = obs, null_mean = mean(null), null_sd = sd(null),
                  p_zmatched = (sum(null >= obs) + 1) / (B + 1))
fwrite(res, file.path(OUT, "jam_prioritisation_zmatched_test.csv"))
print(res)
