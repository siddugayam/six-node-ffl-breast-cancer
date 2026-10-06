#!/usr/bin/env Rscript
# =============================================================================
# go_01_assemble_modules.R
# Assemble EVERY module of size >= 5 recovered by the four v7 module-detection
# runs into one tidy inventory, ready for BiNGO-equivalent ORA.
#
# SELECTION RULE (stated before any enrichment was run):
#   TIER "primary"   = the partition each method itself designated as its
#                      headline / reference result.
#   TIER "secondary" = the sensitivity or control partitions each method
#                      flagged as decision-relevant in its own report.
#   Within each partition, every module with >= 5 NODES is taken.  Node count,
#   not gene count, is the floor, because that is how the methods report size;
#   the gene/TF subset actually testable by GO/KEGG is recorded separately.
#
# miRNA nodes carry no GO BP or KEGG annotation and are dropped before ORA
# (same convention as scripts/12_enrichment_and_metabolism/over_representation/09_enrichment_ora.R).  Modules whose gene+TF
# content falls below 3 after that drop cannot yield a term at the 3-gene floor
# and are flagged untestable rather than silently discarded.
# =============================================================================
suppressPackageStartupMessages({library(data.table)})

REV <- "/path/to/revision"
V7  <- file.path(REV, "results/v7")
LOGF <- file.path(REV, "logs", "go_01_assemble_modules.log")
lg <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOGF, append = TRUE) }
cat("", file = LOGF)

nodes <- fread(file.path(V7, "cm2_nodes.tsv"))
setnames(nodes, c("name", "type"))
TYPE <- setNames(nodes$type, nodes$name)
lg("node types: ", paste(names(table(nodes$type)), table(nodes$type), collapse = " / "))

recs <- list()
add <- function(method, partition, tier, module, members, note = "") {
  members <- trimws(members); members <- members[nzchar(members)]
  members <- unique(members)
  bad <- setdiff(members, names(TYPE))
  if (length(bad)) lg("  WARNING unknown node(s) in ", partition, "/", module, ": ", paste(bad, collapse = ","))
  ty <- TYPE[members]
  recs[[length(recs) + 1L]] <<- data.table(
    method = method, partition = partition, tier = tier, module = as.character(module),
    size = length(members),
    n_gene = sum(ty == "Gene", na.rm = TRUE), n_TF = sum(ty == "TF", na.rm = TRUE),
    n_miRNA = sum(ty == "miRNA", na.rm = TRUE),
    note = note, members = paste(members, collapse = ";"))
}

## ---------------------------------------------------------------- 1. MCODE
mc <- fread(file.path(V7, "mcode_04_clusters_default.csv"))
mc <- mc[run == "default"][order(rank)]
lg("MCODE default: ", nrow(mc), " clusters, sizes ", paste(mc$size, collapse = ","))
for (i in seq_len(nrow(mc)))
  add("MCODE", "MCODE_default", "primary", paste0("C", mc$rank[i]),
      strsplit(mc$members[i], ";")[[1]],
      sprintf("seed=%s score=%.3f", mc$seed[i], mc$score[i]))

ctl <- fread(file.path(V7, "mcode_16_controls.csv"))
for (netname in c("strong_and_no_miRNA_miRNA", "drop_miRNA_miRNA_and_gene_gene")) {
  s <- ctl[network == netname][order(rank)]
  lg("MCODE control ", netname, ": ", nrow(s), " clusters, sizes ", paste(s$size, collapse = ","))
  for (i in seq_len(nrow(s)))
    add("MCODE", paste0("MCODE_", netname), "secondary", paste0("C", s$rank[i]),
        strsplit(s$members[i], ";")[[1]],
        sprintf("seed=%s score=%.3f", s$seed[i], s$score[i]))
}

## --------------------------------------------------------- 2. clusterMaker2
cy <- fread(file.path(V7, "cm2_cytoscape_partitions.tsv"))
stopifnot("louvain_cm2" %in% names(cy), "name" %in% names(cy))
sp <- split(cy$name, cy$louvain_cm2)
lg("cm2 Louvain (clusterMaker2 Java): ", length(sp), " modules, sizes ",
   paste(sort(lengths(sp), decreasing = TRUE), collapse = ","))
for (k in names(sp)) add("clusterMaker2", "cm2_Louvain_clusterMaker2", "primary", paste0("M", k), sp[[k]],
                         "reference partition (battery default)")

ig <- fread(file.path(V7, "cm2_igraph_partitions.csv"))
sp <- split(ig$name, ig[["Leiden_mod_r1.0"]])
lg("cm2 Leiden r1.0 (highest Q in battery): ", length(sp), " modules, sizes ",
   paste(sort(lengths(sp), decreasing = TRUE), collapse = ","))
for (k in names(sp)) add("clusterMaker2", "cm2_Leiden_mod_r1.0", "secondary", paste0("M", k), sp[[k]],
                         "highest-modularity partition in the battery (Q=0.272)")

cons <- fread(file.path(V7, "cm2_consensus_modules.csv"))
lg("cm2 consensus (>=50% co-clustering over 41 partitions): ", nrow(cons),
   " modules, ", sum(cons$size >= 5), " of size >=5")
for (i in seq_len(nrow(cons)))
  add("clusterMaker2", "cm2_consensus", "primary", paste0("K", cons$cluster[i]),
      strsplit(cons$members[i], "[ ]+")[[1]],
      sprintf("consensus density=%.3f", cons$density[i]))

## ----------------------------------------------------------------- 3. BioNet
bn <- fread(file.path(V7, "bionet_module_summary.csv"))
tiermap <- c(primary = "primary", sens_mid = "secondary", sens_vign = "secondary")
for (i in seq_len(nrow(bn))) {
  r <- bn[i]
  add("BioNet", paste0("BioNet_", r$run), tiermap[[r$run]], "M1",
      strsplit(r$members, ";")[[1]],
      sprintf("FDR=%.3g n=%d score=%.1f", r$fdr, r$module_n, r$module_score))
  lg("BioNet ", r$run, ": 1 module, size ", r$module_n)
}

## -------------------------------------------------------- 4. jActiveModules
jf <- fread(file.path(V7, "jam_node_selection_frequency.csv"))
core <- jf[n_seeds_in_top_module == 10L]$node
lg("jAM consensus core (selected in 10/10 seeds): ", length(core), " nodes")
add("jActiveModules", "jAM_consensus_core", "primary", "M1", core,
    "nodes in the top module of all 10 annealing seeds")

ja <- fread(file.path(V7, "jam_modules_all_runs.csv"))
j1 <- ja[variant == "primary" & seed == 1L][order(rank)]
lg("jAM seed 1: ", nrow(j1), " modules, sizes ", paste(j1$size, collapse = ","))
for (i in seq_len(nrow(j1)))
  add("jActiveModules", "jAM_seed1", "secondary", paste0("M", j1$rank[i]),
      strsplit(j1$members[i], ";")[[1]],
      sprintf("score=%.3f", j1$score[i]))

## ---------------------------------------------------------------- write out
inv <- rbindlist(recs)
inv[, n_geneTF := n_gene + n_TF]
inv[, above_size_floor := size >= 5]
inv[, testable_GO := above_size_floor & n_geneTF >= 3]
setorder(inv, method, partition, -size)
fwrite(inv, file.path(V7, "goora_00_module_inventory_all.csv"))
keep <- inv[above_size_floor == TRUE]
fwrite(keep, file.path(V7, "goora_00_module_inventory_size5.csv"))

lg("")
lg("TOTAL modules assembled: ", nrow(inv))
lg("  size >= 5 (to be tested): ", nrow(keep))
lg("  size <  5 (below floor, excluded): ", nrow(inv) - nrow(keep))
lg("  of the size>=5 set, with < 3 gene/TF members (untestable by GO/KEGG): ",
   sum(keep$n_geneTF < 3))
lg("")
print(keep[, .(n_modules = .N, min_size = min(size), max_size = max(size),
               median_geneTF = as.numeric(median(n_geneTF)),
               n_untestable = sum(n_geneTF < 3)), by = .(method, partition, tier)])
