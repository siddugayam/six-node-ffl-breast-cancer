#!/usr/bin/env Rscript
# Test of the reactive-stroma module.
#   Q1  Do the ECM / reactive-stroma genes concentrate in ONE module?
#   Q2  Does that module also hold miR-29a and the TF arm (ETS1, NFKB1, RELA,
#       SP1), whose collagen association tracks compartment composition?
suppressMessages({library(igraph); library(data.table)})
setwd("/path/to/revision")
set.seed(20260912)

g <- readRDS("results/v7/cm2_graph_undirected.rds"); nm <- V(g)$name; N <- length(nm)
P <- readRDS("results/v7/cm2_all_partitions.rds")
s <- fread("results/v7/cm2_run_summary.csv")
nodes <- fread("data/canonical_nodes.tsv"); ntype <- setNames(nodes$type, nodes$name)
prot <- nm[ntype[nm] %in% c("Gene", "TF")]
EMT <- intersect(unique(fread("results/v7/cm2_genesets_msigdb.csv")[
  gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", gene_symbol]), nm)
## secondary, hand-built core-ECM list: collagens, matrix glycoproteins,
## matrix remodellers and myofibroblast markers present in the network.
ECM <- intersect(c("COL1A1","COL3A1","COL7A1","COL18A1","FN1","SPARC","POSTN","LOX",
                   "MMP2","MMP14","THBS1","THBS2","TGFB2","TGFBI","SERPINE1","ACTA2",
                   "TAGLN","VIM","CDH11","FBN1","CTGF","PDGFRB","ITGA5","ITGB1"), nm)
TFarm <- c("ETS1","NFKB1","RELA","SP1")

conc_test <- function(m, set, B = 5000) {
  # concentration of `set` in the single best module, vs. random gene sets of the
  # same size drawn from the protein-coding nodes only (degree ignored).
  tb <- table(m[set]); obs <- max(tb)
  null <- replicate(B, { idx <- sample(prot, length(set)); max(table(m[idx])) })
  list(obs = obs, exp = mean(null), p = (sum(null >= obs) + 1) / (B + 1),
       best_module = as.integer(names(tb)[which.max(tb)]))
}

runs <- s[degenerate == FALSE & k_clusters <= 30, run]   # partitions that give usable modules
out <- rbindlist(lapply(runs, function(k) {
  m <- P[[k]]
  e <- conc_test(m, EMT); c_ <- conc_test(m, ECM)
  colmod <- m["COL1A1"]; mem <- names(m)[m == colmod]
  pm <- intersect(mem, prot)
  data.table(
    run = k, k_clusters = s[run == k, k_clusters], modularity = s[run == k, modularity],
    EMT_max_in_one_module = e$obs, EMT_expected = round(e$exp, 1), EMT_conc_p = e$p,
    EMT_best_module = e$best_module, EMT_module_is_COL1A1_module = e$best_module == colmod,
    ECM_max_in_one_module = c_$obs, ECM_expected = round(c_$exp, 1), ECM_conc_p = c_$p,
    ECM_module_is_COL1A1_module = c_$best_module == colmod,
    COL1A1_module_size = length(mem),
    ECM_in_COL1A1_module = length(intersect(mem, ECM)),
    ECM_fold = round((length(intersect(mem, ECM)) / length(pm)) /
                     (length(ECM) / length(prot)), 2),
    EMT_in_COL1A1_module = length(intersect(mem, EMT)),
    EMT_fold = round((length(intersect(mem, EMT)) / length(pm)) /
                     (length(EMT) / length(prot)), 2),
    miR29a_in = unname(m["hsa-miR-29a"] == colmod),
    ETS1_in = unname(m["ETS1"] == colmod), NFKB1_in = unname(m["NFKB1"] == colmod),
    RELA_in = unname(m["RELA"] == colmod), SP1_in = unname(m["SP1"] == colmod),
    TFarm_all_in = all(m[TFarm] == colmod),
    EZH2_with_miR101 = unname(m["EZH2"] == m["hsa-miR-101"]))
}))
setorder(out, -modularity)
fwrite(out, "results/v7/cm2_stroma_module_test.csv")
cat("== ECM / EMT concentration and the COL1A1 module ==\n")
print(out[, .(run, k = k_clusters, Q = round(modularity, 3),
              EMTmax = EMT_max_in_one_module, EMTexp = EMT_expected, EMTp = EMT_conc_p,
              ECMmax = ECM_max_in_one_module, ECMexp = ECM_expected, ECMp = ECM_conc_p,
              colsize = COL1A1_module_size, ECMfold = ECM_fold, EMTfold = EMT_fold,
              m29 = miR29a_in, TFall = TFarm_all_in)], nrows = 40)

cat("\n== where the TF arm sits relative to COL1A1 ==\n")
print(out[, .(run, ETS1_in, NFKB1_in, RELA_in, SP1_in, TFarm_all_in, EZH2_with_miR101)],
      nrows = 40)
