#!/usr/bin/env Rscript
# v7 / 04 : harmonise every partition (clusterMaker2 2.3.4 via CyREST + igraph),
#           then the four questions the convergent-validity test asks:
#             (a) how many clusters, what modularity, what size distribution
#             (b) which cluster holds COL1A1 / COL3A1 / hsa-miR-29a
#             (c) do the algorithms agree (adjusted Rand index)
#             (d) do the leading modules overlap the 30 prioritised nodes
suppressMessages({library(igraph); library(data.table); library(mclust)})
setwd("/path/to/revision")

g    <- readRDS("results/v7/cm2_graph_undirected.rds")
nm   <- V(g)$name
nodes<- fread("data/canonical_nodes.tsv")
setkey(nodes, name)
ntype <- setNames(nodes$type, nodes$name)

## ---------- 1. collect partitions ----------------------------------------
ig <- fread("results/v7/cm2_igraph_partitions.csv")
setkey(ig, name); ig <- ig[nm]
P <- as.list(ig[, -"name"])

cy <- fread("results/v7/cm2_cytoscape_partitions.tsv", na.strings = c("", "NA", "None"))
setnames(cy, "name", "node")
setkey(cy, node); cy <- cy[nm]
cyc <- setdiff(names(cy), c("SUID", "node", "selected", "shared name"))
for (cc in cyc) {
  v <- cy[[cc]]
  # clusterMaker2 omits singletons from the cluster attribute: NA = its own cluster
  na <- is.na(v)
  v <- as.character(v)
  v[na] <- paste0("singleton_", seq_len(sum(na)))
  P[[paste0("cm2_", cc)]] <- as.integer(factor(v))
}
for (k in names(P)) names(P[[k]]) <- nm
saveRDS(P, "results/v7/cm2_all_partitions.rds")

## ---------- 2. per-run summary -------------------------------------------
q_contrib <- function(g, m) {                       # per-module modularity term
  mm <- ecount(g); A <- as_adjacency_matrix(g, sparse = TRUE)
  dt <- degree(g)
  sapply(sort(unique(m)), function(cl) {
    idx <- which(m == cl)
    e_in <- sum(A[idx, idx]) / 2
    d    <- sum(dt[idx])
    e_in / mm - (d / (2 * mm))^2
  })
}
summ <- rbindlist(lapply(names(P), function(k) {
  m <- P[[k]]; tb <- sort(table(m), decreasing = TRUE)
  data.table(run = k, k_clusters = length(tb), modularity = modularity(g, m),
             largest = tb[1], second = if (length(tb) > 1) tb[2] else 0L,
             third = if (length(tb) > 2) tb[3] else 0L,
             median_size = median(as.numeric(tb)), singletons = sum(tb == 1),
             n_ge5 = sum(tb >= 5), n_ge10 = sum(tb >= 10),
             frac_in_top3 = sum(head(tb, 3)) / length(m))
}))
fwrite(summ, "results/v7/cm2_run_summary.csv")

## ---------- 3. key nodes --------------------------------------------------
key <- c("COL1A1", "COL3A1", "hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c",
         "EZH2", "hsa-miR-101", "ETS1", "NFKB1", "RELA", "SP1", "FN1")
EMT <- unique(fread("results/v7/cm2_genesets_msigdb.csv")[
  gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", gene_symbol])
EMT_net <- intersect(EMT, nm)
prot_universe <- nm[ntype[nm] %in% c("Gene", "TF")]
prior30 <- fread("results/v5/tables/Table4_prioritised_30.csv")$name

kt <- rbindlist(lapply(names(P), function(k) {
  m <- P[[k]]; cl <- m["COL1A1"]; mem <- names(m)[m == cl]
  prot_mem <- intersect(mem, prot_universe)
  hit <- length(intersect(mem, EMT_net))
  p_emt <- phyper(hit - 1, length(EMT_net), length(prot_universe) - length(EMT_net),
                  length(prot_mem), lower.tail = FALSE)
  hitP <- length(intersect(mem, prior30))
  p_pri <- phyper(hitP - 1, length(prior30), length(nm) - length(prior30),
                  length(mem), lower.tail = FALSE)
  data.table(run = k,
             COL1A1_cluster_size = length(mem),
             COL3A1_same = unname(m["COL3A1"] == cl),
             miR29a_same = unname(m["hsa-miR-29a"] == cl),
             miR29b_same = unname(m["hsa-miR-29b"] == cl),
             miR29c_same = unname(m["hsa-miR-29c"] == cl),
             FN1_same = unname(m["FN1"] == cl),
             ETS1_same = unname(m["ETS1"] == cl), NFKB1_same = unname(m["NFKB1"] == cl),
             RELA_same = unname(m["RELA"] == cl), SP1_same = unname(m["SP1"] == cl),
             EZH2_miR101_same = unname(m["EZH2"] == m["hsa-miR-101"]),
             n_gene = sum(ntype[mem] == "Gene"), n_TF = sum(ntype[mem] == "TF"),
             n_miRNA = sum(ntype[mem] == "miRNA"),
             EMT_in_cluster = hit, EMT_expected = length(prot_mem) * length(EMT_net) /
               length(prot_universe), EMT_p = p_emt,
             prior30_in_cluster = hitP, prior30_p = p_pri,
             members_head = paste(head(mem, 40), collapse = " "))
}))
fwrite(kt, "results/v7/cm2_COL1A1_cluster_tracking.csv")

## ---------- 4. pairwise adjusted Rand index -------------------------------
runs <- names(P)
ARI <- matrix(NA_real_, length(runs), length(runs), dimnames = list(runs, runs))
for (i in seq_along(runs)) for (j in seq_along(runs))
  ARI[i, j] <- adjustedRandIndex(P[[i]], P[[j]])
fwrite(data.table(run = runs, as.data.table(ARI)), "results/v7/cm2_ARI_matrix.csv")

## default-parameter battery, one run per algorithm (clusterMaker2 defaults where
## they exist, igraph defaults otherwise)
defaults <- c(MCL = "MCL_I2.5", Louvain = "Louvain_r1.00", Walktrap = "Walktrap_s4",
              Infomap = "Infomap_t10", LabelProp = "LabelProp_best",
              FastGreedy = "FastGreedy",
              cm2_MCL = "cm2_mcl_I25", cm2_Louvain = "cm2_louvain_cm2",
              cm2_FastGreedy = "cm2_fastgreedy_cm2", cm2_Infomap = "cm2_infomap_cm2_t10",
              cm2_LabelProp = "cm2_labelprop_cm2", cm2_GLay = "cm2_glay_cm2")
defaults <- defaults[defaults %in% runs]
A2 <- ARI[defaults, defaults]; dimnames(A2) <- list(names(defaults), names(defaults))
fwrite(data.table(run = names(defaults), as.data.table(A2)),
       "results/v7/cm2_ARI_defaults.csv")

## ---------- 5. leading modules of the best-modularity partitions ----------
best <- summ[order(-modularity)][1:6, run]
mods <- rbindlist(lapply(best, function(k) {
  m <- P[[k]]; qc <- q_contrib(g, m); tb <- table(m)
  dt <- data.table(run = k, cluster = as.integer(names(tb)), size = as.integer(tb),
                   q = qc[match(as.integer(names(tb)), sort(unique(m)))])
  dt[, `:=`(n_gene = sapply(cluster, function(cl) sum(ntype[names(m)[m == cl]] == "Gene")),
            n_TF   = sapply(cluster, function(cl) sum(ntype[names(m)[m == cl]] == "TF")),
            n_miRNA= sapply(cluster, function(cl) sum(ntype[names(m)[m == cl]] == "miRNA")),
            EMT    = sapply(cluster, function(cl) length(intersect(names(m)[m == cl], EMT_net))),
            prior30= sapply(cluster, function(cl) length(intersect(names(m)[m == cl], prior30))),
            has_COL1A1 = sapply(cluster, function(cl) "COL1A1" %in% names(m)[m == cl]),
            has_COL3A1 = sapply(cluster, function(cl) "COL3A1" %in% names(m)[m == cl]),
            has_miR29a = sapply(cluster, function(cl) "hsa-miR-29a" %in% names(m)[m == cl]),
            members = sapply(cluster, function(cl) paste(names(m)[m == cl], collapse = " ")))]
  dt[, EMT_p := mapply(function(cl, sz) {
        mem <- names(m)[m == cl]; pm <- intersect(mem, prot_universe)
        phyper(length(intersect(mem, EMT_net)) - 1, length(EMT_net),
               length(prot_universe) - length(EMT_net), length(pm), lower.tail = FALSE)
      }, cluster, size)]
  dt[, prior30_p := mapply(function(cl, sz)
        phyper(length(intersect(names(m)[m == cl], prior30)) - 1, length(prior30),
               length(nm) - length(prior30), sz, lower.tail = FALSE), cluster, size)]
  dt[order(-size)]
}))
fwrite(mods, "results/v7/cm2_leading_modules.csv")

## ---------- 6. consensus co-clustering over the battery -------------------
sel <- c(defaults, c(Leiden = "Leiden_mod_r1.0", LeadingEigen = "LeadingEigen",
                     Walktrap6 = "Walktrap_s6", MCL2 = "MCL_I2.0"))
sel <- sel[sel %in% runs]
stroma <- intersect(c("COL1A1","COL3A1","COL7A1","COL18A1","FN1","SPARC","POSTN","LOX",
                      "MMP2","MMP14","THBS1","THBS2","TGFB2","TGFBI","SERPINE1","ACTA2",
                      "TAGLN","VIM","CDH11","FBN1","CTGF","PDGFRB","ITGA5","ITGB1",
                      "hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"), nm)
co <- matrix(0, length(stroma), length(stroma), dimnames = list(stroma, stroma))
for (k in sel) { m <- P[[k]][stroma]; co <- co + outer(m, m, "==") }
co <- co / length(sel)
fwrite(data.table(node = stroma, as.data.table(co)),
       "results/v7/cm2_stroma_coclustering.csv")

## whole-network consensus rate for the three headline nodes
trio <- c("COL1A1", "COL3A1", "hsa-miR-29a")
pairs <- t(combn(trio, 2))
trio_dt <- rbindlist(lapply(seq_len(nrow(pairs)), function(i)
  data.table(a = pairs[i,1], b = pairs[i,2],
             n_runs = length(sel),
             n_together = sum(sapply(sel, function(k) P[[k]][pairs[i,1]] == P[[k]][pairs[i,2]])),
             runs_together = paste(names(sel)[sapply(sel, function(k)
               P[[k]][pairs[i,1]] == P[[k]][pairs[i,2]])], collapse = ";"))))
fwrite(trio_dt, "results/v7/cm2_trio_cocluster.csv")

cat("== run summary (sorted by modularity) ==\n"); print(summ[order(-modularity)][1:20])
cat("\n== ARI, default-parameter battery ==\n"); print(round(A2, 3))
cat("\n== COL1A1 cluster ==\n")
print(kt[, .(run, COL1A1_cluster_size, COL3A1_same, miR29a_same, EMT_in_cluster,
             EMT_p = signif(EMT_p, 3), prior30_in_cluster, prior30_p = signif(prior30_p, 3))])
cat("\n== trio co-clustering ==\n"); print(trio_dt[, .(a, b, n_together, n_runs)])
