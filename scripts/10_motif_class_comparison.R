#!/usr/bin/env Rscript
# =============================================================================
# 10_motif_class_comparison.R
#
# Produces results/motif_class_comparison.csv - a long tidy table, one row per
# comparison statistic, with columns:
#   block, metric, class_a, class_b, value, statistic, p_value, p_adjust, n_a, n_b, note
#
# Blocks:
#   GENESET_JACCARD   pairwise Jaccard of the protein-coding gene sets, + permutation p
#   NODESET_JACCARD   pairwise Jaccard of the full node sets (incl. miRNAs)
#   TERMSET_JACCARD   pairwise Jaccard of the BH-significant enriched term sets
#   COMPOSITION       node-type composition per class + chi-square across classes
#   TOPOLOGY_CLASS    per-class node topology, Kruskal-Wallis across classes and
#                     per-class vs rest-of-network Mann-Whitney
#   TOPOLOGY_EXCLUSIVE same on the class-EXCLUSIVE nodes only (non-overlapping)
#   TOPOLOGY_BYTYPE   in/out degree by node type, Kruskal-Wallis + pairwise Wilcoxon
#   COMPARECLUSTER    clusterProfiler::compareCluster summary per class
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(igraph); library(clusterProfiler); library(org.Hs.eg.db)
  library(AnnotationDbi)
})
REV <- "/path/to/revision"
LOGF <- file.path(REV, "logs", "10_motif_class_comparison.log")
lg <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOGF, append = TRUE) }
cat("", file = LOGF); set.seed(1234)

CLS <- c("3-miR", "3-TF", "3-Comp", "4-node", "5-node", "6-node")
OUT <- list(); add <- function(...) OUT[[length(OUT) + 1]] <<- data.table(...)

ns <- fread(file.path(REV, "results/motif_node_sets.tsv"))
nodes_all <- fread(file.path(REV, "data/canonical_nodes.tsv"))
ed <- fread(file.path(REV, "results/graph_used_edges.tsv"))
lg("edges in graph used: ", nrow(ed), "  nodes: ", nrow(nodes_all))

nodeset <- lapply(CLS, function(c) sort(unique(ns[motif_set == c, node])));  names(nodeset) <- CLS
geneset <- lapply(CLS, function(c) sort(unique(ns[motif_set == c & type %in% c("TF","Gene"), node])))
names(geneset) <- CLS

jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))

# ------------------------------------------------------------- 1. set Jaccards
pool_nodes <- nodes_all$name
pool_genes <- nodes_all[type %in% c("TF", "Gene"), name]

perm_p <- function(a, b, pool, obs, B = 2000) {
  # null: two sets of the same sizes drawn uniformly at random from pool
  na <- length(a); nb <- length(b)
  nullj <- replicate(B, jac(sample(pool, na), sample(pool, nb)))
  # two-sided empirical p
  p <- (1 + min(sum(nullj >= obs), sum(nullj <= obs)) * 2) / (B + 1)
  c(p = min(1, p), null_mean = mean(nullj))
}

for (i in 1:(length(CLS) - 1)) for (j in (i + 1):length(CLS)) {
  a <- geneset[[CLS[i]]]; b <- geneset[[CLS[j]]]; o <- jac(a, b)
  pp <- perm_p(a, b, pool_genes, o)
  add(block = "GENESET_JACCARD", metric = "jaccard_protein_coding", class_a = CLS[i],
      class_b = CLS[j], value = o, statistic = pp["null_mean"], p_value = pp["p"],
      n_a = length(a), n_b = length(b),
      note = paste0("shared=", length(intersect(a, b)),
                    "; null mean Jaccard for random sets of same size=",
                    round(pp["null_mean"], 3)))
  a <- nodeset[[CLS[i]]]; b <- nodeset[[CLS[j]]]; o <- jac(a, b)
  pp <- perm_p(a, b, pool_nodes, o)
  add(block = "NODESET_JACCARD", metric = "jaccard_all_nodes", class_a = CLS[i],
      class_b = CLS[j], value = o, statistic = pp["null_mean"], p_value = pp["p"],
      n_a = length(a), n_b = length(b),
      note = paste0("shared=", length(intersect(a, b))))
}

# ------------------------------------------------------------- 2. composition
comp <- dcast(ns[motif_set %in% CLS], motif_set ~ type, fun.aggregate = length, value.var = "node")
lg("node-type composition:"); lg(paste(capture.output(print(comp)), collapse = "\n"))
M <- as.matrix(comp[, -1]); rownames(M) <- comp$motif_set
cs <- suppressWarnings(chisq.test(M))
add(block = "COMPOSITION", metric = "chisq_nodetype_across_classes", class_a = "ALL",
    class_b = NA_character_, value = as.numeric(cs$statistic),
    statistic = as.numeric(cs$parameter), p_value = cs$p.value,
    n_a = sum(M), n_b = NA_integer_,
    note = "chi-square of TF/Gene/miRNA counts across the 6 classes; df in statistic column")
for (k in seq_len(nrow(comp))) {
  add(block = "COMPOSITION", metric = "n_TF", class_a = comp$motif_set[k], class_b = NA_character_,
      value = comp$TF[k], n_a = sum(M[k, ]), note = "count")
  add(block = "COMPOSITION", metric = "n_Gene", class_a = comp$motif_set[k], class_b = NA_character_,
      value = comp$Gene[k], n_a = sum(M[k, ]), note = "count")
  add(block = "COMPOSITION", metric = "n_miRNA", class_a = comp$motif_set[k], class_b = NA_character_,
      value = comp$miRNA[k], n_a = sum(M[k, ]), note = "count")
}

# ------------------------------------------------------------- 3. topology
g <- graph_from_data_frame(ed[, .(source, target)], directed = TRUE,
                           vertices = nodes_all[, .(name, type)])
TOP <- data.table(name = V(g)$name, type = V(g)$type,
                  degree = degree(g, mode = "all"),
                  indegree = degree(g, mode = "in"),
                  outdegree = degree(g, mode = "out"),
                  betweenness = betweenness(g, directed = TRUE, normalized = TRUE),
                  clustering = transitivity(g, type = "local", isolates = "zero"),
                  closeness = closeness(g, mode = "out", normalized = TRUE))
TOP[is.na(clustering), clustering := 0]; TOP[is.na(closeness), closeness := 0]
fwrite(TOP, file.path(REV, "results/network_node_topology.csv"))
lg("wrote network_node_topology.csv rows=", nrow(TOP))

METRICS <- c("degree", "indegree", "outdegree", "betweenness", "clustering", "closeness")

memb <- rbindlist(lapply(CLS, function(c) data.table(motif_set = c, name = nodeset[[c]])))
TT <- merge(memb, TOP, by = "name")

for (m in METRICS) {
  kw <- kruskal.test(TT[[m]] ~ factor(TT$motif_set))
  add(block = "TOPOLOGY_CLASS", metric = paste0("kruskal_", m), class_a = "ALL_6_CLASSES",
      class_b = NA_character_, value = as.numeric(kw$statistic),
      statistic = as.numeric(kw$parameter), p_value = kw$p.value, n_a = nrow(TT),
      note = "Kruskal-Wallis across the 6 class node sets. WARNING: the sets overlap heavily, so observations are not independent; interpret as a descriptive contrast only.")
  for (c in CLS) {
    inn <- TOP[name %in% nodeset[[c]]][[m]]; outn <- TOP[!name %in% nodeset[[c]]][[m]]
    wt <- if (length(outn) > 2) suppressWarnings(wilcox.test(inn, outn)) else NULL
    add(block = "TOPOLOGY_CLASS", metric = paste0("median_", m), class_a = c,
        class_b = "rest_of_network", value = median(inn),
        statistic = if (is.null(wt)) NA_real_ else as.numeric(wt$statistic),
        p_value = if (is.null(wt)) NA_real_ else wt$p.value,
        n_a = length(inn), n_b = length(outn),
        note = paste0("median in class vs median outside = ", round(median(inn), 4), " vs ",
                      ifelse(length(outn) == 0, NA, round(median(outn), 4))))
  }
}
# BH across the per-class Mann-Whitney tests
tmp <- rbindlist(OUT, fill = TRUE)

# ---- exclusive nodes (non-overlapping subsets)
excl <- lapply(CLS, function(c) setdiff(nodeset[[c]], unlist(nodeset[setdiff(CLS, c)])))
names(excl) <- CLS
lg("class-exclusive node counts: ", paste(CLS, sapply(excl, length), sep = "=", collapse = "  "))
ex <- rbindlist(lapply(CLS, function(c) if (length(excl[[c]])) data.table(motif_set = c, name = excl[[c]]) else NULL))
if (nrow(ex) > 0) {
  EX <- merge(ex, TOP, by = "name")
  for (m in METRICS) {
    if (length(unique(EX$motif_set)) > 1 && nrow(EX) > 3) {
      kw <- kruskal.test(EX[[m]] ~ factor(EX$motif_set))
      add(block = "TOPOLOGY_EXCLUSIVE", metric = paste0("kruskal_", m),
          class_a = paste(unique(EX$motif_set), collapse = "+"), class_b = NA_character_,
          value = as.numeric(kw$statistic), statistic = as.numeric(kw$parameter),
          p_value = kw$p.value, n_a = nrow(EX),
          note = "Kruskal-Wallis on CLASS-EXCLUSIVE nodes only (non-overlapping, independent)")
    }
    for (c in unique(EX$motif_set))
      add(block = "TOPOLOGY_EXCLUSIVE", metric = paste0("median_", m), class_a = c,
          class_b = NA_character_, value = median(EX[motif_set == c][[m]]),
          n_a = sum(EX$motif_set == c), note = "class-exclusive nodes only")
  }
}
for (c in CLS)
  add(block = "TOPOLOGY_EXCLUSIVE", metric = "n_exclusive_nodes", class_a = c,
      class_b = NA_character_, value = length(excl[[c]]), n_a = length(nodeset[[c]]),
      note = "nodes present in this class and in NO other class")

# ---- in/out degree by node type
for (m in c("indegree", "outdegree", "degree", "betweenness", "clustering")) {
  kw <- kruskal.test(TOP[[m]] ~ factor(TOP$type))
  add(block = "TOPOLOGY_BYTYPE", metric = paste0("kruskal_", m, "_by_nodetype"),
      class_a = "TF|Gene|miRNA", class_b = NA_character_, value = as.numeric(kw$statistic),
      statistic = as.numeric(kw$parameter), p_value = kw$p.value, n_a = nrow(TOP),
      note = "whole network, node types are disjoint so this test is valid")
  pw <- suppressWarnings(pairwise.wilcox.test(TOP[[m]], TOP$type, p.adjust.method = "BH"))
  pm <- pw$p.value
  for (r in rownames(pm)) for (cc in colnames(pm)) if (!is.na(pm[r, cc]))
    add(block = "TOPOLOGY_BYTYPE", metric = paste0("wilcox_BH_", m), class_a = r, class_b = cc,
        value = NA_real_, p_adjust = pm[r, cc],
        n_a = sum(TOP$type == r), n_b = sum(TOP$type == cc), note = "pairwise Wilcoxon, BH")
  for (ty in unique(TOP$type))
    add(block = "TOPOLOGY_BYTYPE", metric = paste0("median_", m), class_a = ty,
        class_b = NA_character_, value = median(TOP[type == ty][[m]]),
        n_a = sum(TOP$type == ty), note = "whole network")
}

# ------------------------------------------------------------- 4. compareCluster
sym2ent <- function(x) {
  s <- suppressWarnings(suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(x), keytype = "SYMBOL", columns = "ENTREZID")))
  unique(s$ENTREZID[!is.na(s$ENTREZID)])
}
de <- fread(file.path(REV, "results/BRCA_DEX_genes.csv"))
uniA <- sym2ent(unique(c(de$feature[grepl("^[A-Za-z]", de$feature)], pool_genes)))
gl <- lapply(geneset, sym2ent)
lg("compareCluster input sizes: ", paste(names(gl), sapply(gl, length), sep = "=", collapse = " "))

cc_res <- list()
for (dbn in c("GO_BP", "KEGG")) {
  cc <- tryCatch({
    if (dbn == "GO_BP")
      compareCluster(geneClusters = gl, fun = "enrichGO", OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID", ont = "BP", universe = uniA,
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                     minGSSize = 10, maxGSSize = 500, readable = TRUE)
    else
      compareCluster(geneClusters = gl, fun = "enrichKEGG", organism = "hsa",
                     universe = uniA, pAdjustMethod = "BH", pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500)
  }, error = function(e) { lg("compareCluster ", dbn, " FAILED: ", conditionMessage(e)); NULL })
  if (is.null(cc)) next
  d <- as.data.table(cc)
  d[, ontology := dbn]
  cc_res[[dbn]] <- d
  lg("compareCluster ", dbn, ": ", nrow(d), " significant rows over ",
     length(unique(d$Cluster)), " clusters")
}
if (length(cc_res)) {
  CCD <- rbindlist(cc_res, fill = TRUE)
  fwrite(CCD, file.path(REV, "results/compareCluster_motif_classes.csv"))
  lg("wrote compareCluster_motif_classes.csv rows=", nrow(CCD))
  for (dbn in names(cc_res)) {
    d <- cc_res[[dbn]]
    tset <- split(d$ID, d$Cluster)
    for (i in 1:(length(CLS) - 1)) for (j in (i + 1):length(CLS)) {
      a <- tset[[CLS[i]]]; b <- tset[[CLS[j]]]
      a <- if (is.null(a)) character(0) else a; b <- if (is.null(b)) character(0) else b
      v <- if (length(union(a, b)) == 0) NA_real_ else jac(a, b)
      add(block = "TERMSET_JACCARD", metric = paste0("jaccard_sigterms_", dbn),
          class_a = CLS[i], class_b = CLS[j], value = v, n_a = length(a), n_b = length(b),
          note = paste0("BH<0.05 & q<0.2 terms; shared=", length(intersect(a, b))))
    }
    for (c in CLS)
      add(block = "COMPARECLUSTER", metric = paste0("n_sig_terms_", dbn), class_a = c,
          class_b = NA_character_, value = length(tset[[c]]), n_a = length(gl[[c]]),
          note = "significant terms in compareCluster")
  }
}

R <- rbindlist(OUT, fill = TRUE)
for (col in c("p_adjust")) if (!col %in% names(R)) R[[col]] <- NA_real_
R[block %in% c("TOPOLOGY_CLASS") & !is.na(p_value) & grepl("^median_", metric),
  p_adjust := p.adjust(p_value, "BH"), by = metric]
setcolorder(R, c("block", "metric", "class_a", "class_b", "value", "statistic",
                 "p_value", "p_adjust", "n_a", "n_b", "note"))
fwrite(R, file.path(REV, "results/motif_class_comparison.csv"))
lg("wrote motif_class_comparison.csv rows=", nrow(R))
lg("DONE")
