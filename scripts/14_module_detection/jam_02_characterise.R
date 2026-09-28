#!/usr/bin/env Rscript
## jam_02_characterise.R -- what the jActiveModules search recovered:
## reproducibility across seeds, overlap with the paper's prioritised nodes and
## module, and whether the active modules are stromal or proliferative.
source("/path/to/revision/scripts/14_module_detection/jam_common.R")
sourceCpp(file.path(SCR, "jam_engine.cpp"))
suppressPackageStartupMessages(library(msigdbr))

R  <- readRDS(file.path(BASE, "cache/v7/jam/jam_runs.rds"))
g  <- R$g; mods <- R$modules
z  <- R$scores$primary$z; cal <- R$cals$sets
nodes <- V(g)$name; ntype <- setNames(V(g)$type, nodes)
prot  <- nodes[ntype != "miRNA"]
memb  <- function(s) strsplit(s, ";", fixed = TRUE)[[1]]

## ------------------------------------------------------ reference gene sets --
sets <- list()
sets$paper_MS_module   <- readRDS(file.path(BASE, "data/ffl_module_sets.rds"))$MS_MODULE
sets$prioritised_30    <- fread(file.path(BASE, "results/v5/tables/Table4_prioritised_30.csv"))$name
sets$farmer_stroma_50  <- readLines(file.path(BASE, "cache/v5/farmer/farmer2009_stromal_signature_50.txt"))
sets$CAF_scRNA_50      <- readLines(file.path(BASE, "cache/celltype/caf_signature_scrna.txt"))
sets$named_axes        <- c("COL1A1","COL3A1","ETS1","NFKB1","RELA","SP1",
                            "hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","hsa-miR-101","EZH2")
hm <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
hml <- split(unique(hm$gene_symbol), hm$gs_name)[unique(hm$gs_name)]
hml <- split(hm$gene_symbol, hm$gs_name); hml <- lapply(hml, unique)
sets$HALLMARK_EMT        <- hml[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]]
sets$HALLMARK_E2F        <- hml[["HALLMARK_E2F_TARGETS"]]
sets$HALLMARK_G2M        <- hml[["HALLMARK_G2M_CHECKPOINT"]]
sets$HALLMARK_MYC_V1     <- hml[["HALLMARK_MYC_TARGETS_V1"]]
sets$HALLMARK_ANGIOGENESIS <- hml[["HALLMARK_ANGIOGENESIS"]]
sets$HALLMARK_TGF_BETA   <- hml[["HALLMARK_TGF_BETA_SIGNALING"]]
stroma_sets <- c("farmer_stroma_50","CAF_scRNA_50","HALLMARK_EMT","HALLMARK_ANGIOGENESIS","HALLMARK_TGF_BETA")
prolif_sets <- c("HALLMARK_E2F","HALLMARK_G2M","HALLMARK_MYC_V1")

hyper <- function(mod, set, universe) {
  mod <- intersect(mod, universe); set <- intersect(set, universe)
  q <- length(intersect(mod, set))
  data.table(n_module = length(mod), n_set = length(set), n_universe = length(universe),
             overlap = q, expected = length(mod) * length(set) / length(universe),
             fold = (q / length(mod)) / (length(set) / length(universe)),
             p_hyper = phyper(q - 1, length(set), length(universe) - length(set), length(mod),
                              lower.tail = FALSE),
             members = paste(sort(intersect(mod, set)), collapse = ";"))
}

## ------------------------------------------------- module-level annotation --
ann <- rbindlist(lapply(seq_len(nrow(mods)), function(i) {
  m <- memb(mods$members[i]); ty <- ntype[m]
  data.table(mods[i, .(variant, seed, rank, size, score)],
    n_gene = sum(ty == "Gene"), n_TF = sum(ty == "TF"), n_miRNA = sum(ty == "miRNA"),
    mean_z = mean(z[m]), mean_logFC = mean(R$scores$primary$logFC[match(m, nodes)], na.rm = TRUE),
    frac_up = mean(R$scores$primary$logFC[match(m, nodes)] > 0, na.rm = TRUE),
    has_COL1A1 = "COL1A1" %in% m, has_COL3A1 = "COL3A1" %in% m,
    n_MS_module = length(intersect(m, sets$paper_MS_module)),
    n_prior30 = length(intersect(m, sets$prioritised_30)),
    p_prior30 = hyper(m, sets$prioritised_30, nodes)$p_hyper,
    fold_prior30 = hyper(m, sets$prioritised_30, nodes)$fold,
    n_farmer = length(intersect(m, sets$farmer_stroma_50)),
    p_farmer = hyper(m, sets$farmer_stroma_50, prot)$p_hyper,
    n_E2F = length(intersect(m, sets$HALLMARK_E2F)),
    p_E2F = hyper(m, sets$HALLMARK_E2F, prot)$p_hyper)
}))
fwrite(ann, file.path(OUT, "jam_modules_annotated.csv"))

## ------------------------------------------- reproducibility across seeds --
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
rep_tab <- rbindlist(lapply(unique(mods$variant), function(v) {
  rbindlist(lapply(sort(unique(mods[variant == v]$rank)), function(rk) {
    L <- lapply(mods[variant == v & rank == rk]$members, memb)
    if (length(L) < 2) return(NULL)
    J <- outer(seq_along(L), seq_along(L), Vectorize(function(i, j) jac(L[[i]], L[[j]])))
    pj <- J[upper.tri(J)]
    tab <- table(unlist(L)); nseed <- length(L)
    data.table(variant = v, rank = rk, n_seeds = nseed,
               size_min = min(lengths(L)), size_med = median(lengths(L)), size_max = max(lengths(L)),
               jaccard_mean = mean(pj), jaccard_min = min(pj), jaccard_max = max(pj),
               n_union = length(tab), n_core_all_seeds = sum(tab == nseed),
               frac_core = sum(tab == nseed) / length(tab))
  }))
}))
fwrite(rep_tab, file.path(OUT, "jam_reproducibility.csv"))

## node-level selection frequency in the rank-1 module of the 10 primary seeds
L1 <- lapply(mods[variant == "primary" & rank == 1]$members, memb)
freq <- data.table(node = nodes, type = ntype[nodes], z = z[nodes],
                   logFC = R$scores$primary$logFC,
                   n_seeds_in_top_module = vapply(nodes, function(x) sum(vapply(L1, function(l) x %in% l, TRUE)), 0L))
freq[, freq_top_module := n_seeds_in_top_module / length(L1)]
freq[, in_prioritised_30 := node %in% sets$prioritised_30]
freq[, in_paper_MS_module := node %in% sets$paper_MS_module]
freq[, in_farmer_50 := node %in% sets$farmer_stroma_50]
fwrite(freq[order(-freq_top_module, -z)], file.path(OUT, "jam_node_selection_frequency.csv"))

## ---------------------------------------------------- enrichment of modules --
## consensus module = nodes in the rank-1 module of every primary seed
consensus <- freq[n_seeds_in_top_module == length(L1)]$node
targets <- list(consensus_top_module = consensus,
                seed1_top_module = L1[[1]],
                union_top_module = unique(unlist(L1)))
for (rk in 2:5) targets[[paste0("primary_rank", rk, "_seed1")]] <-
  memb(mods[variant == "primary" & rank == rk & seed == 1]$members)

enr <- rbindlist(lapply(names(targets), function(tn) {
  m <- targets[[tn]]
  rbindlist(lapply(names(sets), function(sn) {
    uni <- if (sn %in% c("paper_MS_module", "prioritised_30", "named_axes")) nodes else prot
    cbind(data.table(module = tn, gene_set = sn), hyper(m, sets[[sn]], uni))
  }))
}))
fwrite(enr, file.path(OUT, "jam_module_enrichment.csv"))

## full hallmark ORA for the consensus module (protein-coding universe)
ora <- rbindlist(lapply(names(hml), function(h)
  cbind(data.table(gene_set = h), hyper(consensus, hml[[h]], prot))))
ora[, q := p.adjust(p_hyper, "BH")]
fwrite(ora[order(p_hyper)], file.path(OUT, "jam_hallmark_ORA_consensus.csv"))

## ------------------------------------- how the paper's own sets would score --
ref <- list(paper_MS_module_10 = sets$paper_MS_module,
            prioritised_30 = sets$prioritised_30,
            collagen_axis_5 = c("COL1A1","COL3A1","hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"),
            TF_collagen_arm = c("ETS1","NFKB1","RELA","SP1","COL1A1","COL3A1"),
            mir101_EZH2 = c("hsa-miR-101","EZH2"))
conn_null <- readRDS(file.path(BASE, "cache/v7/jam/jam_runs.rds"))$cals$connected
refsc <- rbindlist(lapply(names(ref), function(nn) {
  m <- intersect(ref[[nn]], nodes); k <- length(m)
  sub <- induced_subgraph(g, m)
  data.table(set = nn, k = k, n_components = components(sub)$no,
             connected = components(sub)$no == 1,
             mean_z = mean(z[m]),
             agg_z = sum(z[m]) / sqrt(k),
             score_sets_null = jam_score_set(m, z, cal),
             score_connected_null = jam_score_set(m, z, conn_null),
             best_jam_score = max(mods[variant == "primary" & rank == 1]$score))
}))
fwrite(refsc, file.path(OUT, "jam_paper_sets_scored.csv"))
print(refsc)
print(rep_tab[variant == "primary"])
print(enr[module == "consensus_top_module"])
jam_log("consensus module size ", length(consensus))

## ------------------------- is the module anything more than a z threshold? --
## If the active module is simply "the k most differentially expressed nodes",
## the search has recovered differential expression, not network structure.
k <- length(consensus)
topk <- names(sort(z, decreasing = TRUE))[1:k]
sub_topk <- induced_subgraph(g, topk)
zt <- sort(z, decreasing = TRUE)[k]
thr <- rbindlist(list(
  data.table(set = "consensus_top_module", k = k,
             min_z = min(z[consensus]), median_z = median(z[consensus]),
             n_components = components(induced_subgraph(g, consensus))$no,
             score = jam_score_set(consensus, z, cal),
             jaccard_with_topk_by_z = jac(consensus, topk)),
  data.table(set = "top_k_nodes_by_z", k = k,
             min_z = min(z[topk]), median_z = median(z[topk]),
             n_components = components(sub_topk)$no,
             score = jam_score_set(topk, z, cal),
             jaccard_with_topk_by_z = 1)))
fwrite(thr, file.path(OUT, "jam_module_vs_z_threshold.csv"))
print(thr)

## AUC of node z for membership of the rank-1 module (union over primary seeds)
u <- unique(unlist(L1)); lab <- nodes %in% u
r <- rank(z[nodes]); n1 <- sum(lab); n0 <- sum(!lab)
auc <- (sum(r[lab]) - n1 * (n1 + 1) / 2) / (n1 * n0)
fwrite(data.table(metric = "AUC_z_predicts_top_module_membership", value = auc,
                  n_in = n1, n_out = n0),
       file.path(OUT, "jam_membership_auc.csv"))
jam_log("AUC of z for top-module membership: ", round(auc, 3))

## ------------------------------- per-variant consensus modules and their -----
## stromal / proliferative character
cons_by_var <- lapply(unique(mods$variant), function(v) {
  L <- lapply(mods[variant == v & rank == 1]$members, memb)
  tb <- table(unlist(L)); names(tb)[tb == length(L)]
})
names(cons_by_var) <- unique(mods$variant)
enr2 <- rbindlist(lapply(names(cons_by_var), function(v) {
  m <- cons_by_var[[v]]
  rbindlist(lapply(names(sets), function(sn) {
    uni <- if (sn %in% c("paper_MS_module", "prioritised_30", "named_axes")) nodes else prot
    cbind(data.table(variant = v, gene_set = sn), hyper(m, sets[[sn]], uni))
  }))
}))
enr2[, class := fifelse(gene_set %in% stroma_sets, "stroma",
                 fifelse(gene_set %in% prolif_sets, "proliferation", "paper"))]
fwrite(enr2, file.path(OUT, "jam_variant_consensus_enrichment.csv"))

## stroma vs proliferation summary per variant: mean log2 fold-enrichment
summ <- enr2[class %in% c("stroma", "proliferation"),
             .(mean_fold = mean(fold, na.rm = TRUE), min_p = min(p_hyper)), by = .(variant, class)]
fwrite(dcast(summ, variant ~ class, value.var = c("mean_fold", "min_p")),
       file.path(OUT, "jam_stroma_vs_proliferation.csv"))
print(dcast(summ, variant ~ class, value.var = c("mean_fold", "min_p")))

## which paper nodes are in / out of the primary consensus module
pn <- data.table(node = union(sets$prioritised_30, sets$paper_MS_module))
pn[, `:=`(type = ntype[node], z = z[node],
          in_consensus = node %in% consensus,
          n_seeds = freq$n_seeds_in_top_module[match(node, freq$node)],
          in_prioritised_30 = node %in% sets$prioritised_30,
          in_paper_MS_module = node %in% sets$paper_MS_module)]
fwrite(pn[order(-z)], file.path(OUT, "jam_paper_nodes_membership.csv"))
print(pn[order(-z)][in_paper_MS_module == TRUE])
jam_log("prioritised-30 in consensus: ", sum(pn$in_prioritised_30 & pn$in_consensus), "/30; ",
        "MS-module in consensus: ", sum(pn$in_paper_MS_module & pn$in_consensus), "/10")
