#!/usr/bin/env Rscript
# v4/04 -- (A) Leading-edge overlap analysis across the ranked lists
suppressPackageStartupMessages({ library(data.table) })
ROOT <- "/path/to/revision"
RES  <- file.path(ROOT, "results", "v4")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

A <- fread(file.path(RES, "gsea_all_results.csv"))
S <- A[padj < 0.05 & collection != "FFL_CLASS"]
msg("significant rows used (excluding the custom FFL_CLASS collection): ", nrow(S))

## explode leading edges
ex <- S[, .(gene = unlist(strsplit(leadingEdge, ";", fixed = TRUE))),
        by = .(ranked_list, collection, collection_label, pathway, NES, padj)]
msg("leading-edge gene-in-set records: ", nrow(ex), " ; distinct genes: ", uniqueN(ex$gene))
ex[, direction := ifelse(NES > 0, "up", "down")]

## ---- per gene x ranked list -------------------------------------------------
per_list <- ex[, .(n_sets = .N,
                   n_collections = uniqueN(collection),
                   n_up = sum(direction == "up"),
                   n_down = sum(direction == "down")),
               by = .(gene, ranked_list)]
fwrite(per_list, file.path(RES, "gsea_leadingedge_gene_by_list.csv"))

## ---- cross-list summary -----------------------------------------------------
glob <- ex[, .(n_sets_total = .N,
               n_ranked_lists = uniqueN(ranked_list),
               n_collections = uniqueN(collection),
               ranked_lists = paste(sort(unique(ranked_list)), collapse = ";"),
               n_up = sum(direction == "up"), n_down = sum(direction == "down")),
           by = gene]
## per-list set counts as wide columns
w <- dcast(per_list, gene ~ ranked_list, value.var = "n_sets", fill = 0L)
setnames(w, setdiff(names(w), "gene"), paste0("nsets_", setdiff(names(w), "gene")))
glob <- merge(glob, w, by = "gene", all.x = TRUE)

## annotate: network node? miR-130a strong target? DE stats?
nodes <- fread(file.path(ROOT, "data", "canonical_nodes.tsv"))
glob[, is_network_node := gene %in% nodes$name]
glob[, node_type := nodes$type[match(gene, nodes$name)]]
mn <- fread(file.path(ROOT, "results", "motif_node_sets.tsv"))
glob[, in_6node_class := gene %in% mn$node[mn$motif_set == "6-node"]]
strong5 <- c("PTEN", "SMAD4", "TGFBR2", "DICER1", "KLF4")
glob[, is_mir130a_strong_TS_target := gene %in% strong5]
tvn <- fread(file.path(RES, "rank_tumour_vs_normal.csv"))
glob[, t_tumour_vs_normal := tvn$stat[match(gene, tvn$feature)]]
glob[, logFC_tumour_vs_normal := tvn$logFC[match(gene, tvn$feature)]]
c130 <- fread(file.path(RES, "rank_mir130a_corr.csv"))
glob[, rho_mir130a := c130$rho[match(gene, c130$feature)]]
c29 <- fread(file.path(RES, "rank_mir29a_corr.csv"))
glob[, rho_mir29a := c29$rho[match(gene, c29$feature)]]
r5 <- fread(file.path(RES, "rank_mir130a_high_vs_low.csv"))
glob[, t_mir130a_high_vs_low := r5$stat[match(gene, r5$feature)]]

setorder(glob, -n_ranked_lists, -n_sets_total)
fwrite(glob, file.path(RES, "gsea_leadingedge_gene_frequency.csv"))
msg("genes in >=1 leading edge: ", nrow(glob),
    " | in all 5 lists: ", sum(glob$n_ranked_lists == 5),
    " | in >=4 lists: ", sum(glob$n_ranked_lists >= 4))
print(head(glob[, .(gene, n_ranked_lists, n_sets_total, is_network_node,
                    t_tumour_vs_normal, rho_mir130a, t_mir130a_high_vs_low)], 25))

## ---- pairwise Jaccard of leading-edge unions across ranked lists ------------
un <- split(ex$gene, ex$ranked_list); un <- lapply(un, unique)
ls_names <- names(un)
J <- data.table(expand.grid(list_a = ls_names, list_b = ls_names, stringsAsFactors = FALSE))
J[, n_a := lengths(un)[list_a]]; J[, n_b := lengths(un)[list_b]]
J[, n_shared := mapply(function(a, b) length(intersect(un[[a]], un[[b]])), list_a, list_b)]
J[, jaccard := n_shared / (n_a + n_b - n_shared)]
fwrite(J, file.path(RES, "gsea_leadingedge_list_overlap.csv"))
print(dcast(J, list_a ~ list_b, value.var = "jaccard"))

## ---- pathway-level overlap: terms significant in multiple lists -------------
pl <- S[, .(n_lists = uniqueN(ranked_list),
            lists = paste(sort(unique(ranked_list)), collapse = ";"),
            NES_str = paste(paste0(ranked_list, "=", round(NES, 3)), collapse = ";")),
        by = .(collection, collection_label, pathway)]
setorder(pl, -n_lists, pathway)
fwrite(pl, file.path(RES, "gsea_pathway_recurrence_across_lists.csv"))
msg("pathways significant in all 5 lists: ", sum(pl$n_lists == 5),
    " | in >=4: ", sum(pl$n_lists >= 4))

## ---- hub-gene focus: leading-edge genes that are network nodes -------------
hub <- glob[is_network_node == TRUE][order(-n_ranked_lists, -n_sets_total)]
fwrite(hub, file.path(RES, "gsea_leadingedge_network_nodes.csv"))
msg("network nodes appearing in a leading edge: ", nrow(hub), " of ",
    sum(nodes$type %in% c("Gene", "TF")), " protein-coding nodes")
msg("DONE 04")
