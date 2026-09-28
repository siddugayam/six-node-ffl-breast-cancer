#!/usr/bin/env Rscript
# v4/08 -- (D) one tidy master table + the C3:MIR internal positive control
#             + the in-vitro readout shortlist from the miR-130a high/low contrast
suppressPackageStartupMessages({ library(data.table) })
ROOT <- "/path/to/revision"
RES  <- file.path(ROOT, "results", "v4")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

A   <- fread(file.path(RES, "gsea_all_results.csv"))
prov<- fread(file.path(RES, "ranklist_provenance.csv"))
inv <- fread(file.path(RES, "gsea_geneset_inventory.csv"))
der <- tryCatch(fread(file.path(RES, "gsea_ffl_derived_sets.csv")), error = function(e) NULL)

## --------------------------------- master tidy table ------------------------
M <- copy(A)
if (!is.null(der)) {
  d <- copy(der)
  d[, `:=`(collection_label = "FFL derived node sets", padj_global = NA_real_,
           minSize_used = 5L, maxSize_used = 600L, nPermSimple = 10000L,
           method = "fgseaMultilevel(eps=0)")]
  M <- rbind(M, d[, names(M), with = FALSE], use.names = TRUE)
}
M <- merge(M, prov[, .(ranked_list, ranking_statistic = statistic,
                       n_features_ranked = n_features, n_samples, ranking_detail = detail)],
           by = "ranked_list", all.x = TRUE)
M[, significant_FDR05 := padj < 0.05]
M[, direction := ifelse(NES > 0, "enriched at the positive/high end",
                        "enriched at the negative/low end")]
setcolorder(M, c("ranked_list", "ranking_statistic", "ranking_detail",
                 "n_features_ranked", "n_samples",
                 "collection", "collection_label", "pathway", "size",
                 "ES", "NES", "pval", "padj", "padj_global", "log2err",
                 "direction", "significant_FDR05", "leadingEdge_size", "leadingEdge",
                 "minSize_used", "maxSize_used", "nPermSimple", "method"))
setorder(M, ranked_list, collection, pval)
fwrite(M, file.path(RES, "gsea_master_table.csv"))
msg("gsea_master_table.csv: ", nrow(M), " rows (", sum(M$significant_FDR05), " FDR<0.05)")

## --------------------------------- headline terms ---------------------------
H <- M[significant_FDR05 == TRUE]
setorder(H, ranked_list, collection, pval)
head_up <- H[NES > 0, head(.SD, 10), by = .(ranked_list, collection)]
head_dn <- H[NES < 0, head(.SD, 10), by = .(ranked_list, collection)]
HD <- rbind(head_up, head_dn); setorder(HD, ranked_list, collection, -NES)
fwrite(HD, file.path(RES, "gsea_headline_terms.csv"))
msg("gsea_headline_terms.csv: ", nrow(HD), " rows")

## --------------------------------- miR positive control ---------------------
ctrl_sets <- c("TTGCACT_MIR130A_MIR301_MIR130B", "MIR130A_3P", "MIR130A_5P",
               "TGGTGCT_MIR29A_MIR29B_MIR29C", "MIR29A_3P", "MIR29A_5P",
               "MIR130B_3P", "MIR301A_3P", "MIR29B_3P", "MIR29C_3P")
C <- A[pathway %in% ctrl_sets]
C[, expectation := fifelse(grepl("130A|130B|301", pathway),
    "miR-130a family: expect NEGATIVE NES where miR-130a is high/positive",
    "miR-29a family: expect NEGATIVE NES where miR-29a is high/positive")]
setorder(C, pathway, ranked_list)
fwrite(C[, .(ranked_list, collection, pathway, size, ES, NES, pval, padj, padj_global,
             leadingEdge_size, leadingEdge, expectation)],
       file.path(RES, "gsea_mir_seed_positive_control.csv"))
cat("\n===== C3:MIR internal positive control =====\n")
print(C[pathway %in% c("TTGCACT_MIR130A_MIR301_MIR130B", "MIR130A_3P",
                       "TGGTGCT_MIR29A_MIR29B_MIR29C", "MIR29A_3P"),
        .(pathway, ranked_list, size, NES, pval, padj)])

## rank of the miR-130a set among all sets in its own collection, per list
rk <- A[collection %in% c("C3_MIR_LEGACY", "C3_MIR_MIRDB")]
rk[, rank_by_NES_ascending := frank(NES), by = .(ranked_list, collection)]
rk[, n_sets_in_collection := .N, by = .(ranked_list, collection)]
fwrite(rk[pathway %in% ctrl_sets,
          .(ranked_list, collection, pathway, NES, pval, padj,
            rank_by_NES_ascending, n_sets_in_collection)],
       file.path(RES, "gsea_mir_seed_control_rankings.csv"))
cat("\n===== rank of the miR-130a target set among all C3:MIR sets (1 = most depleted) =====\n")
print(rk[pathway %in% c("TTGCACT_MIR130A_MIR301_MIR130B", "MIR130A_3P"),
         .(ranked_list, collection, pathway, NES, padj,
           rank_by_NES_ascending, n_sets_in_collection)][order(collection, ranked_list)])

## --------------------------------- in-vitro readout shortlist ---------------
S5 <- A[ranked_list == "mir130a_high_vs_low" & padj < 0.05 &
          collection %in% c("H", "C2_CP_KEGG", "C2_CP_REACTOME", "C5_GO_BP", "C6",
                            "C3_MIR_LEGACY", "C3_MIR_MIRDB")]
ex <- S5[, .(gene = unlist(strsplit(leadingEdge, ";", fixed = TRUE))),
         by = .(collection, pathway, NES, padj)]
ex[, direction := ifelse(NES > 0, "up_in_130a_high", "down_in_130a_high")]
sh <- ex[, .(n_sets = .N, n_sets_down_in_high = sum(direction == "down_in_130a_high"),
             n_sets_up_in_high = sum(direction == "up_in_130a_high"),
             collections = paste(sort(unique(collection)), collapse = ";"),
             example_sets = paste(head(unique(pathway), 4), collapse = ";")), by = gene]
r5  <- fread(file.path(RES, "rank_mir130a_high_vs_low.csv"))
c130<- fread(file.path(RES, "rank_mir130a_corr.csv"))
tvn <- fread(file.path(RES, "rank_tumour_vs_normal.csv"))
sh[, `:=`(t_high_vs_low = r5$stat[match(gene, r5$feature)],
          logFC_high_vs_low = r5$logFC[match(gene, r5$feature)],
          FDR_high_vs_low = r5$adj.P.Val[match(gene, r5$feature)],
          rho_mir130a = c130$rho[match(gene, c130$feature)],
          t_tumour_vs_normal = tvn$stat[match(gene, tvn$feature)],
          logFC_tumour_vs_normal = tvn$logFC[match(gene, tvn$feature)])]
## annotate against the verified miR-130a target evidence
tgt <- tryCatch(fread(file.path(ROOT, "results", "v3", "mir130a_target_set.csv")), error = function(e) NULL)
if (!is.null(tgt)) {
  sh[, mir130a_target_tier := tgt$tier[match(gene, tgt$gene)]]
  sh[, mir130a_strong_evidence := gene %in% tgt$gene[tgt$strong == TRUE]]
} else { sh[, mir130a_target_tier := NA_character_]; sh[, mir130a_strong_evidence := NA] }
sh[, is_TS_anchor := gene %in% c("PTEN", "SMAD4", "TGFBR2", "DICER1", "KLF4")]
nodes <- fread(file.path(ROOT, "data", "canonical_nodes.tsv"))
sh[, is_network_node := gene %in% nodes$name]
sh[, priority := (-t_high_vs_low) * (1 + 2 * as.integer(mir130a_strong_evidence %in% TRUE)) ]
setorder(sh, -n_sets_down_in_high, priority)
fwrite(sh, file.path(RES, "gsea_mir130a_invitro_readouts.csv"))
msg("in-vitro readout shortlist: ", nrow(sh), " leading-edge genes; ",
    sum(sh$mir130a_strong_evidence %in% TRUE), " with strong miR-130a target evidence")
cat("\n===== top readouts DOWN in miR-130a-high tumours (i.e. expected UP after knockdown) =====\n")
print(head(sh[t_high_vs_low < 0][order(t_high_vs_low),
          .(gene, n_sets, t_high_vs_low, logFC_high_vs_low, rho_mir130a,
            mir130a_target_tier, is_TS_anchor, is_network_node)], 30))
cat("\n===== validated miR-130a targets appearing in leading edges =====\n")
print(sh[mir130a_strong_evidence %in% TRUE][order(t_high_vs_low),
        .(gene, n_sets, t_high_vs_low, logFC_high_vs_low, FDR_high_vs_low,
          rho_mir130a, mir130a_target_tier, is_TS_anchor)])
msg("DONE 08")
