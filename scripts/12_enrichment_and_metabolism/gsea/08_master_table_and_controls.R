#!/usr/bin/env Rscript
# v4/08 -- (D) one tidy master table + the C3:MIR internal positive control
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
ctrl_sets <- c("TGGTGCT_MIR29A_MIR29B_MIR29C", "MIR29A_3P", "MIR29A_5P",
               "MIR29B_3P", "MIR29C_3P")
C <- A[pathway %in% ctrl_sets]
C[, expectation := "miR-29a family: expect NEGATIVE NES where miR-29a is high/positive"]
setorder(C, pathway, ranked_list)
fwrite(C[, .(ranked_list, collection, pathway, size, ES, NES, pval, padj, padj_global,
             leadingEdge_size, leadingEdge, expectation)],
       file.path(RES, "gsea_mir_seed_positive_control.csv"))
cat("\n===== C3:MIR internal positive control =====\n")
print(C[pathway %in% c("TGGTGCT_MIR29A_MIR29B_MIR29C", "MIR29A_3P"),
        .(pathway, ranked_list, size, NES, pval, padj)])

## rank of each control set among all sets in its own collection, per list
rk <- A[collection %in% c("C3_MIR_LEGACY", "C3_MIR_MIRDB")]
rk[, rank_by_NES_ascending := frank(NES), by = .(ranked_list, collection)]
rk[, n_sets_in_collection := .N, by = .(ranked_list, collection)]
fwrite(rk[pathway %in% ctrl_sets,
          .(ranked_list, collection, pathway, NES, pval, padj,
            rank_by_NES_ascending, n_sets_in_collection)],
       file.path(RES, "gsea_mir_seed_control_rankings.csv"))
msg("DONE 08")
