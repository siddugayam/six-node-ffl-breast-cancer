#!/usr/bin/env Rscript
# v4/10 -- (D) one consolidated tidy results table across every analysis in v4,
#              plus the figure that sets the FFL classes against both backgrounds.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
ROOT  <- "/path/to/revision"
RES   <- file.path(ROOT, "results", "v4")
FIG   <- file.path(ROOT, "figures", "v4"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
CACHE <- file.path(ROOT, "cache", "v4")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
theme_set(theme_bw(base_size = 9))

## ============================ consolidated tidy table =======================
mt  <- fread(file.path(RES, "gsea_master_table.csv"))       # 3 primary lists x 14 collections + FFL derived
prov<- fread(file.path(RES, "ranklist_provenance.csv"))

mt[, analysis := ifelse(collection == "FFL_DERIVED", "FFL derived node sets",
                 ifelse(collection == "FFL_CLASS", "FFL motif-class node sets (analysis C)",
                        "primary GSEA (analyses 1-3)"))]
ALL <- copy(mt)
setcolorder(ALL, c("analysis", "ranked_list", "ranking_statistic", "ranking_detail",
                   "n_features_ranked", "n_samples", "collection", "collection_label",
                   "pathway", "size", "ES", "NES", "pval", "padj", "padj_global",
                   "log2err", "direction", "significant_FDR05",
                   "leadingEdge_size", "leadingEdge"))
setorder(ALL, analysis, ranked_list, collection, pval)
fwrite(ALL, file.path(RES, "gsea_MASTER_TABLE_ALL.csv"))
msg("gsea_MASTER_TABLE_ALL.csv: ", nrow(ALL), " rows | FDR<0.05: ", sum(ALL$significant_FDR05))

## a compact, human-readable headline table (no leading-edge strings)
HEAD <- ALL[significant_FDR05 == TRUE,
            .(analysis, ranked_list, collection, collection_label, pathway, size,
              NES, pval, padj, leadingEdge_size)]
setorder(HEAD, ranked_list, collection, pval)
top <- rbind(HEAD[NES > 0, head(.SD, 15), by = .(ranked_list, collection)],
             HEAD[NES < 0, head(.SD, 15), by = .(ranked_list, collection)])
setorder(top, ranked_list, collection, -NES)
fwrite(top, file.path(RES, "gsea_MASTER_TABLE_TOP15.csv"))
msg("gsea_MASTER_TABLE_TOP15.csv: ", nrow(top), " rows")

## per ranked-list x collection counts
cnt <- ALL[, .(n_tested = .N, n_sig_FDR05 = sum(significant_FDR05),
               n_sig_pos = sum(significant_FDR05 & NES > 0),
               n_sig_neg = sum(significant_FDR05 & NES < 0),
               top_pos = pathway[which.max(ifelse(significant_FDR05, NES, NA))],
               top_neg = pathway[which.min(ifelse(significant_FDR05, NES, NA))]),
           by = .(analysis, ranked_list, collection, collection_label)]
setorder(cnt, analysis, ranked_list, collection)
fwrite(cnt, file.path(RES, "gsea_MASTER_COUNTS.csv"))
print(cnt[analysis == "primary GSEA (analyses 1-3)"], nrows = 80)

## ============================ figure: FFL class vs both backgrounds =========
nb <- fread(file.path(RES, "gsea_ffl_class_network_background_test.csv"))
tw <- fread(file.path(RES, "gsea_ffl_class_transcriptome_background_test.csv"))
KEEP <- c("3-miR", "3-TF", "3-Comp", "4-node", "5-node", "6-node",
          "4node_not_3node", "5node_not_3node", "6node_not_3node",
          "higherorder_not_3node", "ALL_NETWORK_NODES")
b <- merge(nb[set %in% KEEP, .(set, n_in_bg, network_bg = bg_mean_absT, p_network = p_vs_network_bg)],
           tw[set %in% KEEP, .(set, n, obs_mean_absT, transcriptome_bg = transcriptome_mean_absT,
                               p_transcriptome = p_vs_transcriptome)],
           by = "set")
## ALL_NETWORK_NODES *is* the within-network background, so it has no network-background test
b[set == "ALL_NETWORK_NODES", `:=`(network_bg = NA_real_, p_network = NA_real_)]
b[, set := factor(set, levels = rev(KEEP))]
bl <- melt(b, id.vars = c("set", "obs_mean_absT"),
           measure.vars = c("network_bg", "transcriptome_bg"),
           variable.name = "background", value.name = "bg_mean")
bl <- bl[!is.na(bg_mean)]
bl[, background := factor(background, levels = c("transcriptome_bg", "network_bg"),
                          labels = c("size-matched draw from the transcriptome (20,250 genes)",
                                     "size-matched draw from the 364 network nodes"))]
pb <- ggplot(bl, aes(y = set)) +
  geom_segment(aes(x = bg_mean, xend = obs_mean_absT, yend = set), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(x = bg_mean), colour = "grey45", size = 1.7, shape = 1) +
  geom_point(aes(x = obs_mean_absT), colour = "#B2182B", size = 2.3) +
  facet_wrap(~ background, ncol = 2) +
  labs(x = expression("mean |t| (tumour vs normal)"), y = NULL,
       title = "Are the FFL motif-class node sets more disease-associated than chance?",
       subtitle = "Red = observed set mean |t|; open circle = mean of 10,000 size-matched random draws.\nEvery class beats the transcriptome; only the 3-node classes beat a within-network background.") +
  theme(axis.text.y = element_text(size = 7), strip.text = element_text(size = 7.2),
        plot.title = element_text(face = "bold"))
ggsave(file.path(FIG, "fig_gsea_ffl_class_backgrounds.png"), pb, width = 9.4, height = 4.6, dpi = 300)
ggsave(file.path(FIG, "fig_gsea_ffl_class_backgrounds.pdf"), pb, width = 9.4, height = 4.6)
msg("FFL background figure written")
fwrite(b, file.path(RES, "gsea_ffl_class_background_summary.csv"))

## ============================ file inventory ================================
fs <- data.table(file = list.files(RES, pattern = "\\.(csv|rds)$", full.names = FALSE))
fs[, size_kb := round(file.size(file.path(RES, file)) / 1024, 1)]
fs[, n_rows := vapply(file, function(f) if (grepl("\\.csv$", f))
    tryCatch(nrow(fread(file.path(RES, f), select = 1L)), error = function(e) NA_integer_)
    else NA_integer_, integer(1))]
setorder(fs, file)
fwrite(fs, file.path(RES, "gsea_OUTPUT_INVENTORY.csv"))
print(fs, nrows = 100)
msg("DONE 10")
