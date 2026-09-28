#!/usr/bin/env Rscript
## ============================================================================
## 13_hub_rank_validation.R
## The hypergeometric test in step 10 is a cut-off test and is under-powered
## (only 4 of the 70 ExIR-significant drivers are canonical-network nodes at
## really calls for:
##   (i)  are FFL hubs ranked higher in the ExIR driver ranking than non-hub
##        network nodes?  (Wilcoxon on ExIR Driver rank)
##   (ii) does the ExIR driver score correlate with FFL participation / degree /
##        betweenness across the network nodes?  (Spearman)
##  (iii) what does ExIR call each of the manuscript's top FFL hubs?
## ============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- "/path/to/revision"
LOG  <- file.path(ROOT, "logs", "exir_agent.log")
say <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | 13 | ", paste0(..., collapse = ""))
  cat(msg, "\n", file = LOG, append = TRUE); cat(msg, "\n"); flush.console()
}
say("=== 13 START ===")

cls <- fread(file.path(ROOT, "results", "exir_classification.csv"))
hub <- fread(file.path(ROOT, "results", "ffl_network_hubs.csv"))

drv <- cls[class == "Driver", .(feature, driver_rank = rank, driver_score = score,
                                driver_padj = P_adj, subtype, in_network, logFC)]
med <- cls[class == "nonDE_mediator", .(feature, med_rank = rank, med_score = score,
                                        med_padj = P_adj)]

## ---- (i) Wilcoxon on driver rank: hubs vs non-hub network nodes ------------
h <- merge(hub, drv, by.x = "node", by.y = "feature", all.x = TRUE)
inD <- h[!is.na(driver_rank)]
say("canonical-network nodes present in the ExIR Driver table (i.e. significantly DE): ", nrow(inD),
    " of ", nrow(hub))
say("  of those, FFL hubs (centrality union) = ", sum(inD$hub_centrality),
    " ; FFL-participation hubs = ", sum(inD$hub_ffl))

wt <- function(flagcol) {
  a <- inD[get(flagcol) == TRUE,  driver_rank]
  b <- inD[get(flagcol) == FALSE, driver_rank]
  if (length(a) < 3 || length(b) < 3) { say("  too few for ", flagcol); return(NULL) }
  w <- wilcox.test(a, b)
  say(sprintf("  Wilcoxon ExIR driver RANK, %s: hubs n=%d median=%.0f | non-hubs n=%d median=%.0f | W=%s p=%.3g (lower rank = stronger driver)",
              flagcol, length(a), median(a), length(b), median(b), format(w$statistic), w$p.value))
  data.table(test = paste0("driver_rank ~ ", flagcol), n_hub = length(a), n_nonhub = length(b),
             median_hub = median(a), median_nonhub = median(b),
             W = as.numeric(w$statistic), p = w$p.value)
}
wtab <- rbindlist(list(wt("hub_centrality"), wt("hub_ffl"), wt("hub_degree"), wt("hub_betweenness")))

## ---- (ii) Spearman: ExIR driver score vs network centrality ----------------
sp <- function(x, y, lab) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  say(sprintf("  Spearman %s: rho=%.3f p=%.3g (n=%d)", lab, ct$estimate, ct$p.value, length(x)))
  data.table(test = lab, rho = as.numeric(ct$estimate), p = ct$p.value, n = length(x))
}
say("Spearman correlations across the ", nrow(inD), " DE network nodes:")
stab <- rbindlist(list(
  sp(inD$driver_score, inD$ffl_participation, "ExIR driver score ~ FFL participation"),
  sp(inD$driver_score, inD$degree_all,        "ExIR driver score ~ degree"),
  sp(inD$driver_score, inD$betweenness,       "ExIR driver score ~ betweenness")))

## same for the mediator score
hm <- merge(hub, med, by.x = "node", by.y = "feature")
say("canonical-network nodes present in the ExIR nonDE-mediator table: ", nrow(hm))
if (nrow(hm) > 10) {
  stab <- rbind(stab, rbindlist(list(
    sp(hm$med_score, hm$ffl_participation, "ExIR nonDE-mediator score ~ FFL participation"),
    sp(hm$med_score, hm$degree_all,        "ExIR nonDE-mediator score ~ degree"),
    sp(hm$med_score, hm$betweenness,       "ExIR nonDE-mediator score ~ betweenness"))))
}
fwrite(rbind(wtab, stab, fill = TRUE), file.path(ROOT, "results", "exir_hub_rank_tests.csv"))
say("wrote results/exir_hub_rank_tests.csv")

## ---- (iii) what does ExIR call the top FFL hubs? ---------------------------
top <- hub[order(-ffl_participation)][1:30]
top <- merge(top, drv[, .(feature, driver_rank, driver_score, driver_padj, subtype, logFC)],
             by.x = "node", by.y = "feature", all.x = TRUE)
top <- merge(top, med[, .(feature, med_rank, med_score, med_padj)],
             by.x = "node", by.y = "feature", all.x = TRUE)
top[, exir_call := fifelse(!is.na(driver_padj) & driver_padj < 0.05, "significant Driver",
                    fifelse(!is.na(driver_rank), "DE, driver-ranked but not ExIR-significant",
                    fifelse(!is.na(med_rank), "nonDE MEDIATOR", "not ranked by ExIR")))]
setorder(top, -ffl_participation)
fwrite(top, file.path(ROOT, "results", "exir_top_ffl_hubs_called.csv"))
say("wrote results/exir_top_ffl_hubs_called.csv rows=", nrow(top))
say("ExIR verdict on the top-30 FFL hubs: ",
    paste(sprintf("%s=%s", top$node, top$exir_call), collapse = "; "))
tt <- top[, .N, by = exir_call]
say("summary top-30 FFL hubs: ", paste(sprintf("%s=%d", tt$exir_call, tt$N), collapse = "; "))

## every hub (union definition), what ExIR calls it
allh <- hub[hub_centrality == TRUE | hub_ffl == TRUE]
allh <- merge(allh, drv[, .(feature, driver_rank, driver_padj)], by.x="node", by.y="feature", all.x=TRUE)
allh <- merge(allh, med[, .(feature, med_rank)], by.x="node", by.y="feature", all.x=TRUE)
allh[, exir_call := fifelse(!is.na(driver_padj) & driver_padj < 0.05, "significant Driver",
                     fifelse(!is.na(driver_rank), "DE, not ExIR-significant",
                     fifelse(!is.na(med_rank), "nonDE MEDIATOR", "not ranked by ExIR")))]
at <- allh[, .N, by = exir_call]
say("ALL ", nrow(allh), " FFL hubs (centrality U participation): ",
    paste(sprintf("%s=%d", at$exir_call, at$N), collapse = "; "))
fwrite(allh, file.path(ROOT, "results", "exir_all_ffl_hubs_called.csv"))
say("=== 13 DONE ===")
