#!/usr/bin/env Rscript
## ============================================================================
## 10_exir_classify_overlap.R
##  B) build results/exir_classification.csv
##  C) ExIR-driver  vs  FFL-hub overlap (hypergeometric, BOTH directions)
##  D) mediator |logFC| distribution vs drivers vs biomarkers
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(igraph)
})
ROOT <- "/path/to/revision"
LOG  <- file.path(ROOT, "logs", "exir_agent.log")
say <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | 10 | ", paste0(..., collapse = ""))
  cat(msg, "\n", file = LOG, append = TRUE); cat(msg, "\n"); flush.console()
}
say("=== 10 START ===")

res   <- readRDS(file.path(ROOT, "data", "exir_result.rds"))
hub   <- fread(file.path(ROOT, "results", "ffl_network_hubs.csv"))
nodes <- fread(file.path(ROOT, "data", "canonical_nodes.tsv"))
de    <- fread(file.path(ROOT, "results", "BRCA_DEX_ALL_nodes.csv"))
meta  <- readRDS(file.path(ROOT, "data", "exir_inputs_meta.rds"))
trr   <- fread(file.path(ROOT, "data", "db", "trrust_human.tsv"), header = FALSE)
tf_trrust <- unique(trr$V1)

say("ExIR result tables: ", paste(names(res), collapse = " | "))

## ------------------------------------------------------------------ helpers
gname <- V(res$Graph)$name
say("ExIR association graph: V=", vcount(res$Graph), " E=", ecount(res$Graph))

mirna_ids <- rownames(readRDS(file.path(ROOT, "data", "brca_mirna_expr_canonical.rds")))
netmap <- setNames(nodes$type, nodes$name)

node_type_of <- function(x) {
  ty <- netmap[x]
  out <- ifelse(!is.na(ty), ty,
         ifelse(x %in% mirna_ids, "miRNA_notInNetwork",
         ifelse(x %in% tf_trrust, "TF_notInNetwork", "Gene_notInNetwork")))
  unname(out)
}

de_l <- copy(de); setkey(de_l, Gene)
lfc_of <- function(x) de_l[x, logFC]
adjp_of <- function(x) de_l[x, adj.P.Val]

## ------------------------------------------------- B) CLASSIFICATION TABLE
mk <- function(tab, idcol, cls) {
  if (is.null(tab)) return(NULL)
  d <- as.data.table(tab)
  setnames(d, idcol, "feature")
  d[, class := cls]
  d
}
drv <- mk(res[["Driver table"]],        "Driver",         "Driver")
bio <- mk(res[["Biomarker table"]],     "Biomarker",      "Biomarker")
dem <- mk(res[["DE-mediator table"]],   "DE.mediator",    "DE_mediator")
ndm <- mk(res[["nonDE-mediator table"]],"non.DE.mediator","nonDE_mediator")

cols <- c("feature","class","Rank","Score","Z.score","P.value","P.adj")
sub  <- function(d, extra = NULL) {
  if (is.null(d)) return(NULL)
  d[, c(cols, extra), with = FALSE]
}
drv2 <- sub(drv, "Type"); if (!is.null(drv2)) setnames(drv2, "Type", "subtype")
bio2 <- sub(bio, "Type"); if (!is.null(bio2)) setnames(bio2, "Type", "subtype")
dem2 <- sub(dem); if (!is.null(dem2)) dem2[, subtype := NA_character_]
ndm2 <- sub(ndm); if (!is.null(ndm2)) ndm2[, subtype := NA_character_]

cls <- rbindlist(list(drv2, bio2, dem2, ndm2), use.names = TRUE, fill = TRUE)
setnames(cls, c("Rank","Score","P.adj"), c("rank","score","P_adj"))
cls[, node_type   := node_type_of(feature)]
cls[, in_network  := feature %in% nodes$name]
cls[, logFC       := lfc_of(feature)]
cls[, DE_adj_P    := adjp_of(feature)]
cls[, exir_significant := P_adj < 0.05]
setcolorder(cls, c("feature","class","rank","score","node_type","subtype",
                   "Z.score","P.value","P_adj","exir_significant",
                   "in_network","logFC","DE_adj_P"))
setorder(cls, class, rank)
fwrite(cls, file.path(ROOT, "results", "exir_classification.csv"))
say("wrote results/exir_classification.csv rows=", nrow(cls))

cnt <- cls[, .(n = .N, n_sig = sum(exir_significant, na.rm = TRUE)), by = class]
print(cnt)
say("class counts: ", paste(sprintf("%s=%d(sig %d)", cnt$class, cnt$n, cnt$n_sig), collapse = "; "))
if (!is.null(drv2)) {
  st <- drv2[, .N, by = subtype]
  say("driver subtypes: ", paste(sprintf("%s=%d", st$subtype, st$N), collapse = "; "))
  sts <- drv2[P.adj < 0.05, .N, by = subtype]
  say("driver subtypes (P.adj<0.05): ", paste(sprintf("%s=%d", sts$subtype, sts$N), collapse = "; "))
}
if (!is.null(bio2)) {
  bt <- bio2[, .N, by = subtype]
  say("biomarker types: ", paste(sprintf("%s=%d", bt$subtype, bt$N), collapse = "; "))
}
by_ty <- cls[, .N, by = .(class, node_type)]
say("class x node_type:"); print(dcast(by_ty, class ~ node_type, value.var = "N", fill = 0))

## ------------------------------------------------- C) OVERLAP WITH FFL HUBS
sig_drivers <- drv2[P.adj < 0.05, feature]
top100_drv  <- drv2[order(Rank)][1:min(100, .N), feature]
say("significant ExIR drivers (Driver P.adj<0.05): ", length(sig_drivers))

hub_cent <- hub[hub_centrality == TRUE, node]
hub_fflp <- hub[hub_ffl == TRUE, node]
say("FFL hubs: centrality(top-decile degree|betweenness) = ", length(hub_cent),
    " ; FFL-participation(top decile) = ", length(hub_fflp))

exir_feats <- unique(c(rownames(meta$Diff_data), gname))   # everything ExIR could rank
say("ExIR feature space (Diff_data U association-network nodes): ", length(exir_feats))

hyper <- function(A, B, U, label) {
  A <- intersect(A, U); B <- intersect(B, U)
  ov <- intersect(A, B)
  p  <- phyper(length(ov) - 1, length(B), length(U) - length(B), length(A), lower.tail = FALSE)
  expct <- length(A) * length(B) / length(U)
  say(sprintf("%s : |U|=%d |A(drivers)|=%d |B(hubs)|=%d overlap=%d expected=%.2f fold=%.2f p=%.3g",
              label, length(U), length(A), length(B), length(ov), expct,
              length(ov)/max(expct, 1e-12), p))
  data.table(test = label, universe = length(U), n_drivers = length(A), n_hubs = length(B),
             overlap = length(ov), expected = expct,
             fold_enrichment = length(ov)/max(expct, 1e-12), p_hyper = p,
             overlap_features = paste(sort(ov), collapse = ";"))
}

U_net    <- intersect(nodes$name, exir_feats)      # network nodes ExIR could rank
U_global <- exir_feats
say("universe A (canonical network nodes inside ExIR feature space): ", length(U_net))

ov_tab <- rbindlist(list(
  hyper(sig_drivers, hub_cent, U_net,    "sigDrivers vs centralityHubs | universe = network nodes in ExIR space"),
  hyper(sig_drivers, hub_fflp, U_net,    "sigDrivers vs FFLparticipationHubs | universe = network nodes in ExIR space"),
  hyper(sig_drivers, hub_cent, U_global, "sigDrivers vs centralityHubs | universe = whole ExIR feature space"),
  hyper(sig_drivers, hub_fflp, U_global, "sigDrivers vs FFLparticipationHubs | universe = whole ExIR feature space"),
  hyper(top100_drv,  hub_cent, U_global, "top100Drivers vs centralityHubs | universe = whole ExIR feature space"),
  hyper(top100_drv,  hub_fflp, U_global, "top100Drivers vs FFLparticipationHubs | universe = whole ExIR feature space")
))
## sensitivity: universe restricted to network nodes that are significantly DE
U_netDE <- intersect(U_net, rownames(meta$Diff_data))
ov_tab <- rbind(ov_tab,
  hyper(sig_drivers, hub_cent, U_netDE, "sigDrivers vs centralityHubs | universe = DE-significant network nodes"),
  hyper(sig_drivers, hub_fflp, U_netDE, "sigDrivers vs FFLparticipationHubs | universe = DE-significant network nodes"))
fwrite(ov_tab, file.path(ROOT, "results", "exir_driver_hub_overlap.csv"))
say("wrote results/exir_driver_hub_overlap.csv rows=", nrow(ov_tab))

## ---- both directions, named -------------------------------------------------
drvrank <- drv2[order(Rank)]
drvrank[, is_hub_cent := feature %in% hub_cent]
drvrank[, is_hub_ffl  := feature %in% hub_fflp]
drvrank[, in_network  := feature %in% nodes$name]
drvrank[, node_type   := node_type_of(feature)]

miss <- drvrank[P.adj < 0.05 & !is_hub_cent & !is_hub_ffl][order(Rank)]
say("significant drivers that are NOT FFL hubs (either definition): ", nrow(miss))
say("  of which absent from the 587-node network entirely: ", sum(!miss$in_network))
say("TOP 30 drivers missed by the FFL analysis: ",
    paste(sprintf("%s(%s,rank%d,logFC%.2f)", miss$feature[1:min(30,nrow(miss))],
                  miss$node_type[1:min(30,nrow(miss))],
                  miss$Rank[1:min(30,nrow(miss))],
                  lfc_of(miss$feature)[1:min(30,nrow(miss))]), collapse = "; "))
fwrite(drvrank, file.path(ROOT, "results", "exir_drivers_vs_hubs.csv"))

hit <- drvrank[P.adj < 0.05 & (is_hub_cent | is_hub_ffl)][order(Rank)]
say("significant drivers that ARE FFL hubs: ", nrow(hit), " -> ",
    paste(sprintf("%s(rank%d)", hit$feature, hit$Rank), collapse = "; "))

## hubs that ExIR does NOT call a driver (reverse direction)
hub_not_driver <- setdiff(intersect(union(hub_cent, hub_fflp), U_net), sig_drivers)
say("FFL hubs (either definition) NOT significant ExIR drivers: ", length(hub_not_driver),
    " of ", length(intersect(union(hub_cent, hub_fflp), U_net)))
say("   e.g. ", paste(head(hub_not_driver, 40), collapse = "; "))

## ------------------------------------------------- D) MEDIATOR |logFC|
mediators <- if (!is.null(ndm2)) ndm2$feature else character(0)
say("nonDE mediators: ", length(mediators))
grp <- rbindlist(list(
  data.table(feature = drv2$feature,  grp = "Driver"),
  data.table(feature = bio2$feature,  grp = "Biomarker"),
  data.table(feature = mediators,     grp = "nonDE_mediator")
))
grp[, absLFC := abs(lfc_of(feature))]
grp[, adjP   := adjp_of(feature)]
qs <- grp[!is.na(absLFC), as.list(c(n = .N, summary(absLFC),
                                    pct_absLFC_gt1 = 100*mean(absLFC > 1),
                                    pct_adjP_lt05  = 100*mean(adjP < 0.05, na.rm = TRUE))), by = grp]
print(qs)
fwrite(qs, file.path(ROOT, "results", "exir_mediator_logFC_summary.csv"))
say("|logFC| summary written to results/exir_mediator_logFC_summary.csv")
for (i in seq_len(nrow(qs))) {
  say(sprintf("  %s: n=%d median|logFC|=%.3f mean=%.3f max=%.3f  %%|logFC|>1 = %.1f%%  %%adjP<0.05 = %.1f%%",
      qs$grp[i], qs$n[i], qs$Median[i], qs$Mean[i], qs$Max.[i],
      qs$pct_absLFC_gt1[i], qs$pct_adjP_lt05[i]))
}
w1 <- wilcox.test(grp[grp=="nonDE_mediator", absLFC], grp[grp=="Driver", absLFC])
w2 <- wilcox.test(grp[grp=="nonDE_mediator", absLFC], grp[grp=="Biomarker", absLFC])
say("Wilcoxon |logFC| mediators vs drivers   : W=", format(w1$statistic), " p=", format.pval(w1$p.value))
say("Wilcoxon |logFC| mediators vs biomarkers: W=", format(w2$statistic), " p=", format.pval(w2$p.value))

## top mediators, with type + whether they are FFL hubs
if (length(mediators)) {
  med <- ndm2[order(Rank)]
  med[, node_type := node_type_of(feature)]
  med[, in_network := feature %in% nodes$name]
  med[, is_hub := feature %in% union(hub_cent, hub_fflp)]
  med[, logFC := lfc_of(feature)][, DE_adj_P := adjp_of(feature)]
  fwrite(med, file.path(ROOT, "results", "exir_mediators.csv"))
  say("wrote results/exir_mediators.csv rows=", nrow(med))
  say("TOP 25 mediators: ", paste(sprintf("%s(%s,logFC %.2f,adjP %.2g)",
      med$feature[1:min(25,nrow(med))], med$node_type[1:min(25,nrow(med))],
      med$logFC[1:min(25,nrow(med))], med$DE_adj_P[1:min(25,nrow(med))]), collapse = "; "))
  say("mediators that are canonical-network nodes: ", sum(med$in_network),
      " -> ", paste(med[in_network == TRUE, feature], collapse = "; "))
  say("mediators that are FFL hubs: ", sum(med$is_hub),
      " -> ", paste(med[is_hub == TRUE, feature], collapse = "; "))
  say("mediator node_type breakdown: ",
      paste(sprintf("%s=%d", med[, .N, by=node_type]$node_type, med[, .N, by=node_type]$N), collapse="; "))
}
say("=== 10 DONE ===")
