#!/usr/bin/env Rscript
# 10_survival_cox_hubs.R
# PART 1A: Univariable Cox regression for every node of the canonical FFL network
# in TCGA-BRCA primary tumours. Endpoints: OS and PFI. BH-FDR across all features.
# Output: results/survival_cox_hubs.csv
#
# Design notes (for Methods):
#  * Features are z-scored, so HR = hazard ratio per +1 SD of log2 expression.
#  * ALL 576 expression-matched network nodes are tested (not a cherry-picked hub
#    list) so that the FDR denominator is honest. Hub status is a *label* column.
#  * Hubs are defined on the RECONSTRUCTED directed network (evidence-based layers),
#    which is the network the revision defends, not the authors' deposited one.

suppressPackageStartupMessages({
  library(data.table); library(survival); library(igraph)
})

BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs", "10_survival_cox.log")
logf <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(msg, "\n"); cat(msg, "\n", file = LOG, append = TRUE)
}
cat("", file = LOG)
logf("START 10_survival_cox_hubs.R")

## ---------------------------------------------------------------- network ----
nodes <- fread(file.path(BASE, "data/canonical_nodes.tsv"))
cedge <- fread(file.path(BASE, "data/canonical_edges.tsv"))
logf("canonical nodes=", nrow(nodes), " canonical edges=", nrow(cedge))

# Reconstructed, evidence-based directed network:
#   miRNA_target  <- canonical_edges (no reconstructed layer file exists)
#   TF_target/TF-TF <- layer_TF_target.tsv (TRRUST)
#   TF_miRNA      <- layer_TF_miRNA.tsv (TransmiR)
#   gene_gene     <- layer_gene_gene.tsv (TRRUST directed + STRING undirected)
#   miRNA_miRNA   <- layer_miRNA_miRNA.tsv (polycistronic co-transcription)
mir_t <- cedge[edge_type == "miRNA_target", .(source, target, edge_type = "miRNA_target")]
tf_t  <- fread(file.path(BASE, "data/layer_TF_target.tsv"))[, .(source, target, edge_type = "TF_target")]
tf_m  <- fread(file.path(BASE, "data/layer_TF_miRNA.tsv"))[, .(source, target, edge_type = "TF_miRNA")]
gg    <- fread(file.path(BASE, "data/layer_gene_gene.tsv"))[, .(source, target, edge_type = "gene_gene")]
mm    <- fread(file.path(BASE, "data/layer_miRNA_miRNA.tsv"))[, .(source = miRNA_1, target = miRNA_2,
                                                                 edge_type = "miRNA_miRNA")]
recon <- unique(rbindlist(list(mir_t, tf_t, tf_m, gg, mm)))
recon <- recon[source %in% nodes$name & target %in% nodes$name]
logf("reconstructed network edges=", nrow(recon), " by type: ",
     paste(sprintf("%s=%d", names(table(recon$edge_type)), table(recon$edge_type)), collapse = " "))

g <- graph_from_data_frame(recon[, .(source, target)], directed = TRUE,
                           vertices = data.frame(name = nodes$name))
topo <- data.table(
  name        = V(g)$name,
  degree_tot  = degree(g, mode = "all"),
  degree_out  = degree(g, mode = "out"),
  degree_in   = degree(g, mode = "in"),
  betweenness = betweenness(g, directed = TRUE)
)
topo <- merge(topo, nodes[, .(name, type)], by = "name")
setorder(topo, -degree_tot)
TOPN <- 50
hub_deg <- topo[order(-degree_tot)][1:TOPN, name]
hub_btw <- topo[order(-betweenness)][1:TOPN, name]
named_hubs <- c("NFKB1", "RELA", "SP1", "ETS1", "COL1A1", "COL3A1", "VEGFA", "CCND2",
                "MYC", "E2F1", "TP53", "hsa-miR-130a", "hsa-miR-124", "hsa-miR-101",
                "hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "hsa-let-7b", "hsa-let-7e")
topo[, hub_by_degree      := name %in% hub_deg]
topo[, hub_by_betweenness := name %in% hub_btw]
topo[, hub_named_in_MS    := name %in% named_hubs]
topo[, is_hub := hub_by_degree | hub_by_betweenness | hub_named_in_MS]
logf("hubs: degree-top50=", length(hub_deg), " betw-top50=", length(hub_btw),
     " named=", length(named_hubs), " union=", sum(topo$is_hub))
fwrite(topo, file.path(BASE, "results/network_topology_hubs.csv"))

## ------------------------------------------------------------- expression ----
gex <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
mex <- readRDS(file.path(BASE, "data/brca_mirna_expr_canonical.rds"))
pheno <- readRDS(file.path(BASE, "data/brca_pheno.rds"))
surv  <- as.data.table(readRDS(file.path(BASE, "data/brca_survival.rds")))
setnames(surv, "sample", "sample_id")
logf("gene expr ", paste(dim(gex), collapse = "x"), "; mirna expr ",
     paste(dim(mex), collapse = "x"), "; survival rows=", nrow(surv))

tumour <- pheno$sample[pheno$sample_type == "Primary Tumor"]
gex <- gex[, colnames(gex) %in% tumour, drop = FALSE]
mex <- mex[, colnames(mex) %in% tumour, drop = FALSE]
logf("primary-tumour columns: gene=", ncol(gex), " mirna=", ncol(mex))

# clinical covariates
surv[, age := as.numeric(age_at_initial_pathologic_diagnosis)]
stg <- surv$ajcc_pathologic_tumor_stage
stage_group <- rep(NA_character_, length(stg))
stage_group[grepl("^Stage I($|[AB])",   stg)] <- "I"
stage_group[grepl("^Stage II($|[AB])",  stg)] <- "II"
stage_group[grepl("^Stage III($|[ABC])",stg)] <- "III"
stage_group[grepl("^Stage IV",          stg)] <- "IV"
surv[, stage_group := factor(stage_group, levels = c("I", "II", "III", "IV"))]
logf("stage_group: ", paste(sprintf("%s=%d", names(table(surv$stage_group, useNA = "ifany")),
                                    table(surv$stage_group, useNA = "ifany")), collapse = " "))
saveRDS(surv, file.path(BASE, "data/brca_survival_clean.rds"))

## ------------------------------------------------- feature -> expression map --
# node -> row of gex or mex.  Use the match table produced in the prep phase.
nm <- fread(file.path(BASE, "results/node_DE_match.tsv"))
logf("node_DE_match rows=", nrow(nm), " matched=", sum(nm$matched))

feat_rows <- list()
for (i in seq_len(nrow(nodes))) {
  nn <- nodes$name[i]; tp <- nodes$type[i]
  if (tp == "miRNA") {
    if (nn %in% rownames(mex)) feat_rows[[nn]] <- list(mat = "mirna", row = nn)
  } else {
    if (nn %in% rownames(gex)) feat_rows[[nn]] <- list(mat = "gene", row = nn)
  }
}
logf("features with expression: ", length(feat_rows), " / ", nrow(nodes))

## --------------------------------------------------------------- Cox loop ----
run_cox <- function(x, time, event) {
  ok <- is.finite(x) & is.finite(time) & is.finite(event) & time > 0
  if (sum(ok) < 30 || sum(event[ok] == 1) < 5) return(NULL)
  xs <- as.numeric(scale(x[ok]))
  if (!is.finite(sd(xs)) || sd(xs) == 0) return(NULL)
  fit <- try(coxph(Surv(time[ok], event[ok]) ~ xs), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  s <- summary(fit)
  list(n = sum(ok), n_event = sum(event[ok] == 1),
       HR = unname(s$coefficients[1, "exp(coef)"]),
       lo = unname(s$conf.int[1, "lower .95"]),
       hi = unname(s$conf.int[1, "upper .95"]),
       p  = unname(s$coefficients[1, "Pr(>|z|)"]),
       z  = unname(s$coefficients[1, "z"]),
       cidx = unname(s$concordance[1]))
}

res <- list(); k <- 0L
for (nn in names(feat_rows)) {
  fr <- feat_rows[[nn]]
  vec <- if (fr$mat == "gene") gex[fr$row, ] else mex[fr$row, ]
  d <- data.table(sample_id = names(vec), x = as.numeric(vec))
  d <- merge(d, surv[, .(sample_id, OS, OS.time, PFI, PFI.time, DSS, DSS.time)],
             by = "sample_id")
  for (ep in c("OS", "PFI", "DSS")) {
    tt <- d[[paste0(ep, ".time")]]; ee <- d[[ep]]
    r <- run_cox(d$x, tt, ee)
    if (is.null(r)) next
    k <- k + 1L
    res[[k]] <- data.table(feature = nn, endpoint = ep, n = r$n, n_event = r$n_event,
                           HR_per_SD = r$HR, CI_low = r$lo, CI_high = r$hi,
                           z = r$z, p_value = r$p, C_index = r$cidx)
  }
}
cox <- rbindlist(res)
logf("Cox models fitted: ", nrow(cox), " over ", uniqueN(cox$feature), " features")

# BH-FDR *within endpoint* across all features tested (the honest denominator)
cox[, q_value := p.adjust(p_value, method = "BH"), by = endpoint]
# also a single global FDR across every test performed
cox[, q_value_global := p.adjust(p_value, method = "BH")]

cox <- merge(cox, topo[, .(feature = name, node_type = type, degree_tot, betweenness,
                           hub_by_degree, hub_by_betweenness, hub_named_in_MS, is_hub)],
             by = "feature", all.x = TRUE)

# attach DE stats so the table shows expression change AND prognosis together
de_g <- fread(file.path(BASE, "results/BRCA_DEX_genes.csv"))[, .(feature, logFC_TvsN = logFC,
                                                                 adjP_TvsN = adj.P.Val)]
de_m <- fread(file.path(BASE, "results/BRCA_DEX_mirnas.csv"))[, .(feature, logFC_TvsN = logFC,
                                                                  adjP_TvsN = adj.P.Val)]
de <- unique(rbindlist(list(de_g, de_m)))
cox <- merge(cox, de, by = "feature", all.x = TRUE)

setorder(cox, endpoint, p_value)
setcolorder(cox, c("feature", "node_type", "endpoint", "n", "n_event", "HR_per_SD",
                   "CI_low", "CI_high", "z", "p_value", "q_value", "q_value_global",
                   "C_index", "degree_tot", "betweenness", "hub_by_degree",
                   "hub_by_betweenness", "hub_named_in_MS", "is_hub",
                   "logFC_TvsN", "adjP_TvsN"))
outf <- file.path(BASE, "results/survival_cox_hubs.csv")
fwrite(cox, outf)
logf("WROTE ", outf, " rows=", nrow(cox))

## ------------------------------------------------------------------ report ---
for (ep in c("OS", "PFI", "DSS")) {
  sub <- cox[endpoint == ep]
  logf(sprintf("%s: tested=%d  nominal p<0.05=%d  BH q<0.05=%d  q<0.10=%d",
               ep, nrow(sub), sum(sub$p_value < 0.05), sum(sub$q_value < 0.05),
               sum(sub$q_value < 0.10)))
  hs <- sub[is_hub == TRUE]
  logf(sprintf("   hubs only: tested=%d nominal=%d q<0.05=%d",
               nrow(hs), sum(hs$p_value < 0.05), sum(hs$q_value < 0.05)))
  top <- head(sub[order(p_value)], 15)
  for (i in seq_len(nrow(top)))
    logf(sprintf("   %-14s HR=%.3f (%.3f-%.3f) p=%.3g q=%.3g  hub=%s",
                 top$feature[i], top$HR_per_SD[i], top$CI_low[i], top$CI_high[i],
                 top$p_value[i], top$q_value[i], top$is_hub[i]))
}
# the specifically named hubs
logf("---- named manuscript hubs ----")
for (ep in c("OS", "PFI")) {
  sub <- cox[endpoint == ep & hub_named_in_MS == TRUE][order(p_value)]
  for (i in seq_len(nrow(sub)))
    logf(sprintf("%s %-14s HR=%.3f (%.3f-%.3f) p=%.4g q=%.4g",
                 ep, sub$feature[i], sub$HR_per_SD[i], sub$CI_low[i], sub$CI_high[i],
                 sub$p_value[i], sub$q_value[i]))
}
logf("DONE")
