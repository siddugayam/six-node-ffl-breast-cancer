#!/usr/bin/env Rscript
# 16_external_geo_de.R
# EXTERNAL COHORT 3 (+4): GEO breast tumour-vs-normal series on Affymetrix GPL570.
#   GSE42568 - 17 normal breast + 104 breast tumours
#   GSE45827 - breast tumour subtypes incl. normal breast, if normals are present
# Replicates the DIRECTION of the TCGA tumour-vs-normal differential expression for
# independent tissue-based cohort for the DE claims.
# Output: results/external_GEO_DE.csv

suppressPackageStartupMessages({
  library(data.table); library(limma); library(hgu133plus2.db); library(AnnotationDbi)
})
BASE <- "/path/to/revision"
EXT  <- file.path(BASE, "cache/external")
LOG  <- file.path(BASE, "logs", "16_external_geo_de.log")
logf <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                        cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); logf("START")

## ---------------------------------------------------- series-matrix reader ---
read_series <- function(path) {
  ln <- readLines(path)
  gsm <- grep("^!Sample_geo_accession", ln, value = TRUE)[1]
  gsm <- gsub('"', "", strsplit(gsm, "\t")[[1]][-1])
  ch  <- grep("^!Sample_characteristics_ch1", ln, value = TRUE)
  ch  <- lapply(ch, function(x) gsub('"', "", strsplit(x, "\t")[[1]][-1]))
  ti  <- grep("^!Sample_title", ln, value = TRUE)[1]
  ti  <- gsub('"', "", strsplit(ti, "\t")[[1]][-1])
  b <- grep("^!series_matrix_table_begin", ln); e <- grep("^!series_matrix_table_end", ln)
  tab <- data.table::fread(text = paste(ln[(b + 1):(e - 1)], collapse = "\n"))
  setnames(tab, 1, "ID_REF")
  list(gsm = gsm, ch = ch, title = ti, tab = tab)
}

probe2sym <- function(ids) {
  m <- suppressMessages(AnnotationDbi::select(hgu133plus2.db, keys = ids,
                                              keytype = "PROBEID", columns = "SYMBOL"))
  m <- as.data.table(m)[!is.na(SYMBOL)]
  m[, .SD[1], by = PROBEID]
}

de_one <- function(gse, path, group_fun) {
  logf("==== ", gse, " ====")
  s <- read_series(path)
  logf(gse, ": ", length(s$gsm), " samples, expression table ",
       nrow(s$tab), " probes x ", ncol(s$tab) - 1, " cols")
  grp <- group_fun(s)
  logf(gse, " group assignment: ",
       paste(sprintf("%s=%d", names(table(grp, useNA = "ifany")),
                     table(grp, useNA = "ifany")), collapse = " "))
  keep <- !is.na(grp)
  X <- as.matrix(s$tab[, -1]); rownames(X) <- s$tab$ID_REF
  X <- X[, keep, drop = FALSE]; grp <- droplevels(factor(grp[keep], levels = c("normal", "tumour")))
  X <- X[stats::complete.cases(X), , drop = FALSE]
  logf(gse, ": probes with complete data = ", nrow(X))
  if (max(X, na.rm = TRUE) > 50) { X <- log2(X + 1); logf(gse, ": log2-transformed (values were on linear scale)") }
  logf(gse, ": value range ", paste(round(range(X), 2), collapse = " to "))

  map <- probe2sym(rownames(X))
  X <- X[map$PROBEID, , drop = FALSE]
  sym <- map$SYMBOL
  # collapse probes to genes by the probe with the highest mean expression
  o <- order(rowMeans(X), decreasing = TRUE)
  X <- X[o, , drop = FALSE]; sym <- sym[o]
  X <- X[!duplicated(sym), , drop = FALSE]; rownames(X) <- sym[!duplicated(sym)]
  logf(gse, ": gene-level matrix ", paste(dim(X), collapse = " x "))

  des <- model.matrix(~ grp)
  fit <- eBayes(lmFit(X, des))
  tt <- as.data.table(topTable(fit, coef = 2, number = Inf), keep.rownames = "feature")
  tt[, gse := gse]
  logf(gse, ": DE done; adj.P<0.05 in ", sum(tt$adj.P.Val < 0.05), " / ", nrow(tt), " genes")
  tt[, .(gse, feature, logFC, P.Value, adj.P.Val, n_normal = sum(grp == "normal"),
         n_tumour = sum(grp == "tumour"))]
}

## ------------------------------------------------------------- GSE42568 -----
g1 <- de_one("GSE42568", file.path(EXT, "GSE42568_series_matrix.txt.gz"), function(s) {
  tis <- s$ch[[1]]
  ifelse(grepl("normal", tis, ignore.case = TRUE), "normal",
         ifelse(grepl("breast cancer|tumor|tumour|cancer", tis, ignore.case = TRUE),
                "tumour", NA_character_))
})

## ------------------------------------------------------------- GSE45827 -----
## The tumor-subtype characteristics field is mis-aligned across samples in this
## series (normals appear as "N/A" and cell lines as "cell origin: ..."), so the
## group is taken from Sample_title, which carries it unambiguously.
g2 <- try({
  de_one("GSE45827", file.path(EXT, "GSE45827_series_matrix.txt.gz"), function(s) {
    ti <- s$title
    logf("GSE45827 title prefixes: ",
         paste(sprintf("%s=%d", names(table(sub("[0-9].*$", "", ti))),
                       table(sub("[0-9].*$", "", ti))), collapse = " | "))
    ifelse(grepl("^normal", ti, ignore.case = TRUE), "normal",
           ifelse(grepl("^(basal|her2|luminal)", ti, ignore.case = TRUE),
                  "tumour", NA_character_))   # cell lines excluded
  })
}, silent = TRUE)
if (inherits(g2, "try-error")) { logf("GSE45827 FAILED: ", as.character(g2)); g2 <- NULL }

## ------------------------------------------------------------- GSE10780 -----
g3 <- try({
  de_one("GSE10780", file.path(EXT, "GSE10780_series_matrix.txt.gz"), function(s) {
    tis <- s$ch[[1]]
    logf("GSE10780 characteristics: ",
         paste(sprintf("%s=%d", names(table(tis)), table(tis)), collapse = " | "))
    ifelse(grepl("normal", tis, ignore.case = TRUE), "normal",
           ifelse(grepl("IDC|carcinoma|tumor|tumour|malignant", tis, ignore.case = TRUE),
                  "tumour", NA_character_))
  })
}, silent = TRUE)
if (inherits(g3, "try-error")) { logf("GSE10780 FAILED: ", as.character(g3)); g3 <- NULL }

geo <- rbindlist(list(g1, g2, g3), use.names = TRUE, fill = TRUE)
fwrite(geo, file.path(BASE, "results/external_GEO_DE.csv"))
logf("WROTE results/external_GEO_DE.csv rows=", nrow(geo))

## ------------------------------------- replication vs TCGA DE (direction) ----
tcga <- fread(file.path(BASE, "results/BRCA_DEX_genes.csv"))
nodes <- fread(file.path(BASE, "data/canonical_nodes.tsv"))
prot <- nodes[type != "miRNA", name]

rep_rows <- list(); k <- 0L
for (this_gse in unique(geo$gse)) {
  gg <- geo[gse == this_gse][, .(feature, logFC_ext = logFC, p_ext = P.Value,
                                 q_ext = adj.P.Val, n_normal, n_tumour)]
  gse <- this_gse
  m <- merge(tcga[, .(feature, logFC_tcga = logFC, q_tcga = adj.P.Val)], gg, by = "feature")
  m[, is_network_node := feature %in% prot]
  m[, gse := this_gse]
  # replication is judged on nodes that were SIGNIFICANT in TCGA
  sig <- m[is_network_node == TRUE & q_tcga < 0.05]
  same <- sig[sign(logFC_tcga) == sign(logFC_ext)]
  same_sig <- same[q_ext < 0.05]
  logf(sprintf("%s: network protein-coding nodes shared with TCGA = %d ; TCGA-significant = %d",
               gse, sum(m$is_network_node), nrow(sig)))
  logf(sprintf("%s: same DIRECTION as TCGA = %d/%d (%.1f%%) ; same direction AND q_ext<0.05 = %d/%d (%.1f%%)",
               gse, nrow(same), nrow(sig), 100 * nrow(same) / nrow(sig),
               nrow(same_sig), nrow(sig), 100 * nrow(same_sig) / nrow(sig)))
  ct <- stats::cor.test(sig$logFC_tcga, sig$logFC_ext, method = "spearman", exact = FALSE)
  logf(sprintf("%s: Spearman rho(logFC TCGA, logFC %s) over TCGA-significant nodes = %.3f (p=%.3g)",
               gse, gse, ct$estimate, ct$p.value))
  k <- k + 1L
  rep_rows[[k]] <- m
}
repl <- rbindlist(rep_rows)
fwrite(repl, file.path(BASE, "results/external_GEO_DE_vs_TCGA.csv"))
logf("WROTE results/external_GEO_DE_vs_TCGA.csv rows=", nrow(repl))

## ---- consensus across the GEO series: agreement in ALL cohorts tested -------
cons <- repl[is_network_node == TRUE & q_tcga < 0.05,
             .(n_cohorts = .N,
               n_same_dir = sum(sign(logFC_tcga) == sign(logFC_ext)),
               n_same_dir_sig = sum(sign(logFC_tcga) == sign(logFC_ext) & q_ext < 0.05),
               logFC_tcga = logFC_tcga[1]), by = feature]
cons[, replicates_all_cohorts := n_same_dir == n_cohorts]
cons[, replicates_all_cohorts_sig := n_same_dir_sig == n_cohorts]
fwrite(cons, file.path(BASE, "results/external_GEO_DE_consensus.csv"))
logf("CONSENSUS over ", uniqueN(repl$gse), " GEO series, TCGA-significant network nodes n=",
     nrow(cons))
logf(sprintf("  same direction in ALL series: %d/%d (%.1f%%)",
             sum(cons$replicates_all_cohorts), nrow(cons),
             100 * mean(cons$replicates_all_cohorts)))
logf(sprintf("  same direction AND q<0.05 in ALL series: %d/%d (%.1f%%)",
             sum(cons$replicates_all_cohorts_sig), nrow(cons),
             100 * mean(cons$replicates_all_cohorts_sig)))

named <- c("NFKB1", "RELA", "SP1", "ETS1", "COL1A1", "COL3A1", "VEGFA", "CCND2",
           "MYC", "E2F1", "TP53")
for (this_gse in unique(repl$gse)) {
  logf("---- named hubs in ", this_gse, " (logFC tumour vs normal) ----")
  for (gn in named) {
    x <- repl[gse == this_gse & feature == gn]
    if (nrow(x) == 0) { logf(sprintf("%-8s not on array", gn)); next }
    logf(sprintf("%-8s TCGA logFC=%+.3f (q=%.3g) | %s logFC=%+.3f (q=%.3g) | same direction=%s",
                 gn, x$logFC_tcga, x$q_tcga, this_gse, x$logFC_ext, x$q_ext,
                 sign(x$logFC_tcga) == sign(x$logFC_ext)))
  }
}
logf("DONE")
