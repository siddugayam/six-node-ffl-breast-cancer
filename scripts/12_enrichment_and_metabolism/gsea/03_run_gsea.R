#!/usr/bin/env Rscript
# v4/03 -- GSEA of every ranked list against every collection (fgsea, BH-adjusted)
suppressPackageStartupMessages({
  library(fgsea); library(data.table); library(BiocParallel)
})
ROOT  <- "/path/to/revision"
RES   <- file.path(ROOT, "results", "v4")
CACHE <- file.path(ROOT, "cache", "v4")
set.seed(20260908)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
NPERM <- 10000L   # nPermSimple for fgseaMultilevel
BPP   <- MulticoreParam(workers = 8, RNGseed = 20260908)

sets <- readRDS(file.path(CACHE, "genesets.rds"))
inv  <- fread(file.path(RES, "gsea_geneset_inventory.csv"))
lab  <- setNames(inv$label, inv$collection)

LISTS <- c("tumour_vs_normal", "mir130a_corr", "mir29a_corr",
           "ffl_signed_influence", "mir130a_high_vs_low")

read_rank <- function(tag) {
  d <- fread(file.path(RES, paste0("rank_", tag, ".csv")))
  stopifnot(all(c("feature", "stat") %in% names(d)))
  if (any(grepl("^[0-9]+$", d$feature)) && mean(grepl("^[0-9]+$", d$feature)) > 0.5)
    stop("FATAL: ranked-list feature column is not gene symbols for ", tag)
  d <- d[!is.na(stat)]
  d <- d[!duplicated(feature)]
  s <- setNames(d$stat, d$feature)
  sort(s, decreasing = TRUE)
}

all_res <- list()
for (L in LISTS) {
  st <- read_rank(L)
  nties <- sum(duplicated(st))
  msg("=== ranked list: ", L, " | n=", length(st),
      " | range ", paste(signif(range(st), 4), collapse = " .. "),
      " | tied values: ", nties)
  ## small network-node list: relax minSize, cap maxSize at list length
  small <- length(st) < 1000
  minS  <- if (small) 5L  else 10L
  maxS  <- if (small) length(st) else 500L
  for (cl in names(sets)) {
    ## FFL_CLASS sets can exceed 500 genes -> allow up to full network
    mx <- if (cl == "FFL_CLASS") max(maxS, 400L) else maxS
    t0 <- Sys.time()
    r <- suppressWarnings(fgsea::fgseaMultilevel(
      pathways = sets[[cl]], stats = st,
      minSize = minS, maxSize = mx,
      nPermSimple = NPERM, eps = 0, nproc = 0, BPPARAM = BPP))
    if (nrow(r) == 0) { msg("   ", cl, ": 0 testable sets"); next }
    r <- as.data.table(r)
    r[, padj := p.adjust(pval, method = "BH")]   # BH within list x collection
    r[, `:=`(ranked_list = L, collection = cl, collection_label = lab[[cl]],
             leadingEdge_size = lengths(leadingEdge),
             leadingEdge = vapply(leadingEdge, paste, "", collapse = ";"),
             minSize_used = minS, maxSize_used = mx,
             nPermSimple = NPERM, method = "fgseaMultilevel(eps=0)")]
    all_res[[paste(L, cl)]] <- r
    msg("   ", cl, ": ", nrow(r), " sets tested, padj<0.05: ", sum(r$padj < 0.05),
        " (", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s)")
  }
}

A <- rbindlist(all_res, use.names = TRUE, fill = TRUE)
A[, padj_global := p.adjust(pval, method = "BH")]
setcolorder(A, c("ranked_list", "collection", "collection_label", "pathway", "size",
                 "ES", "NES", "pval", "padj", "padj_global", "log2err",
                 "leadingEdge_size", "leadingEdge",
                 "minSize_used", "maxSize_used", "nPermSimple", "method"))
setorder(A, ranked_list, collection, pval)
fwrite(A, file.path(RES, "gsea_all_results.csv"))
msg("wrote gsea_all_results.csv: ", nrow(A), " rows")

## significant subset (BH within list x collection)
S <- A[padj < 0.05]
fwrite(S, file.path(RES, "gsea_significant_all.csv"))
msg("significant (padj<0.05): ", nrow(S), " rows")

## per-list files
for (L in LISTS) {
  fwrite(A[ranked_list == L], file.path(RES, paste0("gsea_", L, "_all.csv")))
  fwrite(S[ranked_list == L][order(pval)], file.path(RES, paste0("gsea_", L, "_significant.csv")))
}

## summary counts
smy <- A[, .(n_tested = .N, n_sig = sum(padj < 0.05),
             n_sig_up = sum(padj < 0.05 & NES > 0),
             n_sig_dn = sum(padj < 0.05 & NES < 0),
             n_sig_global = sum(padj_global < 0.05)),
         by = .(ranked_list, collection, collection_label)]
fwrite(smy, file.path(RES, "gsea_summary_counts.csv"))
print(smy)
msg("DONE 03")
