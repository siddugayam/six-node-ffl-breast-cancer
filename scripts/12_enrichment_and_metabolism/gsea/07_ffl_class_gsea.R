#!/usr/bin/env Rscript
# v4/07 -- (C) GSEA comparison of the FFL motif classes against the disease signature
suppressPackageStartupMessages({
  library(fgsea); library(data.table); library(BiocParallel)
})
ROOT  <- "/path/to/revision"
RES   <- file.path(ROOT, "results", "v4")
CACHE <- file.path(ROOT, "cache", "v4")
set.seed(20260908)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
BPP <- MulticoreParam(workers = 8, RNGseed = 20260908)

sets  <- readRDS(file.path(CACHE, "genesets.rds"))
ffl   <- sets$FFL_CLASS
A     <- fread(file.path(RES, "gsea_all_results.csv"))
CLS   <- c("3-miR", "3-TF", "3-Comp", "4-node", "5-node", "6-node")

## ---------------- overlap between class node sets (honesty check) -----------
J <- CJ(a = CLS, b = CLS, sorted = FALSE)
J[, n_a := lengths(ffl)[a]]; J[, n_b := lengths(ffl)[b]]
J[, n_shared := mapply(function(x, y) length(intersect(ffl[[x]], ffl[[y]])), a, b)]
J[, jaccard := n_shared / (n_a + n_b - n_shared)]
fwrite(J, file.path(RES, "gsea_ffl_class_setoverlap.csv"))
msg("class-set Jaccard (off-diagonal) range: ",
    paste(signif(range(J[a != b]$jaccard), 4), collapse = " .. "))
print(dcast(J, a ~ b, value.var = "jaccard"))

## ---------------- derived, non-overlapping sets -----------------------------
u3 <- unique(c(ffl[["3-miR"]], ffl[["3-TF"]], ffl[["3-Comp"]]))
extra <- list(
  "3node_union"            = u3,
  "4node_not_3node"        = setdiff(ffl[["4-node"]], u3),
  "5node_not_3node"        = setdiff(ffl[["5-node"]], u3),
  "6node_not_3node"        = setdiff(ffl[["6-node"]], u3),
  "higherorder_union"      = unique(c(ffl[["4-node"]], ffl[["5-node"]], ffl[["6-node"]])),
  "higherorder_not_3node"  = setdiff(unique(c(ffl[["4-node"]], ffl[["5-node"]], ffl[["6-node"]])), u3),
  "3miR_only"              = setdiff(ffl[["3-miR"]], unique(c(ffl[["3-TF"]], ffl[["3-Comp"]]))),
  "network_not_in_any_FFL" = setdiff(ffl[["ALL_NETWORK_NODES"]],
                                     unique(c(u3, ffl[["4-node"]], ffl[["5-node"]], ffl[["6-node"]]))),
  "ALL_NETWORK_NODES"      = ffl[["ALL_NETWORK_NODES"]]
)
extra <- extra[lengths(extra) >= 5]
msg("derived sets: ", paste(names(extra), lengths(extra), sep = "=", collapse = "  "))

LISTS <- c("tumour_vs_normal", "mir29a_corr")
read_rank <- function(tag) {
  d <- fread(file.path(RES, paste0("rank_", tag, ".csv")))
  d <- d[!is.na(stat)][!duplicated(feature)]
  sort(setNames(d$stat, d$feature), decreasing = TRUE)
}
ex_res <- rbindlist(lapply(LISTS, function(L) {
  st <- read_rank(L)
  r <- suppressWarnings(fgseaMultilevel(extra, st, minSize = 5, maxSize = 600,
                                        nPermSimple = 10000, eps = 0, nproc = 0, BPPARAM = BPP))
  r <- as.data.table(r)
  r[, padj := p.adjust(pval, "BH")]
  r[, `:=`(ranked_list = L, collection = "FFL_DERIVED",
           leadingEdge_size = lengths(leadingEdge),
           leadingEdge = vapply(leadingEdge, paste, "", collapse = ";"))]
  r
}))
fwrite(ex_res, file.path(RES, "gsea_ffl_derived_sets.csv"))
msg("derived-set GSEA rows: ", nrow(ex_res))

## ---------------- assemble the class comparison table -----------------------
cls <- A[collection == "FFL_CLASS"]
comp <- rbind(
  cls[, .(ranked_list, set = pathway, set_kind = "class", size, ES, NES, pval, padj,
          padj_global, log2err, leadingEdge_size, leadingEdge)],
  ex_res[, .(ranked_list, set = pathway, set_kind = "derived", size, ES, NES, pval, padj,
             padj_global = NA_real_, log2err, leadingEdge_size, leadingEdge)]
)

## ---------------- network-background permutation ---------------------------
## For each set: is its mean tumour-vs-normal |t| more extreme than size-matched
## random draws from the 364 protein-coding network nodes?
tvn <- fread(file.path(RES, "rank_tumour_vs_normal.csv"))
tvec <- setNames(tvn$stat, tvn$feature)
bg   <- intersect(ffl[["ALL_NETWORK_NODES"]], names(tvec))
NB   <- 10000L
allsets <- c(ffl[CLS], extra)
nbres <- rbindlist(lapply(names(allsets), function(nm) {
  g <- intersect(allsets[[nm]], names(tvec)); n <- length(g)
  if (n < 5 || n >= length(bg)) return(data.table(set = nm, n_in_bg = n,
      obs_mean_absT = NA_real_, bg_mean_absT = NA_real_, p_vs_network_bg = NA_real_,
      obs_mean_T = NA_real_, p_signed_vs_network_bg = NA_real_))
  obs  <- mean(abs(tvec[g])); obss <- mean(tvec[g])
  draw <- replicate(NB, { s <- sample(bg, n); c(mean(abs(tvec[s])), mean(tvec[s])) })
  data.table(set = nm, n_in_bg = n, obs_mean_absT = obs,
             bg_mean_absT = mean(draw[1, ]),
             p_vs_network_bg = (1 + sum(draw[1, ] >= obs)) / (NB + 1),
             obs_mean_T = obss,
             p_signed_vs_network_bg = (1 + sum(abs(draw[2, ]) >= abs(obss))) / (NB + 1))
}))
fwrite(nbres, file.path(RES, "gsea_ffl_class_network_background_test.csv"))
print(nbres)

## also: transcriptome-wide size-matched null on mean |t| (context)
allg <- names(tvec)
twres <- rbindlist(lapply(names(allsets), function(nm) {
  g <- intersect(allsets[[nm]], allg); n <- length(g)
  if (n < 5) return(NULL)
  obs <- mean(abs(tvec[g]))
  draw <- replicate(NB, mean(abs(tvec[sample(allg, n)])))
  data.table(set = nm, n = n, obs_mean_absT = obs, transcriptome_mean_absT = mean(draw),
             p_vs_transcriptome = (1 + sum(draw >= obs)) / (NB + 1))
}))
fwrite(twres, file.path(RES, "gsea_ffl_class_transcriptome_background_test.csv"))
print(twres)

comp <- merge(comp, nbres[, .(set, n_in_network_bg = n_in_bg, obs_mean_absT,
                              network_bg_mean_absT = bg_mean_absT, p_vs_network_bg)],
              by = "set", all.x = TRUE)
comp <- merge(comp, twres[, .(set, transcriptome_mean_absT, p_vs_transcriptome)],
              by = "set", all.x = TRUE)
setorder(comp, ranked_list, -NES)
fwrite(comp, file.path(RES, "gsea_ffl_class_comparison.csv"))
msg("wrote gsea_ffl_class_comparison.csv (", nrow(comp), " rows)")

cat("\n===== FFL classes vs the tumour-vs-normal disease signature =====\n")
print(comp[ranked_list == "tumour_vs_normal" & set %in% c(CLS, names(extra)),
           .(set, set_kind, size, NES, pval, padj, leadingEdge_size,
             obs_mean_absT, network_bg_mean_absT, p_vs_network_bg,
             transcriptome_mean_absT, p_vs_transcriptome)][order(-NES)])
msg("DONE 07")
