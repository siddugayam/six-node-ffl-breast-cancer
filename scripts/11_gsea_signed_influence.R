#!/usr/bin/env Rscript
# =============================================================================
# 11_gsea_signed_influence.R
# The manuscript (Methods 5-6, MS.md lines 230-256) states that a Signed Random Walk
# with Restart was run from the top 60 hubs (composite centrality) to give every node a
# "Signed Influence Score", and that GSEA was then run on the COMPLETE ranked list of
# network nodes sorted by that score.  This script re-implements that and reports
# NES, p.adjust and leading-edge genes.
#
# TWO SIGN SCHEMES ARE RUN, because the sign assignment is the single largest
# methodological error in the deposited pipeline:
#   as_published : miRNA-sourced edge = -1, every other edge = +1
#                  (exactly MS.md lines 234-236 and Scripts/03.1 line 508)
#   curated      : sign taken from TRRUST / TransmiR curation
#                  (Activation +1, Repression -1); miRNA_target = -1;
#                  edges with no curated sign, and the non-regulatory
#                  gene_gene (STRING association) and miRNA_miRNA (genomic
#                  co-transcription) layers, are DROPPED rather than assumed +1.
#
# TWO FIXES to the deposited implementation, both disclosed:
#   (i)  the deposited code propagates with A %*% v, i.e. from targets back to sources.
#        Influence must flow source -> target, so we use t(W_norm) %*% v.
#   (ii) column normalisation is by out-strength of the SOURCE, not of the target.
#
# Output: results/gsea_results.csv  and  results/signed_influence_scores.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(igraph); library(clusterProfiler); library(org.Hs.eg.db)
  library(AnnotationDbi)
})
REV <- "/path/to/revision"
LOGF <- file.path(REV, "logs", "11_gsea.log")
lg <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOGF, append = TRUE) }
cat("", file = LOGF); set.seed(1234)

ALPHA <- 0.85; MAXIT <- 100; TOL <- 1e-9; N_SEED <- 60

nodes <- fread(file.path(REV, "data/canonical_nodes.tsv"))
ced <- fread(file.path(REV, "data/canonical_edges.tsv"))
ltf <- fread(file.path(REV, "data/layer_TF_target.tsv"))
ltm <- fread(file.path(REV, "data/layer_TF_miRNA.tsv"))
ed <- fread(file.path(REV, "results/graph_used_edges.tsv"))
lg("graph_used_edges rows: ", nrow(ed))

# ---- sign lookup from curation
cur <- rbind(ltf[, .(source, target, mode)], ltm[, .(source, target, mode)])
cur[, sgn := fifelse(mode == "Activation", 1, fifelse(mode == "Repression", -1, NA_real_))]
setkey(cur, source, target)

E <- copy(ed)
E <- merge(E, unique(cur[, .(source, target, sgn)]), by = c("source", "target"), all.x = TRUE)
ty <- setNames(nodes$type, nodes$name)
E[, src_type := ty[source]]

E[, w_pub := fifelse(src_type == "miRNA", -1, 1)]
E[, w_cur := NA_real_]
E[edge_type == "miRNA_target", w_cur := -1]
E[!is.na(sgn), w_cur := sgn]
E[edge_type %in% c("gene_gene", "miRNA_miRNA"), w_cur := NA_real_]
lg("edges: total=", nrow(E), "  as_published signed=", sum(!is.na(E$w_pub)),
   " (+1=", sum(E$w_pub == 1), ", -1=", sum(E$w_pub == -1), ")")
lg("       curated signed=", sum(!is.na(E$w_cur)), " (+1=", sum(E$w_cur == 1, na.rm = TRUE),
   ", -1=", sum(E$w_cur == -1, na.rm = TRUE), "), dropped as unsignable/non-regulatory=",
   sum(is.na(E$w_cur)))

# ---- composite centrality on the unsigned directed graph (as in the deposited script)
g <- graph_from_data_frame(E[, .(source, target)], directed = TRUE,
                           vertices = nodes[, .(name, type)])
z <- function(x) as.numeric(scale(x))
CEN <- data.table(name = V(g)$name, type = V(g)$type,
                  Degree = degree(g, normalized = TRUE),
                  Betweenness = betweenness(g, directed = TRUE, normalized = TRUE),
                  Closeness = closeness(g, mode = "out", normalized = TRUE))
CEN[is.na(Closeness), Closeness := 0]
CEN[, Composite_Score := (z(Degree) + z(Betweenness) + z(Closeness)) / 3]
setorder(CEN, -Composite_Score)
seeds <- CEN$name[seq_len(min(N_SEED, nrow(CEN)))]
lg("top ", length(seeds), " composite-centrality hubs used as RWR seeds; top 10: ",
   paste(head(seeds, 10), collapse = ", "))

rwr <- function(Edt, wcol) {
  d <- Edt[!is.na(get(wcol))]
  nm <- sort(unique(c(d$source, d$target, nodes$name)))
  idx <- setNames(seq_along(nm), nm)
  W <- Matrix::sparseMatrix(i = idx[d$source], j = idx[d$target], x = d[[wcol]],
                            dims = c(length(nm), length(nm)), dimnames = list(nm, nm))
  outs <- Matrix::rowSums(abs(W)); outs[outs == 0] <- 1
  Wn <- W / outs                       # row (= source) normalised
  p0 <- setNames(rep(0, length(nm)), nm); p0[intersect(seeds, nm)] <- 1
  v <- p0
  for (i in seq_len(MAXIT)) {
    vn <- ALPHA * as.numeric(Matrix::crossprod(Wn, v)) + (1 - ALPHA) * p0
    names(vn) <- nm
    if (sum(abs(vn - v)) < TOL) { lg("   RWR (", wcol, ") converged at iter ", i); break }
    v <- vn
  }
  v
}
suppressPackageStartupMessages(library(Matrix))
sc_pub <- rwr(E, "w_pub"); sc_cur <- rwr(E, "w_cur")

SC <- data.table(name = names(sc_pub), type = ty[names(sc_pub)],
                 influence_as_published = as.numeric(sc_pub),
                 influence_curated = as.numeric(sc_cur[names(sc_pub)]))
SC[, is_seed := name %in% seeds]
setorder(SC, -influence_as_published)
fwrite(SC, file.path(REV, "results/signed_influence_scores.csv"))
lg("wrote signed_influence_scores.csv rows=", nrow(SC))
lg("Spearman correlation between the two influence rankings: ",
   round(cor(SC$influence_as_published, SC$influence_curated, method = "spearman"), 3))
lg("sign agreement (both >0 or both <0): ",
   round(100 * mean(sign(SC$influence_as_published) == sign(SC$influence_curated)), 1), "%")

# ---- GSEA on the ranked protein-coding node list
sym2ent <- function(x) {
  s <- suppressWarnings(suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(x), keytype = "SYMBOL", columns = "ENTREZID")))
  s <- s[!is.na(s$ENTREZID), ]; s[!duplicated(s$SYMBOL), ]
}
res <- list(); k <- 0
for (scheme in c("as_published", "curated")) {
  col <- if (scheme == "as_published") "influence_as_published" else "influence_curated"
  d <- SC[type %in% c("TF", "Gene")]
  m <- as.data.table(sym2ent(d$name))
  d <- merge(d, m, by.x = "name", by.y = "SYMBOL")
  d <- d[!duplicated(ENTREZID)]
  gl <- setNames(d[[col]], d$ENTREZID)
  gl <- sort(gl, decreasing = TRUE)
  lg("== GSEA scheme=", scheme, ": ranked list length=", length(gl),
     " range=[", signif(min(gl), 3), ", ", signif(max(gl), 3), "]  n_nonzero=",
     sum(gl != 0), "  n_negative=", sum(gl < 0))
  if (sum(gl != 0) < 10) { lg("   SKIPPED: fewer than 10 non-zero scores"); next }
  for (ont in c("BP", "MF", "CC")) {
    gg <- tryCatch(gseGO(geneList = gl, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
                         minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1,
                         pAdjustMethod = "BH", eps = 0, seed = TRUE, verbose = FALSE),
                   error = function(e) { lg("   gseGO ", ont, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(gg) && nrow(as.data.frame(gg)) > 0) {
      gg <- setReadable(gg, org.Hs.eg.db, "ENTREZID")
      k <- k + 1; res[[k]] <- data.table(sign_scheme = scheme, ontology = paste0("GO_", ont),
                                         as.data.frame(gg))
    }
  }
  kk <- tryCatch(gseKEGG(geneList = gl, organism = "hsa", keyType = "kegg", minGSSize = 10,
                         maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH",
                         eps = 0, seed = TRUE, verbose = FALSE),
                 error = function(e) { lg("   gseKEGG failed: ", conditionMessage(e)); NULL })
  if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
    kk <- tryCatch(setReadable(kk, org.Hs.eg.db, "ENTREZID"), error = function(e) kk)
    k <- k + 1; res[[k]] <- data.table(sign_scheme = scheme, ontology = "KEGG",
                                       as.data.frame(kk))
  }
}
if (length(res) == 0) { lg("NO GSEA RESULTS AT ALL"); quit(status = 0) }
G <- rbindlist(res, fill = TRUE)
if ("ONTOLOGY" %in% names(G)) G[, ONTOLOGY := NULL]
keep <- c("sign_scheme", "ontology", "ID", "Description", "setSize", "enrichmentScore", "NES",
          "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment")
G <- G[, intersect(keep, names(G)), with = FALSE]
G[, significant_BH_0.05 := p.adjust < 0.05]
G[, n_leading_edge := sapply(strsplit(core_enrichment, "/"), length)]
setorder(G, sign_scheme, ontology, pvalue)
fwrite(G, file.path(REV, "results/gsea_results.csv"))
lg("wrote gsea_results.csv rows=", nrow(G),
   "  BH-significant=", sum(G$significant_BH_0.05, na.rm = TRUE))
lg(paste(capture.output(print(G[, .(n_tested = .N, n_sig = sum(significant_BH_0.05)),
                                by = .(sign_scheme, ontology)])), collapse = "\n"))
lg("DONE")
