#!/usr/bin/env Rscript
## ============================================================================
## 11_exir_mediator_enrichment.R
## EXPLICIT BACKGROUND UNIVERSE = every gene tested in the DE analysis
## (not the default whole-genome background).
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
})
ROOT <- "/path/to/revision"
LOG  <- file.path(ROOT, "logs", "exir_agent.log")
say <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | 11 | ", paste0(..., collapse = ""))
  cat(msg, "\n", file = LOG, append = TRUE); cat(msg, "\n"); flush.console()
}
say("=== 11 START === clusterProfiler ", as.character(packageVersion("clusterProfiler")),
    " org.Hs.eg.db ", as.character(packageVersion("org.Hs.eg.db")))

de  <- fread(file.path(ROOT, "results", "BRCA_DEX_ALL_nodes.csv"))
med <- fread(file.path(ROOT, "results", "exir_mediators.csv"))

genes_tested <- de[class == "gene", Gene]
say("gene features tested in DE: ", length(genes_tested))
de_sig <- de[class == "gene" & abs(logFC) > 1 & adj.P.Val < 0.05, Gene]
say("DE-significant genes (|logFC|>1 & adj.P<0.05): ", length(de_sig))

## mediators are non-DE by construction; keep only the gene-level ones
med_genes <- med[!grepl("^hsa-", feature), feature]
say("nonDE mediators total ", nrow(med), " ; gene-level (non-miRNA) mediators: ", length(med_genes))
med_top   <- med[!grepl("^hsa-", feature)][order(Rank)][1:min(100, .N), feature]
med_top500<- med[!grepl("^hsa-", feature)][order(Rank)][1:min(500, .N), feature]
say("top-100 gene-level mediators by ExIR rank: ", length(med_top),
    " ; top-500: ", length(med_top500))

## ------------------------------------------------------------- ENTREZ mapping
map_entrez <- function(sym) {
  s <- unique(sym[!grepl("^[0-9]+$", sym)])          # drop unmapped numeric TCGA ids
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = s,
        keytype = "SYMBOL", columns = "ENTREZID"))
  unique(na.omit(m$ENTREZID))
}
UNIV <- map_entrez(genes_tested)
say("background universe mapped to ENTREZ: ", length(UNIV), " of ", length(genes_tested))
S_de   <- intersect(map_entrez(de_sig), UNIV)
S_med  <- intersect(map_entrez(med_genes), UNIV)
S_medT <- intersect(map_entrez(med_top),  UNIV)
S_medT5<- intersect(map_entrez(med_top500), UNIV)
S_both     <- union(S_de, S_med)
S_bothTop  <- union(S_de, S_medT)
S_bothTop5 <- union(S_de, S_medT5)
say("ENTREZ sets: DE=", length(S_de), " mediators(all genes)=", length(S_med),
    " mediators(top100)=", length(S_medT), " mediators(top500)=", length(S_medT5),
    " DE+mediators=", length(S_both), " DE+top100mediators=", length(S_bothTop),
    " DE+top500mediators=", length(S_bothTop5))

runGO <- function(g, tag) {
  e <- enrichGO(gene = g, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", universe = UNIV, pAdjustMethod = "BH",
                pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(data.table())
  d <- as.data.table(as.data.frame(e)); d[, `:=`(ontology = "GO_BP", gene_set = tag)]; d
}
runKEGG <- function(g, tag) {
  e <- tryCatch(enrichKEGG(gene = g, organism = "hsa", universe = UNIV,
                           pAdjustMethod = "BH", pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05),
                error = function(err) { say("KEGG failed for ", tag, ": ", conditionMessage(err)); NULL })
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(data.table())
  d <- as.data.table(as.data.frame(e)); d[, `:=`(ontology = "KEGG", gene_set = tag)]; d
}

res <- rbindlist(list(
  runGO(S_de,        "DE_only"),
  runGO(S_both,      "DE_plus_all_mediators"),
  runGO(S_bothTop,   "DE_plus_top100_mediators"),
  runGO(S_bothTop5,  "DE_plus_top500_mediators"),
  runKEGG(S_de,      "DE_only"),
  runKEGG(S_both,    "DE_plus_all_mediators"),
  runKEGG(S_bothTop, "DE_plus_top100_mediators"),
  runKEGG(S_bothTop5,"DE_plus_top500_mediators")
), use.names = TRUE, fill = TRUE)
say("total enrichment rows: ", nrow(res))

for (o in unique(res$ontology)) {
  for (s in unique(res$gene_set)) {
    n <- res[ontology == o & gene_set == s, .N]
    say("  ", o, " / ", s, " : ", n, " significant terms (p.adjust<0.05 & q<0.05)")
  }
}

## ------------------------------------------------------------- GAINED TERMS
res[, gained_vs_DE_only := FALSE]
for (o in unique(res$ontology)) {
  base_ids <- res[ontology == o & gene_set == "DE_only", ID]
  res[ontology == o & gene_set != "DE_only" & !(ID %in% base_ids), gained_vs_DE_only := TRUE]
  for (s in setdiff(unique(res$gene_set), "DE_only")) {
    ng <- res[ontology == o & gene_set == s & gained_vs_DE_only == TRUE, .N]
    nl <- length(setdiff(base_ids, res[ontology == o & gene_set == s, ID]))
    say("  ", o, " / ", s, " : GAINED ", ng, " terms vs DE_only ; LOST ", nl, " terms")
  }
}
setorder(res, ontology, gene_set, p.adjust)
fwrite(res, file.path(ROOT, "results", "exir_mediator_enrichment.csv"))
say("wrote results/exir_mediator_enrichment.csv rows=", nrow(res))

for (o in unique(res$ontology)) {
  for (s2 in c("DE_plus_top100_mediators","DE_plus_top500_mediators","DE_plus_all_mediators")) {
    g <- res[ontology == o & gene_set == s2 & gained_vs_DE_only == TRUE][order(p.adjust)]
    say("TOP GAINED ", o, " terms (", s2, "): ",
        if (nrow(g)) paste(sprintf("%s|%s|padj=%.3g", g$ID[1:min(20,nrow(g))],
            g$Description[1:min(20,nrow(g))], g$p.adjust[1:min(20,nrow(g))]), collapse = " ;; ") else "NONE")
  }
}
say("=== 11 DONE ===")
