#!/usr/bin/env Rscript
# =============================================================================
# v7_go_01_ora.R   BiNGO-equivalent GO-BP + KEGG over-representation for every
#                  module of size >= 5 recovered by the four v7 methods.
#
# CONVENTIONS - identical to scripts/09_enrichment_ora.R, the manuscript's own:
#   pAdjustMethod = "BH"
#   pvalueCutoff = 1, qvalueCutoff = 1 at the enrich* call, so the FULL tested
#     term table is returned; significance (p.adjust < 0.05) is applied
#     downstream and recorded per row.  This is what lets a NEGATIVE result be
#     reported with its actual p rather than as an empty table.
#   minGSSize = 10, maxGSSize = 500
#   Reporting floor: Count >= 3 genes in the term (LOW_SUPPORT in script 09).
#
# BACKGROUND (universe_A, as pre-specified): the 20,250 genes tested for differential
#   expression in TCGA-BRCA (results/v2/v2_DE_genes.csv), NOT the genome.
# SENSITIVITY (universe_B): the 364 protein-coding nodes of the canonical
#   network.  Every gene in it was already hand-picked as breast-cancer
#   related, so this is the curation-matched control: it asks whether a module
#   is enriched RELATIVE TO THE NETWORK IT WAS CUT OUT OF, which is the only
#   version of the question that is not answered by the curation itself.
#
# miRNA members carry no GO BP / KEGG annotation and are dropped before ORA.
#
# KEGG: clusterProfiler::download_KEGG("hsa") is called ONCE, cached to
#   cache/v7/kegg_hsa_<date>.rds, and the identical hypergeometric test is run
#   through enricher() with that TERM2GENE.  This is arithmetically the same as
#   enrichKEGG() but avoids 128 live REST calls and pins the KEGG release.
# =============================================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
  library(data.table)
})

REV <- "/path/to/revision"
V7  <- file.path(REV, "results/v7")
CACHE <- file.path(REV, "cache/v7"); dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
LOGF <- file.path(REV, "logs", "v7_go_01_ora.log")
lg <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOGF, append = TRUE) }
cat("", file = LOGF)
set.seed(1234)

P_ADJUST <- "BH"; MINGS <- 10; MAXGS <- 500; SIG_CUT <- 0.05; GENE_FLOOR <- 3

## --------------------------------------------------------------- universes
de <- fread(file.path(REV, "results/v2/v2_DE_genes.csv"))
uniA_sym <- unique(de$feature)
lg("DE table rows: ", nrow(de), "  unique symbols: ", length(uniA_sym))

map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = uniA_sym,
        keytype = "SYMBOL", columns = "ENTREZID"))
map <- as.data.table(map)[!is.na(ENTREZID)]
map <- map[!duplicated(SYMBOL)]
UNIV_A <- unique(map$ENTREZID)
SYM2EG <- setNames(map$ENTREZID, map$SYMBOL)
lg("universe_A (DE-tested genes) : ", length(uniA_sym), " symbols -> ",
   length(UNIV_A), " Entrez IDs (", length(uniA_sym) - length(UNIV_A), " unmapped)")

nodes <- fread(file.path(V7, "cm2_nodes.tsv")); setnames(nodes, c("name", "type"))
netpc <- nodes[type %in% c("Gene", "TF")]$name
# network protein-coding nodes not in the DE table still belong in universe_B
map2 <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = netpc,
         keytype = "SYMBOL", columns = "ENTREZID"))
map2 <- as.data.table(map2)[!is.na(ENTREZID)][!duplicated(SYMBOL)]
UNIV_B <- unique(map2$ENTREZID)
SYM2EG <- c(SYM2EG, setNames(map2$ENTREZID, map2$SYMBOL)[setdiff(map2$SYMBOL, names(SYM2EG))])
lg("universe_B (network protein-coding): ", length(netpc), " nodes -> ", length(UNIV_B), " Entrez IDs")
lg("network PC nodes absent from universe_A: ",
   paste(setdiff(map2$ENTREZID, UNIV_A), collapse = ","), " (added to universe_B only)")

EG2SYM <- setNames(names(SYM2EG), SYM2EG)

## ------------------------------------------------------------------- KEGG
kf <- file.path(CACHE, "kegg_hsa_download.rds")
if (file.exists(kf)) { KG <- readRDS(kf); lg("KEGG: loaded cache ", kf) } else {
  KG <- clusterProfiler::download_KEGG("hsa", keggType = "KEGG", keyType = "kegg")
  attr(KG, "downloaded") <- Sys.time(); saveRDS(KG, kf); lg("KEGG: downloaded and cached ", kf)
}
K2G <- as.data.frame(KG$KEGGPATHID2EXTID); names(K2G) <- c("term", "gene")
K2N <- as.data.frame(KG$KEGGPATHID2NAME);  names(K2N) <- c("term", "name")
lg("KEGG: ", length(unique(K2G$term)), " hsa pathways, ", nrow(K2G), " pathway-gene pairs")

## ---------------------------------------------------------------- modules
inv <- fread(file.path(V7, "goora_00_module_inventory_size5.csv"))
lg("modules of size >= 5: ", nrow(inv))

run_one <- function(eg, univ, univ_name) {
  out <- list()
  go <- tryCatch(enrichGO(gene = eg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                          universe = univ, pAdjustMethod = P_ADJUST,
                          pvalueCutoff = 1, qvalueCutoff = 1,
                          minGSSize = MINGS, maxGSSize = MAXGS, readable = TRUE),
                 error = function(e) NULL)
  if (!is.null(go) && nrow(as.data.frame(go))) {
    d <- as.data.table(as.data.frame(go)); d[, ontology := "GO:BP"]; out[[1]] <- d
  }
  kk <- tryCatch(enricher(gene = eg, TERM2GENE = K2G, TERM2NAME = K2N, universe = univ,
                          pAdjustMethod = P_ADJUST, pvalueCutoff = 1, qvalueCutoff = 1,
                          minGSSize = MINGS, maxGSSize = MAXGS),
                 error = function(e) NULL)
  if (!is.null(kk) && nrow(as.data.frame(kk))) {
    d <- as.data.table(as.data.frame(kk)); d[, ontology := "KEGG"]
    d[, geneID := vapply(strsplit(geneID, "/"), function(v)
        paste(ifelse(is.na(EG2SYM[v]), v, EG2SYM[v]), collapse = "/"), character(1))]
    out[[2]] <- d
  }
  if (!length(out)) return(NULL)
  r <- rbindlist(out, fill = TRUE); r[, universe := univ_name]; r[]
}

res <- list(); diag <- list()
for (i in seq_len(nrow(inv))) {
  m <- inv[i]
  mem <- strsplit(m$members, ";")[[1]]
  pc  <- mem[nodes[match(mem, name)]$type %in% c("Gene", "TF")]
  eg  <- unique(na.omit(SYM2EG[pc]))
  egA <- intersect(eg, UNIV_A); egB <- intersect(eg, UNIV_B)
  diag[[i]] <- data.table(method = m$method, partition = m$partition, tier = m$tier,
                          module = m$module, size = m$size, n_geneTF = m$n_geneTF,
                          n_mapped_entrez = length(eg),
                          n_in_universeA = length(egA), n_in_universeB = length(egB),
                          unmapped_symbols = paste(setdiff(pc, names(SYM2EG)), collapse = ";"))
  if (length(egA) >= 3) {
    r <- run_one(egA, UNIV_A, "A_DEtested_20250")
    if (!is.null(r)) { r[, `:=`(method = m$method, partition = m$partition, tier = m$tier,
                                module = m$module, module_size = m$size,
                                n_geneTF = m$n_geneTF, n_tested = length(egA))]
                       res[[length(res) + 1L]] <- r }
  }
  if (length(egB) >= 3) {
    r <- run_one(egB, UNIV_B, "B_network364")
    if (!is.null(r)) { r[, `:=`(method = m$method, partition = m$partition, tier = m$tier,
                                module = m$module, module_size = m$size,
                                n_geneTF = m$n_geneTF, n_tested = length(egB))]
                       res[[length(res) + 1L]] <- r }
  }
  if (i %% 10 == 0) lg("  ...", i, "/", nrow(inv))
}

ALL <- rbindlist(res, fill = TRUE)
D   <- rbindlist(diag)
ALL[, sig := p.adjust < SIG_CUT & Count >= GENE_FLOOR]
setcolorder(ALL, c("method", "partition", "tier", "module", "module_size", "n_geneTF",
                   "n_tested", "universe", "ontology", "ID", "Description",
                   "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count",
                   "geneID", "sig"))
fwrite(ALL, file.path(V7, "goora_01_ORA_full_term_table.csv"))
fwrite(ALL[sig == TRUE][order(method, partition, module, universe, ontology, p.adjust)],
       file.path(V7, "goora_01_ORA_significant_terms.csv"))
fwrite(D, file.path(V7, "goora_01_module_gene_mapping.csv"))

lg("")
lg("term rows tested (all modules x both universes x GO+KEGG): ", nrow(ALL))
lg("significant (BH < 0.05 AND Count >= 3): ", sum(ALL$sig))
lg("")
s <- ALL[universe == "A_DEtested_20250", .(n_terms_tested = .N,
        n_sig = sum(sig), min_padj = min(p.adjust)),
        by = .(method, partition, module, module_size, n_tested)]
setorder(s, method, partition, -module_size)
print(s, nrows = 200)
