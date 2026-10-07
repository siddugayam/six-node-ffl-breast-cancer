#!/usr/bin/env Rscript
# =============================================================================
# 09b_claimed_terms_hypergeometric.R
# A targeted, filter-free test of every named functional term (from the Results and
# Discussion), for every motif set and both background universes.
#
# Why this exists in addition to 09_enrichment_ora.R: clusterProfiler's standard
# minGSSize = 10 / maxGSSize = 500 filters silently REMOVE terms from the tested set.
# "supramolecular fiber organization" (GO:0097435) annotates 871 genes and is therefore
# never tested under default settings, so its absence from a clusterProfiler table is
# not evidence of non-enrichment.  Here every named term is tested directly with the
# one-sided hypergeometric (Fisher) test that ORA is built on, with no size filter:
#     p = phyper(k - 1, n, N - n, K, lower.tail = FALSE)
#     k = genes shared by the motif set and the term (within the universe)
#     n = term size within the universe
#     K = motif set size within the universe
#     N = universe size
# BH adjustment is applied across the named terms within each motif_set x universe.
# Output: results/enrichment_claimed_terms_hypergeometric.csv
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(org.Hs.eg.db); library(GO.db); library(AnnotationDbi)
  library(clusterProfiler)
})
REV <- "/path/to/revision"
LOGF <- file.path(REV, "logs", "09b_claimed_terms.log")
lg <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOGF, append = TRUE) }
cat("", file = LOGF)

GO_TERMS <- c(
  "GO:0005583" = "fibrillar collagen trimer",
  "GO:0098643" = "banded collagen fibril",
  "GO:0062023" = "collagen-containing extracellular matrix",
  "GO:0043062" = "extracellular structure organization",
  "GO:0097435" = "supramolecular fiber organization",
  "GO:0030199" = "collagen fibril organization",
  "GO:0031047" = "regulatory ncRNA-mediated gene silencing",
  "GO:0016441" = "post-transcriptional gene silencing",
  "GO:0035195" = "miRNA-mediated post-transcriptional gene silencing",
  "GO:0019222" = "regulation of metabolic process",
  "GO:0051726" = "regulation of cell cycle",
  "GO:0010468" = "regulation of gene expression",
  "GO:0009059" = "macromolecule biosynthetic process",
  "GO:0003700" = "DNA-binding transcription factor activity",
  "GO:0140110" = "transcription regulator activity")
KEGG_TERMS <- c(
  "hsa05205" = "Proteoglycans in cancer", "hsa01523" = "Antifolate resistance",
  "hsa04926" = "Relaxin signaling pathway",
  "hsa04933" = "AGE-RAGE signaling pathway in diabetic complications",
  "hsa05206" = "MicroRNAs in cancer", "hsa04218" = "Cellular senescence",
  "hsa04110" = "Cell cycle", "hsa04115" = "p53 signaling pathway",
  "hsa05224" = "Breast cancer", "hsa05202" = "Transcriptional misregulation in cancer",
  "hsa04151" = "PI3K-Akt signaling pathway", "hsa05415" = "Diabetic cardiomyopathy",
  "hsa05418" = "Fluid shear stress and atherosclerosis")
REACTOME_TERMS <- c("R-HSA-9755511" = "KEAP1-NFE2L2 pathway",
                    "R-HSA-9759194" = "Nuclear events mediated by NFE2L2")

# ---------------------------------------------------------------- gene sets
sym2ent <- function(x) {
  s <- suppressWarnings(suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(x), keytype = "SYMBOL", columns = "ENTREZID")))
  s <- s[!is.na(s$ENTREZID), ]; s[!duplicated(s$SYMBOL), ]
}
ent2sym <- function(x) {
  s <- suppressWarnings(suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(x), keytype = "ENTREZID", columns = "SYMBOL")))
  setNames(s$SYMBOL, s$ENTREZID)
}
go_genes <- function(id) {
  g <- tryCatch(suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = id, keytype = "GOALL", columns = "ENTREZID")$ENTREZID),
    error = function(e) character(0))
  unique(g[!is.na(g)])
}
kegg <- tryCatch(clusterProfiler::download_KEGG("hsa"), error = function(e) NULL)
kegg_map <- if (!is.null(kegg)) as.data.table(kegg$KEGGPATHID2EXTID) else NULL
if (!is.null(kegg_map)) setnames(kegg_map, c("pathway", "gene"))
lg("KEGG pathway-gene rows downloaded: ", if (is.null(kegg_map)) 0 else nrow(kegg_map))
rp <- tryCatch({ suppressPackageStartupMessages(library(reactome.db))
  as.data.table(as.data.frame(reactomePATHID2EXTID)) }, error = function(e) NULL)
if (!is.null(rp)) lg("Reactome pathway-gene rows: ", nrow(rp))

term_genes <- function(ont, id) {
  if (ont == "GO") go_genes(id)
  else if (ont == "KEGG") { if (is.null(kegg_map)) character(0) else unique(kegg_map[pathway == id, gene]) }
  else { if (is.null(rp)) character(0) else unique(rp[DB_ID == id, gene_id]) }
}

ns <- fread(file.path(REV, "results/motif_node_sets.tsv"))
sets_pc <- split(ns[type %in% c("TF", "Gene"), node], ns[type %in% c("TF", "Gene"), motif_set])
ORDER <- c("3-miR", "3-TF", "3-Comp", "4-node", "5-node", "6-node",
           "exemplar_4node", "exemplar_5node", "exemplar_6node")
sets_pc <- sets_pc[ORDER]

nodes_all <- fread(file.path(REV, "data/canonical_nodes.tsv"))
de <- fread(file.path(REV, "results/BRCA_DEX_genes.csv"))
net_pc <- sort(unique(nodes_all[type %in% c("TF", "Gene"), name]))
uniA <- unique(sym2ent(unique(c(de$feature[grepl("^[A-Za-z]", de$feature)], net_pc)))$ENTREZID)
uniB <- unique(sym2ent(net_pc)$ENTREZID)
UNIS <- list(TCGA_tested_plus_network = uniA, network_protein_coding = uniB)
lg("universe sizes: A=", length(uniA), " B=", length(uniB))

TERMS <- rbindlist(list(
  data.table(ontology = "GO", ID = names(GO_TERMS), Description = unname(GO_TERMS)),
  data.table(ontology = "KEGG", ID = names(KEGG_TERMS), Description = unname(KEGG_TERMS)),
  data.table(ontology = "Reactome", ID = names(REACTOME_TERMS), Description = unname(REACTOME_TERMS))))
TERMS[, genes := lapply(seq_len(.N), function(i) term_genes(ontology[i], ID[i]))]
TERMS[, n_genes_total := sapply(genes, length)]
lg("term sizes:"); for (i in seq_len(nrow(TERMS)))
  lg("   ", TERMS$ID[i], "  ", TERMS$Description[i], "  n=", TERMS$n_genes_total[i])

out <- list(); q <- 0
for (sn in names(sets_pc)) {
  ent <- unique(sym2ent(sets_pc[[sn]])$ENTREZID)
  for (un in names(UNIS)) {
    U <- UNIS[[un]]; K <- intersect(ent, U); N <- length(U)
    for (i in seq_len(nrow(TERMS))) {
      tg <- intersect(TERMS$genes[[i]], U); n <- length(tg)
      sh <- intersect(K, tg); k <- length(sh)
      p <- if (n == 0 || length(K) == 0) NA_real_ else
        stats::phyper(k - 1, n, N - n, length(K), lower.tail = FALSE)
      sym <- if (k > 0) paste(sort(unname(ent2sym(sh)[sh])), collapse = "/") else ""
      q <- q + 1
      out[[q]] <- data.table(
        motif_set = sn, universe = un, ontology = TERMS$ontology[i], ID = TERMS$ID[i],
        Description = TERMS$Description[i], Count = k,
        GeneRatio = paste0(k, "/", length(K)), BgRatio = paste0(n, "/", N),
        term_size_in_universe = n, term_size_total = TERMS$n_genes_total[i],
        pvalue = p, geneID = sym,
        tested_under_default_filters = (TERMS$n_genes_total[i] >= 10 &
                                          TERMS$n_genes_total[i] <= 500))
    }
  }
}
R <- rbindlist(out)
R[, p.adjust := p.adjust(pvalue, "BH"), by = .(motif_set, universe)]
R[, significant_BH_0.05 := !is.na(p.adjust) & p.adjust < 0.05]
R[, low_support_lt3_genes := Count < 3]
setorder(R, Description, universe, motif_set)
fwrite(R, file.path(REV, "results/enrichment_claimed_terms_hypergeometric.csv"))
lg("wrote enrichment_claimed_terms_hypergeometric.csv rows=", nrow(R))
lg("DONE")
