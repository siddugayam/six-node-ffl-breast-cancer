#!/usr/bin/env Rscript
# =============================================================================
# 09_enrichment_ora.R
# Over-representation analysis (ORA) for every FFL motif class, with an EXPLICIT
# background universe, full statistics, and a low-support flag.
#
# Answers:
#                         "'Proteoglycans in cancer' does not appear in the bar plot"
#
# PARAMETERS
#   pAdjustMethod = "BH"
#   pvalueCutoff  = 1, qvalueCutoff = 1  at the enrich* call, so that the FULL tested
#                   term table is returned; significance is then applied downstream and
#                   recorded per row.  This is required to be able to report the actual
#                   p / p.adjust of terms the manuscript claims but which are NOT enriched.
#   minGSSize = 10, maxGSSize = 500   (term size limits; clusterProfiler defaults)
#   readable  = TRUE                  (Entrez -> HGNC symbols in geneID)
#
# BACKGROUND UNIVERSES
#   universe_A "TCGA_tested_plus_network" - the 20,250 genes tested for differential
#      expression in TCGA-BRCA (results/BRCA_DEX_genes.csv) UNION the protein-coding
#      nodes of the canonical network.  This is the universe the task specifies.
#      The 1,057-gene breast-cancer universe quoted in the manuscript is NOT recoverable:
#      it is not deposited in the authors' GitHub repository (only the 6 node-attribute
#      and 6 SIF files are), and the manuscript itself gives two contradictory
#      compositions for it (line 140: "1,057 genes including 233 TFs"; line 286:
#      "1,057 genes ... among which 665 were identified as transcription factors").
#   universe_B "network_protein_coding" - the 364 protein-coding nodes of the canonical
#      network only.  This is the curation-matched control: every gene in it was already
#      hand-picked as breast-cancer related, so enrichment for cancer terms against it is
#      NOT inflated by the curation of the input list.  Reported as a sensitivity analysis.
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(data.table); library(AnnotationDbi)
})

REV <- "/path/to/revision"
LOGF <- file.path(REV, "logs", "09_enrichment_ora.log")
lg <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOGF, append = TRUE) }
cat("", file = LOGF)
set.seed(1234)

P_ADJUST <- "BH"; MINGS <- 10; MAXGS <- 500; SIG_CUT <- 0.05; LOW_SUPPORT <- 3

# ---------------------------------------------------------------- input sets
ns <- fread(file.path(REV, "results/motif_node_sets.tsv"))
lg("motif_node_sets.tsv rows: ", nrow(ns))
sets <- split(ns$node, ns$motif_set)
ORDER <- c("3-miR", "3-TF", "3-Comp", "4-node", "5-node", "6-node",
           "exemplar_4node", "exemplar_5node", "exemplar_6node")
sets <- sets[ORDER[ORDER %in% names(sets)]]

# protein-coding only: GO/KEGG have no annotation for miRNA nodes
pc <- ns[type %in% c("TF", "Gene")]
sets_pc <- split(pc$node, pc$motif_set)[names(sets)]

# ---------------------------------------------------------------- universes
nodes_all <- fread(file.path(REV, "data/canonical_nodes.tsv"))
net_pc <- sort(unique(nodes_all[type %in% c("TF", "Gene"), name]))
de <- fread(file.path(REV, "results/BRCA_DEX_genes.csv"))
lg("BRCA_DEX_genes.csv rows: ", nrow(de))
tcga_tested <- unique(de$feature)
tcga_tested <- tcga_tested[grepl("^[A-Za-z]", tcga_tested)]   # drop 29 numeric-only ids

uniA_sym <- sort(unique(c(tcga_tested, net_pc)))
uniB_sym <- net_pc
lg("universe_A symbols: ", length(uniA_sym), "   universe_B symbols: ", length(uniB_sym))

sym2ent <- function(x) {
  s <- suppressWarnings(suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = unique(x), keytype = "SYMBOL",
                          columns = "ENTREZID")))
  s <- s[!is.na(s$ENTREZID), ]
  s[!duplicated(s$SYMBOL), ]
}
mapA <- sym2ent(uniA_sym); mapB <- sym2ent(uniB_sym)
uniA <- unique(mapA$ENTREZID); uniB <- unique(mapB$ENTREZID)
lg("universe_A entrez: ", length(uniA), " (mapped ", nrow(mapA), "/", length(uniA_sym), ")")
lg("universe_B entrez: ", length(uniB), " (mapped ", nrow(mapB), "/", length(uniB_sym), ")")

UNIS <- list(TCGA_tested_plus_network = uniA, network_protein_coding = uniB)

# ---------------------------------------------------------------- run ORA
grab <- function(obj) {
  if (is.null(obj)) return(NULL)
  d <- as.data.frame(obj)
  if (nrow(d) == 0) return(NULL)
  d
}

res <- list(); i <- 0
map_input <- list()
for (sn in names(sets_pc)) {
  syms <- unique(sets_pc[[sn]])
  m <- sym2ent(syms)
  ent <- unique(m$ENTREZID)
  map_input[[sn]] <- data.table(motif_set = sn, n_input_symbols = length(syms),
                                n_mapped_entrez = length(ent))
  lg("== ", sn, ": ", length(syms), " protein-coding nodes -> ", length(ent), " Entrez IDs")
  for (un in names(UNIS)) {
    U <- UNIS[[un]]
    ent_u <- ent[ent %in% U]
    lg("   universe ", un, " (", length(U), "): input genes inside universe = ", length(ent_u))
    if (length(ent_u) < 2) next
    for (ont in c("BP", "MF", "CC")) {
      e <- tryCatch(enrichGO(gene = ent_u, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                             ont = ont, universe = U, pAdjustMethod = P_ADJUST,
                             pvalueCutoff = 1, qvalueCutoff = 1,
                             minGSSize = MINGS, maxGSSize = MAXGS, readable = TRUE),
                    error = function(x) { lg("   ERROR GO ", ont, ": ", conditionMessage(x)); NULL })
      d <- grab(e)
      if (!is.null(d)) { i <- i + 1
        res[[i]] <- data.table(motif_set = sn, universe = un, ontology = paste0("GO_", ont), d) }
    }
    ek <- tryCatch(enrichKEGG(gene = ent_u, organism = "hsa", keyType = "kegg",
                              universe = U, pAdjustMethod = P_ADJUST,
                              pvalueCutoff = 1, qvalueCutoff = 1,
                              minGSSize = MINGS, maxGSSize = MAXGS),
                   error = function(x) { lg("   ERROR KEGG: ", conditionMessage(x)); NULL })
    if (!is.null(ek)) ek <- tryCatch(setReadable(ek, org.Hs.eg.db, "ENTREZID"),
                                     error = function(x) ek)
    d <- grab(ek)
    if (!is.null(d)) { i <- i + 1
      res[[i]] <- data.table(motif_set = sn, universe = un, ontology = "KEGG", d) }

    if (requireNamespace("ReactomePA", quietly = TRUE)) {
      er <- tryCatch(ReactomePA::enrichPathway(gene = ent_u, organism = "human",
                                               universe = U, pAdjustMethod = P_ADJUST,
                                               pvalueCutoff = 1, qvalueCutoff = 1,
                                               minGSSize = MINGS, maxGSSize = MAXGS,
                                               readable = TRUE),
                     error = function(x) { lg("   ERROR Reactome: ", conditionMessage(x)); NULL })
      d <- grab(er)
      if (!is.null(d)) { i <- i + 1
        res[[i]] <- data.table(motif_set = sn, universe = un, ontology = "Reactome", d) }
    }
  }
}

ALL <- rbindlist(res, fill = TRUE)
lg("total raw enrichment rows across all sets/universes/ontologies: ", nrow(ALL))
if ("ONTOLOGY" %in% names(ALL)) ALL[, ONTOLOGY := NULL]
keep <- c("motif_set", "universe", "ontology", "ID", "Description", "GeneRatio", "BgRatio",
          "pvalue", "p.adjust", "qvalue", "Count", "geneID")
miss <- setdiff(keep, names(ALL)); if (length(miss)) lg("MISSING COLS: ", paste(miss, collapse = ","))
ALL <- ALL[, intersect(keep, names(ALL)), with = FALSE]
ALL[, significant_BH_0.05 := p.adjust < SIG_CUT]
ALL[, low_support_lt3_genes := Count < LOW_SUPPORT]

fwrite(ALL, file.path(REV, "results/enrichment_all_motifs_FULL.csv"))
lg("wrote enrichment_all_motifs_FULL.csv rows=", nrow(ALL))

# ---------------------------------------------------------- manuscript-claimed terms
CLAIMED <- c(
  "fibrillar collagen trimer", "banded collagen fibril",
  "Proteoglycans in cancer", "Antifolate resistance", "Relaxin signaling pathway",
  "AGE-RAGE signaling pathway in diabetic complications", "KEAP1-NFE2L2 pathway",
  "MicroRNAs in cancer", "Cellular senescence", "Cell cycle", "p53 signaling pathway",
  "Breast cancer", "Transcriptional misregulation in cancer", "PI3K-Akt signaling pathway",
  "Diabetic cardiomyopathy", "Fluid shear stress and atherosclerosis",
  "collagen-containing extracellular matrix", "extracellular structure organization",
  "supramolecular fiber organization", "collagen fibril organization",
  # NOTE: "gene silencing" (GO:0016458), the term printed in the manuscript's Figure 7B,
  # is OBSOLETE in the current GO release and cannot be tested. Its surviving relatives
  # are included instead.
  "gene silencing", "regulatory ncRNA-mediated gene silencing",
  "post-transcriptional gene silencing", "miRNA-mediated post-transcriptional gene silencing",
  "regulation of metabolic process", "regulation of cell cycle",
  "regulation of gene expression", "macromolecule biosynthetic process",
  "DNA-binding transcription factor activity", "transcription regulator activity",
  "Nuclear events mediated by NFE2L2")
ALL[, manuscript_claimed_term := Description %in% CLAIMED]

# main deliverable: significant terms + every manuscript-claimed term (so that terms the
# manuscript names but which are NOT enriched are still reported, with their real p)
MAIN <- ALL[significant_BH_0.05 == TRUE | manuscript_claimed_term == TRUE]
setorder(MAIN, motif_set, universe, ontology, pvalue)
fwrite(MAIN, file.path(REV, "results/enrichment_all_motifs.csv"))
lg("wrote enrichment_all_motifs.csv rows=", nrow(MAIN))

CLAIM <- ALL[manuscript_claimed_term == TRUE]
setorder(CLAIM, Description, motif_set, universe)
fwrite(CLAIM, file.path(REV, "results/enrichment_manuscript_claimed_terms.csv"))
lg("wrote enrichment_manuscript_claimed_terms.csv rows=", nrow(CLAIM))
lg("claimed terms NEVER found in any tested term table: ",
   paste(setdiff(CLAIMED, unique(ALL$Description)), collapse = " | "))

# ---------------------------------------------------------- low-support summary
SUMM <- ALL[significant_BH_0.05 == TRUE,
            .(n_sig_terms = .N,
              n_sig_lt3_genes = sum(Count < LOW_SUPPORT),
              pct_sig_lt3_genes = round(100 * mean(Count < LOW_SUPPORT), 1),
              n_sig_1gene = sum(Count == 1),
              median_Count = as.numeric(median(Count)),
              min_padj = min(p.adjust)),
            by = .(motif_set, universe, ontology)]
setorder(SUMM, motif_set, universe, ontology)
fwrite(SUMM, file.path(REV, "results/enrichment_low_support_summary.csv"))
lg("wrote enrichment_low_support_summary.csv rows=", nrow(SUMM))

fwrite(rbindlist(map_input), file.path(REV, "results/enrichment_input_sizes.csv"))
lg("DONE")
