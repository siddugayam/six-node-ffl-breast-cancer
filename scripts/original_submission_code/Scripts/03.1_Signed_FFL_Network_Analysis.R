# ==============================================================================
# SCRIPT 02: Signed Feed-Forward Loop (FFL) Network Analysis — Main Pipeline
# ==============================================================================
#
# PURPOSE:
#   Performs multi-layered transcriptional regulatory network analysis for
#   breast cancer. Constructs 3- to 6-node Feed-Forward Loops (FFLs),
#   computes centrality-based hub rankings, runs Signed Random Walk with
#   Restart (RWR) propagation, performs functional enrichment, and validates
#   predicted hubs against TCGA-BRCA experimental expression data.
#
# WORKFLOW POSITION:
#   [01_BRCA_miRNA_DESeq2_Preprocessing.R]
#            |  (produces BRCA_miRNA_DEX_results.csv)
#            v
#   [02_Signed_FFL_Network_Analysis.R]      <-- YOU ARE HERE
#
# INPUTS (all must be in working directory):
#   {PREFIX}.sif        — Network edge list (whitespace-separated, 3 columns:
#                         source, interaction_type, target)
#   {PREFIX}.csv        — Node attribute table (must have 'name' and 'type' columns)
#   BRCA_DEX_ALL_nodes.csv  — Combined gene + miRNA DEX table from Script 01
#                             (columns: Gene, logFC, adj.P.Val)
#   interactions.tsv    — [Optional] DGIdb drug-gene interaction file
#   miRWalk_Pathways_mature.gmt — [Optional] miRNA pathway GMT for unified GSEA
#
# OUTPUTS (saved to {PREFIX}_FFL_Integrated_RESULTS/):
#   01_Top_Hubs.csv                  — Top N hubs ranked by Composite Centrality
#   02_Propagation_Scores.csv        — Signed Influence Scores for all nodes
#   03_Network_Modules.csv           — Louvain community assignments
#   04_Drug_Targets.csv              — Druggable hub–drug pairs (if DGIdb present)
#   {PREFIX}_FFL_triplets_*.csv      — Detected FFL motifs (miRNA–TF–Gene)
#   Enr_Global_Reactome_ORA.csv/.png — Reactome ORA plot and table
#   Enr_Global_GO_Dotplot.png        — GO ORA plot (BP/MF/CC)
#   Enr_Global_KEGG_Dotplot.png      — KEGG ORA plot
#   Enr_Global_GO_Ridge.png          — GO GSEA ridgeplot
#   Enr_Global_Reactome_Ridge.png    — Reactome GSEA ridgeplot
#   Enr_Unified_Dotplot.png          — Unified GSEA (genes + miRNAs via msigdbr)
#   Enr_FFL_Reactome_*.csv/.png      — FFL-subset Reactome enrichment
#   Final_Validation_Barplot.png     — Hub LogFC barplot (TCGA-BRCA)
#   Final_Validation_Table.csv       — Merged hub + DEX data (sig. filtered)
#
# NETWORK TYPE CONFIGURATION:
#   Set NETWORK_PREFIX to match your SIF/CSV filename stem:
#     "3-miR"  → 3-node miRNA-FFL (miRNA regulates TF and Gene)
#     "3-TF"   → 3-node TF-FFL   (TF regulates miRNA and Gene)
#     "3-Comp" → 3-node Composite-FFL (TF ↔ miRNA, both target Gene)
#     "4-TF"   → 4-node TF-FFL   (adds gene-gene interaction)
#     "5-TF"   → 5-node TF-FFL   (adds miRNA-miRNA interaction)
#     "6TF"    → 6-node TF-FFL   (adds TF-TF interaction)
#
# ANALYSIS OVERVIEW:
#   1. Network construction from SIF + node attribute files
#   2. FFL motif detection (type-aware: miRNA-FFL / TF-FFL / Composite)
#   3. Topological centrality (degree, betweenness, closeness, composite score)
#   4. Louvain community detection
#   5. Signed RWR propagation from hub seeds (α = 0.85)
#   6. Functional enrichment: GO ORA, KEGG ORA, Reactome ORA, GO GSEA,
#      Reactome GSEA, and unified GSEA (genes + miRNAs)
#   7. Drug-gene interaction screening (DGIdb, optional)
#   8. TCGA-BRCA validation: hub expression vs. differential expression data
#
# KEY PARAMETERS (see SECTION 0 below):
#   NUM_TOP_GENES_CENTRALITY     — Number of top hubs to report (default: 60)
#   NUM_HUBS_TO_TEST_PROPAGATION — Seeds for RWR (default: 60)
#   NUM_HUBS_TO_REPORT_PER_TYPE  — Hubs per node type for validation (default: 20)
#   RWR_ALPHA                    — RWR restart probability (default: 0.85)
#   ORA_TOP_N                    — Genes submitted to ORA (default: 200)
#   LOGFC_THRESHOLD              — |logFC| cutoff for validation (default: 1.0)
#   PVAL_THRESHOLD               — adj.p cutoff for validation (default: 0.05)
#
# DEPENDENCIES:
#   Bioconductor: clusterProfiler, org.Hs.eg.db, ReactomePA, DOSE, enrichplot
#   CRAN:         tidyverse, igraph, ggraph, patchwork, tidygraph, scales,
#                 ggrepel, data.table, msigdbr
#
# AUTHORS:  Gayam Prasanna Kumar Reddy, Jesil Mathew A, Fayaz Shaik Mahammad
# VERSION:  3.0 (cleaned — bug fixes: ORA top-N filter, DEX significance
#           thresholds, FFL-aware seed selection, validation significance markers)
# REFERENCE: Bhat et al. (2023) EurekaSelect; TCGA-BRCA via UCSC Xena
# ==============================================================================


# ==============================================================================
# SECTION 0: User Configuration  ← EDIT THESE VALUES
# ==============================================================================

# Network prefix — must match your .sif and .csv file names exactly
NETWORK_PREFIX <- "5-TF"    # e.g. "3-miR", "3-TF", "3-Comp", "4-TF", "5-TF", "6TF"

# External input files
DEX_FILE_NAME  <- "BRCA_DEX_ALL_nodes.csv"         # Combined gene + miRNA DEX (Script 01 output)
DGIDB_FILE     <- "interactions.tsv"               # [Optional] DGIdb drug-gene interactions
MIRNA_GMT_FILE <- "miRWalk_Pathways_mature.gmt"    # [Optional] miRNA pathway GMT for unified GSEA

# Centrality & propagation parameters
NUM_TOP_GENES_CENTRALITY     <- 60    # Top N hubs to write to 01_Top_Hubs.csv
NUM_HUBS_TO_TEST_PROPAGATION <- 60    # Number of seeds for RWR
NUM_HUBS_TO_REPORT_PER_TYPE  <- 20    # Hubs per node type used in validation
RWR_ALPHA                    <- 0.85  # Restart probability for RWR (0 = pure random walk)

# Enrichment parameters
ORA_TOP_N     <- 200    # Top N genes (by |Influence Score|) submitted to ORA
                        # Note: paper uses top 100; set to 200 for broader coverage

# Validation significance thresholds (paper: |logFC| > 1.0, padj < 0.05)
LOGFC_THRESHOLD <- 1.0
PVAL_THRESHOLD  <- 0.05


# ==============================================================================
# SECTION 1: Package Installation and Loading
# ==============================================================================

cat("=== RUNNING INTEGRATED PIPELINE FOR:", NETWORK_PREFIX, "===\n\n")
cat("--- [1/10] Installing and loading packages ---\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "DOSE", "enrichplot")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

cran_packages <- c("tidyverse", "igraph", "ggraph", "patchwork",
                   "tidygraph", "scales", "ggrepel", "data.table", "msigdbr")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(patchwork)
  library(tidygraph)
  library(scales)
  library(ggrepel)
  library(data.table)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  library(ReactomePA)
  library(org.Hs.eg.db)
})

# Create output directory
out_dir <- paste0(NETWORK_PREFIX, "_FFL_Integrated_RESULTS")
dir.create(out_dir, showWarnings = FALSE)
cat("Output directory:", out_dir, "\n\n")


# ==============================================================================
# SECTION 2: Core Functions
# ==============================================================================

cat("--- [2/10] Defining core functions ---\n")

# ------------------------------------------------------------------------------
# signed_influence_propagation()
#   Implements Signed Random Walk with Restart (RWR).
#   Propagates regulatory signals through the signed, directed network from
#   a set of seed nodes. Edges from miRNAs carry weight -1 (inhibition),
#   edges from TFs carry weight +1 (activation).
#
# Args:
#   graph          — igraph object with a 'weight_signed' edge attribute
#   initial_vector — Named numeric vector of seed scores (1 = hub seed, 0 = other)
#   alpha          — Restart probability (default 0.85; higher = stays near seeds)
#   max_iter       — Maximum iterations before forced stop
#   tolerance      — Convergence threshold (L1 norm of score change)
#
# Returns:
#   Named numeric vector of Signed Influence Scores for all nodes.
#   Positive score → predicted net activation; negative → net repression.
# ------------------------------------------------------------------------------
signed_influence_propagation <- function(graph,
                                         initial_vector,
                                         alpha     = 0.85,
                                         max_iter  = 100,
                                         tolerance = 1e-6) {
  A_sparse <- igraph::as_adjacency_matrix(
    graph, attr = "weight_signed", type = "both", sparse = TRUE
  )
  A <- as.matrix(A_sparse)

  # Column-normalize by absolute column sum (handles negative weights)
  D_inv <- diag(1 / colSums(abs(A)))
  D_inv[!is.finite(D_inv)] <- 0
  A_norm <- A %*% D_inv

  V_curr <- initial_vector
  V_init <- initial_vector

  for (i in seq_len(max_iter)) {
    V_next <- alpha * (A_norm %*% V_curr) + (1 - alpha) * V_init
    if (sum(abs(V_next - V_curr)) < tolerance) {
      cat("   RWR converged at iteration", i, "\n")
      break
    }
    V_curr <- V_next
  }

  result_vec        <- as.vector(V_curr)
  names(result_vec) <- rownames(A)   # preserve node names through matrix math
  return(result_vec)
}


# ------------------------------------------------------------------------------
# symbol_to_entrez_df()
#   Maps HGNC gene symbols to Entrez IDs using org.Hs.eg.db.
#   Silently drops symbols that cannot be mapped (e.g., non-coding, aliases).
#
# Args:
#   symbols — Character vector of gene symbols
#
# Returns:
#   Data frame with columns SYMBOL and ENTREZID.
# ------------------------------------------------------------------------------
symbol_to_entrez_df <- function(symbols) {
  ids <- tryCatch(
    bitr(unique(symbols),
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = org.Hs.eg.db),
    error = function(e) NULL
  )
  if (is.null(ids)) return(data.frame())
  ids
}


# ------------------------------------------------------------------------------
# clean_name_for_match()
#   Standardises node names for matching between the network and DEX tables.
#   miRNA names (hsa-miR-XX or hsa-let-XX) are upper-cased and truncated to
#   the first three dash-delimited fields; gene/TF names are simply upper-cased.
#
# Args:
#   x — Single character string (node name)
#
# Returns:
#   Cleaned, upper-cased character string.
# ------------------------------------------------------------------------------
clean_name_for_match <- function(x) {
  parts <- str_split(x, "-")[[1]]
  if (length(parts) >= 3 &&
      grepl("^(hsa|HSA)", parts[1], ignore.case = TRUE)) {
    return(toupper(paste(parts[1], parts[2], parts[3], sep = "-")))
  }
  return(toupper(x))
}


# ==============================================================================
# SECTION 3: Data Import and Graph Construction
# ==============================================================================

cat("--- [3/10] Loading network files ---\n")

sif_file <- paste0(NETWORK_PREFIX, ".sif")
att_file <- paste0(NETWORK_PREFIX, ".csv")

if (!file.exists(sif_file)) stop(paste("SIF file not found:", sif_file))
if (!file.exists(att_file)) stop(paste("Attribute file not found:", att_file))

# --- Parse SIF edge list ---
# Format: source <whitespace> interaction_type <whitespace> target
raw_lines <- readLines(sif_file)
raw_lines <- raw_lines[nzchar(trimws(raw_lines))]   # drop blank lines
sif_list  <- strsplit(raw_lines, "\\s+")

edges <- data.frame(
  from = sapply(sif_list, `[`, 1),
  to   = sapply(sif_list, function(x) tail(x, 1)),
  stringsAsFactors = FALSE
) %>%
  dplyr::distinct()

# --- Parse node attribute table ---
nodes_raw <- read.csv(att_file, check.names = FALSE)
cols      <- colnames(nodes_raw)

# Accept columns named 'name', 'shared name', or 'shared.name' (Cytoscape variants)
name_col <- grep("^(shared[ ._]?)?name$", cols, ignore.case = TRUE, value = TRUE)[1]
type_col <- grep("^(shared[ ._]?)?type$", cols, ignore.case = TRUE, value = TRUE)[1]

if (is.na(name_col) || is.na(type_col)) {
  stop("Could not find 'name' or 'type' columns in the attribute CSV. ",
       "Found columns: ", paste(cols, collapse = ", "))
}

nodes_clean <- nodes_raw %>%
  dplyr::select(name = all_of(name_col), Type = all_of(type_col)) %>%
  dplyr::mutate(
    name = trimws(as.character(name)),
    Type = trimws(as.character(Type))
  ) %>%
  dplyr::distinct(name, .keep_all = TRUE)

# Add any edge endpoints missing from the attribute file as 'Unknown'
missing_nodes <- setdiff(unique(c(edges$from, edges$to)), nodes_clean$name)
if (length(missing_nodes) > 0) {
  cat(paste("   Adding", length(missing_nodes), "nodes missing from attribute file as 'Unknown'\n"))
  nodes_clean <- bind_rows(
    nodes_clean,
    tibble(name = missing_nodes, Type = "Unknown")
  )
}

# --- Build undirected graph; retain only the largest connected component ---
g_full     <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes_clean)
components <- igraph::decompose(g_full)
g_comp     <- components[[which.max(sapply(components, vcount))]]

cat(paste("Network built — Nodes:", vcount(g_comp),
          "| Edges:", ecount(g_comp), "\n\n"))


# ==============================================================================
# SECTION 4: FFL Motif Detection
# ==============================================================================

cat("--- [4/10] Detecting FFL motifs ---\n")

# Helper: return neighbour names of a named vertex within a graph
get_neighbors_by_name <- function(graph, vname) {
  if (!(vname %in% V(graph)$name)) return(character(0))
  V(graph)[neighbors(graph, V(graph)[name == vname])]$name
}

# Classify nodes by type within the largest component
comp_nodes <- tibble(name = V(g_comp)$name) %>%
  left_join(nodes_clean, by = "name")

mirnas <- comp_nodes %>% filter(grepl("mirna",  Type, ignore.case = TRUE)) %>% pull(name)
tfs    <- comp_nodes %>% filter(grepl("^tf$",   Type, ignore.case = TRUE)) %>% pull(name)
genes  <- comp_nodes %>% filter(grepl("^gene$", Type, ignore.case = TRUE)) %>% pull(name)

cat("Node counts — miRNA:", length(mirnas),
    "| TF:", length(tfs),
    "| Gene:", length(genes), "\n")

# Determine FFL detection mode from prefix
ffl_mode <- dplyr::case_when(
  grepl("3-miR",  NETWORK_PREFIX, ignore.case = TRUE) ~ "miRNA",
  grepl("3-Comp", NETWORK_PREFIX, ignore.case = TRUE) ~ "Comp",
  grepl("TF",     NETWORK_PREFIX, ignore.case = TRUE) ~ "TF",
  TRUE                                                 ~ "miRNA"
)
cat("FFL detection mode:", ffl_mode, "\n")

# ------------------------------------------------------------------------------
# detect_miRNA_FFL():
#   Finds miRNA-FFLs: miRNA → TF, miRNA → Gene, TF → Gene
#   The miRNA is the master regulator that acts on both a TF and a shared Gene.
# ------------------------------------------------------------------------------
detect_miRNA_FFL <- function() {
  res <- list()
  for (mi in mirnas) {
    mi_nei  <- get_neighbors_by_name(g_comp, mi)
    mi_tf   <- intersect(mi_nei, tfs)
    mi_gene <- intersect(mi_nei, genes)
    if (length(mi_tf) == 0 || length(mi_gene) == 0) next

    for (tf in mi_tf) {
      tf_nei       <- get_neighbors_by_name(g_comp, tf)
      common_genes <- intersect(mi_gene, intersect(tf_nei, genes))
      for (gn in common_genes) {
        res[[length(res) + 1]] <- data.frame(miRNA = mi, TF = tf, Gene = gn,
                                              stringsAsFactors = FALSE)
      }
    }
  }
  bind_rows(res)
}

# ------------------------------------------------------------------------------
# detect_TF_FFL():
#   Finds TF-FFLs: TF → miRNA, TF → Gene, miRNA → Gene
#   The TF is the master regulator that drives both a miRNA and a shared Gene.
# ------------------------------------------------------------------------------
detect_TF_FFL <- function() {
  res <- list()
  for (tf in tfs) {
    tf_nei  <- get_neighbors_by_name(g_comp, tf)
    tf_mi   <- intersect(tf_nei, mirnas)
    tf_gene <- intersect(tf_nei, genes)
    if (length(tf_mi) == 0 || length(tf_gene) == 0) next

    for (mi in tf_mi) {
      mi_nei       <- get_neighbors_by_name(g_comp, mi)
      common_genes <- intersect(tf_gene, intersect(mi_nei, genes))
      for (gn in common_genes) {
        res[[length(res) + 1]] <- data.frame(miRNA = mi, TF = tf, Gene = gn,
                                              stringsAsFactors = FALSE)
      }
    }
  }
  bind_rows(res)
}

# Run detection based on mode
triplets <- switch(ffl_mode,
  miRNA = detect_miRNA_FFL(),
  TF    = detect_TF_FFL(),
  Comp  = bind_rows(detect_miRNA_FFL(), detect_TF_FFL()),
          detect_miRNA_FFL()   # default fallback
) %>%
  dplyr::distinct()

cat("FFL core motifs detected (miRNA–TF–Gene triples):", nrow(triplets), "\n")

if (nrow(triplets) == 0) {
  warning("No FFL triplets detected. Continuing with global network analysis only.")
} else {
  write.csv(
    triplets,
    file.path(out_dir, paste0(NETWORK_PREFIX, "_FFL_triplets_miRNA_TF_Gene.csv")),
    row.names = FALSE
  )
}

FFL_genes     <- sort(unique(triplets$Gene))
FFL_nodes_all <- sort(unique(c(triplets$miRNA, triplets$TF, triplets$Gene)))
V(g_comp)$in_FFL <- V(g_comp)$name %in% FFL_nodes_all


# ==============================================================================
# SECTION 5: Topological Centrality and Module Detection
# ==============================================================================

cat("--- [5/10] Computing centrality and Louvain modules ---\n")

# Three standard centrality measures; all normalized to [0, 1]
V(g_comp)$degree      <- igraph::degree(g_comp,      normalized = TRUE)
V(g_comp)$betweenness <- igraph::betweenness(g_comp, normalized = TRUE)
V(g_comp)$closeness   <- igraph::closeness(g_comp,   normalized = TRUE)

# Composite Score = mean of z-scored centrality measures
# (equal weighting; nodes with consistently high centrality across all three
#  measures rank highest)
results_df <- tibble(
  name        = V(g_comp)$name,
  Degree      = V(g_comp)$degree,
  Betweenness = V(g_comp)$betweenness,
  Closeness   = V(g_comp)$closeness,
  Type        = comp_nodes$Type,
  in_FFL      = V(g_comp)$in_FFL
) %>%
  dplyr::mutate(
    Composite_Score = (scale(Degree) + scale(Betweenness) + scale(Closeness)) / 3
  ) %>%
  dplyr::arrange(desc(Composite_Score))

write.csv(
  head(results_df, NUM_TOP_GENES_CENTRALITY),
  file.path(out_dir, "01_Top_Hubs.csv"),
  row.names = FALSE
)
cat("Top hubs written to 01_Top_Hubs.csv\n")

# Louvain community detection (unsupervised; resolution not fixed)
cl <- cluster_louvain(g_comp)
V(g_comp)$module <- membership(cl)

module_stats <- tibble(
  name   = V(g_comp)$name,
  Module = V(g_comp)$module
) %>%
  dplyr::group_by(Module) %>%
  dplyr::mutate(Size = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(Size))

write.csv(
  module_stats,
  file.path(out_dir, "03_Network_Modules.csv"),
  row.names = FALSE
)
cat("Modules written to 03_Network_Modules.csv\n\n")


# ==============================================================================
# SECTION 6: Signed Random Walk with Restart (RWR) Propagation
# ==============================================================================

cat("--- [6/10] Running Signed RWR propagation ---\n")

# Assign signed edge weights:
#   miRNA source → -1  (miRNA-mediated inhibition / post-transcriptional silencing)
#   TF or Gene source → +1  (transcriptional activation)
edges_signed <- edges %>%
  left_join(nodes_clean, by = c("from" = "name"))

if (!"Type" %in% colnames(edges_signed)) {
  stop("'Type' column missing after joining nodes to edges. Check attribute CSV column names.")
}

edges_signed <- edges_signed %>%
  dplyr::rename(Type_Source = Type) %>%
  dplyr::mutate(
    weight_signed = ifelse(grepl("mirna", Type_Source, ignore.case = TRUE), -1, 1)
  )

g_signed     <- igraph::graph_from_data_frame(edges_signed, directed = TRUE, vertices = nodes_clean)
g_signed_sub <- igraph::induced_subgraph(g_signed, V(g_comp)$name)

# Seed selection: prefer FFL hubs (they have confirmed regulatory topology);
# fall back to overall top hubs if fewer than 5 FFL hubs are in the top set
top_hubs          <- head(results_df, NUM_HUBS_TO_TEST_PROPAGATION)$name
ffl_hub_candidates <- top_hubs[top_hubs %in% FFL_nodes_all]

if (length(ffl_hub_candidates) >= 5) {
  seeds <- ffl_hub_candidates
  cat("Seed strategy: FFL-aware (n =", length(seeds), "FFL hubs)\n")
} else {
  seeds <- top_hubs
  cat("Seed strategy: top hubs (insufficient FFL hubs in top set)\n")
}

# Initialise propagation vector (1 for seeds, 0 for all others)
p_vec <- setNames(rep(0, vcount(g_signed_sub)), V(g_signed_sub)$name)
p_vec[seeds[seeds %in% names(p_vec)]] <- 1

scores <- signed_influence_propagation(
  g_signed_sub,
  initial_vector = p_vec,
  alpha          = RWR_ALPHA
)

prop_results <- tibble(
  name  = names(scores),
  score = scores
) %>%
  dplyr::arrange(desc(abs(score)))

write.csv(prop_results, file.path(out_dir, "02_Propagation_Scores.csv"), row.names = FALSE)
cat("Propagation scores written to 02_Propagation_Scores.csv\n\n")


# ==============================================================================
# SECTION 7: Functional Enrichment — ORA and GSEA
# ==============================================================================

cat("--- [7/10] Running functional enrichment ---\n")

# Exclude miRNAs from gene-centric enrichment databases
gene_nodes <- prop_results %>%
  dplyr::filter(!grepl("mir|let-", name, ignore.case = TRUE))

# --- 7A. Over-Representation Analysis (ORA) on top-N influenced genes ---
# Uses the top ORA_TOP_N genes by absolute Influence Score (per paper methods)
top_genes_global <- gene_nodes %>%
  dplyr::arrange(desc(abs(score))) %>%
  dplyr::slice_head(n = ORA_TOP_N) %>%
  dplyr::pull(name)

ids_ora <- symbol_to_entrez_df(top_genes_global)
cat("Genes mapped to Entrez IDs for ORA:", nrow(ids_ora), "\n")

if (nrow(ids_ora) > 5) {
  entrez_vec <- ids_ora$ENTREZID

  # Reactome ORA
  tryCatch({
    er <- ReactomePA::enrichPathway(
      gene         = entrez_vec,
      pvalueCutoff = 0.2,
      readable     = TRUE
    )
    if (!is.null(er) && nrow(as.data.frame(er)) > 0) {
      write.csv(as.data.frame(er),
                file.path(out_dir, "Enr_Global_Reactome_ORA.csv"), row.names = FALSE)
      ggsave(file.path(out_dir, "Enr_Global_Reactome_Dotplot.png"),
             dotplot(er, showCategory = 15),
             width = 10, height = 8, dpi = 300)
      cat("   Reactome ORA saved\n")
    }
  }, error = function(e) cat("   WARNING: Reactome ORA failed:", conditionMessage(e), "\n"))

  # GO ORA (all ontologies: BP, MF, CC)
  tryCatch({
    ego <- clusterProfiler::enrichGO(
      gene         = entrez_vec,
      OrgDb        = org.Hs.eg.db,
      ont          = "ALL",
      pvalueCutoff = 0.1,
      readable     = TRUE
    )
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      p_go <- dotplot(ego, showCategory = 10, split = "ONTOLOGY") +
        facet_grid(ONTOLOGY ~ ., scale = "free")
      ggsave(file.path(out_dir, "Enr_Global_GO_Dotplot.png"),
             p_go, width = 10, height = 12, dpi = 300)
      cat("   GO ORA saved\n")
    }
  }, error = function(e) cat("   WARNING: GO ORA failed:", conditionMessage(e), "\n"))

  # KEGG ORA (requires internet connection)
  tryCatch({
    ekegg <- clusterProfiler::enrichKEGG(
      gene         = entrez_vec,
      organism     = "hsa",
      pvalueCutoff = 0.1
    )
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      ggsave(file.path(out_dir, "Enr_Global_KEGG_Dotplot.png"),
             dotplot(ekegg, showCategory = 15),
             width = 10, height = 8, dpi = 300)
      cat("   KEGG ORA saved\n")
    }
  }, error = function(e) cat("   WARNING: KEGG ORA failed (check internet):", conditionMessage(e), "\n"))
}

# --- 7A-bis. FFL-subset ORA ---
# Runs Reactome ORA on FFL-member genes only, for FFL-specific pathway insight
if (length(FFL_genes) > 0) {
  ffl_gene_inputs <- gene_nodes %>%
    dplyr::filter(name %in% FFL_genes) %>%
    dplyr::arrange(desc(abs(score))) %>%
    dplyr::slice_head(n = ORA_TOP_N) %>%
    dplyr::pull(name)

  ids_ffl <- symbol_to_entrez_df(ffl_gene_inputs)

  if (nrow(ids_ffl) > 5) {
    tryCatch({
      er_ffl <- ReactomePA::enrichPathway(
        gene         = ids_ffl$ENTREZID,
        pvalueCutoff = 0.3,
        readable     = TRUE
      )
      if (!is.null(er_ffl) && nrow(as.data.frame(er_ffl)) > 0) {
        write.csv(as.data.frame(er_ffl),
                  file.path(out_dir, "Enr_FFL_Reactome_ORA.csv"), row.names = FALSE)
        ggsave(file.path(out_dir, "Enr_FFL_Reactome_Dotplot.png"),
               dotplot(er_ffl, showCategory = 15),
               width = 10, height = 8, dpi = 300)
        cat("   FFL-subset Reactome ORA saved\n")
      }
    }, error = function(e) NULL)
  }
}

# --- 7B. GSEA on ranked Influence Scores (GO BP + Reactome) ---
cat("   Running GSEA (GO BP + Reactome)...\n")

gsea_input <- gene_nodes %>%
  inner_join(symbol_to_entrez_df(gene_nodes$name), by = c("name" = "SYMBOL")) %>%
  dplyr::arrange(desc(score)) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

gene_list_gsea <- setNames(gsea_input$score, gsea_input$ENTREZID)
gene_list_gsea <- sort(gene_list_gsea, decreasing = TRUE)

if (length(gene_list_gsea) > 10) {
  # GO BP GSEA
  tryCatch({
    gse_go <- clusterProfiler::gseGO(
      geneList     = gene_list_gsea,
      OrgDb        = org.Hs.eg.db,
      ont          = "BP",
      pvalueCutoff = 1.0,    # permissive; filter visually if needed
      minGSSize    = 5,
      seed         = 123
    )
    if (!is.null(gse_go) && nrow(gse_go@result) > 0) {
      gse_go@result <- gse_go@result %>% dplyr::arrange(pvalue)
      p <- ridgeplot(gse_go, showCategory = 15) +
        ggtitle("GO Biological Process — GSEA Trends")
      ggsave(file.path(out_dir, "Enr_Global_GO_Ridge.png"),
             p, width = 10, height = 8, dpi = 300)
      cat("   GO GSEA ridgeplot saved\n")
    }
  }, error = function(e) cat("   WARNING: GO GSEA failed:", conditionMessage(e), "\n"))

  # Reactome GSEA
  tryCatch({
    gse_react <- ReactomePA::gsePathway(
      geneList     = gene_list_gsea,
      pvalueCutoff = 1.0,
      minGSSize    = 5,
      seed         = 123
    )
    if (!is.null(gse_react) && nrow(gse_react@result) > 0) {
      gse_react@result <- gse_react@result %>% dplyr::arrange(pvalue)
      p <- ridgeplot(gse_react, showCategory = 15) +
        ggtitle("Reactome Pathways — GSEA Trends")
      ggsave(file.path(out_dir, "Enr_Global_Reactome_Ridge.png"),
             p, width = 10, height = 8, dpi = 300)
      cat("   Reactome GSEA ridgeplot saved\n")
    }
  }, error = function(e) cat("   WARNING: Reactome GSEA failed:", conditionMessage(e), "\n"))
}


# ==============================================================================
# SECTION 8: Unified GSEA (Genes + miRNAs via msigdbr + optional miRWalk GMT)
# ==============================================================================

cat("--- [8/10] Unified GSEA (Genes + miRNAs) ---\n")

# Build a combined ranked list (upper-case names, strip strand suffixes -3p/-5p
# so miRNA names can match pathway GMT gene symbols)
ranked_list_df <- prop_results %>%
  dplyr::mutate(
    MatchName = toupper(name),
    MatchName = str_remove_all(MatchName, "-(3|5)P")
  ) %>%
  dplyr::group_by(MatchName) %>%
  dplyr::summarise(score = max(score), .groups = "drop") %>%
  dplyr::mutate(score = score + rnorm(n(), 0, 1e-10)) %>%   # break ties
  dplyr::arrange(desc(score))

gene_list_uni     <- setNames(ranked_list_df$score, ranked_list_df$MatchName)
unified_term2gene <- NULL

# Load MSigDB KEGG (C2:CP:KEGG) and GO BP (C5:GO:BP) gene sets
tryCatch({
  h_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::mutate(gene_symbol = toupper(gene_symbol))

  h_go <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::mutate(gene_symbol = toupper(gene_symbol))

  unified_term2gene <- rbind(h_kegg, h_go)
  cat("   MSigDB KEGG + GO sets loaded\n")
}, error = function(e) cat("   WARNING: MSigDB load failed:", conditionMessage(e), "\n"))

# Optionally append miRNA-specific pathways from miRWalk GMT
if (file.exists(MIRNA_GMT_FILE)) {
  mirna_gmt <- clusterProfiler::read.gmt(MIRNA_GMT_FILE) %>%
    dplyr::mutate(
      gene = toupper(gene),
      gene = str_remove_all(gene, "-(3|5)P")
    )
  colnames(mirna_gmt) <- c("gs_name", "gene_symbol")
  unified_term2gene   <- if (is.null(unified_term2gene)) mirna_gmt else rbind(unified_term2gene, mirna_gmt)
  cat("   miRWalk GMT appended\n")
}

if (!is.null(unified_term2gene)) {
  tryCatch({
    uni_gsea <- clusterProfiler::GSEA(
      geneList     = gene_list_uni,
      TERM2GENE    = unified_term2gene,
      pvalueCutoff = 1.0,
      minGSSize    = 3,
      seed         = 123
    )
    if (!is.null(uni_gsea) && nrow(uni_gsea@result) > 0) {
      uni_plot         <- uni_gsea
      uni_plot@result  <- uni_gsea@result %>% dplyr::arrange(pvalue) %>% head(15)
      p_dot <- dotplot(uni_plot, showCategory = 15) +
        ggtitle("Top Unified Trends (Genes + miRNAs)")
      ggsave(file.path(out_dir, "Enr_Unified_Dotplot.png"),
             p_dot, width = 12, height = 10, dpi = 300)
      cat("   Unified GSEA dotplot saved\n")
    }
  }, error = function(e) cat("   WARNING: Unified GSEA failed:", conditionMessage(e), "\n"))
}


# ==============================================================================
# SECTION 9: Drug Target Identification (DGIdb)
# ==============================================================================

cat("--- [9/10] Drug-gene interaction screening (DGIdb) ---\n")

if (file.exists(DGIDB_FILE)) {
  header_preview <- read.delim(DGIDB_FILE, nrows = 2, header = TRUE)
  all_cols <- colnames(header_preview)
  gene_col <- all_cols[grep("gene", all_cols, ignore.case = TRUE)][1]
  drug_col <- all_cols[grep("drug", all_cols, ignore.case = TRUE)][1]

  if (!is.na(gene_col) && !is.na(drug_col)) {
    dgidb_full  <- read.delim(DGIDB_FILE, sep = "\t", header = TRUE,
                               stringsAsFactors = FALSE, quote = "", check.names = FALSE)
    top_hub_names <- read.csv(file.path(out_dir, "01_Top_Hubs.csv")) %>% dplyr::pull(name)

    drug_clean <- dgidb_full %>%
      dplyr::filter(!!sym(gene_col) %in% top_hub_names) %>%
      dplyr::select(Gene = !!sym(gene_col), Drug = !!sym(drug_col)) %>%
      dplyr::distinct() %>%
      dplyr::filter(nzchar(Drug))

    if (nrow(drug_clean) > 0) {
      write.csv(drug_clean, file.path(out_dir, "04_Drug_Targets.csv"), row.names = FALSE)

      drug_counts <- drug_clean %>%
        dplyr::group_by(Gene) %>%
        dplyr::summarise(Count = n(), .groups = "drop") %>%
        dplyr::arrange(desc(Count))

      p_drug <- ggplot(head(drug_counts, 15), aes(x = reorder(Gene, Count), y = Count)) +
        geom_col(fill = "#E74C3C") +
        coord_flip() +
        theme_minimal(base_size = 12) +
        labs(title = "Druggable Network Hubs", y = "Number of Drug Interactions", x = "Hub Gene")
      ggsave(file.path(out_dir, "04_Drug_Counts.png"),
             p_drug, width = 6, height = 6, dpi = 300)

      cat("   Drug targets saved (", nrow(drug_clean), "interactions )\n")
    }
  }
} else {
  cat("   DGIdb file not found — skipping drug screening\n")
}


# ==============================================================================
# SECTION 10: Experimental Validation (TCGA-BRCA)
# ==============================================================================

cat("--- [10/10] TCGA-BRCA experimental validation ---\n")

# Select top hubs per node type for comparison against DEX data
final_targets <- results_df %>%
  dplyr::filter(Type != "Unknown") %>%
  dplyr::group_by(Type) %>%
  dplyr::slice_max(order_by = Composite_Score,
                   n         = NUM_HUBS_TO_REPORT_PER_TYPE,
                   with_ties = FALSE) %>%
  dplyr::ungroup()

if (file.exists(DEX_FILE_NAME)) {
  dex_raw <- read.csv(DEX_FILE_NAME)

  # Clean and deduplicate DEX table
  dex_clean <- dex_raw %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Node_Name_Match = clean_name_for_match(Gene)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(adj.P.Val = ifelse(is.na(adj.P.Val), 1, adj.P.Val)) %>%
    dplyr::group_by(Node_Name_Match) %>%
    dplyr::slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      Node_Name_Match,
      Experimental_LogFC   = logFC,
      Experimental_PValue  = adj.P.Val
    )

  # Join network hubs with DEX data
  merged_data <- final_targets %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Node_Name_Match = clean_name_for_match(name)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(dex_clean, by = "Node_Name_Match") %>%
    dplyr::filter(!is.na(Experimental_LogFC))

  # *** BUG FIX: Apply paper significance thresholds before plotting ***
  # Paper: "significantly dysregulated if |Log2FC| > 1.0 AND adjusted p < 0.05"
  merged_data <- merged_data %>%
    dplyr::mutate(
      Significant = abs(Experimental_LogFC) > LOGFC_THRESHOLD &
                    Experimental_PValue      < PVAL_THRESHOLD
    )

  n_sig   <- sum(merged_data$Significant,  na.rm = TRUE)
  n_total <- nrow(merged_data)
  cat("   Matched hubs:", n_total, "| Significant (|logFC|>", LOGFC_THRESHOLD,
      "& padj<", PVAL_THRESHOLD, "):", n_sig, "\n")

  if (n_total > 0) {
    # Barplot — colour by node type; black outline marks significance
    # (replicates Figure 6 design: "black outlines denote statistical significance")
    p_bar <- ggplot(
      merged_data,
      aes(x    = reorder(name, Experimental_LogFC),
          y    = Experimental_LogFC,
          fill = Type,
          colour = Significant,
          linewidth = Significant)
    ) +
      geom_col() +
      scale_colour_manual(
        values = c("TRUE" = "black", "FALSE" = NA),
        labels = c("TRUE" = "Significant", "FALSE" = ""),
        name   = paste0("|logFC|>", LOGFC_THRESHOLD, " & padj<", PVAL_THRESHOLD)
      ) +
      scale_linewidth_manual(values = c("TRUE" = 0.6, "FALSE" = 0)) +
      scale_fill_brewer(palette = "Set2") +
      coord_flip() +
      facet_wrap(~Type, scales = "free_y") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom") +
      labs(
        title = paste0("Hub Expression in TCGA-BRCA: ", NETWORK_PREFIX),
        x     = NULL,
        y     = "Log2 Fold Change (Tumor vs. Normal)"
      )

    ggsave(file.path(out_dir, "Final_Validation_Barplot.png"),
           p_bar, width = 12, height = 10, dpi = 300)

    write.csv(merged_data,
              file.path(out_dir, "Final_Validation_Table.csv"), row.names = FALSE)

    cat("   Validation barplot and table saved\n")
  } else {
    cat("   WARNING: No matched hubs found between network and DEX data\n")
  }
} else {
  cat("   DEX file not found — skipping validation\n",
      "   Run 01_BRCA_miRNA_DESeq2_Preprocessing.R first and ensure",
      DEX_FILE_NAME, "is in the working directory.\n")
}

cat("\n=== PIPELINE COMPLETE ===\n")
cat("All results saved to:", out_dir, "\n")
