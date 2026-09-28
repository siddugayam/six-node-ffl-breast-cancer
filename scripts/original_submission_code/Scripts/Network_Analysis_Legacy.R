# ==============================================================================
# SCRIPT 03: One-Click Signed Network Analysis — Legacy Pipeline
# ==============================================================================
#
# STATUS:  LEGACY / SIMPLIFIED
#   This is an earlier, single-network version of the analysis pipeline.
#   For the full analysis (multi-network, FFL detection, Reactome, drug
#   screening), use Script 02 (02_Signed_FFL_Network_Analysis.R) instead.
#   This script is retained as a self-contained quick-run option for a
#   single network without optional dependencies (ReactomePA, msigdbr, DOSE).
#
# DIFFERENCES FROM SCRIPT 02:
#   - No FFL topology detection
#   - No Louvain module detection
#   - No Reactome / Disease Ontology enrichment
#   - No drug screening
#   - No unified miRNA + gene GSEA
#   - Single hardcoded network (configure below); cannot loop over prefixes
#
# PURPOSE:
#   Performs centrality analysis, Signed RWR propagation, GO + KEGG enrichment,
#   and TCGA-BRCA validation on a single miRNA-TF-gene regulatory network.
#
# INPUTS (must be in working directory):
#   {NODE_ATTR_FILE}   — Node attribute CSV (columns: name, type)
#   {SIF_FILE}         — Edge list in SIF format (whitespace-separated)
#   {DEX_FILE_NAME}    — DEX results CSV (columns: Gene, logFC, adj.P.Val)
#                        Produced by Script 01 (01_BRCA_miRNA_DESeq2_Preprocessing.R)
#
# OUTPUTS (saved to {prefix}_RESULTS/):
#   01_Top_Hubs.csv              — Top N hubs by Composite Centrality
#   02_Propagation_Scores.csv    — Signed Influence Scores for all nodes
#   Enr_01_GO_ORA.png            — GO ORA barplot
#   Enr_02_KEGG_ORA.png/.csv     — KEGG ORA dotplot and table
#   Enr_04_GSEA_Ridge_Clean.png  — GO GSEA ridgeplot (high-confidence pathways)
#   Final_Validation_Barplot.png — Hub LogFC barplot (sig. outlined)
#   Final_Validation_Stats.png   — ANOVA violin/box comparison by node type
#   Final_Validation_Table.csv   — Merged hub + DEX data (significance-filtered)
#
# DEPENDENCIES:
#   Bioconductor: clusterProfiler, org.Hs.eg.db, enrichplot, DOSE
#   CRAN:         tidyverse, igraph, ggraph, patchwork, tidygraph, scales,
#                 ggrepel, forcats, ggpubr, stringr, data.table
#
# AUTHORS:  Gayam Prasanna Kumar Reddy, Jesil Mathew A, Fayaz Shaik Mahammad
# VERSION:  21.1 (cleaned — bugs fixed: ORA top-N filter applied, GSEA KEGG
#           added, DEX significance thresholds enforced, significance markers
#           on barplot, section numbering corrected)
# ==============================================================================


# ==============================================================================
# SECTION 0: User Configuration  ← EDIT THESE VALUES
# ==============================================================================

# Network file names (stems must match)
NODE_ATTR_FILE <- "3-Comp.csv"    # node attribute CSV
SIF_FILE       <- "3-Comp.sif"   # edge list SIF file

# DEX validation file (produced by Script 01)
DEX_FILE_NAME  <- "BRCA_DEX_ALL_nodes.csv"

# Analysis parameters
NUM_TOP_GENES_CENTRALITY     <- 60    # top hubs to write (01_Top_Hubs.csv)
NUM_HUBS_TO_TEST_PROPAGATION <- 60    # seeds for RWR
NUM_HUBS_TO_REPORT_PER_TYPE  <- 20    # hubs per node type in validation
NUM_GENES_FOR_ORA_ENRICHMENT <- 100   # top N genes (by |influence|) for ORA
RWR_ALPHA                    <- 0.85  # RWR restart probability

# Validation significance thresholds (paper: |logFC| > 1.0, padj < 0.05)
LOGFC_THRESHOLD <- 1.0
PVAL_THRESHOLD  <- 0.05


# ==============================================================================
# SECTION 1: Package Installation and Loading
# ==============================================================================

cat("--- [1/7] Loading packages ---\n")

cran_packages <- c("tidyverse", "igraph", "ggraph", "patchwork", "tidygraph",
                   "scales", "ggrepel", "forcats", "ggpubr", "stringr", "data.table")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}


# ==============================================================================
# SECTION 2: Core Functions
# ==============================================================================

cat("--- [2/7] Defining core functions ---\n")

# ------------------------------------------------------------------------------
# signed_influence_propagation()
#   Signed Random Walk with Restart. Propagates regulatory signals through a
#   signed directed network from a set of seed hubs.
#   miRNA-sourced edges: weight = -1 (inhibition)
#   TF/gene-sourced edges: weight = +1 (activation)
#
# Args:
#   graph          — igraph object with 'weight_signed' edge attribute
#   initial_vector — named numeric seed vector (1 = seed, 0 = other)
#   alpha          — restart probability (default 0.85)
#   max_iter       — iteration cap
#   tolerance      — L1 convergence threshold
#
# Returns: named numeric vector of Signed Influence Scores.
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

  # Preserve row names for re-attachment after matrix multiplication
  node_names <- rownames(A)

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
  names(result_vec) <- node_names    # force names back after matrix ops
  return(result_vec)
}


# ------------------------------------------------------------------------------
# clean_name_for_match()
#   Standardises node name for DEX table matching.
#   miRNAs (hsa-*) are upper-cased and truncated to 3 fields; others upper-cased.
# ------------------------------------------------------------------------------
clean_name_for_match <- function(x) {
  parts <- str_split(x, "-")[[1]]
  if (length(parts) >= 3 && grepl("^(hsa|HSA)", parts[1], ignore.case = TRUE)) {
    return(toupper(paste(parts[1], parts[2], parts[3], sep = "-")))
  }
  return(toupper(x))
}


# ==============================================================================
# SECTION 3: Data Import and Graph Construction
# ==============================================================================

cat("--- [3/7] Loading network and building graph ---\n")

# --- Parse SIF edge list ---
raw_lines <- readLines(SIF_FILE)
raw_lines <- raw_lines[nzchar(trimws(raw_lines))]
sif_list  <- strsplit(raw_lines, "\\s+")

edges <- data.frame(
  from = sapply(sif_list, `[`, 1),
  to   = sapply(sif_list, function(x) tail(x, 1)),
  stringsAsFactors = FALSE
) %>%
  dplyr::distinct()

# --- Parse node attribute table ---
nodes_raw <- read.csv(NODE_ATTR_FILE, check.names = FALSE)
cols      <- colnames(nodes_raw)

name_idx <- grep("^(shared[ ._]?)?name$", cols, ignore.case = TRUE)[1]
type_idx <- grep("^(shared[ ._]?)?type$", cols, ignore.case = TRUE)[1]

if (is.na(name_idx) || is.na(type_idx)) {
  stop("Could not find 'name' or 'type' columns. Found: ", paste(cols, collapse = ", "))
}

nodes_clean <- nodes_raw %>%
  dplyr::select(name = all_of(name_idx), Type = all_of(type_idx)) %>%
  dplyr::mutate(
    name = trimws(as.character(name)),
    Type = trimws(as.character(Type))
  ) %>%
  dplyr::distinct(name, .keep_all = TRUE)

# Add missing nodes as 'Unknown'
missing_nodes <- setdiff(unique(c(edges$from, edges$to)), nodes_clean$name)
if (length(missing_nodes) > 0) {
  cat(paste("   Adding", length(missing_nodes), "missing nodes as 'Unknown'\n"))
  nodes_clean <- bind_rows(nodes_clean, tibble(name = missing_nodes, Type = "Unknown"))
}

# Build undirected graph; extract largest connected component
g_full <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes_clean)
g_comp <- igraph::decompose(g_full)[[which.max(sapply(igraph::decompose(g_full), vcount))]]
cat(paste("Network built — Nodes:", vcount(g_comp), "| Edges:", ecount(g_comp), "\n"))

# Create output directory named after the network file
out_dir <- paste0(tools::file_path_sans_ext(basename(NODE_ATTR_FILE)), "_RESULTS")
dir.create(out_dir, showWarnings = FALSE)


# ==============================================================================
# SECTION 4: Centrality Analysis
# ==============================================================================

cat("--- [4/7] Computing centrality ---\n")

V(g_comp)$degree      <- igraph::degree(g_comp,      normalized = TRUE)
V(g_comp)$betweenness <- igraph::betweenness(g_comp, normalized = TRUE)
V(g_comp)$closeness   <- igraph::closeness(g_comp,   normalized = TRUE)

# Composite Score = mean z-score across three centrality measures
results_df <- tibble(
  name        = V(g_comp)$name,
  Degree      = V(g_comp)$degree,
  Betweenness = V(g_comp)$betweenness,
  Closeness   = V(g_comp)$closeness,
  Type        = V(g_comp)$Type
) %>%
  dplyr::mutate(
    Composite_Score = (scale(Degree) + scale(Betweenness) + scale(Closeness)) / 3
  ) %>%
  dplyr::arrange(desc(Composite_Score))

write.csv(head(results_df, NUM_TOP_GENES_CENTRALITY),
          file.path(out_dir, "01_Top_Hubs.csv"), row.names = FALSE)
cat("   Top hubs written to 01_Top_Hubs.csv\n")


# ==============================================================================
# SECTION 5: Signed RWR Propagation
# ==============================================================================

cat("--- [5/7] Running Signed RWR propagation ---\n")

# Assign signed weights:
#   miRNA source → -1 (post-transcriptional silencing)
#   TF/gene source → +1 (transcriptional activation)
edges_signed <- edges %>%
  left_join(nodes_clean, by = c("from" = "name"))

if (!"Type" %in% colnames(edges_signed)) {
  stop("'Type' column not found after joining edges to nodes. Check CSV column names.")
}

edges_signed <- edges_signed %>%
  dplyr::rename(Type_Source = Type) %>%
  dplyr::mutate(
    weight_signed = ifelse(grepl("miRNA", Type_Source, ignore.case = TRUE), -1, 1)
  )

g_signed     <- igraph::graph_from_data_frame(edges_signed, directed = TRUE, vertices = nodes_clean)
g_signed_sub <- igraph::induced_subgraph(g_signed, V(g_comp)$name)

# Seed initialisation: top-N hubs get score 1; all others 0
seeds <- head(results_df, NUM_HUBS_TO_TEST_PROPAGATION)$name
p_vec <- setNames(rep(0, vcount(g_signed_sub)), V(g_signed_sub)$name)
p_vec[seeds[seeds %in% names(p_vec)]] <- 1

scores <- signed_influence_propagation(g_signed_sub, p_vec, alpha = RWR_ALPHA)

if (is.null(names(scores))) stop("Propagation scores lost their names. Check the function.")

prop_results <- tibble(name = names(scores), score = scores) %>%
  dplyr::arrange(desc(abs(score)))

write.csv(prop_results, file.path(out_dir, "02_Propagation_Scores.csv"), row.names = FALSE)
cat("   Propagation scores written to 02_Propagation_Scores.csv\n")


# ==============================================================================
# SECTION 6: Functional Enrichment (GO ORA, KEGG ORA, KEGG GSEA, GO GSEA)
# ==============================================================================

cat("--- [6/7] Running functional enrichment ---\n")

# Exclude miRNAs — standard enrichment databases operate on gene symbols only
gene_nodes_all <- nodes_clean %>%
  dplyr::filter(!grepl("miRNA", Type, ignore.case = TRUE)) %>%
  dplyr::pull(name)

# *** BUG FIX: apply NUM_GENES_FOR_ORA_ENRICHMENT filter (was previously unused) ***
# Select top N gene nodes by absolute Signed Influence Score for ORA
enrich_input <- prop_results %>%
  dplyr::filter(name %in% gene_nodes_all) %>%
  dplyr::arrange(desc(abs(score))) %>%
  dplyr::slice_head(n = NUM_GENES_FOR_ORA_ENRICHMENT)

cat("   Genes submitted to ORA:", nrow(enrich_input),
    "(top", NUM_GENES_FOR_ORA_ENRICHMENT, "by |Influence Score|)\n")

# Map gene symbols to Entrez IDs
ids <- tryCatch(
  bitr(enrich_input$name,
       fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
  error = function(e) NULL
)

if (!is.null(ids) && nrow(ids) > 0) {

  # --- GO ORA (all ontologies) ---
  cat("   Running GO ORA...\n")
  tryCatch({
    ora_go <- clusterProfiler::enrichGO(
      ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05
    )
    if (!is.null(ora_go)) {
      ggsave(file.path(out_dir, "Enr_01_GO_ORA.png"),
             barplot(ora_go, showCategory = 20), width = 10, height = 8)
      cat("   GO ORA saved\n")
    }
  }, error = function(e) cat("   WARNING: GO ORA failed:", conditionMessage(e), "\n"))

  # --- KEGG ORA (requires internet) ---
  cat("   Running KEGG ORA...\n")
  tryCatch({
    kk <- clusterProfiler::enrichKEGG(
      gene = ids$ENTREZID, organism = "hsa", pvalueCutoff = 0.05
    )
    if (!is.null(kk)) {
      ggsave(file.path(out_dir, "Enr_02_KEGG_ORA.png"),
             dotplot(kk, showCategory = 20) + ggtitle("KEGG Pathway Enrichment"),
             width = 10, height = 8)
      write.csv(as.data.frame(kk), file.path(out_dir, "Enr_02_KEGG_Table.csv"))
      cat("   KEGG ORA saved\n")
    }
  }, error = function(e) cat("   WARNING: KEGG ORA failed (check internet):", conditionMessage(e), "\n"))

  # --- KEGG GSEA on full ranked gene list ---
  # *** BUG FIX: added gseKEGG (was missing; paper states GSEA on KEGG + GO) ***
  cat("   Running KEGG GSEA...\n")
  tryCatch({
    gsea_kegg_input <- prop_results %>%
      dplyr::filter(name %in% gene_nodes_all) %>%
      left_join(ids, by = c("name" = "SYMBOL")) %>%
      dplyr::filter(!is.na(ENTREZID)) %>%
      dplyr::mutate(score = score + rnorm(n(), 0, 1e-10)) %>%
      dplyr::arrange(desc(score)) %>%
      dplyr::distinct(ENTREZID, .keep_all = TRUE)

    gene_list_kegg <- setNames(gsea_kegg_input$score, gsea_kegg_input$ENTREZID)
    gene_list_kegg <- sort(gene_list_kegg, decreasing = TRUE)

    if (length(gene_list_kegg) > 10) {
      gsea_kk <- clusterProfiler::gseKEGG(
        geneList     = gene_list_kegg,
        organism     = "hsa",
        pvalueCutoff = 0.05,
        verbose      = FALSE
      )
      if (!is.null(gsea_kk) && nrow(gsea_kk@result) > 0) {
        ggsave(file.path(out_dir, "Enr_03_KEGG_GSEA_Ridge.png"),
               ridgeplot(gsea_kk, showCategory = 20) + ggtitle("KEGG GSEA"),
               width = 10, height = 8)
        write.csv(as.data.frame(gsea_kk),
                  file.path(out_dir, "Enr_03_KEGG_GSEA_Table.csv"))
        cat("   KEGG GSEA saved\n")
      }
    }
  }, error = function(e) cat("   WARNING: KEGG GSEA failed:", conditionMessage(e), "\n"))

  # --- GO GSEA on full ranked gene list ---
  cat("   Running GO GSEA...\n")
  tryCatch({
    gsea_go_input <- prop_results %>%
      dplyr::filter(name %in% gene_nodes_all) %>%
      left_join(ids, by = c("name" = "SYMBOL")) %>%
      dplyr::filter(!is.na(ENTREZID)) %>%
      dplyr::mutate(score = score + rnorm(n(), 0, 1e-10)) %>%
      dplyr::arrange(desc(score)) %>%
      dplyr::distinct(ENTREZID, .keep_all = TRUE)

    gene_list_gsea <- setNames(gsea_go_input$score, gsea_go_input$ENTREZID)
    gene_list_gsea <- sort(gene_list_gsea, decreasing = TRUE)

    if (length(gene_list_gsea) > 10) {
      gsea_go <- clusterProfiler::gseGO(
        gene_list_gsea, OrgDb = org.Hs.eg.db, ont = "ALL",
        pvalueCutoff = 0.05, verbose = FALSE
      )

      if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
        # Keep only high-confidence pathways (setSize >= 10) for clean plot
        gsea_clean         <- gsea_go
        gsea_clean@result  <- gsea_go@result %>%
          dplyr::filter(setSize >= 10) %>%
          dplyr::arrange(p.adjust)

        p_ridge <- ridgeplot(gsea_clean, showCategory = 20, fill = "p.adjust") +
          labs(title = "Top 20 Enriched GO Pathways (GSEA — High Confidence)") +
          theme(axis.text.y = element_text(size = 10))

        ggsave(file.path(out_dir, "Enr_04_GSEA_Ridge_Clean.png"),
               p_ridge, width = 10, height = 8)
        cat("   GO GSEA ridgeplot saved\n")
      }
    }
  }, error = function(e) cat("   WARNING: GO GSEA failed:", conditionMessage(e), "\n"))
}


# ==============================================================================
# SECTION 7: Experimental Validation (TCGA-BRCA)
# ==============================================================================

cat("--- [7/7] TCGA-BRCA experimental validation ---\n")

# Select top hubs per node type for comparison
final_targets <- results_df %>%
  dplyr::filter(Type != "Unknown") %>%
  dplyr::group_by(Type) %>%
  dplyr::slice_max(order_by = Composite_Score,
                   n         = NUM_HUBS_TO_REPORT_PER_TYPE,
                   with_ties = FALSE) %>%
  dplyr::ungroup()

if (file.exists(DEX_FILE_NAME)) {
  dex_raw <- read.csv(DEX_FILE_NAME)

  # Standardise DEX table; deduplicate by best adjusted p-value
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
      Experimental_LogFC  = logFC,
      Experimental_PValue = adj.P.Val
    )

  # Merge and filter by DEX significance thresholds
  # *** BUG FIX: DEX thresholds now enforced; significance column added for plot ***
  merged_data <- final_targets %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Node_Name_Match = clean_name_for_match(name)) %>%
    dplyr::ungroup() %>%
    left_join(dex_clean, by = "Node_Name_Match") %>%
    dplyr::filter(!is.na(Experimental_LogFC)) %>%
    dplyr::mutate(
      Significant = abs(Experimental_LogFC) > LOGFC_THRESHOLD &
                    Experimental_PValue      < PVAL_THRESHOLD
    )

  cat("   Matched hubs:", nrow(merged_data),
      "| Significant (|logFC|>", LOGFC_THRESHOLD, "& padj<", PVAL_THRESHOLD, "):",
      sum(merged_data$Significant, na.rm = TRUE), "\n")

  if (nrow(merged_data) > 0) {
    # *** BUG FIX: black outline on significant bars (per Figure 6 description) ***
    p_bar <- ggplot(
      merged_data,
      aes(x         = reorder(name, Experimental_LogFC),
          y         = Experimental_LogFC,
          fill      = Type,
          colour    = Significant,
          linewidth = Significant)
    ) +
      geom_col() +
      scale_colour_manual(
        values = c("TRUE" = "black", "FALSE" = NA),
        labels = c("TRUE" = "Significant", "FALSE" = ""),
        name   = paste0("|logFC|>", LOGFC_THRESHOLD, " & padj<", PVAL_THRESHOLD)
      ) +
      scale_linewidth_manual(values = c("TRUE" = 0.6, "FALSE" = 0)) +
      coord_flip() +
      facet_wrap(~Type, scales = "free_y") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom") +
      labs(
        title = "Hub Expression in TCGA-BRCA",
        x     = NULL,
        y     = "Log2 Fold Change (Tumor vs. Normal)"
      )

    ggsave(file.path(out_dir, "Final_Validation_Barplot.png"),
           p_bar, width = 10, height = 10, dpi = 300)

    # Violin + boxplot ANOVA comparison across node types
    p_stats <- ggplot(merged_data,
                      aes(x = Type, y = Experimental_LogFC, fill = Type)) +
      geom_violin(alpha = 0.6, trim = FALSE) +
      geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
      ggpubr::stat_compare_means(method = "anova", label.y.npc = "top") +
      tryCatch(
        ggpubr::stat_compare_means(label     = "p.signif",
                                   method    = "t.test",
                                   ref.group = ".all."),
        error = function(e) NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "none") +
      labs(title = "LogFC Distribution by Node Type (ANOVA)",
           y = "Log2 Fold Change", x = "Node Type")

    ggsave(file.path(out_dir, "Final_Validation_Stats.png"),
           p_stats, width = 8, height = 6, dpi = 300)

    write.csv(merged_data,
              file.path(out_dir, "Final_Validation_Table.csv"), row.names = FALSE)

    cat("   Validation plots and table saved\n")
    cat("--- SUCCESS ---\n")
  } else {
    cat("   WARNING: No hubs matched between network and DEX data\n")
  }
} else {
  cat("   DEX file not found:", DEX_FILE_NAME, "— skipping validation\n")
}

cat("\n--- DONE. Results saved to:", out_dir, "---\n")
