# ==============================================================================
# MASTER SCRIPT 03: Network Analysis — Topology, Hub Identification,
#                   Functional Enrichment, Propagation, and Perturbation
# ==============================================================================
#
# CONSOLIDATES (11 original scripts):
#   - Network Topology Analysis Script (v17).R
#       → Section 3: Subnetwork topology tables + community detection + plots
#   - Network-Based Gene Prioritization Script- Hub Gene Identification.R
#       → Section 4: Centrality-based hub ranking
#   - consolidated_network.R
#       → Section 4: Hub-target table generation
#   - Comprehensive Multi-Analysis Script (v10 - FINAL).R
#       → Section 5: K-means clustering + node metric comparison + enrichment
#   - # Script 1 (Upgraded) - Functional Enrichment of Communities.R
#       → Section 6: Community-level enrichment for all molecule types
#   - SIF Network Functional Enrichment Analysis in R .R
#       → Section 6: 4-part enrichment (all genes, miRNA targets, TF-only, non-TF)
#   - GSEA with All Visualizations.R
#       → Section 7: GSEA with ridge, dot, concept network, UpSet, heatmap plots
#   - Directed Signed Network Propagation Script.R
#       → Section 8: Signed RWR propagation + network visualisation
#   - Combined Structural Perturbation & Network Propagation Script.R
#       → Section 9: Node knockout perturbation + criticality ranking
#   - Targtetted subnetworknalaysis .R
#       → Section 10: Per-hub neighborhood enrichment (GO, KEGG, Reactome, DO)
#   - Compare Signed Propagation Scores by Node Type.R  [EMPTY — skipped]
#
# PURPOSE:
#   Comprehensive network analysis pipeline for miRNA–TF–Gene regulatory
#   networks exported from Cytoscape as SIF + node-attribute CSV files.
#   Runs topology, hub identification, enrichment, propagation, perturbation,
#   and per-hub subnetwork analysis in sequence.
#
# RUN MODES (set RUN_SECTIONS in SECTION 0):
#   Each section can be toggled independently with TRUE/FALSE, allowing you
#   to run only the parts you need without re-running the full pipeline.
#
# REQUIRED INPUTS:
#   {NETWORK_PREFIX}.sif     — Edge list (whitespace-separated: source type target)
#   {NETWORK_PREFIX}.csv     — Node attributes (must have 'name' and 'Type' columns)
#   {RANKED_LIST_FILE}       — [Section 7 only] Pre-ranked gene list (.rnk or .csv)
#                              Two columns: gene_symbol, rank_score (no header needed)
#
# OPTIONAL INPUTS:
#   {HUB_SEEDS_FILE}         — CSV of seed genes for propagation (Section 8)
#                              If absent, top hubs from Section 4 are used
#
# OUTPUTS (all written to {NETWORK_PREFIX}_Network_Analysis/):
#   01_topology/             — Global summary + per-subnetwork node metrics
#   02_network_plots/        — ggraph visualisations per subnetwork
#   03_hubs/                 — Hub rankings, hub-target table
#   04_clustering/           — K-means elbow plot, cluster summary
#   05_enrichment/           — GO, KEGG, Reactome, DO for all node subsets
#   06_community_enrichment/ — Enrichment per Walktrap community
#   07_gsea/                 — GSEA plots (ridge, dot, concept, UpSet, heatmap)
#   08_propagation/          — Signed RWR scores + network plot
#   09_perturbation/         — Structural perturbation + criticality barplot
#   10_subnetwork/           — Per-hub neighborhood enrichment folders
#
# KEY PARAMETERS (edit in SECTION 0):
#   NETWORK_PREFIX           — Stem of your .sif and .csv files (e.g. "3-TF")
#   NUM_TOP_HUBS             — Top hubs for enrichment and propagation seeds
#   NUM_HUBS_TO_PERTURB      — Number of top-degree nodes to knock out
#   NUM_HUBS_TO_DISPLAY      — Hubs per type in the hub-target output table
#   RWR_ALPHA                — Restart probability for Random Walk (default 0.85)
#   KMEANS_K                 — Number of clusters (0 = auto-detect via elbow)
#   RUN_SECTIONS             — Named logical vector: toggle individual sections
#
# DEPENDENCIES:
#   CRAN:         tidyverse, igraph, ggraph, patchwork, tidygraph, scales,
#                 ggrepel, ggpubr, factoextra, ggupset, data.table
#   Bioconductor: clusterProfiler, org.Hs.eg.db, enrichplot, ReactomePA, DOSE,
#                 AnnotationDbi, multiMiR
#
# AUTHORS:  Gayam Prasanna Kumar Reddy, Jesil Mathew A, Fayaz Shaik Mahammad
# VERSION:  1.0 (consolidated)
# ==============================================================================


# ==============================================================================
# SECTION 0: Configuration  ← EDIT THESE
# ==============================================================================

NETWORK_PREFIX   <- "3-TF"          # file stem (no extension)
RANKED_LIST_FILE <- "BRCA_ranked_list.rnk"  # for Section 7 GSEA (optional)
HUB_SEEDS_FILE   <- ""              # CSV of seed genes for propagation; "" = use top hubs

NUM_TOP_HUBS        <- 50           # top hubs to report and use as seeds
NUM_HUBS_TO_PERTURB <- 20           # nodes to knock out in perturbation
NUM_HUBS_TO_DISPLAY <- 10           # hubs per type in hub-target table
RWR_ALPHA           <- 0.85         # RWR restart probability
KMEANS_K            <- 0            # 0 = auto-detect; set e.g. 4 to fix

# Toggle individual analysis sections (TRUE = run, FALSE = skip)
RUN_SECTIONS <- c(
  topology      = TRUE,    # Section 3: topology tables + community plots
  hubs          = TRUE,    # Section 4: centrality + hub-target table
  clustering    = TRUE,    # Section 5: K-means + node metric comparison
  enrichment    = TRUE,    # Section 6: GO/KEGG/Reactome/DO enrichment
  gsea          = TRUE,    # Section 7: GSEA visualisations
  propagation   = TRUE,    # Section 8: signed RWR
  perturbation  = TRUE,    # Section 9: structural perturbation
  subnetwork    = TRUE     # Section 10: per-hub neighborhood enrichment
)


# ==============================================================================
# SECTION 1: Package Installation and Loading
# ==============================================================================

cat("=== NETWORK ANALYSIS PIPELINE ===\n")
cat("Network prefix:", NETWORK_PREFIX, "\n\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("tidyverse", "igraph", "ggraph", "patchwork", "tidygraph", "scales",
               "ggrepel", "ggpubr", "factoextra", "data.table")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
# ggupset is optional (for UpSet plots in Section 7)
if (!requireNamespace("ggupset", quietly = TRUE)) {
  tryCatch(install.packages("ggupset"), error = function(e) NULL)
}

bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot",
               "ReactomePA", "DOSE", "AnnotationDbi")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
}
# multiMiR is optional (for miRNA validated target lookup)
if (!requireNamespace("multiMiR", quietly = TRUE)) {
  tryCatch(BiocManager::install("multiMiR", ask = FALSE), error = function(e) NULL)
}

suppressPackageStartupMessages({
  library(tidyverse); library(igraph); library(ggraph); library(patchwork)
  library(tidygraph); library(scales); library(ggrepel); library(ggpubr)
  library(factoextra); library(data.table)
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(ReactomePA); library(DOSE); library(AnnotationDbi)
})

out_dir <- paste0(NETWORK_PREFIX, "_Network_Analysis")
dir.create(out_dir, showWarnings = FALSE)
cat("Output directory:", out_dir, "\n\n")


# ==============================================================================
# SECTION 2: Core Helper Functions
# ==============================================================================

# ------------------------------------------------------------------------------
# signed_influence_propagation()
#   Signed Random Walk with Restart (RWR) on a directed, weighted graph.
#   Propagates activating (+1) and inhibitory (-1) signals from seed nodes.
#   Converges when L1 score change < tolerance.
#
# Args:
#   graph          — igraph with 'weight_signed' edge attribute
#   initial_vector — named numeric seed vector
#   alpha          — restart probability
# Returns: named numeric Signed Influence Score vector
# ------------------------------------------------------------------------------
signed_influence_propagation <- function(graph, initial_vector,
                                         alpha = 0.85, max_iter = 100, tolerance = 1e-6) {
  A     <- as.matrix(igraph::as_adjacency_matrix(graph, attr = "weight_signed",
                                                  type = "both", sparse = TRUE))
  D_inv <- diag(1 / colSums(abs(A))); D_inv[!is.finite(D_inv)] <- 0
  A_norm <- A %*% D_inv
  V_curr <- initial_vector; V_init <- initial_vector
  for (i in seq_len(max_iter)) {
    V_next <- alpha * (A_norm %*% V_curr) + (1 - alpha) * V_init
    if (sum(abs(V_next - V_curr)) < tolerance) { cat("   RWR converged at iter", i, "\n"); break }
    V_curr <- V_next
  }
  result <- as.vector(V_curr); names(result) <- rownames(A)
  return(result)
}

# ------------------------------------------------------------------------------
# run_full_enrichment()
#   Runs GO, KEGG, Reactome, and Disease Ontology ORA on a gene symbol vector.
#   Saves CSV and dotplot PNG for each database that returns results.
#
# Args:
#   gene_syms  — character vector of gene symbols
#   prefix     — output filename prefix
#   save_dir   — directory to write outputs
# ------------------------------------------------------------------------------
run_full_enrichment <- function(gene_syms, prefix, save_dir) {
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  ids <- tryCatch(
    clusterProfiler::bitr(gene_syms, fromType = "SYMBOL", toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  if (is.null(ids) || nrow(ids) < 5) {
    message("  -> Too few Entrez IDs for ", prefix, " — skipping enrichment")
    return(invisible(NULL))
  }
  ev <- ids$ENTREZID

  save_enrich <- function(res, tag) {
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      write.csv(as.data.frame(res), file.path(save_dir, paste0(prefix, "_", tag, ".csv")),
                row.names = FALSE)
      tryCatch({
        p <- dotplot(res, showCategory = 20, title = paste("Top 20:", gsub("_", " ", tag)))
        ggplot2::ggsave(file.path(save_dir, paste0(prefix, "_", tag, ".png")),
                        p, width = 12, height = 10, dpi = 300)
        message("    Saved: ", prefix, "_", tag)
      }, error = function(e) NULL)
    }
  }

  tryCatch(save_enrich(
    clusterProfiler::enrichGO(ev, OrgDb = org.Hs.eg.db, ont = "ALL",
                               pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE),
    "GO"), error = function(e) NULL)

  tryCatch(save_enrich(
    clusterProfiler::enrichKEGG(ev, organism = "hsa", pvalueCutoff = 0.05),
    "KEGG"), error = function(e) NULL)

  tryCatch(save_enrich(
    ReactomePA::enrichPathway(ev, pvalueCutoff = 0.05, readable = TRUE),
    "Reactome"), error = function(e) NULL)

  tryCatch(save_enrich(
    DOSE::enrichDO(ev, pvalueCutoff = 0.05),
    "DiseaseOntology"), error = function(e) NULL)
}


# ==============================================================================
# SECTION 3: Data Import and Graph Construction
# ==============================================================================

cat("--- [2] Loading network files ---\n")

sif_file <- paste0(NETWORK_PREFIX, ".sif")
att_file <- paste0(NETWORK_PREFIX, ".csv")

if (!file.exists(sif_file)) stop("SIF file not found: ", sif_file)
if (!file.exists(att_file)) stop("Attribute file not found: ", att_file)

# Parse SIF edge list
raw_lines <- readLines(sif_file)
raw_lines <- raw_lines[nzchar(trimws(raw_lines))]
edges <- do.call(rbind, lapply(raw_lines, function(l) {
  tok <- strsplit(trimws(l), "\\s+")[[1]]
  if (length(tok) >= 2) data.frame(from = toupper(trimws(tok[1])),
                                    to   = toupper(trimws(tok[length(tok)])),
                                    stringsAsFactors = FALSE)
})) %>% dplyr::distinct() %>% dplyr::filter(from != to)

# Parse node attributes
nodes_raw  <- read.csv(att_file, check.names = FALSE)
cols       <- colnames(nodes_raw)
name_col   <- grep("^(shared[ ._]?)?name$", cols, ignore.case = TRUE, value = TRUE)[1]
type_col   <- grep("^(shared[ ._]?)?type$", cols, ignore.case = TRUE, value = TRUE)[1]

if (is.na(name_col) || is.na(type_col))
  stop("Cannot find 'name' or 'type' columns. Available: ", paste(cols, collapse = ", "))

nodes_clean <- nodes_raw %>%
  dplyr::select(name = all_of(name_col), Type = all_of(type_col)) %>%
  dplyr::mutate(name = toupper(trimws(as.character(name))),
                Type = trimws(as.character(Type))) %>%
  dplyr::distinct(name, .keep_all = TRUE)

# Add missing nodes
missing <- setdiff(unique(c(edges$from, edges$to)), nodes_clean$name)
if (length(missing) > 0) {
  nodes_clean <- dplyr::bind_rows(nodes_clean, tibble::tibble(name = missing, Type = "Unknown"))
  cat("Added", length(missing), "missing nodes as 'Unknown'\n")
}

# Build undirected graph; extract largest component
g_full <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes_clean)
subs   <- igraph::decompose(g_full)
g_main <- subs[[which.max(sapply(subs, igraph::vcount))]]

cat("Network built — Nodes:", igraph::vcount(g_main),
    "| Edges:", igraph::ecount(g_main),
    "| Components:", length(subs), "\n\n")


# ==============================================================================
# SECTION 3: Topology Analysis and Network Visualisation
# ==============================================================================

if (RUN_SECTIONS["topology"]) {
  cat("--- [3] Topology Analysis ---\n")

  topo_dir <- file.path(out_dir, "01_topology")
  plot_dir <- file.path(out_dir, "02_network_plots")
  dir.create(topo_dir, showWarnings = FALSE)
  dir.create(plot_dir, showWarnings = FALSE)

  global_summaries <- list()

  for (i in seq_along(subs)) {
    sg <- subs[[i]]
    if (igraph::vcount(sg) <= 5) next

    cat("  Subnetwork", i, "— Nodes:", igraph::vcount(sg),
        "| Edges:", igraph::ecount(sg), "\n")

    # Global topology metrics
    global_summaries[[i]] <- tibble::tibble(
      Subnetwork_ID          = paste0("Subnetwork_", i),
      Nodes                  = igraph::vcount(sg),
      Edges                  = igraph::ecount(sg),
      Density                = igraph::edge_density(sg),
      Avg_Path_Length        = igraph::mean_distance(sg),
      Diameter               = igraph::diameter(sg),
      Clustering_Coefficient = igraph::transitivity(sg, type = "global"),
      Degree_Centralization  = igraph::centr_degree(sg)$centralization
    )

    # Per-node metrics
    comms <- igraph::cluster_walktrap(sg)
    igraph::V(sg)$community <- comms$membership

    node_metrics <- tibble::tibble(
      node_name              = igraph::V(sg)$name,
      degree                 = igraph::degree(sg),
      betweenness_centrality = igraph::betweenness(sg, normalized = TRUE),
      closeness_centrality   = igraph::closeness(sg, normalized = TRUE),
      eigenvector_centrality = igraph::eigen_centrality(sg)$vector,
      community              = igraph::V(sg)$community
    )
    write.csv(node_metrics,
              file.path(topo_dir, paste0("subnetwork_", i, "_node_metrics.csv")),
              row.names = FALSE)

    # Network visualisation (label top-10% degree nodes)
    if (!"Degree" %in% igraph::vertex_attr_names(sg)) igraph::V(sg)$Degree <- igraph::degree(sg)
    if (!"Type"   %in% igraph::vertex_attr_names(sg)) igraph::V(sg)$Type   <- "Unknown"

    set.seed(42)
    p_net <- ggraph::ggraph(sg, layout = "fr") +
      ggraph::geom_edge_fan(alpha = 0.2, width = 0.4) +
      ggraph::geom_node_point(ggplot2::aes(size = Degree, colour = as.factor(Type),
                                           shape = as.factor(community))) +
      ggraph::geom_node_text(
        ggplot2::aes(label = ifelse(Degree > stats::quantile(Degree, 0.90), name, "")),
        repel = TRUE, size = 3
      ) +
      ggraph::theme_graph(base_family = "sans") +
      ggplot2::labs(title    = paste("Subnetwork", i),
                    subtitle = paste(igraph::vcount(sg), "nodes |",
                                     igraph::ecount(sg), "edges |",
                                     length(comms), "communities"),
                    size = "Degree", colour = "Node Type", shape = "Community")

    ggplot2::ggsave(file.path(plot_dir, paste0("subnetwork_", i, "_plot.png")),
                    p_net, width = 12, height = 10, dpi = 300)
  }

  write.csv(dplyr::bind_rows(global_summaries),
            file.path(topo_dir, "global_topology_summary.csv"), row.names = FALSE)
  cat("  Topology results saved\n\n")
}


# ==============================================================================
# SECTION 4: Centrality-Based Hub Identification and Hub-Target Table
# ==============================================================================

if (RUN_SECTIONS["hubs"]) {
  cat("--- [4] Hub Identification ---\n")

  hub_dir <- file.path(out_dir, "03_hubs")
  dir.create(hub_dir, showWarnings = FALSE)

  igraph::V(g_main)$degree      <- igraph::degree(g_main,      normalized = TRUE)
  igraph::V(g_main)$betweenness <- igraph::betweenness(g_main, normalized = TRUE)
  igraph::V(g_main)$closeness   <- igraph::closeness(g_main,   normalized = TRUE)

  centrality_df <- tibble::tibble(
    name        = igraph::V(g_main)$name,
    Type        = igraph::V(g_main)$Type,
    Degree      = igraph::V(g_main)$degree,
    Betweenness = igraph::V(g_main)$betweenness,
    Closeness   = igraph::V(g_main)$closeness
  ) %>%
    dplyr::mutate(
      Composite_Score = (scale(Degree) + scale(Betweenness) + scale(Closeness)) / 3
    ) %>%
    dplyr::arrange(dplyr::desc(Composite_Score))

  write.csv(head(centrality_df, NUM_TOP_HUBS),
            file.path(hub_dir, "01_top_hubs_composite_score.csv"), row.names = FALSE)

  # Top hubs by individual measures
  for (metric in c("Degree", "Betweenness", "Closeness")) {
    write.csv(
      head(dplyr::arrange(centrality_df, dplyr::desc(.data[[metric]])), NUM_TOP_HUBS),
      file.path(hub_dir, paste0("01_top_hubs_by_", metric, ".csv")),
      row.names = FALSE
    )
  }

  # Hub-target table: top N hubs per node type with all their direct targets
  node_type_map <- setNames(nodes_clean$Type, nodes_clean$name)

  top_per_type <- centrality_df %>%
    dplyr::filter(Type != "Unknown") %>%
    dplyr::group_by(Type) %>%
    dplyr::slice_max(order_by = Composite_Score, n = NUM_HUBS_TO_DISPLAY, with_ties = FALSE) %>%
    dplyr::ungroup()

  target_edges <- edges %>%
    dplyr::mutate(Target_Type = node_type_map[to]) %>%
    dplyr::filter(!is.na(Target_Type))

  hub_target_table <- top_per_type %>%
    dplyr::select(Source_Node = name, Source_Type = Type) %>%
    dplyr::left_join(
      dplyr::rename(target_edges, Source_Node = from, Target_Node = to),
      by = "Source_Node"
    ) %>%
    dplyr::arrange(Source_Type, Source_Node, Target_Node)

  write.csv(hub_target_table,
            file.path(hub_dir, paste0("02_top", NUM_HUBS_TO_DISPLAY, "_hubs_per_type_targets.csv")),
            row.names = FALSE)

  cat("  Hub tables saved\n\n")
  hub_genes <- centrality_df$name[1:min(NUM_TOP_HUBS, nrow(centrality_df))]
}


# ==============================================================================
# SECTION 5: K-Means Clustering and Node Metric Comparison
# ==============================================================================

if (RUN_SECTIONS["clustering"] && exists("hub_genes")) {
  cat("--- [5] K-Means Clustering + Node Metric Comparison ---\n")

  clust_dir <- file.path(out_dir, "04_clustering")
  dir.create(clust_dir, showWarnings = FALSE)

  net_csv   <- read.csv(att_file, check.names = FALSE)
  num_cols  <- net_csv %>%
    dplyr::select(dplyr::where(is.numeric)) %>%
    dplyr::select(dplyr::where(~ stats::sd(.x, na.rm = TRUE) > 0))

  if (ncol(num_cols) >= 2) {
    scaled_data <- scale(num_cols)
    n_distinct  <- nrow(unique(scaled_data))

    if (n_distinct >= 3) {
      max_k    <- min(10, n_distinct - 1)
      p_elbow  <- factoextra::fviz_nbclust(scaled_data, stats::kmeans,
                                            method = "wss", k.max = max_k) +
        ggplot2::labs(title = "Elbow Method for Optimal k")
      ggplot2::ggsave(file.path(clust_dir, "kmeans_elbow.png"),
                      p_elbow, width = 8, height = 6)

      k_use <- if (KMEANS_K > 0) KMEANS_K else min(4, max_k)
      cat("  Using k =", k_use, "\n")
      set.seed(123)
      km     <- stats::kmeans(scaled_data, centers = k_use, nstart = 25)
      net_csv$cluster <- km$cluster

      clust_summary <- net_csv %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(dplyr::across(dplyr::where(is.numeric), ~ mean(.x, na.rm = TRUE)),
                         count = dplyr::n(), .groups = "drop")
      write.csv(clust_summary, file.path(clust_dir, "cluster_summary.csv"), row.names = FALSE)

      # Node type metric comparison (violin + Kruskal-Wallis)
      type_col_name <- grep("^type$", colnames(net_csv), ignore.case = TRUE, value = TRUE)[1]
      if (!is.na(type_col_name)) {
        metrics <- c("Degree", "BetweennessCentrality", "ClosenessCentrality", "ClusteringCoefficient")
        avail   <- intersect(metrics, colnames(net_csv))
        if (length(avail) > 0) {
          comp_data <- net_csv %>%
            dplyr::select(Type = all_of(type_col_name), dplyr::all_of(avail)) %>%
            tidyr::pivot_longer(-Type, names_to = "Metric", values_to = "Value")

          p_comp <- ggplot2::ggplot(comp_data,
                                    ggplot2::aes(x = Type, y = Value, fill = Type)) +
            ggplot2::geom_violin() +
            ggplot2::geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
            ggplot2::facet_wrap(~Metric, scales = "free_y") +
            ggpubr::stat_compare_means(label = "p.format", method = "kruskal.test",
                                       label.y.npc = 0.9) +
            ggplot2::labs(title = "Node Metric Comparison by Type") +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                           legend.position = "none")
          ggplot2::ggsave(file.path(clust_dir, "node_metrics_by_type.png"),
                          p_comp, width = 10, height = 8)
        }
      }
    }
  }
  cat("  Clustering saved\n\n")
}


# ==============================================================================
# SECTION 6: Functional Enrichment (GO, KEGG, Reactome, DO)
#   Four sub-analyses: all nodes, miRNA targets, TF-only, non-TF genes
# ==============================================================================

if (RUN_SECTIONS["enrichment"]) {
  cat("--- [6] Functional Enrichment ---\n")

  enr_dir <- file.path(out_dir, "05_enrichment")
  dir.create(enr_dir, showWarnings = FALSE)

  all_node_names <- toupper(igraph::V(g_main)$name)
  mirna_nodes    <- all_node_names[grepl("^HSA-MIR-|^HSA-LET-", all_node_names)]
  non_mirna      <- setdiff(all_node_names, mirna_nodes)

  # --- 6A. All gene/TF nodes ---
  cat("  Enrichment: all gene/TF nodes\n")
  run_full_enrichment(non_mirna, "AllGenesAndTFs", file.path(enr_dir, "all_genes_TFs"))

  # --- 6B. TF-only nodes (identified via GO DNA-binding annotation) ---
  cat("  Enrichment: TF nodes\n")
  tryCatch({
    go_tf_genes <- AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0003700",
                                          columns = "SYMBOL", keytype = "GOALL")
    tf_node_names <- intersect(non_mirna, toupper(unique(go_tf_genes$SYMBOL)))
    if (length(tf_node_names) >= 5)
      run_full_enrichment(tf_node_names, "TF_Nodes_Only", file.path(enr_dir, "tf_nodes"))
  }, error = function(e) cat("  WARNING: TF enrichment failed:", conditionMessage(e), "\n"))

  # --- 6C. Non-TF gene nodes ---
  cat("  Enrichment: non-TF gene nodes\n")
  tryCatch({
    go_tf_genes <- AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0003700",
                                          columns = "SYMBOL", keytype = "GOALL")
    nontf_names <- setdiff(non_mirna, toupper(unique(go_tf_genes$SYMBOL)))
    if (length(nontf_names) >= 5)
      run_full_enrichment(nontf_names, "NonTF_Gene_Nodes", file.path(enr_dir, "non_tf_genes"))
  }, error = function(e) NULL)

  # --- 6D. Validated targets of miRNA nodes (requires multiMiR) ---
  if (length(mirna_nodes) > 0 && requireNamespace("multiMiR", quietly = TRUE)) {
    cat("  Enrichment: miRNA validated targets\n")
    tryCatch({
      val_targets <- multiMiR::get_multimir(mirna = mirna_nodes, summary = TRUE,
                                             table = "validated")
      targ_syms   <- unique(val_targets@data$target_symbol)
      if (length(targ_syms) >= 5)
        run_full_enrichment(targ_syms, "miRNA_Validated_Targets",
                            file.path(enr_dir, "mirna_targets"))
    }, error = function(e) cat("  WARNING: miRNA target enrichment failed:", conditionMessage(e), "\n"))
  }

  # --- 6E. Per-community enrichment (Walktrap communities) ---
  cat("  Enrichment: per-community\n")
  comm_dir <- file.path(out_dir, "06_community_enrichment")
  dir.create(comm_dir, showWarnings = FALSE)

  comms  <- igraph::cluster_walktrap(g_main)
  comm_vec <- igraph::membership(comms)

  for (cid in sort(unique(comm_vec))) {
    c_nodes <- names(comm_vec[comm_vec == cid])
    c_genes <- setdiff(c_nodes, mirna_nodes)
    if (length(c_genes) >= 5) {
      run_full_enrichment(c_genes, paste0("Comm", cid),
                          file.path(comm_dir, paste0("community_", cid)))
    }
  }

  cat("  Enrichment analysis complete\n\n")
}


# ==============================================================================
# SECTION 7: GSEA with Full Visualisation Suite
# ==============================================================================

if (RUN_SECTIONS["gsea"] && file.exists(RANKED_LIST_FILE)) {
  cat("--- [7] GSEA Visualisations ---\n")

  gsea_dir <- file.path(out_dir, "07_gsea")
  dir.create(gsea_dir, showWarnings = FALSE)

  # Load ranked list
  rl_raw <- tryCatch(
    read.table(RANKED_LIST_FILE, header = FALSE, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) read.csv(RANKED_LIST_FILE, header = FALSE, stringsAsFactors = FALSE)
  )
  gene_list_gsea <- setNames(as.numeric(rl_raw[[2]]), rl_raw[[1]])
  gene_list_gsea <- sort(gene_list_gsea[!is.na(gene_list_gsea)], decreasing = TRUE)

  cat("  Ranked list loaded:", length(gene_list_gsea), "genes\n")

  # Entrez ID conversion
  ids_gsea <- tryCatch(
    clusterProfiler::bitr(names(gene_list_gsea), fromType = "SYMBOL",
                           toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )

  if (!is.null(ids_gsea) && nrow(ids_gsea) > 10) {
    gene_list_entrez <- gene_list_gsea[ids_gsea$SYMBOL]
    names(gene_list_entrez) <- ids_gsea$ENTREZID
    gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

    tryCatch({
      # GO GSEA
      gsea_go <- clusterProfiler::gseGO(
        gene_list_entrez, OrgDb = org.Hs.eg.db, ont = "ALL",
        pvalueCutoff = 0.05, verbose = FALSE
      )

      if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
        gsea_clean         <- gsea_go
        gsea_clean@result  <- gsea_go@result %>%
          dplyr::filter(setSize >= 10) %>%
          dplyr::arrange(p.adjust)

        # Ridge plot
        p_ridge <- enrichplot::ridgeplot(gsea_clean, showCategory = 20) +
          ggplot2::labs(title = "Top Enriched GO Pathways — GSEA Ridgeplot")
        ggplot2::ggsave(file.path(gsea_dir, "GSEA_GO_Ridgeplot.png"),
                        p_ridge, width = 10, height = 8, dpi = 300)

        # Dot plot
        p_dot <- enrichplot::dotplot(gsea_clean, showCategory = 20,
                                      split = "ONTOLOGY") +
          ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free_y") +
          ggplot2::labs(title = "GO GSEA Dot Plot")
        ggplot2::ggsave(file.path(gsea_dir, "GSEA_GO_Dotplot.png"),
                        p_dot, width = 10, height = 12, dpi = 300)

        # Concept network (cnetplot) for top 5 terms
        p_cnet <- enrichplot::cnetplot(gsea_clean, showCategory = 5,
                                        foldChange = gene_list_gsea)
        ggplot2::ggsave(file.path(gsea_dir, "GSEA_GO_ConceptNetwork.png"),
                        p_cnet, width = 12, height = 10, dpi = 300)

        # UpSet plot
        if (requireNamespace("ggupset", quietly = TRUE)) {
          p_up <- enrichplot::upsetplot(gsea_clean, n = 10)
          ggplot2::ggsave(file.path(gsea_dir, "GSEA_GO_UpSet.png"),
                          p_up, width = 10, height = 6, dpi = 300)
        }

        # Heatmap plot
        p_heat <- enrichplot::heatplot(gsea_clean, showCategory = 10,
                                        foldChange = gene_list_gsea)
        ggplot2::ggsave(file.path(gsea_dir, "GSEA_GO_Heatmap.png"),
                        p_heat, width = 14, height = 8, dpi = 300)

        cat("  GO GSEA plots saved\n")
      }
    }, error = function(e) cat("  WARNING: GO GSEA failed:", conditionMessage(e), "\n"))
  }
  cat("  GSEA section complete\n\n")
} else if (RUN_SECTIONS["gsea"]) {
  cat("  Skipping GSEA (ranked list file not found:", RANKED_LIST_FILE, ")\n\n")
}


# ==============================================================================
# SECTION 8: Directed Signed Network Propagation (RWR)
# ==============================================================================

if (RUN_SECTIONS["propagation"]) {
  cat("--- [8] Signed RWR Propagation ---\n")

  prop_dir <- file.path(out_dir, "08_propagation")
  dir.create(prop_dir, showWarnings = FALSE)

  # Build signed directed graph
  # Signed rules:
  #   TF  → Gene  : +1 (transcriptional activation)
  #   miRNA → Gene: -1 (post-transcriptional inhibition)
  #   TF  → miRNA : +1 (transcriptional induction)
  #   Gene → miRNA: +1 (indirect activation)
  #   all others  :  0 (excluded)
  node_type_map <- setNames(nodes_clean$Type, nodes_clean$name)

  edges_signed <- edges %>%
    dplyr::mutate(
      type_from = node_type_map[from],
      type_to   = node_type_map[to],
      weight_signed = dplyr::case_when(
        grepl("TF",    type_from, ignore.case = TRUE) &
          grepl("Gene",  type_to,   ignore.case = TRUE) ~ 1,
        grepl("miRNA", type_from, ignore.case = TRUE) &
          grepl("Gene",  type_to,   ignore.case = TRUE) ~ -1,
        grepl("Gene",  type_from, ignore.case = TRUE) &
          grepl("miRNA", type_to,   ignore.case = TRUE) ~ 1,
        grepl("TF",    type_from, ignore.case = TRUE) &
          grepl("miRNA", type_to,   ignore.case = TRUE) ~ 1,
        TRUE ~ 0
      )
    ) %>%
    dplyr::filter(weight_signed != 0)

  if (nrow(edges_signed) == 0)
    stop("No signed edges found. Check that node Types contain 'TF', 'miRNA', 'Gene'.")

  g_signed <- igraph::graph_from_data_frame(
    edges_signed %>% dplyr::select(from, to, weight_signed),
    directed = TRUE,
    vertices = nodes_clean
  )
  subs_signed <- igraph::decompose(g_signed, mode = "weak")
  g_signed_main <- subs_signed[[which.max(sapply(subs_signed, igraph::vcount))]]

  # Seed selection: use user-provided file OR top hubs from Section 4
  if (nzchar(HUB_SEEDS_FILE) && file.exists(HUB_SEEDS_FILE)) {
    seed_genes <- read.csv(HUB_SEEDS_FILE, header = FALSE, stringsAsFactors = FALSE)[[1]]
    cat("  Seeds from file:", length(seed_genes), "\n")
  } else if (exists("hub_genes")) {
    seed_genes <- hub_genes[hub_genes %in% igraph::V(g_signed_main)$name]
    cat("  Seeds from top hubs:", length(seed_genes), "\n")
  } else {
    # Fall back to highest-degree nodes
    deg_vec    <- igraph::degree(g_signed_main, mode = "all")
    seed_genes <- names(sort(deg_vec, decreasing = TRUE))[1:min(20, length(deg_vec))]
    cat("  Seeds: top-degree nodes (", length(seed_genes), ")\n")
  }

  valid_seeds <- intersect(seed_genes, igraph::V(g_signed_main)$name)
  p_vec <- setNames(rep(0, igraph::vcount(g_signed_main)), igraph::V(g_signed_main)$name)
  p_vec[valid_seeds] <- 1

  scores <- signed_influence_propagation(g_signed_main, p_vec, alpha = RWR_ALPHA)

  prop_df <- tibble::tibble(name = names(scores), score = scores) %>%
    dplyr::arrange(dplyr::desc(abs(score)))
  write.csv(prop_df, file.path(prop_dir, "signed_propagation_scores.csv"), row.names = FALSE)

  # Visualise: diverging colour scale (blue=inhibited, red=activated), seeds highlighted
  igraph::V(g_signed_main)$prop_score <- scores[igraph::V(g_signed_main)$name]
  igraph::V(g_signed_main)$node_deg   <- igraph::degree(g_signed_main)
  max_abs <- max(abs(scores))

  set.seed(42)
  p_prop <- ggraph::ggraph(g_signed_main, layout = "fr") +
    ggraph::geom_edge_fan(ggplot2::aes(colour = factor(weight_signed), alpha = 0.3),
                          width = 0.5) +
    ggraph::geom_node_point(ggplot2::aes(size = node_deg, colour = prop_score)) +
    ggraph::geom_node_point(
      data = . %>% dplyr::filter(name %in% valid_seeds),
      ggplot2::aes(size = node_deg),
      colour = "black", shape = 23, fill = "yellow", stroke = 1.5
    ) +
    ggplot2::scale_colour_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = c(-max_abs, max_abs), name = "Influence Score"
    ) +
    ggraph::scale_edge_colour_manual(
      values = c("-1" = "blue", "1" = "red"), name = "Interaction"
    ) +
    ggraph::theme_graph() +
    ggplot2::labs(
      title    = paste("Signed RWR from", length(valid_seeds), "hub seeds"),
      subtitle = "Blue = inhibited | Red = activated | Yellow diamond = seed"
    )

  ggplot2::ggsave(file.path(prop_dir, "signed_propagation_network.png"),
                  p_prop, width = 14, height = 12, dpi = 300)
  cat("  Propagation complete\n\n")
}


# ==============================================================================
# SECTION 9: Structural Perturbation (Node Knockout Analysis)
# ==============================================================================

if (RUN_SECTIONS["perturbation"]) {
  cat("--- [9] Structural Perturbation ---\n")

  pert_dir <- file.path(out_dir, "09_perturbation")
  dir.create(pert_dir, showWarnings = FALSE)

  # ------------------------------------------------------------------------------
  # calc_largest_component()
  #   Returns the number of nodes in the largest connected component of a graph.
  # Args: g — igraph object
  # Returns: integer
  # ------------------------------------------------------------------------------
  calc_largest_component <- function(g) {
    if (igraph::vcount(g) == 0) return(0)
    g_u <- igraph::as.undirected(g, mode = "collapse")
    max(igraph::components(g_u)$csize)
  }

  # Rank by raw degree for perturbation (undirected)
  deg_raw   <- igraph::degree(g_main)
  top_nodes <- names(sort(deg_raw, decreasing = TRUE))[1:min(NUM_HUBS_TO_PERTURB,
                                                              igraph::vcount(g_main))]
  orig_size <- calc_largest_component(g_main)

  pert_results <- dplyr::bind_rows(lapply(top_nodes, function(nd) {
    g_pert <- igraph::delete_vertices(g_main, igraph::V(g_main)[name == nd])
    new_size <- calc_largest_component(g_pert)
    tibble::tibble(
      Removed_Node                 = nd,
      Remaining_Component_Size     = new_size,
      Nodes_Lost                   = orig_size - new_size - 1,
      Pct_Main_Component_Lost      = ((orig_size - new_size - 1) / (orig_size - 1)) * 100
    )
  })) %>%
    dplyr::arrange(dplyr::desc(Pct_Main_Component_Lost))

  write.csv(pert_results, file.path(pert_dir, "perturbation_results.csv"), row.names = FALSE)

  p_pert <- ggplot2::ggplot(pert_results,
                             ggplot2::aes(x = reorder(Removed_Node, Pct_Main_Component_Lost),
                                          y = Pct_Main_Component_Lost)) +
    ggplot2::geom_col(fill = "firebrick") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(title = "Node Criticality (Structural Perturbation)",
                  x = "Removed Hub", y = "% Main Component Lost")

  ggplot2::ggsave(file.path(pert_dir, "node_criticality_barplot.png"),
                  p_pert, width = 9, height = 7, dpi = 300)
  cat("  Perturbation complete\n\n")
}


# ==============================================================================
# SECTION 10: Per-Hub Neighborhood Enrichment
#   For each hub gene, extract its direct neighbors and run GO/KEGG/Reactome/DO
# ==============================================================================

if (RUN_SECTIONS["subnetwork"] && exists("hub_genes")) {
  cat("--- [10] Per-Hub Neighborhood Enrichment ---\n")

  sub_dir <- file.path(out_dir, "10_subnetwork")
  dir.create(sub_dir, showWarnings = FALSE)

  # Limit to top N hubs to avoid excessive runtime
  hubs_for_sub <- head(hub_genes, min(NUM_TOP_HUBS, 20))

  for (hub in hubs_for_sub) {
    if (!(hub %in% igraph::V(g_main)$name)) next
    cat("  Subnetwork enrichment:", hub, "\n")

    hub_nbrs    <- igraph::neighbors(g_main, igraph::V(g_main)[name == hub])
    nbr_names   <- igraph::V(g_main)[hub_nbrs]$name
    nbr_types   <- igraph::V(g_main)[hub_nbrs]$Type
    gene_nbrs   <- nbr_names[!grepl("miRNA", nbr_types, ignore.case = TRUE)]

    if (length(gene_nbrs) < 5) {
      cat("    Skipping (fewer than 5 gene neighbors)\n"); next
    }

    hub_out_dir <- file.path(sub_dir, paste0("hub_", gsub("[^[:alnum:]_]", "_", hub)))
    run_full_enrichment(gene_nbrs, paste0(hub, "_Neighborhood"), hub_out_dir)
  }

  cat("  Per-hub enrichment complete\n\n")
}

cat("=== NETWORK ANALYSIS PIPELINE COMPLETE ===\n")
cat("All results saved to:", out_dir, "\n")
