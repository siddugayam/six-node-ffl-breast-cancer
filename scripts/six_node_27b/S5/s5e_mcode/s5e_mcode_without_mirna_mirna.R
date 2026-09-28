## S5(e): MCODE on the six-node composite FFL network of INBOX_2026-09-26c (Bhat et al. 2024 Fig. 1d,
## nolegacy graph) with and without its miRNA-miRNA (10 kb co-transcription) links.  Same pipeline and
## parameters as INBOX_2026-09-26c/bhat_networks/mcode_six_node_composite.R: scripts/v7/mcode.R
## (Cytoscape-faithful port), degree cutoff 2, k-core 2, max depth 100, node score cutoff 0.2, haircut,
## no fluff; undirected, multi-edges collapsed; node order = first appearance in the ORIGINAL SIF (the
## tie-break order of the 26c run), plus 100 random node orders.  Read-only on the project.
REV  <- "/path/to/revision"
SRC  <- file.path(REV, "INBOX_2026-09-26c/bhat_networks")
OUT  <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(REV, "scripts/v7/mcode.R"))
suppressPackageStartupMessages(library(igraph))
set.seed(20260927)
sif0 <- read.delim(file.path(SRC, "sif", "nolegacy_6node_Composite_FFL.sif"), header = FALSE,
                   col.names = c("s", "lab", "t"), stringsAsFactors = FALSE)
attr <- read.delim(file.path(SRC, "sif", "nolegacy_node_attributes.tsv"), stringsAsFactors = FALSE)
ord0 <- unique(as.vector(t(as.matrix(sif0[, c("s", "t")]))))
for (variant in c("with_miRNA_miRNA_links", "without_miRNA_miRNA_links")) {
  sif <- if (variant == "with_miRNA_miRNA_links") sif0 else sif0[sif0$lab != "miRNA-miRNA", ]
  ord <- ord0[ord0 %in% unique(c(sif$s, sif$t))]
  g <- graph_from_data_frame(sif[, c("s", "t")], directed = FALSE, vertices = data.frame(name = ord))
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  V(g)$type <- attr$type[match(V(g)$name, attr$node)]
  cat(sprintf("\n== %s: %d SIF rows (%d miRNA-miRNA), %d nodes, %d undirected edges\n", variant, nrow(sif),
              sum(sif$lab == "miRNA-miRNA"), vcount(g), ecount(g)))
  info <- mcode_score_graph(g)
  cl <- mcode_find_clusters(g, info, degree_cutoff = 2L, k_core = 2L, max_depth = 100L,
                            node_score_cutoff = 0.2, haircut = TRUE, fluff = FALSE)
  nm <- V(g)$name; ty <- V(g)$type
  memb <- function(v, t) paste(sort(nm[v][ty[v] == t]), collapse = ";")
  tab <- do.call(rbind, lapply(cl, function(x) data.frame(
    rank = x$rank, score = round(x$score, 4), nodes = x$size, edges = x$edges, seed = nm[x$seed],
    n_TF = sum(ty[x$nodes] == "TF"), n_miRNA = sum(ty[x$nodes] == "miRNA"), n_Gene = sum(ty[x$nodes] == "Gene"),
    TF = memb(x$nodes, "TF"), miRNA = memb(x$nodes, "miRNA"), Gene = memb(x$nodes, "Gene"))))
  write.csv(tab, file.path(OUT, paste0("s5e_mcode_clusters_", variant, ".csv")), row.names = FALSE)
  print(head(tab, 6), row.names = FALSE)
  # where do MYC, E2F1 and the miR-17~92 members end up?
  key <- c("MYC", "E2F1", "hsa-miR-17", "hsa-miR-18a", "hsa-miR-19a", "hsa-miR-19b", "hsa-miR-20a", "hsa-miR-92a")
  loc <- sapply(key, function(k) { if (!k %in% nm) return("not in network")
    r <- which(sapply(cl, function(x) k %in% nm[x$nodes])); if (length(r)) paste(sapply(cl[r], `[[`, "rank"), collapse = ",") else "no cluster" })
  cat("   cluster rank of key nodes:", paste(names(loc), loc, sep = "=", collapse = "; "), "\n")
  same <- 0L; together <- 0L
  for (i in 1:100) {
    ci <- mcode_find_clusters(g, info, seed_order = sample(vcount(g)))
    same <- same + as.integer(setequal(ci[[1]]$nodes, cl[[1]]$nodes))
    tog <- any(sapply(ci, function(x) all(c("MYC", "E2F1") %in% nm[x$nodes]) &&
                                      any(key[3:8] %in% nm[x$nodes])))
    together <- together + as.integer(tog)
  }
  cat(sprintf("   top cluster identical in %d of 100 random node orders; MYC, E2F1 and >=1 miR-17~92 member in one cluster in %d of 100\n",
              same, together))
  writeLines(c(paste("variant", variant), paste("cluster_rank_of_key_nodes", paste(names(loc), loc, sep = "=", collapse = "; ")),
               paste("top_cluster_identical_in_random_orders", same, "of 100"),
               paste("MYC_E2F1_and_miR17_92_member_together_in_random_orders", together, "of 100")),
             file.path(OUT, paste0("s5e_summary_", variant, ".txt")))
}
