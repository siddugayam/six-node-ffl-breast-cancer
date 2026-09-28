#!/usr/bin/env Rscript
## 43_farmer_module_overlap.R  -- THE KEY TEST
## How much of the higher-order FFL module IS the Farmer reactive-stroma programme?
##  (1) gene-level overlap, hypergeometric p, two universes
##  (2) is the higher-order module more stromal than the rest of the network? (Fisher, network universe)
##  (3) score-level correlation, including after making the two sets disjoint
##  (4) per-gene correlation of module members with the Farmer score vs matched random null
##  (5) CAF-compartment enrichment of module members
##  (6) robustness across the 49 Farmer leave-one-out variant metagenes
suppressPackageStartupMessages({library(data.table); library(matrixStats)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v5/43_farmer_module_overlap.log")
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); set.seed(2024)

sets <- readRDS(file.path(BASE, "cache/v5/farmer_signature_sets.rds"))
sc   <- readRDS(file.path(BASE, "cache/v5/farmer_tcga_scores.rds"))
gex  <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
nodes<- fread(file.path(BASE, "data/canonical_nodes.tsv"))
S    <- sc$samples
X    <- gex[, S, drop = FALSE]
X    <- X[rowSums(is.na(X)) == 0 & rowSds(X) > 0, , drop = FALSE]
UNIV <- rownames(X)                       # 20,248 measurable genes
prot <- intersect(nodes[type != "miRNA", name], UNIV)
logf("universe A (TCGA-measurable genes): ", length(UNIV))
logf("universe B (protein-coding network nodes measurable in TCGA): ", length(prot))

restr <- function(v) intersect(v, UNIV)
MODS <- c("FFL_3node", "FFL_4node", "FFL_5node", "FFL_6node_all", "FFL_higher_order_only",
          "MS_exemplar_module")
SIGS <- c("FARMER_STROMAL", "WEST_DTF_FIBROMATOSIS", "WEST_SFT", "WINSLOW_STROMAL_SIG1",
          "FINAK_SDPP", "FINAK_REF_STROMA_I", "FARMER2005_STROMAL_CLUSTER4",
          "CHANG_CSR_LATE_UP", "CHANG_CSR_EARLY_UP", "TRIULZI_ECM", "NABA_CORE_MATRISOME",
          "NABA_COLLAGENS", "HALLMARK_EMT", "HALLMARK_TGF_BETA", "HALLMARK_ANGIOGENESIS",
          "ESTIMATE_STROMAL", "ESTIMATE_STROMAL_noCOL", "CAF_scRNA_50", "ESTIMATE_IMMUNE",
          "FINAK_REF_EPITHELIUM_E", "MAMMOSPHERE")

## ---------------------------------------------------------------- (1) overlap tests ----
ov <- list(); k <- 0L
for (m in MODS) for (s in SIGS) {
  A <- restr(sets[[m]]); B <- restr(sets[[s]])
  int <- intersect(A, B); N <- length(UNIV)
  p_hyper <- phyper(length(int) - 1L, length(B), N - length(B), length(A), lower.tail = FALSE)
  exp_ov  <- length(A) * length(B) / N
  k <- k + 1L
  ov[[k]] <- data.table(module = m, signature = s, n_module = length(A), n_signature = length(B),
    n_overlap = length(int), expected = round(exp_ov, 2),
    fold_enrichment = round(length(int) / exp_ov, 2),
    frac_of_module = round(length(int) / length(A), 4),
    frac_of_signature = round(length(int) / length(B), 4),
    jaccard = round(length(int) / length(union(A, B)), 4),
    p_hypergeometric = p_hyper, universe = "TCGA_measurable_genes",
    overlap_genes = paste(sort(int), collapse = ";"))
}
ovt <- rbindlist(ov)
ovt[, q_BH := p.adjust(p_hypergeometric, "BH")]
setorder(ovt, module, p_hypergeometric)
fwrite(ovt, file.path(BASE, "results/v5/farmer_module_overlap.csv"))
logf("WROTE results/v5/farmer_module_overlap.csv  rows=", nrow(ovt))
logf("\n=== overlap with the FARMER 50-gene stroma-related signature ===")
for (m in MODS) {
  r <- ovt[module == m & signature == "FARMER_STROMAL"]
  logf(sprintf("%-24s module n=%3d | overlap %2d/%d (%.1f%% of module, %.1f%% of Farmer) | exp %.2f | fold %.1f | hyperg p=%.3g",
    m, r$n_module, r$n_overlap, r$n_signature, 100*r$frac_of_module, 100*r$frac_of_signature,
    r$expected, r$fold_enrichment, r$p_hypergeometric))
  logf("     genes: ", r$overlap_genes)
}

## ---------------- (2) network-conditioned: higher-order vs the rest of the network ------
logf("\n=== is the higher-order module more stromal than the rest of the network? ===")
net_rows <- list(); k2 <- 0L
for (m in c("FFL_higher_order_only", "FFL_3node", "FFL_6node_all")) for (s in SIGS) {
  A <- intersect(sets[[m]], prot); B <- restr(sets[[s]])
  inA_inB <- length(intersect(A, B)); inA_out <- length(A) - inA_inB
  rest <- setdiff(prot, A)
  out_inB <- length(intersect(rest, B)); out_out <- length(rest) - out_inB
  ft <- fisher.test(matrix(c(inA_inB, inA_out, out_inB, out_out), 2, 2), alternative = "greater")
  k2 <- k2 + 1L
  net_rows[[k2]] <- data.table(module = m, signature = s, universe = "network_protein_nodes",
    n_universe = length(prot), n_module = length(A), n_module_in_sig = inA_inB,
    frac_module_in_sig = round(inA_inB / length(A), 4),
    n_rest = length(rest), n_rest_in_sig = out_inB,
    frac_rest_in_sig = round(out_inB / length(rest), 4),
    odds_ratio = unname(ft$estimate), p_fisher_greater = ft$p.value)
}
nett <- rbindlist(net_rows)
nett[, q_BH := p.adjust(p_fisher_greater, "BH"), by = module]
fwrite(nett, file.path(BASE, "results/v5/farmer_module_vs_network_background.csv"))
for (m in unique(nett$module)) {
  r <- nett[module == m & signature == "FARMER_STROMAL"]
  logf(sprintf("%-24s %d/%d (%.1f%%) of module members are Farmer genes vs %d/%d (%.1f%%) of the rest of the network | OR=%.2f p=%.3g",
    m, r$n_module_in_sig, r$n_module, 100*r$frac_module_in_sig,
    r$n_rest_in_sig, r$n_rest, 100*r$frac_rest_in_sig, r$odds_ratio, r$p_fisher_greater))
}

## ------------------------------------------------- (3) score-level, incl. disjointified --
logf("\n=== score-level relationship (ssGSEA and mean-z) ===")
ssg <- sc$ssgsea; mzz <- sc$meanz
zmean <- function(g) { g <- intersect(g, rownames(X)); if (length(g) < 3) return(NULL)
  Z <- (X[g, , drop = FALSE] - rowMeans(X[g, , drop = FALSE])) / rowSds(X[g, , drop = FALSE])
  colMeans(Z) }
sr <- list(); k3 <- 0L
for (m in MODS) {
  A <- restr(sets[[m]]); B <- restr(sets$FARMER_STROMAL)
  A_only <- setdiff(A, B); B_only <- setdiff(B, A)
  for (meth in c("ssGSEA", "meanZ")) {
    M <- if (meth == "ssGSEA") ssg else mzz
    x <- as.numeric(M[m, ]); y <- as.numeric(M["FARMER_STROMAL", ])
    h  <- cor.test(x, y, method = "spearman", exact = FALSE)
    hp <- cor.test(x, y, method = "pearson")
    ## disjoint versions (mean-z on non-shared members only)
    dz <- NA_real_; dp <- NA_real_
    if (length(A_only) >= 3 && length(B_only) >= 3) {
      xa <- zmean(A_only); yb <- zmean(B_only)
      dz <- cor(xa, yb, method = "spearman"); dp <- cor(xa, yb)
    }
    k3 <- k3 + 1L
    sr[[k3]] <- data.table(module = m, method = meth, n = length(S),
      n_module = length(A), n_shared = length(intersect(A, B)),
      n_module_excl = length(A_only), n_farmer_excl = length(B_only),
      spearman_rho = unname(h$estimate), spearman_p = h$p.value,
      pearson_r = unname(hp$estimate), r2 = unname(hp$estimate)^2,
      spearman_rho_disjoint = dz, pearson_r_disjoint = dp, r2_disjoint = dp^2)
  }
}
srt <- rbindlist(sr)
fwrite(srt, file.path(BASE, "results/v5/farmer_module_score_correlation.csv"))
for (i in seq_len(nrow(srt))) {
  r <- srt[i]
  logf(sprintf("%-24s %-7s rho=%+.3f (p=%.3g) r2=%.3f | disjoint (module-only %d genes vs Farmer-only %d genes) rho=%+.3f r2=%.3f",
    r$module, r$method, r$spearman_rho, r$spearman_p, r$r2,
    r$n_module_excl, r$n_farmer_excl, r$spearman_rho_disjoint, r$r2_disjoint))
}

## ------------------------------- (4) per-gene correlation with the Farmer score + null ---
logf("\n=== per-gene: how stromal is each module member? ===")
fscore <- as.numeric(ssg["FARMER_STROMAL", ])
rho_all <- apply(X, 1, function(v) cor(v, fscore, method = "spearman"))
gene_rows <- data.table(gene = names(rho_all), rho_vs_farmer = as.numeric(rho_all))
gene_rows[, in_farmer := gene %in% sets$FARMER_STROMAL]
for (m in MODS) gene_rows[[m]] <- gene_rows$gene %in% restr(sets[[m]])
fwrite(gene_rows, file.path(BASE, "results/v5/farmer_pergene_rho.csv"))
pg <- list(); k4 <- 0L
B <- 2000L
for (m in MODS) {
  A <- restr(sets[[m]]); r <- rho_all[A]
  obs_med <- median(r); obs_frac <- mean(r > 0.5)
  nullmed <- nullfrac <- numeric(B)
  for (b in seq_len(B)) { s <- sample(UNIV, length(A)); rr <- rho_all[s]
    nullmed[b] <- median(rr); nullfrac[b] <- mean(rr > 0.5) }
  k4 <- k4 + 1L
  pg[[k4]] <- data.table(module = m, n = length(A), median_rho = obs_med,
    null_median_rho = median(nullmed), p_perm_median = (1 + sum(nullmed >= obs_med)) / (B + 1),
    frac_rho_gt_0.5 = obs_frac, null_frac = median(nullfrac),
    p_perm_frac = (1 + sum(nullfrac >= obs_frac)) / (B + 1))
}
pgt <- rbindlist(pg)
fwrite(pgt, file.path(BASE, "results/v5/farmer_pergene_module_null.csv"))
for (i in seq_len(nrow(pgt))) { r <- pgt[i]
  logf(sprintf("%-24s n=%3d median rho vs Farmer score = %+.3f (null %+.3f, perm p=%.4g) | %.1f%% of members rho>0.5 (null %.1f%%, p=%.4g)",
    r$module, r$n, r$median_rho, r$null_median_rho, r$p_perm_median,
    100*r$frac_rho_gt_0.5, 100*r$null_frac, r$p_perm_frac)) }

## --------------------------------------------- (5) CAF compartment of module members ----
ctf <- file.path(BASE, "results/multiomics/celltype_source.csv")
if (file.exists(ctf)) {
  ct <- fread(ctf)
  ct <- ct[!is.na(scrna_log2_CAF_over_CancerEpi)]
  caf_rows <- list(); k5 <- 0L
  for (m in c(MODS, "FARMER_STROMAL")) {
    A <- restr(sets[[m]]); sub <- ct[gene %in% A]
    k5 <- k5 + 1L
    caf_rows[[k5]] <- data.table(set = m, n_with_celltype_data = nrow(sub),
      median_log2_CAF_over_CancerEpi = median(sub$scrna_log2_CAF_over_CancerEpi),
      frac_CAF_enriched_gt1 = mean(sub$scrna_log2_CAF_over_CancerEpi > 1),
      frac_top_celltype_CAF = mean(sub$scrna_top_celltype == "CAFs"))
  }
  caft <- rbindlist(caf_rows)
  fwrite(caft, file.path(BASE, "results/v5/farmer_module_caf_compartment.csv"))
  logf("\n=== CAF compartment (GSE176078 pseudobulk; only network nodes have these data) ===")
  for (i in seq_len(nrow(caft))) { r <- caft[i]
    logf(sprintf("%-24s n=%3d median log2(CAF/cancer-epi)=%+.2f | %.0f%% CAF-enriched (>2x) | %.0f%% CAF-dominant",
      r$set, r$n_with_celltype_data, r$median_log2_CAF_over_CancerEpi,
      100*r$frac_CAF_enriched_gt1, 100*r$frac_top_celltype_CAF)) }
}

## ------------------------------------- (6) robustness over 49 Farmer variant metagenes ---
vf <- list.files(file.path(BASE, "cache/v5/farmer"), pattern = "^farmer2009_variant_", full.names = TRUE)
logf("\n=== robustness: overlap with the 49 leave-one-out Farmer metagenes (Suppl Table 4 cols 2-50) ===")
vr <- list(); k6 <- 0L
sym_fix <- c(CSPG2 = "VCAN", PRSS11 = "HTRA1", MGC3047 = "MXRA8", LOC83468 = NA,
             KIAA0992 = "PLEC", PLEKHC1 = "FERMT2", LEPRE1 = "P3H1", SART2 = NA)
for (f in vf) {
  g <- readLines(f); g <- ifelse(g %in% names(sym_fix), sym_fix[g], g); g <- g[!is.na(g)]
  g <- restr(unique(g))
  for (m in c("FFL_higher_order_only", "FFL_3node", "FFL_6node_all")) {
    A <- restr(sets[[m]]); int <- intersect(A, g)
    p <- phyper(length(int) - 1L, length(g), length(UNIV) - length(g), length(A), lower.tail = FALSE)
    k6 <- k6 + 1L
    vr[[k6]] <- data.table(variant = sub(".*variant_", "", basename(f)), module = m,
      n_variant_in_TCGA = length(g), n_overlap = length(int),
      frac_of_module = length(int) / length(A), p_hypergeometric = p)
  }
}
vrt <- rbindlist(vr)
fwrite(vrt, file.path(BASE, "results/v5/farmer_variant_metagene_overlap.csv"))
for (m in unique(vrt$module)) { s <- vrt[module == m]
  logf(sprintf("%-24s across 49 variants: overlap median=%d (range %d-%d); %% of module median=%.1f%%; max p=%.3g",
    m, median(s$n_overlap), min(s$n_overlap), max(s$n_overlap),
    100*median(s$frac_of_module), max(s$p_hypergeometric))) }
logf("DONE 43")
