#!/usr/bin/env Rscript
# 17_external_cohort_validation.R
# Consolidates every external-cohort comparison into one table and computes the
# replication rates: of the features significant in TCGA-BRCA, how many
# replicate in the same direction in an independent cohort.
# Output: results/external_cohort_validation.csv
#         results/external_cohort_replication_summary.csv

suppressPackageStartupMessages({ library(data.table) })
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs", "17_external_validation.log")
logf <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                        cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); logf("START")

nodes <- fread(file.path(BASE, "data/canonical_nodes.tsv"))
ntype <- setNames(nodes$type, nodes$name)
topo  <- fread(file.path(BASE, "results/network_topology_hubs.csv"))
hubset <- topo[is_hub == TRUE, name]
logf("network hubs (top-50 degree U top-50 betweenness U 19 named) = ", length(hubset))

rows <- list(); r <- 0L

## ======================= 1. GEO tumour-vs-normal DE ========================
geo <- fread(file.path(BASE, "results/external_GEO_DE_vs_TCGA.csv"))
geo_n <- fread(file.path(BASE, "results/external_GEO_DE.csv"))[, .(n_normal = n_normal[1],
                                                                   n_tumour = n_tumour[1]),
                                                               by = gse]
logf("GEO series loaded: ", paste(geo_n$gse, collapse = ", "))
for (i in seq_len(nrow(geo_n))) {
  gs <- geo_n$gse[i]
  d <- geo[gse == gs & is_network_node == TRUE]
  r <- r + 1L
  rows[[r]] <- data.table(
    cohort = gs, accession = gs, platform = "Affymetrix HG-U133 Plus 2.0 (GPL570)",
    n_samples = geo_n$n_normal[i] + geo_n$n_tumour[i],
    n_detail = sprintf("%d normal / %d tumour", geo_n$n_normal[i], geo_n$n_tumour[i]),
    analysis = "DE_tumour_vs_normal",
    feature = d$feature, feature_type = ntype[d$feature],
    is_hub = d$feature %in% hubset,
    effect_TCGA = d$logFC_tcga, sig_TCGA = d$q_tcga,
    effect_external = d$logFC_ext, sig_external = d$q_ext,
    same_direction = sign(d$logFC_tcga) == sign(d$logFC_ext),
    effect_metric = "logFC (tumour vs normal)")
}

## ======================= 2. METABRIC univariable Cox ========================
mbf <- file.path(BASE, "results/external_METABRIC_cox.csv")
if (file.exists(mbf)) {
  mb <- fread(mbf)
  tc <- fread(file.path(BASE, "results/survival_cox_hubs.csv"))
  # TCGA PFI is the closest analogue of METABRIC RFS
  pair <- list(c("OS", "OS"), c("PFI", "RFS"))
  for (p in pair) {
    a <- tc[endpoint == p[1], .(feature, HR_tcga = HR_per_SD, p_tcga = p_value, q_tcga = q_value)]
    b <- mb[endpoint == p[2], .(feature, HR_mb = HR_per_SD, p_mb = p_value, q_mb = q_value,
                                n = n[1], n_event = n_event[1])]
    m <- merge(a, b, by = "feature")
    r <- r + 1L
    rows[[r]] <- data.table(
      cohort = "METABRIC", accession = "cBioPortal brca_metabric",
      platform = "Illumina HT-12 v3 microarray",
      n_samples = m$n[1],
      n_detail = sprintf("%d tumours, %d %s events", m$n[1], m$n_event[1], p[2]),
      analysis = sprintf("Cox_%s_TCGA%s", p[2], p[1]),
      feature = m$feature, feature_type = ntype[m$feature],
      is_hub = m$feature %in% hubset,
      effect_TCGA = m$HR_tcga, sig_TCGA = m$q_tcga,
      effect_external = m$HR_mb, sig_external = m$q_mb,
      same_direction = sign(log(m$HR_tcga)) == sign(log(m$HR_mb)),
      effect_metric = "Cox HR per +1 SD")
    # keep the nominal p-values too, for the nominal-significance replication rate
    rows[[r]][, p_TCGA_nominal := m$p_tcga][, p_external_nominal := m$p_mb]
  }
} else logf("WARNING: METABRIC Cox results not found, skipping")

## ======================= 3. MET500 metastasis vs primary ====================
mt <- fread(file.path(BASE, "results/external_MET500_breast.csv"))
tde <- fread(file.path(BASE, "results/BRCA_DEX_genes.csv"))[, .(feature, logFC_tcga = logFC,
                                                                q_tcga = adj.P.Val)]
m3 <- merge(tde, mt[, .(feature = gene, delta = delta_pctrank_met_minus_prim,
                        p_ext = wilcox_p, q_ext = wilcox_q_BH, n_met = n_met[1],
                        n_prim = n_prim[1])], by = "feature")
r <- r + 1L
rows[[r]] <- data.table(
  cohort = "MET500", accession = "MET500 (breast cohort)",
  platform = "RNA-seq (poly-A / exome-capture), log2 expression",
  n_samples = m3$n_met[1],
  n_detail = sprintf("%d metastatic breast tumours vs %d TCGA primary", m3$n_met[1], m3$n_prim[1]),
  analysis = "metastasis_vs_primary",
  feature = m3$feature, feature_type = ntype[m3$feature],
  is_hub = m3$feature %in% hubset,
  effect_TCGA = m3$logFC_tcga, sig_TCGA = m3$q_tcga,
  effect_external = m3$delta, sig_external = m3$q_ext,
  same_direction = sign(m3$logFC_tcga) == sign(m3$delta),
  effect_metric = "delta within-sample percentile rank (met - primary)")

val <- rbindlist(rows, use.names = TRUE, fill = TRUE)
fwrite(val, file.path(BASE, "results/external_cohort_validation.csv"))
logf("WROTE results/external_cohort_validation.csv rows=", nrow(val))

## =========================== replication rates ==============================
summ <- list(); s <- 0L
for (co in unique(val$cohort)) for (an in unique(val[cohort == co]$analysis)) {
  d <- val[cohort == co & analysis == an]
  # "significant in TCGA" defined two ways
  strict <- d[sig_TCGA < 0.05]
  loose  <- if ("p_TCGA_nominal" %in% names(d) && any(!is.na(d$p_TCGA_nominal)))
              d[p_TCGA_nominal < 0.05] else d[sig_TCGA < 0.05]
  hubs   <- strict[is_hub == TRUE]
  s <- s + 1L
  summ[[s]] <- data.table(
    cohort = co, accession = d$accession[1], platform = d$platform[1],
    n_samples = d$n_samples[1], n_detail = d$n_detail[1], analysis = an,
    n_features_compared = nrow(d),
    n_TCGA_sig_q05 = nrow(strict),
    n_replicate_direction = sum(strict$same_direction),
    pct_replicate_direction = 100 * mean(strict$same_direction),
    n_replicate_direction_and_sig = sum(strict$same_direction & strict$sig_external < 0.05),
    pct_replicate_direction_and_sig = 100 * mean(strict$same_direction & strict$sig_external < 0.05),
    n_TCGA_sig_nominal = nrow(loose),
    pct_replicate_direction_nominalset = 100 * mean(loose$same_direction),
    n_replicate_direction_nominalset = sum(loose$same_direction),
    binom_p_nominalset = if (nrow(loose) > 0)
      stats::binom.test(sum(loose$same_direction), nrow(loose), 0.5)$p.value else NA_real_,
    # reverse direction: features significant in the EXTERNAL cohort, checked in TCGA
    n_external_sig_q05 = sum(d$sig_external < 0.05, na.rm = TRUE),
    pct_external_sig_same_direction = 100 * mean(d[sig_external < 0.05]$same_direction),
    binom_p_reverse = if (sum(d$sig_external < 0.05, na.rm = TRUE) > 0)
      stats::binom.test(sum(d[sig_external < 0.05]$same_direction),
                        sum(d$sig_external < 0.05, na.rm = TRUE), 0.5)$p.value else NA_real_,
    spearman_rho_effect_sizes = if (grepl("^Cox", an))
      suppressWarnings(stats::cor(log(d$effect_TCGA), log(d$effect_external),
                                  method = "spearman", use = "complete.obs"))
      else suppressWarnings(stats::cor(d$effect_TCGA, d$effect_external,
                                       method = "spearman", use = "complete.obs")),
    n_hubs_TCGA_sig = nrow(hubs),
    pct_hubs_replicate_direction = if (nrow(hubs)) 100 * mean(hubs$same_direction) else NA_real_,
    binom_p_vs_chance = if (nrow(strict) > 0)
      stats::binom.test(sum(strict$same_direction), nrow(strict), 0.5)$p.value else NA_real_)
}
summ <- rbindlist(summ)
fwrite(summ, file.path(BASE, "results/external_cohort_replication_summary.csv"))
logf("WROTE results/external_cohort_replication_summary.csv rows=", nrow(summ))
for (i in seq_len(nrow(summ))) {
  x <- summ[i]
  logf(sprintf("%-10s %-22s n=%s (%s)", x$cohort, x$analysis, x$n_samples, x$n_detail))
  logf(sprintf("   TCGA-significant features compared: %d ; same direction %d (%.1f%%) [binom p=%.3g] ; same direction AND q_ext<0.05 %d (%.1f%%)",
               x$n_TCGA_sig_q05, x$n_replicate_direction, x$pct_replicate_direction,
               x$binom_p_vs_chance, x$n_replicate_direction_and_sig,
               x$pct_replicate_direction_and_sig))
  logf(sprintf("   nominal-significant-in-TCGA set: %d features, %d same direction (%.1f%%) [binom p=%.3g]",
               x$n_TCGA_sig_nominal, x$n_replicate_direction_nominalset,
               x$pct_replicate_direction_nominalset, x$binom_p_nominalset))
  logf(sprintf("   REVERSE (external-significant, direction in TCGA): %d features, %.1f%% same direction [binom p=%.3g]",
               x$n_external_sig_q05, x$pct_external_sig_same_direction, x$binom_p_reverse))
  logf(sprintf("   Spearman rho of effect sizes across all compared features = %.3f",
               x$spearman_rho_effect_sizes))
  logf(sprintf("   hubs among them: %d ; %.1f%% same direction",
               x$n_hubs_TCGA_sig, x$pct_hubs_replicate_direction))
}
logf("DONE")
