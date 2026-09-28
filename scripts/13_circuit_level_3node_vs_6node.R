#!/usr/bin/env Rscript
# 13_circuit_level_3node_vs_6node.R
#
# The aggregate-score test (script 12) asks whether the *union* of higher-order
# nodes carries information. This script asks the sharper question, circuit by
# circuit: score EVERY individual 3-node FFL and a large random sample of
# individual 6-node FFLs, fit Cox for each, and compare the two distributions of
# prognostic strength. Also fits, for every 6-node circuit, a nested model
# testing whether its 3 extension members add anything over its own 3-node core.
# Output: results/circuit_level_3node_vs_6node.csv (per-circuit)
#         results/circuit_level_summary.csv       (the verdict)

suppressPackageStartupMessages({ library(data.table); library(survival) })
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs", "13_circuit_level.log")
logf <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                        cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); logf("START")
set.seed(1234)

sets  <- readRDS(file.path(BASE, "data/ffl_module_sets.rds"))
ffl3  <- as.data.table(sets$ffl3)
nodes <- fread(file.path(BASE, "data/canonical_nodes.tsv")); ntype <- setNames(nodes$type, nodes$name)

gex <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
mex <- readRDS(file.path(BASE, "data/brca_mirna_expr_canonical.rds"))
pheno <- readRDS(file.path(BASE, "data/brca_pheno.rds"))
surv  <- as.data.table(readRDS(file.path(BASE, "data/brca_survival_clean.rds")))
tumour <- pheno$sample[pheno$sample_type == "Primary Tumor"]
common <- Reduce(intersect, list(colnames(gex), colnames(mex), tumour, surv$sample_id))
X <- rbind(gex[, common, drop = FALSE], mex[, common, drop = FALSE])
v <- apply(X, 1, function(z) stats::sd(z, na.rm = TRUE))
X <- X[is.finite(v) & v > 0, , drop = FALSE]
Z <- t(scale(t(X)))                    # per-feature z-score across samples
logf("scoring matrix ", paste(dim(Z), collapse = " x "), " ; samples=", length(common))

cl <- surv[sample_id %in% common]; setkey(cl, sample_id); cl <- cl[common]
stopifnot(identical(cl$sample_id, common))

## ------------------------------------------------- extension edge lists -----
tf_all <- fread(file.path(BASE, "data/layer_TF_target.tsv"))
gg <- fread(file.path(BASE, "data/layer_gene_gene.tsv"))
mm <- fread(file.path(BASE, "data/layer_miRNA_miRNA.tsv"))
tf_tf <- unique(tf_all[ntype[source] == "TF" & ntype[target] == "TF" & source != target,
                       .(source, target)])
tftfu <- unique(rbind(tf_tf[, .(a = source, b = target)], tf_tf[, .(a = target, b = source)]))
ggu <- unique(rbind(gg[, .(a = source, b = target)], gg[, .(a = target, b = source)]))
ggu <- ggu[ntype[a] == "Gene" & ntype[b] == "Gene" & a != b]
mmu <- unique(rbind(mm[, .(a = miRNA_1, b = miRNA_2)], mm[, .(a = miRNA_2, b = miRNA_1)]))
mmu <- mmu[a != b]

have <- rownames(Z)
ffl3 <- ffl3[TF %in% have & miR %in% have & Gene %in% have]
logf("3-node FFLs with full expression coverage: ", nrow(ffl3), " / ", nrow(sets$ffl3))

## ------------------------------------- build the 6-node circuit population ---
# 6-node = (TF1, miR1, G1) core + G2 (gene-gene w/ G1) + miR2 (miR-miR w/ miR1) + TF2 (TF-TF w/ TF1)
g_map <- split(ggu$b,   ggu$a)
m_map <- split(mmu$b,   mmu$a)
t_map <- split(tftfu$b, tftfu$a)
cnt <- ffl3[, .(nG2 = length(intersect(g_map[[Gene]], have)),
                nM2 = length(intersect(m_map[[miR]],  have)),
                nT2 = length(intersect(t_map[[TF]],   have))), by = seq_len(nrow(ffl3))]
ffl3[, `:=`(nG2 = cnt$nG2, nM2 = cnt$nM2, nT2 = cnt$nT2)]
ffl3[, n6 := as.numeric(nG2) * as.numeric(nM2) * as.numeric(nT2)]
TOT6 <- sum(ffl3$n6)
logf("6-node circuits with full expression coverage (exact count) = ",
     format(TOT6, big.mark = ","), " over ", sum(ffl3$n6 > 0), " extendable 3-node cores")

N6_SAMPLE <- 2000L
extendable <- ffl3[n6 > 0]
# sample cores proportional to how many 6-node circuits they generate (uniform over circuits)
idx <- sample.int(nrow(extendable), size = N6_SAMPLE, replace = TRUE,
                  prob = extendable$n6 / sum(extendable$n6))
six <- rbindlist(lapply(idx, function(i) {
  row <- extendable[i]
  data.table(TF = row$TF, miR = row$miR, Gene = row$Gene,
             G2  = sample(intersect(g_map[[row$Gene]], have), 1),
             M2  = sample(intersect(m_map[[row$miR]],  have), 1),
             T2  = sample(intersect(t_map[[row$TF]],   have), 1))
}))
six <- unique(six)
logf("distinct 6-node circuits sampled uniformly at random: ", nrow(six))

## ---------------------------------------------------------------- Cox util --
cox_score <- function(members, tt, ee, Zm) {
  sc <- colMeans(Zm[members, , drop = FALSE])
  scz <- as.numeric(scale(sc))
  f <- try(coxph(Surv(tt, ee) ~ scz), silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  s <- summary(f)
  list(HR = s$coefficients[1, "exp(coef)"], p = s$coefficients[1, "Pr(>|z|)"],
       z = s$coefficients[1, "z"], C = s$concordance[1], sc = scz)
}

out <- list(); k <- 0L
for (ep in c("OS", "PFI")) {
  tt <- cl[[paste0(ep, ".time")]]; ee <- cl[[ep]]
  ok <- is.finite(tt) & tt > 0 & is.finite(ee)
  t2 <- tt[ok]; e2 <- ee[ok]; Zs <- Z[, ok, drop = FALSE]
  logf(sprintf("endpoint %s: n=%d events=%d", ep, sum(ok), sum(e2 == 1)))

  ## --- every 3-node circuit ---
  for (i in seq_len(nrow(ffl3))) {
    mem <- c(ffl3$TF[i], ffl3$miR[i], ffl3$Gene[i])
    r <- cox_score(mem, t2, e2, Zs); if (is.null(r)) next
    k <- k + 1L
    out[[k]] <- data.table(endpoint = ep, circuit_class = "3-node", circuit_size = 3L,
                           members = paste(mem, collapse = ";"),
                           core = paste(mem, collapse = ";"),
                           HR = r$HR, p = r$p, absz = abs(r$z), C_index = r$C,
                           LRT_p_extension_adds = NA_real_, deltaC_vs_core = NA_real_)
  }

  ## --- sampled 6-node circuits, each compared with its own 3-node core ---
  for (i in seq_len(nrow(six))) {
    core <- c(six$TF[i], six$miR[i], six$Gene[i])
    ext  <- c(six$T2[i], six$M2[i], six$G2[i])
    mem  <- unique(c(core, ext))
    if (length(mem) < 6) next
    r6 <- cox_score(mem,  t2, e2, Zs); if (is.null(r6)) next
    r3 <- cox_score(core, t2, e2, Zs); if (is.null(r3)) next
    rx <- cox_score(ext,  t2, e2, Zs); if (is.null(rx)) next
    m3 <- coxph(Surv(t2, e2) ~ r3$sc)
    mb <- coxph(Surv(t2, e2) ~ r3$sc + rx$sc)
    a  <- anova(m3, mb)
    k <- k + 1L
    out[[k]] <- data.table(endpoint = ep, circuit_class = "6-node", circuit_size = 6L,
                           members = paste(mem, collapse = ";"),
                           core = paste(core, collapse = ";"),
                           HR = r6$HR, p = r6$p, absz = abs(r6$z), C_index = r6$C,
                           LRT_p_extension_adds = a[2, "Pr(>|Chi|)"],
                           deltaC_vs_core = r6$C - r3$C)
  }
}
res <- rbindlist(out)
res[, q_BH := p.adjust(p, method = "BH"), by = .(endpoint, circuit_class)]
fwrite(res, file.path(BASE, "results/circuit_level_3node_vs_6node.csv"))
logf("WROTE results/circuit_level_3node_vs_6node.csv rows=", nrow(res))

## ------------------------------------------------------------------ verdict --
summ <- list(); s <- 0L
for (ep in unique(res$endpoint)) {
  a <- res[endpoint == ep & circuit_class == "3-node"]
  b <- res[endpoint == ep & circuit_class == "6-node"]
  wt <- stats::wilcox.test(b$absz, a$absz)
  wc <- stats::wilcox.test(b$C_index, a$C_index)
  ext_fdr <- p.adjust(b$LRT_p_extension_adds, method = "BH")
  s <- s + 1L
  summ[[s]] <- data.table(
    endpoint = ep,
    n_3node = nrow(a), n_6node = nrow(b),
    median_absZ_3node = stats::median(a$absz), median_absZ_6node = stats::median(b$absz),
    wilcox_p_absZ = wt$p.value,
    median_C_3node = stats::median(a$C_index), median_C_6node = stats::median(b$C_index),
    wilcox_p_C = wc$p.value,
    n_3node_q05 = sum(a$q_BH < 0.05), pct_3node_q05 = 100 * mean(a$q_BH < 0.05),
    n_6node_q05 = sum(b$q_BH < 0.05), pct_6node_q05 = 100 * mean(b$q_BH < 0.05),
    best_C_3node = max(a$C_index), best_C_6node = max(b$C_index),
    n_6node_ext_adds_nominal = sum(b$LRT_p_extension_adds < 0.05, na.rm = TRUE),
    pct_6node_ext_adds_nominal = 100 * mean(b$LRT_p_extension_adds < 0.05, na.rm = TRUE),
    n_6node_ext_adds_FDR05 = sum(ext_fdr < 0.05, na.rm = TRUE),
    pct_6node_ext_adds_FDR05 = 100 * mean(ext_fdr < 0.05, na.rm = TRUE),
    median_deltaC_6node_minus_core = stats::median(b$deltaC_vs_core),
    pct_6node_better_than_own_core = 100 * mean(b$deltaC_vs_core > 0))
}
summ <- rbindlist(summ)
fwrite(summ, file.path(BASE, "results/circuit_level_summary.csv"))
logf("WROTE results/circuit_level_summary.csv")
for (i in seq_len(nrow(summ))) {
  x <- summ[i]
  logf(sprintf("== %s ==", x$endpoint))
  logf(sprintf("  n circuits: 3-node=%d  6-node=%d", x$n_3node, x$n_6node))
  logf(sprintf("  median |z|: 3-node=%.3f  6-node=%.3f  (Wilcoxon p=%.3g)",
               x$median_absZ_3node, x$median_absZ_6node, x$wilcox_p_absZ))
  logf(sprintf("  median C  : 3-node=%.4f 6-node=%.4f (Wilcoxon p=%.3g); best C 3=%.4f 6=%.4f",
               x$median_C_3node, x$median_C_6node, x$wilcox_p_C, x$best_C_3node, x$best_C_6node))
  logf(sprintf("  circuits passing BH q<0.05 within class: 3-node %d/%d (%.1f%%)  6-node %d/%d (%.1f%%)",
               x$n_3node_q05, x$n_3node, x$pct_3node_q05,
               x$n_6node_q05, x$n_6node, x$pct_6node_q05))
  logf(sprintf("  6-node extension adds over its OWN 3-node core: nominal p<0.05 %d (%.1f%%), BH q<0.05 %d (%.1f%%)",
               x$n_6node_ext_adds_nominal, x$pct_6node_ext_adds_nominal,
               x$n_6node_ext_adds_FDR05, x$pct_6node_ext_adds_FDR05))
  logf(sprintf("  median deltaC (6-node minus its own core) = %+.4f ; %.1f%% of 6-node circuits beat their core",
               x$median_deltaC_6node_minus_core, x$pct_6node_better_than_own_core))
}
logf("DONE")
