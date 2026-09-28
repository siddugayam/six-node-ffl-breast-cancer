#!/usr/bin/env Rscript
# 09_anova_ffl_named_axes.R
# D) quantitative demonstration that the published TF-vs-miRNA logFC ANOVA is invalid
# E) FFL-level (3-node core) sign-coherence validation
# F) the manuscript's named regulatory axes, measured
# + sensitivity analysis: partial Spearman correlation adjusting for the top 5 PCs
#   (global tumour-composition / purity confounding)

REV <- "/path/to/revision"
con <- file(file.path(REV, "logs", "expr_validation_part2.log"), open = "wt")
say <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n"); cat(msg, "\n", file = con); flush(con)
}
set.seed(1234)
say("=== 09 START ===")

WS <- readRDS(file.path(REV, "data", "edge_corr_workspace.rds"))
E <- WS$E; NL <- WS$NL; sm <- WS$samples; gene_ok <- WS$gene_ok; mirna_ok <- WS$mirna_ok
N_PAIRED <- length(sm)
nodes <- read.delim(file.path(REV, "data", "canonical_nodes.tsv"), stringsAsFactors = FALSE)
edges <- read.delim(file.path(REV, "data", "canonical_edges.tsv"), stringsAsFactors = FALSE)
node_type <- setNames(nodes$type, nodes$name)
say("workspace: ", nrow(E), " edges, ", nrow(NL), " null pairs, n=", N_PAIRED, " paired tumours")

##############################################################################
## D) WHY THE PUBLISHED ANOVA WAS WRONG
##############################################################################
say("=== D) ANOVA CHECK ===")
dg <- read.csv(file.path(REV, "results", "BRCA_DEX_genes.csv"), stringsAsFactors = FALSE)
dm <- read.csv(file.path(REV, "results", "BRCA_DEX_mirnas.csv"), stringsAsFactors = FALSE)
say("DE tables: genes ", nrow(dg), " rows, miRNAs ", nrow(dm), " rows")

DE <- rbind(data.frame(feature = dg$feature, logFC = dg$logFC, AveExpr = dg$AveExpr,
                       adjP = dg$adj.P.Val, platform = "mRNA (log2 RSEM)", stringsAsFactors = FALSE),
            data.frame(feature = dm$feature, logFC = dm$logFC, AveExpr = dm$AveExpr,
                       adjP = dm$adj.P.Val, platform = "miRNA (log2 RPM)", stringsAsFactors = FALSE))
DE$class <- ifelse(!is.na(node_type[DE$feature]), node_type[DE$feature], "not_in_network")
DE$class[DE$platform == "miRNA (log2 RPM)" & DE$class == "not_in_network"] <- "not_in_network_miRNA"

deg_tot <- table(factor(c(edges$source, edges$target), levels = nodes$name))
deg_tot <- setNames(as.integer(deg_tot), names(deg_tot))
top20 <- unlist(lapply(c("TF", "Gene", "miRNA"), function(ty) {
  nn <- nodes$name[nodes$type == ty]; nn <- nn[nn %in% DE$feature]
  head(nn[order(-deg_tot[nn])], 20)
}))
DE$is_top20hub <- DE$feature %in% top20

summ <- function(d, lab) data.frame(
  set = lab, n = nrow(d),
  mean_logFC = mean(d$logFC), sd_logFC = sd(d$logFC), var_logFC = var(d$logFC),
  IQR_logFC = IQR(d$logFC), min_logFC = min(d$logFC), max_logFC = max(d$logFC),
  range_logFC = diff(range(d$logFC)),
  mean_AveExpr = mean(d$AveExpr), sd_AveExpr = sd(d$AveExpr),
  min_AveExpr = min(d$AveExpr), max_AveExpr = max(d$AveExpr),
  range_AveExpr = diff(range(d$AveExpr)), stringsAsFactors = FALSE)

rows <- list()
for (ty in c("TF", "Gene", "miRNA")) {
  d <- DE[DE$class == ty, ]; rows[[paste0("network_", ty)]] <- summ(d, paste0("network ", ty, " nodes"))
  d2 <- d[d$is_top20hub, ]; rows[[paste0("top20_", ty)]] <- summ(d2, paste0("top-20 degree hub ", ty))
}
rows[["all_mRNA"]]  <- summ(DE[DE$platform == "mRNA (log2 RSEM)", ], "ALL tested mRNA features")
rows[["all_miRNA"]] <- summ(DE[DE$platform == "miRNA (log2 RPM)", ], "ALL tested miRNA features")
ANOVATAB <- do.call(rbind, rows)
for (i in seq_len(nrow(ANOVATAB))) say(sprintf("%-26s n=%5d  logFC mean=%+.3f sd=%.3f range=%.2f | AveExpr mean=%.2f sd=%.2f range=%.2f",
  ANOVATAB$set[i], ANOVATAB$n[i], ANOVATAB$mean_logFC[i], ANOVATAB$sd_logFC[i],
  ANOVATAB$range_logFC[i], ANOVATAB$mean_AveExpr[i], ANOVATAB$sd_AveExpr[i], ANOVATAB$range_AveExpr[i]))

## the authors' test, reproduced on the top-20 hubs of each class
top <- DE[DE$is_top20hub, ]
top$class <- factor(top$class, levels = c("TF", "Gene", "miRNA"))
aovfit <- aov(logFC ~ class, data = top); asum <- summary(aovfit)[[1]]
say("REPRODUCED ANOVA on top-20 hubs per class: F(", asum$Df[1], ",", asum$Df[2], ") = ",
    signif(asum$`F value`[1], 4), " , p = ", signif(asum$`Pr(>F)`[1], 4))
bt <- bartlett.test(logFC ~ class, data = top)
say("Bartlett test of equal variance across the three classes: K2 = ", signif(unname(bt$statistic), 4),
    " df=", unname(bt$parameter), " p = ", signif(bt$p.value, 4),
    "  -> homoscedasticity assumption ", ifelse(bt$p.value < 0.05, "VIOLATED", "not rejected"))
## Levene (Brown-Forsythe, median-centred) without car
lev <- function(y, g) {
  g <- factor(g); z <- abs(y - ave(y, g, FUN = median))
  a <- summary(aov(z ~ g))[[1]]; c(F = a$`F value`[1], p = a$`Pr(>F)`[1], df1 = a$Df[1], df2 = a$Df[2])
}
lv <- lev(top$logFC, top$class)
say("Brown-Forsythe/Levene test: F(", lv["df1"], ",", lv["df2"], ") = ", signif(lv["F"], 4),
    " p = ", signif(lv["p"], 4))
sw <- by(top$logFC, top$class, function(x) shapiro.test(x)$p.value)
say("Shapiro-Wilk normality p by class: TF=", signif(sw[["TF"]], 3), " Gene=", signif(sw[["Gene"]], 3),
    " miRNA=", signif(sw[["miRNA"]], 3))
kw <- kruskal.test(logFC ~ class, data = top)
say("Kruskal-Wallis (non-parametric equivalent): chi2 = ", signif(unname(kw$statistic), 4),
    " df=", unname(kw$parameter), " p = ", signif(kw$p.value, 4))

## THE DECISIVE DEMONSTRATION: the same ANOVA on RANDOM features of each class
gpool <- DE[DE$platform == "mRNA (log2 RSEM)", ]
mpool <- DE[DE$platform == "miRNA (log2 RPM)", ]
B <- 2000; pvals <- numeric(B); dvals <- numeric(B)
for (b in seq_len(B)) {
  a1 <- gpool[sample.int(nrow(gpool), 20), ]; a2 <- gpool[sample.int(nrow(gpool), 20), ]
  a3 <- mpool[sample.int(nrow(mpool), 20), ]
  dd <- data.frame(logFC = c(a1$logFC, a2$logFC, a3$logFC),
                   class = rep(c("TF", "Gene", "miRNA"), each = 20))
  s <- summary(aov(logFC ~ class, data = dd))[[1]]
  pvals[b] <- s$`Pr(>F)`[1]
  dvals[b] <- sd(c(a1$logFC, a2$logFC)) / sd(a3$logFC)
}
say("RANDOM-FEATURE CONTROL (", B, " draws of 20 random mRNAs / 20 random mRNAs / 20 random miRNAs):")
say("  fraction of random draws with ANOVA p < 0.05 = ", signif(mean(pvals < 0.05), 4),
    " (median p = ", signif(median(pvals), 4), ")")
say("  median ratio SD(logFC mRNA) / SD(logFC miRNA) in random draws = ", signif(median(dvals), 4))
## variance-ratio F test, whole platforms
vt <- var.test(gpool$logFC, mpool$logFC)
say("VARIANCE RATIO whole platforms: var(logFC mRNA)/var(logFC miRNA) = ",
    signif(unname(vt$estimate), 4), " , F test p = ", signif(vt$p.value, 4),
    " ; SD mRNA = ", signif(sd(gpool$logFC), 4), " SD miRNA = ", signif(sd(mpool$logFC), 4))

extra <- data.frame(
  metric = c("reproduced_ANOVA_F", "reproduced_ANOVA_p", "reproduced_ANOVA_df1", "reproduced_ANOVA_df2",
             "Bartlett_K2", "Bartlett_p", "Levene_BF_F", "Levene_BF_p",
             "ShapiroWilk_p_TF", "ShapiroWilk_p_Gene", "ShapiroWilk_p_miRNA",
             "KruskalWallis_chi2", "KruskalWallis_p",
             "random_draw_frac_ANOVA_p_lt_0.05", "random_draw_median_ANOVA_p",
             "random_draw_median_SDratio_mRNA_over_miRNA",
             "platform_var_ratio_mRNA_over_miRNA", "platform_var_ratio_F_p",
             "SD_logFC_all_mRNA", "SD_logFC_all_miRNA",
             "mean_AveExpr_all_mRNA", "mean_AveExpr_all_miRNA"),
  value = c(asum$`F value`[1], asum$`Pr(>F)`[1], asum$Df[1], asum$Df[2],
            unname(bt$statistic), bt$p.value, unname(lv["F"]), unname(lv["p"]),
            sw[["TF"]], sw[["Gene"]], sw[["miRNA"]],
            unname(kw$statistic), kw$p.value,
            mean(pvals < 0.05), median(pvals), median(dvals),
            unname(vt$estimate), vt$p.value,
            sd(gpool$logFC), sd(mpool$logFC), mean(gpool$AveExpr), mean(mpool$AveExpr)),
  stringsAsFactors = FALSE)
write.csv(ANOVATAB, file.path(REV, "results", "anova_named_axes_class_stats.csv"), row.names = FALSE)
write.csv(extra,    file.path(REV, "results", "anova_named_axes_tests.csv"), row.names = FALSE)
say("WROTE results/anova_named_axes_class_stats.csv rows=", nrow(ANOVATAB),
    " and results/anova_named_axes_tests.csv rows=", nrow(extra))

##############################################################################
## E) FFL-LEVEL COHERENCE
##############################################################################
say("=== E) FFL 3-NODE CORE COHERENCE ===")
cf <- file.path(REV, "results", "ffl_3node_cores.csv")
if (!file.exists(cf)) cf <- file.path(REV, "results", "ffl_3node_cores_crosscheck.csv")
say("using 3-node core file: ", basename(cf))
CO <- read.csv(cf, stringsAsFactors = FALSE)
say("cores read: ", nrow(CO), " ; classes: ", paste(names(table(CO$class)), table(CO$class), sep = "=", collapse = " "))

key <- function(s, t) paste(s, t, sep = "\r")
ei <- setNames(seq_len(nrow(E)), key(E$source, E$target))
gi <- function(s, t) unname(ei[key(s, t)])
CO$i_RM <- gi(CO$regulator, CO$intermediate)
CO$i_MT <- gi(CO$intermediate, CO$target)
CO$i_RT <- gi(CO$regulator, CO$target)
say("arms resolved to canonical edges: RM ", sum(!is.na(CO$i_RM)), " MT ", sum(!is.na(CO$i_MT)),
    " RT ", sum(!is.na(CO$i_RT)), " of ", nrow(CO))

grab <- function(idx, col) ifelse(is.na(idx), NA, E[[col]][idx])
for (arm in c("RM", "MT", "RT")) {
  ii <- CO[[paste0("i_", arm)]]
  CO[[paste0("pred_", arm)]] <- grab(ii, "predicted_sign")
  CO[[paste0("rho_", arm)]]  <- grab(ii, "rho")
  CO[[paste0("p_", arm)]]    <- grab(ii, "p")
  CO[[paste0("fdr_", arm)]]  <- grab(ii, "fdr")
  CO[[paste0("conc_", arm)]] <- grab(ii, "concordant")
  CO[[paste0("etype_", arm)]] <- grab(ii, "edge_type")
}
CO$all3_measurable <- !is.na(CO$rho_RM) & !is.na(CO$rho_MT) & !is.na(CO$rho_RT)
CO$all3_signed <- !is.na(CO$pred_RM) & !is.na(CO$pred_MT) & !is.na(CO$pred_RT)
CO$all3_testable <- CO$all3_measurable & CO$all3_signed
CO$n_arms_concordant <- rowSums(cbind(CO$conc_RM, CO$conc_MT, CO$conc_RT), na.rm = TRUE)
CO$all3_concordant <- ifelse(CO$all3_testable, CO$conc_RM & CO$conc_MT & CO$conc_RT, NA)

## network-predicted vs observed coherence class (sign logic of the loop)
CO$predicted_coherence <- ifelse(CO$all3_signed,
  ifelse(CO$pred_RM * CO$pred_MT == CO$pred_RT, "coherent", "incoherent"), NA)
CO$observed_coherence <- ifelse(CO$all3_measurable,
  ifelse(sign(CO$rho_RM) * sign(CO$rho_MT) == sign(CO$rho_RT), "coherent", "incoherent"), NA)
CO$coherence_match <- ifelse(!is.na(CO$predicted_coherence) & !is.na(CO$observed_coherence),
                             CO$predicted_coherence == CO$observed_coherence, NA)

T3 <- CO[CO$all3_testable, ]
say("cores with all three arms measurable AND sign-annotated: ", nrow(T3), " of ", nrow(CO))
if (nrow(T3) > 0) {
  p_RM <- mean(T3$conc_RM); p_MT <- mean(T3$conc_MT); p_RT <- mean(T3$conc_RT)
  obs <- mean(T3$all3_concordant)
  expect_ind <- p_RM * p_MT * p_RT
  say(sprintf("per-arm concordance in these cores: R->M %.4f  M->T %.4f  R->T %.4f", p_RM, p_MT, p_RT))
  say(sprintf("all-three-arms concordant: observed %.4f (%d/%d) ; expected under independence %.4f",
              obs, sum(T3$all3_concordant), nrow(T3), expect_ind))
  ## permutation: independently permute each arm's concordance across cores (preserves marginals)
  B2 <- 10000; ge <- 0L; nsim <- numeric(B2)
  a <- T3$conc_RM; b <- T3$conc_MT; cc <- T3$conc_RT; K <- nrow(T3)
  for (s in seq_len(B2)) {
    v <- mean(a[sample.int(K)] & b[sample.int(K)] & cc[sample.int(K)])
    nsim[s] <- v; if (v >= obs) ge <- ge + 1L
  }
  perm_p <- (ge + 1) / (B2 + 1)
  say(sprintf("permutation (10000x, arm labels shuffled across cores): mean null %.4f  sd %.4f  p = %.4g",
              mean(nsim), sd(nsim), perm_p))
  bino <- binom.test(sum(T3$all3_concordant), nrow(T3), p = expect_ind, alternative = "greater")
  say("binomial vs independence expectation: p = ", signif(bino$p.value, 4))
  ## coherence-class agreement
  cm <- T3[!is.na(T3$coherence_match), ]
  say("network-predicted vs expression-observed coherence class: agreement ",
      sum(cm$coherence_match), "/", nrow(cm), " = ", signif(mean(cm$coherence_match), 4))
  say("  predicted coherent ", sum(cm$predicted_coherence == "coherent"),
      " ; observed coherent ", sum(cm$observed_coherence == "coherent"))
  bino2 <- binom.test(sum(cm$coherence_match), nrow(cm), p = 0.5, alternative = "two.sided")
  say("  binomial vs 0.5: p = ", signif(bino2$p.value, 4))
  ffl_summary <- data.frame(
    metric = c("n_cores_total", "n_cores_all3_measurable", "n_cores_all3_testable",
               "arm_conc_R_to_M", "arm_conc_M_to_T", "arm_conc_R_to_T",
               "observed_all3_concordant_n", "observed_all3_concordant_rate",
               "expected_all3_under_independence", "permutation_mean_null",
               "permutation_sd_null", "permutation_p", "binomial_p_vs_independence",
               "coherence_class_agreement_n", "coherence_class_agreement_total",
               "coherence_class_agreement_rate", "coherence_binom_p_vs_0.5"),
    value = c(nrow(CO), sum(CO$all3_measurable), nrow(T3), p_RM, p_MT, p_RT,
              sum(T3$all3_concordant), obs, expect_ind, mean(nsim), sd(nsim), perm_p,
              bino$p.value, sum(cm$coherence_match), nrow(cm), mean(cm$coherence_match),
              bino2$p.value), stringsAsFactors = FALSE)
  write.csv(ffl_summary, file.path(REV, "results", "ffl_coherence_summary.csv"), row.names = FALSE)
  say("WROTE results/ffl_coherence_summary.csv rows=", nrow(ffl_summary))
}
outcols <- c("regulator", "intermediate", "target", "class", "type_R", "type_M", "type_T",
             "etype_RM", "etype_MT", "etype_RT",
             "pred_RM", "pred_MT", "pred_RT", "rho_RM", "rho_MT", "rho_RT",
             "p_RM", "p_MT", "p_RT", "fdr_RM", "fdr_MT", "fdr_RT",
             "conc_RM", "conc_MT", "conc_RT", "n_arms_concordant",
             "all3_measurable", "all3_signed", "all3_testable", "all3_concordant",
             "predicted_coherence", "observed_coherence", "coherence_match")
outcols <- intersect(outcols, colnames(CO))
write.csv(CO[, outcols], file.path(REV, "results", "ffl_coherence_validation.csv"), row.names = FALSE)
say("WROTE results/ffl_coherence_validation.csv rows=", nrow(CO))

##############################################################################
## F) NAMED AXES FROM THE DISCUSSION
##############################################################################
say("=== F) NAMED AXES ===")
G <- readRDS(file.path(REV, "data", "brca_gene_expr.rds"))[, sm, drop = FALSE]
M <- readRDS(file.path(REV, "data", "brca_mirna_expr_canonical.rds"))[, sm, drop = FALSE]
sdG <- apply(G, 1, sd); sdM <- apply(M, 1, sd)
zrow <- function(X) { mu <- rowMeans(X); s <- sqrt(rowSums((X - mu)^2)); (X - mu) / s }
allfeat <- rbind(zrow(t(apply(G[sdG > 0, , drop = FALSE], 1, rank))),
                 zrow(t(apply(M[sdM > 0, , drop = FALSE], 1, rank))))
say("full z-rank matrix for named axes: ", nrow(allfeat), " features")
named <- data.frame(rbind(
  c("NFKB1", "COL1A1"), c("SP1", "COL1A1"), c("RELA", "COL1A1"), c("ETS1", "COL1A1"),
  c("COL1A1", "COL3A1"),
  c("hsa-miR-29a", "COL1A1"), c("hsa-miR-29b", "COL1A1"), c("hsa-miR-29c", "COL1A1"),
  c("hsa-miR-29a", "COL3A1"), c("hsa-miR-29b", "COL3A1"), c("hsa-miR-29c", "COL3A1"),
  c("hsa-let-7b", "COL3A1"), c("hsa-let-7e", "COL3A1"),
  c("hsa-miR-130a", "VEGFA"), c("hsa-miR-130a", "MMP2"),
  c("hsa-miR-124", "STAT3"), c("hsa-miR-101", "EZH2"), c("hsa-let-7b", "HK2")),
  stringsAsFactors = FALSE)
colnames(named) <- c("source", "target")
eset <- setNames(E$edge_type, key(E$source, E$target))
esign <- setNames(E$predicted_sign, key(E$source, E$target))
eann <- setNames(E$annotation, key(E$source, E$target))
etier <- setNames(E$evidence_tier, key(E$source, E$target))
kk <- key(named$source, named$target)
named$in_canonical_network <- kk %in% names(eset)
named$edge_type <- unname(eset[kk])
named$predicted_sign <- unname(esign[kk])
named$db_annotation <- unname(eann[kk])
named$evidence_tier <- unname(etier[kk])
named$rho <- NA_real_; named$p <- NA_real_; named$n_samples <- NA_integer_
for (i in seq_len(nrow(named))) {
  a <- named$source[i]; b <- named$target[i]
  if (a %in% rownames(allfeat) && b %in% rownames(allfeat)) {
    r <- sum(allfeat[a, ] * allfeat[b, ]); r <- max(min(r, 1), -1)
    named$rho[i] <- r
    named$p[i] <- 2 * pt(-abs(r * sqrt((N_PAIRED - 2) / max(1 - r^2, 1e-300))), df = N_PAIRED - 2)
    named$n_samples[i] <- N_PAIRED
  }
}
named$fdr_within_named <- p.adjust(named$p, method = "BH")
named$direction_supports_claim <- ifelse(is.na(named$rho), NA,
  ifelse(named$source %in% rownames(M), named$rho < 0, named$rho > 0))
named$concordant_with_predicted_sign <- ifelse(is.na(named$rho) | is.na(named$predicted_sign), NA,
                                               sign(named$rho) == named$predicted_sign)
for (i in seq_len(nrow(named))) say(sprintf("%-14s -> %-8s rho=%+.4f p=%.3g fdr=%.3g inNet=%s type=%s ann=%s tier=%s claimDir=%s",
  named$source[i], named$target[i], named$rho[i], named$p[i], named$fdr_within_named[i],
  named$in_canonical_network[i], ifelse(is.na(named$edge_type[i]), "-", named$edge_type[i]),
  ifelse(is.na(named$db_annotation[i]), "-", named$db_annotation[i]),
  ifelse(is.na(named$evidence_tier[i]), "-", named$evidence_tier[i]),
  named$direction_supports_claim[i]))
write.csv(named, file.path(REV, "results", "named_axis_correlations.csv"), row.names = FALSE)
say("WROTE results/named_axis_correlations.csv rows=", nrow(named))

##############################################################################
## SENSITIVITY: partial Spearman adjusting for the top 5 PCs (composition/purity)
##############################################################################
say("=== SENSITIVITY: PC-adjusted partial Spearman ===")
Zn <- allfeat[c(intersect(gene_ok, rownames(allfeat)), intersect(mirna_ok, rownames(allfeat))), ]
sv <- svd(Zn, nu = 0, nv = 10)
K <- 5L
V <- sv$v[, seq_len(K), drop = FALSE]
varexp <- sv$d^2 / sum(sv$d^2)
say("top-5 PC variance explained of the network-feature rank matrix: ",
    paste(signif(varexp[1:5], 3), collapse = ", "), " (cumulative ", signif(sum(varexp[1:5]), 3), ")")
Rres <- Zn - (Zn %*% V) %*% t(V)
Rres <- Rres / sqrt(rowSums(Rres^2))
dfa <- N_PAIRED - 2 - K
rho_adj <- function(a, b) {
  ok <- a %in% rownames(Rres) & b %in% rownames(Rres)
  out <- rep(NA_real_, length(a))
  if (any(ok)) out[ok] <- rowSums(Rres[a[ok], , drop = FALSE] * Rres[b[ok], , drop = FALSE])
  pmax(pmin(out, 1), -1)
}
E$rho_pcadj <- rho_adj(E$source, E$target)
E$p_pcadj <- 2 * pt(-abs(E$rho_pcadj * sqrt(dfa / pmax(1 - E$rho_pcadj^2, 1e-300))), df = dfa)
E$conc_pcadj <- ifelse(is.na(E$rho_pcadj) | is.na(E$predicted_sign), NA,
                       sign(E$rho_pcadj) == E$predicted_sign)
NL$rho_pcadj <- rho_adj(NL$rnd_source, NL$rnd_target)
NL$conc_pcadj <- ifelse(is.na(NL$predicted_sign), NA, sign(NL$rho_pcadj) == NL$predicted_sign)

ann <- E$annotation
strata <- list(
  "ALL_sign_annotated"    = which(!is.na(E$rho_pcadj) & !is.na(E$predicted_sign)),
  "miRNA_target_ALL"      = which(!is.na(E$rho_pcadj) & E$edge_type == "miRNA_target"),
  "miRNA_target_strong"   = which(!is.na(E$rho_pcadj) & E$edge_type == "miRNA_target" & E$evidence_tier == "strong"),
  "miRNA_target_weak"     = which(!is.na(E$rho_pcadj) & E$edge_type == "miRNA_target" & E$evidence_tier == "weak"),
  "miRNA_target_predonly" = which(!is.na(E$rho_pcadj) & E$edge_type == "miRNA_target" & E$evidence_tier == "predicted_only"),
  "TF_target_Activation"  = which(!is.na(E$rho_pcadj) & E$edge_type == "TF_target" & ann == "Activation"),
  "TF_target_Repression"  = which(!is.na(E$rho_pcadj) & E$edge_type == "TF_target" & ann == "Repression"),
  "TF_miRNA_Activation"   = which(!is.na(E$rho_pcadj) & E$edge_type == "TF_miRNA" & ann == "Activation"),
  "TF_miRNA_Repression"   = which(!is.na(E$rho_pcadj) & E$edge_type == "TF_miRNA" & ann == "Repression"))
sres <- list()
for (nm in names(strata)) {
  rows <- strata[[nm]]; rows <- rows[!is.na(E$conc_pcadj[rows])]
  if (length(rows) < 2) next
  cc <- as.numeric(E$conc_pcadj[rows])
  nn <- NL$conc_pcadj[NL$edge_row %in% rows]; nn <- nn[!is.na(nn)]
  p0 <- mean(nn); k <- sum(cc); n <- length(cc)
  bt <- binom.test(k, n, p = p0, alternative = "greater")
  sres[[nm]] <- data.frame(stratum = nm, n_edges = n, concordance_rate_pcadj = k / n,
    null_concordance_rate_pcadj = p0, excess = k / n - p0, binom_p = bt$p.value,
    mean_rho_pcadj_real = mean(E$rho_pcadj[rows]),
    mean_rho_pcadj_null = mean(NL$rho_pcadj[NL$edge_row %in% rows], na.rm = TRUE),
    stringsAsFactors = FALSE)
  say(sprintf("PCadj %-22s n=%5d conc=%.4f null=%.4f excess=%+.4f binomP=%.4g meanRho real=%+.4f null=%+.4f",
      nm, n, k / n, p0, k / n - p0, bt$p.value, mean(E$rho_pcadj[rows]),
      mean(NL$rho_pcadj[NL$edge_row %in% rows], na.rm = TRUE)))
}
SENS <- do.call(rbind, sres)
write.csv(SENS, file.path(REV, "results", "sign_concordance_pcadjusted.csv"), row.names = FALSE)
say("WROTE results/sign_concordance_pcadjusted.csv rows=", nrow(SENS))

## append the PC-adjusted columns to the main edge table
main <- read.csv(file.path(REV, "results", "edge_correlation.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(main) == nrow(E), all(main$source == E$source), all(main$target == E$target))
main$rho_pcadj <- E$rho_pcadj; main$p_pcadj <- E$p_pcadj; main$concordant_pcadj <- E$conc_pcadj
write.csv(main, file.path(REV, "results", "edge_correlation.csv"), row.names = FALSE)
say("edge_correlation.csv updated with PC-adjusted columns; rows=", nrow(main))
say("=== 09 DONE ===")
close(con)
