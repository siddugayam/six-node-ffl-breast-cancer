#!/usr/bin/env Rscript
# 10_directional_evidence_and_anova_control.R
# Supplementary to 08/09:
#  (i)  per-stratum counts of edges whose correlation is SIGNIFICANT and in the predicted
#       direction vs significant and in the OPPOSITE direction, real edges vs matched nulls
#       -> results/directional_evidence_counts.csv  and the expression-supported sub-network
#  (ii) a sharper, two-group control for the published ANOVA: random mRNAs vs random miRNAs
#       (no hub selection at all) tested exactly as the authors tested their hubs.

REV <- "/path/to/revision"
con <- file(file.path(REV, "logs", "expr_validation_part3.log"), open = "wt")
say <- function(...) { msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n"); cat(msg, "\n", file = con); flush(con) }
set.seed(1234)
say("=== 10 START ===")

E  <- read.csv(file.path(REV, "results", "edge_correlation.csv"), stringsAsFactors = FALSE)
NL <- read.csv(file.path(REV, "results", "edge_correlation_null.csv"), stringsAsFactors = FALSE)
say("edge_correlation.csv rows=", nrow(E), " ; edge_correlation_null.csv rows=", nrow(NL))

## null pairs need their own FDR to be comparable
NL$fdr <- p.adjust(NL$p, method = "BH")

ann <- E$annotation
strata <- list(
  "ALL_sign_annotated"    = which(!is.na(E$rho) & !is.na(E$predicted_sign)),
  "miRNA_target_ALL"      = which(!is.na(E$rho) & E$edge_type == "miRNA_target"),
  "miRNA_target_strong"   = which(!is.na(E$rho) & E$edge_type == "miRNA_target" & E$evidence_tier == "strong"),
  "miRNA_target_weak"     = which(!is.na(E$rho) & E$edge_type == "miRNA_target" & E$evidence_tier == "weak"),
  "miRNA_target_predonly" = which(!is.na(E$rho) & E$edge_type == "miRNA_target" & E$evidence_tier == "predicted_only"),
  "TF_target_Activation"  = which(!is.na(E$rho) & E$edge_type == "TF_target" & ann == "Activation"),
  "TF_target_Repression"  = which(!is.na(E$rho) & E$edge_type == "TF_target" & ann == "Repression"),
  "TF_target_unannotated" = which(!is.na(E$rho) & E$edge_type == "TF_target" & !ann %in% c("Activation","Repression")),
  "TF_miRNA_Activation"   = which(!is.na(E$rho) & E$edge_type == "TF_miRNA" & ann == "Activation"),
  "TF_miRNA_Repression"   = which(!is.na(E$rho) & E$edge_type == "TF_miRNA" & ann == "Repression"),
  "TF_miRNA_unannotated"  = which(!is.na(E$rho) & E$edge_type == "TF_miRNA" & !ann %in% c("Activation","Repression")),
  "miRNA_miRNA_cluster"   = which(!is.na(E$rho) & E$edge_type == "miRNA_miRNA"),
  "gene_gene"             = which(!is.na(E$rho) & E$edge_type == "gene_gene"))

rows <- list()
for (nm in names(strata)) {
  ii <- strata[[nm]]; if (!length(ii)) next
  ps <- E$predicted_sign[ii]
  sig <- !is.na(E$fdr[ii]) & E$fdr[ii] < 0.05
  agree <- if (all(!is.na(ps))) sign(E$rho[ii]) == ps else NA
  nn <- NL[NL$edge_row %in% ii, ]
  nsig <- !is.na(nn$fdr) & nn$fdr < 0.05
  nagree <- if (all(!is.na(nn$predicted_sign))) sign(nn$rho) == nn$predicted_sign else NA
  r <- data.frame(stratum = nm, n_edges = length(ii),
    n_sig_FDR05 = sum(sig), frac_sig_FDR05 = mean(sig),
    n_sig_predicted_direction = if (all(!is.na(ps))) sum(sig & agree) else NA_integer_,
    n_sig_opposite_direction  = if (all(!is.na(ps))) sum(sig & !agree) else NA_integer_,
    frac_sig_predicted_direction = if (all(!is.na(ps))) mean(sig & agree) else NA_real_,
    ratio_supporting_to_contradicting = if (all(!is.na(ps)) && sum(sig & !agree) > 0)
      sum(sig & agree) / sum(sig & !agree) else NA_real_,
    null_n_pairs = nrow(nn), null_frac_sig_FDR05 = mean(nsig),
    null_frac_sig_predicted_direction = if (all(!is.na(nn$predicted_sign))) mean(nsig & nagree) else NA_real_,
    stringsAsFactors = FALSE)
  rows[[nm]] <- r
  say(sprintf("%-22s n=%5d  sig(FDR<.05)=%4d (%.3f)  sig&predicted-dir=%s  sig&opposite=%s  | null sig&pred-dir frac=%s",
      nm, length(ii), sum(sig), mean(sig),
      ifelse(is.na(r$n_sig_predicted_direction), "NA", r$n_sig_predicted_direction),
      ifelse(is.na(r$n_sig_opposite_direction), "NA", r$n_sig_opposite_direction),
      ifelse(is.na(r$null_frac_sig_predicted_direction), "NA", sprintf("%.3f", r$null_frac_sig_predicted_direction))))
}
DIR <- do.call(rbind, rows)
write.csv(DIR, file.path(REV, "results", "directional_evidence_counts.csv"), row.names = FALSE)
say("WROTE results/directional_evidence_counts.csv rows=", nrow(DIR))

## expression-supported sub-network: sign-annotated edges, FDR<0.05, correct direction
keep <- !is.na(E$predicted_sign) & !is.na(E$fdr) & E$fdr < 0.05 & sign(E$rho) == E$predicted_sign
SUB <- E[keep, ]
write.csv(SUB, file.path(REV, "results", "expression_supported_edges.csv"), row.names = FALSE)
say("EXPRESSION-SUPPORTED SUB-NETWORK: ", nrow(SUB), " of ", sum(!is.na(E$predicted_sign) & !is.na(E$rho)),
    " sign-annotated measurable edges (", signif(100 * nrow(SUB) / sum(!is.na(E$predicted_sign) & !is.na(E$rho)), 4), "%)")
say("  by edge type: ", paste(names(table(SUB$edge_type)), table(SUB$edge_type), sep = "=", collapse = " "))
say("  distinct nodes retained: ", length(unique(c(SUB$source, SUB$target))))
say("WROTE results/expression_supported_edges.csv rows=", nrow(SUB))

## --------------------------------------------------------------------------
## (ii) SHARPER ANOVA CONTROL
## --------------------------------------------------------------------------
say("=== SHARPER ANOVA CONTROL: random mRNAs vs random miRNAs, no hub selection ===")
dg <- read.csv(file.path(REV, "results", "BRCA_DEX_genes.csv"), stringsAsFactors = FALSE)
dm <- read.csv(file.path(REV, "results", "BRCA_DEX_mirnas.csv"), stringsAsFactors = FALSE)
say("mean |logFC| all mRNA = ", signif(mean(abs(dg$logFC)), 4),
    " ; all miRNA = ", signif(mean(abs(dm$logFC)), 4),
    " ; ratio = ", signif(mean(abs(dg$logFC)) / mean(abs(dm$logFC)), 4))
wt <- wilcox.test(abs(dg$logFC), abs(dm$logFC))
say("Wilcoxon |logFC| mRNA vs miRNA over ALL tested features: W=", format(unname(wt$statistic), scientific = TRUE),
    " p = ", signif(wt$p.value, 4))
tt <- t.test(abs(dg$logFC), abs(dm$logFC))
say("Welch t on |logFC| over ALL tested features: t=", signif(unname(tt$statistic), 4),
    " df=", signif(unname(tt$parameter), 5), " p=", signif(tt$p.value, 4))

B <- 5000
p_abs <- numeric(B); p_raw <- numeric(B); p_var <- numeric(B)
for (b in seq_len(B)) {
  a <- dg[sample.int(nrow(dg), 20), ]; c2 <- dm[sample.int(nrow(dm), 20), ]
  p_abs[b] <- t.test(abs(a$logFC), abs(c2$logFC))$p.value
  p_raw[b] <- t.test(a$logFC, c2$logFC)$p.value
  p_var[b] <- var.test(a$logFC, c2$logFC)$p.value
}
say("random 20 mRNAs vs 20 miRNAs, NO hub selection, ", B, " draws:")
say("  |logFC| t-test  p<0.05 in ", signif(mean(p_abs < 0.05), 4), " of draws (median p ", signif(median(p_abs), 4), ")")
say("  logFC  t-test   p<0.05 in ", signif(mean(p_raw < 0.05), 4), " of draws (median p ", signif(median(p_raw), 4), ")")
say("  variance F-test p<0.05 in ", signif(mean(p_var < 0.05), 4), " of draws (median p ", signif(median(p_var), 4), ")")

## same, using ALL network nodes rather than 20 (higher power, the honest comparison)
nodes <- read.delim(file.path(REV, "data", "canonical_nodes.tsv"), stringsAsFactors = FALSE)
nt <- setNames(nodes$type, nodes$name)
dg$class <- unname(nt[dg$feature]); dm$class <- unname(nt[dm$feature])
netTF <- dg[which(dg$class == "TF"), ]; netG <- dg[which(dg$class == "Gene"), ]
netMi <- dm[which(dm$class == "miRNA"), ]
say("network node counts with DE stats: TF ", nrow(netTF), " Gene ", nrow(netG), " miRNA ", nrow(netMi))
allc <- data.frame(logFC = c(netTF$logFC, netG$logFC, netMi$logFC),
                   AveExpr = c(netTF$AveExpr, netG$AveExpr, netMi$AveExpr),
                   class = factor(rep(c("TF", "Gene", "miRNA"),
                                      c(nrow(netTF), nrow(netG), nrow(netMi))),
                                  levels = c("TF", "Gene", "miRNA")))
av <- summary(aov(logFC ~ class, data = allc))[[1]]
say("ANOVA on ALL network nodes (n=", nrow(allc), "): F(", av$Df[1], ",", av$Df[2], ")=",
    signif(av$`F value`[1], 4), " p=", signif(av$`Pr(>F)`[1], 4))
bt <- bartlett.test(logFC ~ class, data = allc)
lev <- function(y, g) { g <- factor(g); z <- abs(y - ave(y, g, FUN = median))
  a <- summary(aov(z ~ g))[[1]]; c(F = a$`F value`[1], p = a$`Pr(>F)`[1], df1 = a$Df[1], df2 = a$Df[2]) }
lv <- lev(allc$logFC, allc$class)
say("Bartlett on ALL network nodes: K2=", signif(unname(bt$statistic), 4), " p=", signif(bt$p.value, 4),
    "  -> equal variance ", ifelse(bt$p.value < 0.05, "REJECTED", "not rejected"))
say("Brown-Forsythe on ALL network nodes: F(", lv["df1"], ",", lv["df2"], ")=", signif(lv["F"], 4),
    " p=", signif(lv["p"], 4), "  -> equal variance ", ifelse(lv["p"] < 0.05, "REJECTED", "not rejected"))
sw <- by(allc$logFC, allc$class, function(x) shapiro.test(x)$p.value)
say("Shapiro-Wilk normality p (ALL network nodes): TF=", signif(sw[["TF"]], 4),
    " Gene=", signif(sw[["Gene"]], 4), " miRNA=", signif(sw[["miRNA"]], 4),
    "  -> normality ", ifelse(min(unlist(sw)) < 0.05, "REJECTED for at least one class", "not rejected"))
## AveExpr difference = the dynamic-range confound, quantified
tt2 <- t.test(c(netTF$AveExpr, netG$AveExpr), netMi$AveExpr)
say("mean AveExpr network mRNA (TF+Gene) = ", signif(mean(c(netTF$AveExpr, netG$AveExpr)), 4),
    " vs network miRNA = ", signif(mean(netMi$AveExpr), 4),
    " ; Welch t p = ", signif(tt2$p.value, 4))
vr <- var.test(c(netTF$logFC, netG$logFC), netMi$logFC)
say("variance ratio logFC network mRNA / network miRNA = ", signif(unname(vr$estimate), 4),
    " F-test p = ", signif(vr$p.value, 4))
## correlation between AveExpr and |logFC| -- shows fold change is a function of expression level
ct1 <- cor.test(dg$AveExpr, abs(dg$logFC), method = "spearman", exact = FALSE)
ct2 <- cor.test(dm$AveExpr, abs(dm$logFC), method = "spearman", exact = FALSE)
say("Spearman AveExpr vs |logFC|: mRNA rho=", signif(unname(ct1$estimate), 4), " p=", signif(ct1$p.value, 4),
    " ; miRNA rho=", signif(unname(ct2$estimate), 4), " p=", signif(ct2$p.value, 4))

CTRL <- data.frame(metric = c(
  "mean_abs_logFC_all_mRNA", "mean_abs_logFC_all_miRNA", "ratio_mean_abs_logFC",
  "wilcox_p_abs_logFC_all_features", "welch_t_p_abs_logFC_all_features",
  "rand20_frac_p_lt05_abs_logFC_ttest", "rand20_median_p_abs_logFC_ttest",
  "rand20_frac_p_lt05_logFC_ttest", "rand20_frac_p_lt05_variance_Ftest",
  "ANOVA_all_network_nodes_F", "ANOVA_all_network_nodes_p",
  "Bartlett_all_network_nodes_K2", "Bartlett_all_network_nodes_p",
  "BrownForsythe_all_network_nodes_F", "BrownForsythe_all_network_nodes_p",
  "Shapiro_p_TF", "Shapiro_p_Gene", "Shapiro_p_miRNA",
  "mean_AveExpr_network_mRNA", "mean_AveExpr_network_miRNA", "welch_p_AveExpr_mRNA_vs_miRNA",
  "var_ratio_logFC_network_mRNA_over_miRNA", "var_ratio_F_p",
  "spearman_AveExpr_vs_absLogFC_mRNA", "spearman_AveExpr_vs_absLogFC_miRNA"),
  value = c(mean(abs(dg$logFC)), mean(abs(dm$logFC)), mean(abs(dg$logFC)) / mean(abs(dm$logFC)),
    wt$p.value, tt$p.value,
    mean(p_abs < 0.05), median(p_abs), mean(p_raw < 0.05), mean(p_var < 0.05),
    av$`F value`[1], av$`Pr(>F)`[1], unname(bt$statistic), bt$p.value,
    unname(lv["F"]), unname(lv["p"]), sw[["TF"]], sw[["Gene"]], sw[["miRNA"]],
    mean(c(netTF$AveExpr, netG$AveExpr)), mean(netMi$AveExpr), tt2$p.value,
    unname(vr$estimate), vr$p.value, unname(ct1$estimate), unname(ct2$estimate)),
  stringsAsFactors = FALSE)
write.csv(CTRL, file.path(REV, "results", "anova_named_axes_control.csv"), row.names = FALSE)
say("WROTE results/anova_named_axes_control.csv rows=", nrow(CTRL))
say("=== 10 DONE ===")
close(con)
