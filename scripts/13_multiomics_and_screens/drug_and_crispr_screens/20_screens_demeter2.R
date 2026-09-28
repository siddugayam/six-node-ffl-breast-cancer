#!/usr/bin/env Rscript
## PART 1A) DepMap RNAi (DEMETER2 v6 combined: Achilles + DRIVE + Marcotte) as an
## orthogonal loss-of-function modality to the CRISPR/Chronos analysis already done.
suppressPackageStartupMessages({library(data.table)})
set.seed(20260909)
REV <- "/path/to/revision"
OUT <- file.path(REV,"results/v3"); CA <- file.path(REV,"cache/v7")
RAW <- "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data"
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

## ---------------- load DEMETER2 ----------------
D <- fread(file.path(CA,"D2_combined_gene_dep_scores.csv"))
setnames(D, 1, "gene_raw")
D[, gene := sub(" \\(\\d+\\)$", "", gene_raw)]
msg("DEMETER2 combined matrix:", nrow(D), "genes x", ncol(D)-2, "cell lines")
stopifnot(!any(duplicated(D$gene)))
cls <- setdiff(names(D), c("gene_raw","gene"))
bcl <- grep("_BREAST$", cls, value=TRUE)
msg("breast cell lines in DEMETER2:", length(bcl))
si <- fread(file.path(CA,"D2_sample_info.csv"))
msg("of those, in Achilles:", sum(si[CCLE_ID %in% bcl, in_Achilles]),
    "| DRIVE:", sum(si[CCLE_ID %in% bcl, in_DRIVE]),
    "| Marcotte:", sum(si[CCLE_ID %in% bcl, in_Marcotte]))

M  <- as.matrix(D[, ..bcl]);  rownames(M)  <- D$gene
ocl <- setdiff(cls, bcl)
Mo <- as.matrix(D[, ..ocl]); rownames(Mo) <- D$gene
d2_breast     <- rowMeans(M,  na.rm=TRUE)
d2_other      <- rowMeans(Mo, na.rm=TRUE)
d2_frac_ess   <- rowMeans(M < -0.5, na.rm=TRUE)
d2_n_measured <- rowSums(!is.na(M))
msg("mean D2 across all genes in breast lines:", round(mean(d2_breast, na.rm=TRUE),4),
    "| genes with mean D2 < -0.5:", sum(d2_breast < -0.5, na.rm=TRUE),
    "(", round(100*mean(d2_breast < -0.5, na.rm=TRUE),2), "% )")

## sanity: known common-essential vs known non-essential controls
for (g in c("RPL7","PSMA1","POLR2A","EIF3B","KRAS","MYC","GAPDH","PTEN","CDKN2A")) {
  if (g %in% rownames(M)) msg("  control", g, ": mean D2 breast =", round(d2_breast[g],3),
                              "| n lines =", d2_n_measured[g])
}

## ---------------- CRISPR (Chronos) for the same genes ----------------
mod <- fread(file.path(REV,"data/depmap/Model_24Q4.csv"))
breast_ids <- mod[OncotreeLineage=="Breast", ModelID]
ce <- fread(file.path(RAW,"CRISPRGeneEffect.csv")); setnames(ce,1,"ModelID")
gsC <- sub(" \\(\\d+\\)$","", setdiff(names(ce),"ModelID")); stopifnot(!any(duplicated(gsC)))
bidC <- intersect(ce$ModelID, breast_ids)
msg("DepMap CRISPR breast lines:", length(bidC), "| genes:", length(gsC))
CM <- as.matrix(ce[ModelID %in% bidC, -1]); colnames(CM) <- gsC
chronos_breast <- colMeans(CM, na.rm=TRUE); chronos_frac_ess <- colMeans(CM < -0.5, na.rm=TRUE)

common <- intersect(rownames(M), gsC)
cc <- cor(d2_breast[common], chronos_breast[common], method="spearman", use="complete.obs")
ccp <- cor(d2_breast[common], chronos_breast[common], method="pearson", use="complete.obs")
msg("CRISPR-vs-RNAi concordance over", length(common), "shared genes (breast-line means): Spearman",
    round(cc,4), "| Pearson", round(ccp,4))

## ---------------- focus genes ----------------
nodes <- fread(file.path(REV,"data/canonical_nodes.tsv"))
ctrl  <- fread(file.path(OUT,"systems_controllability_nodes.csv"))
hubs  <- ctrl[is_FFL_hub==TRUE, node]
msg("FFL hubs from the manuscript's criterion:", length(hubs),
    "| protein-coding of those:", sum(hubs %in% rownames(M)))
focus <- c("NFKB1","RELA","SP1","ETS1","COL1A1","COL3A1","MYC")
FOC <- data.table(gene=focus)
FOC[, d2_breast := d2_breast[gene]][, d2_other_lineages := d2_other[gene]]
FOC[, d2_frac_lines_essential := d2_frac_ess[gene]][, d2_n_lines := d2_n_measured[gene]]
FOC[, chronos_breast := chronos_breast[gene]][, chronos_frac_lines_essential := chronos_frac_ess[gene]]
## percentile of each gene within the genome-wide breast-line distribution
## (vectors kept out of the data.table scope so column names cannot mask them)
.d2b <- d2_breast; .chb <- chronos_breast
FOC[, d2_percentile_breast := sapply(focus, function(g)
      if (!g %in% names(.d2b) || is.na(.d2b[[g]])) NA_real_ else mean(.d2b <= .d2b[[g]], na.rm=TRUE))]
FOC[, chronos_percentile_breast := sapply(focus, function(g)
      if (!g %in% names(.chb) || is.na(.chb[[g]])) NA_real_ else mean(.chb <= .chb[[g]], na.rm=TRUE))]
FOC[, call_RNAi := fifelse(is.na(d2_breast), "not screened",
                    fifelse(d2_breast < -0.5, "essential",
                     fifelse(d2_breast < -0.25, "weakly depleting", "non-essential")))]
FOC[, call_CRISPR := fifelse(is.na(chronos_breast), "not screened",
                      fifelse(chronos_breast < -0.5, "essential",
                       fifelse(chronos_breast < -0.25, "weakly depleting", "non-essential")))]
msg("FOCUS GENES (DEMETER2 RNAi vs DepMap CRISPR, breast lines):")
print(FOC, digits=3)
fwrite(FOC, file.path(OUT,"screens_demeter2_focus_genes.csv"))

## all network nodes.  Three network symbols were renamed by HGNC after the network was
## built and are absent from the current DepMap releases under the old name; map them.
alias <- c(MKL1="MRTFA", CTGF="CCN2", CYR61="CCN1")
## prefer the network's own symbol; fall back to the current HGNC symbol only if absent
lookup <- function(v, nm) {
  key <- nm
  swap <- !(nm %in% names(v)) & (nm %in% names(alias))
  key[swap] <- alias[nm[swap]]
  out <- v[key]; names(out) <- nm; out }
ND <- data.table(node=nodes$name, type=nodes$type)
ND[, is_FFL_hub := node %in% hubs]
ND[, d2_breast := lookup(d2_breast, node)][, d2_other_lineages := lookup(d2_other, node)]
ND[, d2_frac_lines_essential := lookup(d2_frac_ess, node)][, d2_n_lines := lookup(d2_n_measured, node)]
ND[, chronos_breast := lookup(chronos_breast, node)][, chronos_frac_lines_essential := lookup(chronos_frac_ess, node)]
msg("alias-resolved nodes:", paste(names(alias), "->", alias, collapse="; "))
fwrite(ND, file.path(OUT,"screens_demeter2_network_nodes.csv"))
msg("network nodes with a DEMETER2 score:", sum(!is.na(ND$d2_breast)), "of", nrow(ND),
    "| with a CRISPR score:", sum(!is.na(ND$chronos_breast)))
msg("FFL hubs: mean D2 breast =", round(mean(ND[is_FFL_hub==TRUE, d2_breast], na.rm=TRUE),4),
    "vs non-hub protein-coding nodes =",
    round(mean(ND[is_FFL_hub==FALSE & !is.na(d2_breast), d2_breast], na.rm=TRUE),4))
wt <- wilcox.test(ND[is_FFL_hub==TRUE, d2_breast], ND[is_FFL_hub==FALSE, d2_breast])
msg("  Wilcoxon hub vs non-hub D2: p =", signif(wt$p.value,3))

## ---------------- set-level test with the project's decile-matched null ----------------
A130 <- fread(file.path(OUT,"mir130a_target_set.csv"))
R130 <- fread(file.path(OUT,"mir130a_tcga_target_correlations.csv"))
A29  <- fread(file.path(OUT,"screens_mir29_target_set.csv"))
R29  <- fread(file.path(OUT,"screens_mir29_tcga_correlations.csv"))
edges <- fread(file.path(REV,"data/canonical_edges.tsv"))
universe <- R130$gene
gs <- rownames(M)

sets <- list(
  mir130a_STRONG              = A130[tier=="STRONG_lowthroughput", gene],
  mir130a_VALIDATED_anticorr  = R130[!is.na(fdr_VALIDATED_any) & fdr_VALIDATED_any<0.05 & rho_3p<0, gene],
  mir130a_TARGETSCAN          = A130[in_targetscan==TRUE, gene],
  mir130a_AUTHOR_NETWORK      = edges[source=="hsa-miR-130a" & edge_type=="miRNA_target", target],
  mir29_STRONG                = A29[tier=="STRONG_lowthroughput", gene],
  mir29_anticorrelated        = fread(file.path(OUT,"screens_mir29_anticorrelated_genes.csv"))$gene,
  mir29_TARGETSCAN            = A29[in_targetscan==TRUE, gene],
  FFL_hubs                    = hubs,
  network_TFs                 = nodes[type=="TF", name],
  network_genes               = nodes[type=="Gene", name]
)
sets <- lapply(sets, function(g) intersect(intersect(unique(g), gs), universe))
for (nm in names(sets)) msg("set", nm, ": screened in DEMETER2 & TCGA-expressed:", length(sets[[nm]]))

mexp <- setNames(R130$mean_expr, R130$gene)
pooluniv <- intersect(gs, universe)
brk <- unique(quantile(mexp[pooluniv], probs=seq(0,1,0.1)))
dec <- setNames(as.integer(cut(mexp[pooluniv], breaks=brk, include.lowest=TRUE)), pooluniv)
B <- 2000L
res <- rbindlist(lapply(names(sets), function(nm){
  g <- sets[[nm]]; if(length(g) < 5) return(NULL)
  ## null pool excludes anything in the set itself and any annotated target of the same miRNA
  excl <- if (grepl("^mir130a", nm)) A130$gene else if (grepl("^mir29", nm)) A29$gene else g
  pool <- setdiff(pooluniv, excl)
  pbd <- split(pool, dec[pool])
  obs_mean <- mean(d2_breast[g], na.rm=TRUE); obs_ess <- mean(d2_frac_ess[g] > 0.5, na.rm=TRUE)
  nmx <- matrix(unlist(lapply(dec[g], function(d) sample(pbd[[as.character(d)]], B, replace=TRUE))),
                nrow=length(g), byrow=TRUE)
  nm_mean <- colMeans(matrix(d2_breast[nmx], nrow=length(g)), na.rm=TRUE)
  nm_ess  <- colMeans(matrix(d2_frac_ess[nmx] > 0.5, nrow=length(g)), na.rm=TRUE)
  data.table(set=nm, n=length(g),
             obs_mean_D2_breast=obs_mean, null_mean_D2=mean(nm_mean), null_sd_D2=sd(nm_mean),
             emp_p_more_essential=(1+sum(nm_mean<=obs_mean))/(1+B),
             obs_frac_essential_most_lines=obs_ess, null_frac_essential=mean(nm_ess),
             emp_p_frac_essential=(1+sum(nm_ess>=obs_ess))/(1+B))
}))
msg("DEMETER2 set-level essentiality vs expression-decile-matched null (B =", B, "):")
print(res, digits=3)
fwrite(res, file.path(OUT,"screens_demeter2_setlevel.csv"))

## miR-130a strong targets, gene by gene
st <- data.table(gene=sets$mir130a_STRONG)
st[, d2_breast := d2_breast[gene]][, d2_frac_lines_essential := d2_frac_ess[gene]]
st[, d2_other_lineages := d2_other[gene]][, chronos_breast := chronos_breast[gene]]
st <- merge(st, R130[, .(gene, rho_3p)], by="gene", all.x=TRUE)
setorder(st, d2_breast)
fwrite(st, file.path(OUT,"screens_demeter2_mir130a_strong_targets.csv"))
msg("miR-130a STRONG targets essential by RNAi (mean D2 < -0.5):",
    sum(st$d2_breast < -0.5, na.rm=TRUE), "of", nrow(st),
    "| by CRISPR (Chronos < -0.5):", sum(st$chronos_breast < -0.5, na.rm=TRUE))
msg("miR-130a STRONG targets with POSITIVE mean D2 (knockdown increases fitness):",
    sum(st$d2_breast > 0, na.rm=TRUE), "of", nrow(st))
print(head(st, 12), digits=3)

## miR-29 strong targets
st29 <- data.table(gene=sets$mir29_STRONG)
st29[, d2_breast := d2_breast[gene]][, chronos_breast := chronos_breast[gene]]
st29 <- merge(st29, R29[, .(gene, rho_mir29)], by="gene", all.x=TRUE)
setorder(st29, d2_breast)
fwrite(st29, file.path(OUT,"screens_demeter2_mir29_strong_targets.csv"))
msg("miR-29 STRONG targets essential by RNAi:", sum(st29$d2_breast < -0.5, na.rm=TRUE), "of", nrow(st29))
msg("DONE 20")
