#!/usr/bin/env Rscript
## PART 1E) One integrated essentiality matrix over the 587 network nodes:
##   CRISPR (DepMap 24Q4 Chronos, breast lines)
##   RNAi   (DEMETER2 v6 combined, breast lines)
##   SCORE  (Sanger Project Score 2, breast models)
##   drug-sensitivity association (best |rho| across PRISM and GDSC2 module-score tests)
##   phenotypic-screen hit rate (BioGRID ORCS breast screens)
## with a consensus call per node.
suppressPackageStartupMessages({library(data.table)})
REV <- "/path/to/revision"; OUT <- file.path(REV,"results/v3")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

N <- fread(file.path(OUT,"screens_demeter2_network_nodes.csv"))   # node, type, is_FFL_hub, D2, Chronos
msg("network nodes:", nrow(N))


## Project SCORE ------------------------------------------------------------
sf <- file.path(OUT,"screens_projectscore_network_nodes.csv")
if (file.exists(sf)) {
  S <- fread(sf); N <- merge(N, S[, .(node, score_breast=score_breast, score_frac_lines_significant=frac_significant, score_n_lines=n_lines)],
                             by="node", all.x=TRUE)
  ## Project Score's scaled Bayesian factor is signed the opposite way to Chronos/DEMETER2
  msg("Project SCORE merged: nodes with a score =", sum(!is.na(N$score_breast)))
} else { N[, `:=`(score_breast=NA_real_, score_frac_lines_significant=NA_real_, score_n_lines=NA_integer_)]
  msg("Project SCORE table not present -- SCORE column left NA") }

## drug association ---------------------------------------------------------
## per node we cannot test a drug association directly (the scores are set-level),
## so we record whether the node belongs to a set that had ANY FDR<0.05 drug association.
P <- fread(file.path(OUT,"screens_prism_associations.csv"))
G <- if (file.exists(file.path(OUT,"screens_gdsc2_associations.csv"))) fread(file.path(OUT,"screens_gdsc2_associations.csv")) else NULL
sig_sets <- unique(c(P[fdr<0.05, score], if(!is.null(G)) G[fdr<0.05, score]))
msg("module scores with any FDR<0.05 drug association:", if(length(sig_sets)) paste(sig_sets, collapse=",") else "NONE")
setmem <- list(
  COLLAGEN = c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL11A1"),
  MIR29_TARGET = fread(file.path(OUT,"screens_mir29_anticorrelated_genes.csv"))$gene,
  MIR130A_TARGET = fread(file.path(OUT,"mir130a_anticorrelated_validated_genes.csv"))$gene,
  MIR29_STRONG = fread(file.path(OUT,"screens_mir29_target_set.csv"))[tier=="STRONG_lowthroughput", gene],
  MIR130A_STRONG = fread(file.path(OUT,"mir130a_target_set.csv"))[tier=="STRONG_lowthroughput", gene])
N[, in_drug_associated_set := node %in% unlist(setmem[sig_sets])]
for (nm in names(setmem)) N[[paste0("in_", nm)]] <- N$node %in% setmem[[nm]]

## ORCS phenotypic screens --------------------------------------------------
of <- file.path(OUT,"screens_orcs_gene_level.csv")
if (file.exists(of)) {
  O <- fread(of)[group != "migration_nonbreast"]
  Og <- O[, .(orcs_breast_screens_tested=.N, orcs_breast_screens_hit=sum(hit),
              orcs_hit_rate=mean(hit)), by=.(node=gene)]
  N <- merge(N, Og, by="node", all.x=TRUE)
  msg("ORCS gene-level merged: nodes with >=1 breast screen =", sum(!is.na(N$orcs_breast_screens_tested)))
}

## consensus ----------------------------------------------------------------
ess <- function(x, thr=-0.5) fifelse(is.na(x), NA, x < thr)
N[, crispr_essential := ess(chronos_breast)]
N[, rnai_essential   := ess(d2_breast)]
## Project Score: a gene is called essential when it is a significant fitness gene in more
## than half of the breast models (the release's own binary significance matrix).
N[, score_essential  := fifelse(is.na(score_frac_lines_significant), NA, score_frac_lines_significant > 0.5)]
N[, n_modalities := (!is.na(chronos_breast)) + (!is.na(d2_breast)) + (!is.na(score_breast))]
N[, n_essential  := rowSums(cbind(crispr_essential, rnai_essential, score_essential), na.rm=TRUE)]
N[, consensus := fifelse(n_modalities==0, "not screened",
              fifelse(n_essential==0, "non-essential in every modality",
               fifelse(n_essential==n_modalities, "essential in every modality tested",
                "essential in some modalities only")))]
setcolorder(N, c("node","type","is_FFL_hub","chronos_breast","d2_breast","score_breast",
                 "crispr_essential","rnai_essential","score_essential","n_modalities","n_essential","consensus"))
setorder(N, -n_essential, chronos_breast)
fwrite(N, file.path(OUT,"screens_integrated_essentiality_matrix.csv"))
msg("wrote screens_integrated_essentiality_matrix.csv:", nrow(N), "rows x", ncol(N), "cols")

msg("\n=== consensus call across the 587 network nodes ===")
print(N[, .N, by=consensus][order(-N)])
msg("\n=== by node type ===")
print(dcast(N[, .N, by=.(type, consensus)], type ~ consensus, value.var="N", fill=0))
msg("\n=== nodes essential in EVERY modality tested (>=2 modalities) ===")
print(N[consensus=="essential in every modality tested" & n_modalities>=2,
        .(node, type, is_FFL_hub, chronos_breast, d2_breast, score_breast, n_modalities)], nrows=60)
msg("\n=== the focus genes ===")
print(N[node %in% c("NFKB1","RELA","SP1","ETS1","COL1A1","COL3A1","MYC","TP53","STAT3","HIF1A","MKL1","VEGFA","ESR1","XIAP","PTEN"),
        .(node, type, is_FFL_hub, chronos_breast, d2_breast, score_breast, n_modalities, n_essential, consensus,
          orcs_hit_rate)], nrows=20)
msg("\n=== FFL hubs ===")
print(N[is_FFL_hub==TRUE, .(node, type, chronos_breast, d2_breast, score_breast, n_essential, consensus)][order(chronos_breast)], nrows=65)
msg("\nagreement between modalities among nodes measured in all three:")
sub <- N[n_modalities==3]
if (nrow(sub) > 5) {
  msg("  n =", nrow(sub))
  msg("  CRISPR vs RNAi  Spearman:", round(cor(sub$chronos_breast, sub$d2_breast, method='spearman'),3))
  msg("  CRISPR vs SCORE Spearman:", round(cor(sub$chronos_breast, sub$score_breast, method='spearman'),3))
  msg("  RNAi   vs SCORE Spearman:", round(cor(sub$d2_breast, sub$score_breast, method='spearman'),3))
}
msg("DONE 24")
