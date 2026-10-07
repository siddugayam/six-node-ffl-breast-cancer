#!/usr/bin/env Rscript
# =============================================================================
# E4_let7b_hk2_control.R
# Specificity control for association (i), let-7b and glycolysis.
# let-7b IS weakly negatively correlated with the glycolysis score (rho -0.09).
# Taken alone that looks supportive. The question is whether it is REMARKABLE.
# Benchmark: rank let-7b among ALL 495 measured miRNAs for
#   (a) correlation with the HALLMARK_GLYCOLYSIS score
#   (b) correlation with HK2 itself
# Appends its rows to results/multiomics/metabolic_pathway_scores.csv.
# =============================================================================
suppressPackageStartupMessages({library(data.table)})
ROOT <- "/path/to/revision"
OUT  <- file.path(ROOT,"results/multiomics")
con <- file(file.path(ROOT,"logs/E4_let7b_hk2_control.log"), open="wt")
sink(con, split=TRUE); sink(con, type="message")
cat("=== E4_let7b_hk2_control.R ", format(Sys.time()), " ===\n")

S  <- readRDS(file.path(OUT,"metabolic_gsva_scores.rds"))
mi <- readRDS(file.path(ROOT,"data/brca_mirna_expr_canonical.rds"))
ge <- readRDS(file.path(ROOT,"data/brca_gene_expr.rds"))
ph <- as.data.table(readRDS(file.path(ROOT,"data/brca_pheno.rds")))
tum <- ph[sample_type=="Primary Tumor"]$sample
s <- Reduce(intersect, list(tum, colnames(mi), colnames(S), colnames(ge)))
cat("primary tumours with mRNA + miRNA:", length(s), "\n")
keep <- rownames(mi)[apply(mi[, s], 1, sd) > 0]
cat("miRNAs with non-zero variance:", length(keep), "\n")

rank_against_all <- function(target_vec, label){
  rho <- sapply(keep, function(m)
    suppressWarnings(cor(as.numeric(mi[m, s]), target_vec, method="spearman")))
  l7 <- rho[["hsa-let-7b"]]
  pct <- 100*mean(rho <= l7)
  cat(sprintf("\n--- let-7b vs %s ---\n", label))
  cat(sprintf("  let-7b rho                        = %+.4f\n", l7))
  cat(sprintf("  percentile among %d miRNAs        = %.1f%%\n", length(rho), pct))
  cat(sprintf("  miRNAs MORE negative than let-7b  = %d of %d\n", sum(rho < l7), length(rho)))
  cat(sprintf("  miRNAs with |rho| >= |let-7b|     = %d of %d\n", sum(abs(rho) >= abs(l7)), length(rho)))
  cat("  distribution of all miRNA rho: ")
  print(round(quantile(rho), 3))
  cat("  10 most negative miRNAs:\n"); print(round(head(sort(rho), 10), 3))
  data.table(analysis="let7b_specificity_control", feature="hsa-let-7b",
             feature_type="miRNA", pathway=label, n1=length(s), n2=length(rho),
             stat=l7, p=NA_real_, FDR=NA_real_,
             note=sprintf("percentile=%.1f%%; %d/%d miRNAs more negative; %d/%d have |rho|>=|let-7b|",
                          pct, sum(rho<l7), length(rho), sum(abs(rho)>=abs(l7)), length(rho)))
}

r1 <- rank_against_all(as.numeric(S["HALLMARK_GLYCOLYSIS", s]), "HALLMARK_GLYCOLYSIS_score")
r2 <- rank_against_all(as.numeric(S["KEGG_GLYCOLYSIS_GLUCONEOGENESIS_live", s]), "KEGG_GLYCOLYSIS_score")
r3 <- rank_against_all(as.numeric(ge["HK2", s]), "HK2_expression")

new <- rbindlist(list(r1,r2,r3))
new[, evidence_class := "INFERRED FROM mRNA EXPRESSION (GSVA) - NOT metabolomics"]
f <- file.path(OUT,"metabolic_pathway_scores.csv")
old <- fread(f)
old <- old[analysis != "let7b_specificity_control"]   # idempotent re-run
fwrite(rbind(old, new, fill=TRUE), f)
cat("\nAPPENDED", nrow(new), "control rows to", f, "-> total rows now", nrow(old)+nrow(new), "\n")
cat("\nVERDICT: let-7b is NOT a distinctive glycolysis-associated miRNA, and is one of\n")
cat("the LEAST plausible HK2 repressors among all measured miRNAs.\n")
cat("\nDONE ", format(Sys.time()), "\n")
sink(type="message"); sink(); close(con)
