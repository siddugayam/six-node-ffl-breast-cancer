## 23_metabolic_associations.R --------------------------------------------
## Part C: correlate every metabolic score with hsa-miR-130a-3p, the miR-29
## family, the let-7 family and the FFL 3-node / higher-order module scores.
## BH-FDR across ALL pathway x feature pairs.
suppressPackageStartupMessages({library(data.table)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
set.seed(1)
ss <- readRDS(file.path(RES,"metabolic_ssgsea_scores.rds"))
gv <- readRDS(file.path(RES,"metabolic_gsva_scores.rds"))
inv <- fread(file.path(RES,"metabolic_geneset_inventory.csv"))
ph  <- readRDS(file.path(ROOT,"data/brca_pheno.rds"))
MI  <- readRDS(file.path(ROOT,"data/brca_mirna_expr.rds"))

tum <- ph$sample[ph$sample_type == "Primary Tumor"]
smp <- intersect(intersect(colnames(ss), tum), colnames(MI))
cat("paired primary tumours with mRNA + miRNA:", length(smp), "\n")
SS <- ss[, smp, drop=FALSE]; GV <- gv[, smp, drop=FALSE]; MM <- MI[, smp, drop=FALSE]

mir_feats <- c("hsa-miR-130a-3p","hsa-miR-130a-5p","hsa-miR-130b-3p","hsa-miR-301a-3p","hsa-miR-301b",
               "hsa-miR-29a-3p","hsa-miR-29b-3p","hsa-miR-29c-3p","hsa-miR-29a-5p","hsa-miR-29c-5p",
               grep("^hsa-let-7", rownames(MM), value=TRUE))
mir_feats <- intersect(mir_feats, rownames(MM))
cat("miRNA features:", length(mir_feats), "\n")
## family aggregate (mean of log2 arms) for miR-130a/301 seed family and miR-29 family
fam <- list(
  FAMILY_miR130_301_seed = c("hsa-miR-130a-3p","hsa-miR-130b-3p","hsa-miR-301a-3p","hsa-miR-301b"),
  FAMILY_miR29_3p        = c("hsa-miR-29a-3p","hsa-miR-29b-3p","hsa-miR-29c-3p"),
  FAMILY_let7_5p         = grep("^hsa-let-7.*-5p$", rownames(MM), value=TRUE))
FM <- do.call(rbind, lapply(fam, function(v) colMeans(MM[intersect(v, rownames(MM)), , drop=FALSE])))
featM <- rbind(MM[mir_feats, , drop=FALSE], FM)

mod_ids <- grep("^FFLCLASS\\|", rownames(SS), value=TRUE)
ctl_ids <- grep("^CONTROL\\|", rownames(SS), value=TRUE)
modM <- SS[c(mod_ids, ctl_ids), , drop=FALSE]
FEAT <- rbind(featM, modM)
cat("total features:", nrow(FEAT), "( miRNA", nrow(featM), "+ module/control", nrow(modM), ")\n")

met_ids <- setdiff(rownames(SS), c(mod_ids, ctl_ids))
corr_block <- function(SCORE, tag) {
  rs <- t(apply(SCORE[met_ids,,drop=FALSE], 1, rank)); rf <- t(apply(FEAT, 1, rank))
  R <- cor(t(rs), t(rf))                      # Spearman via ranks
  n <- length(smp)
  tt <- R*sqrt((n-2)/(1-R^2)); P <- 2*pt(-abs(tt), df=n-2)
  dt <- data.table(set_id=rep(rownames(R), times=ncol(R)),
                   feature=rep(colnames(R), each=nrow(R)),
                   score=tag, rho=as.vector(R), p=as.vector(P), n=n)
  dt
}
out <- rbind(corr_block(SS,"ssGSEA"), corr_block(GV,"GSVA"))
out[, FDR := p.adjust(p, "BH"), by=score]       # BH across all pathway x feature pairs
out[, feature_type := fifelse(grepl("^hsa-|^FAMILY_", feature), "miRNA",
                       fifelse(grepl("^FFLCLASS", feature), "FFL_module", "context_control"))]
out <- merge(out, inv[, .(set_id, collection, n_genes_in_universe)], by="set_id", all.x=TRUE)
setorder(out, score, feature, p)
fwrite(out, file.path(RES,"metabolic_feature_correlations.csv"))
cat("pathway x feature pairs per score matrix:", out[score=="ssGSEA", .N], "\n")

show <- function(f, k=15) {
  d <- out[score=="ssGSEA" & feature==f]
  cat("\n=====", f, " (", nrow(d), "pathways; FDR<0.05:", sum(d$FDR<0.05), ")\n")
  cat("-- most POSITIVE --\n")
  print(d[order(-rho)][1:k, .(set_id, collection, rho=round(rho,3), FDR=signif(FDR,3))])
  cat("-- most NEGATIVE --\n")
  print(d[order(rho)][1:k, .(set_id, collection, rho=round(rho,3), FDR=signif(FDR,3))])
}
show("hsa-miR-130a-3p"); show("hsa-miR-29a-3p", 10); show("FAMILY_let7_5p", 10)
show("FFLCLASS|3node_union", 10); show("FFLCLASS|higherorder_not_3node", 10)

## GSVA sensitivity for miR-130a
a <- out[score=="ssGSEA" & feature=="hsa-miR-130a-3p", .(set_id, rho_ss=rho)]
b <- out[score=="GSVA"   & feature=="hsa-miR-130a-3p", .(set_id, rho_gv=rho)]
m <- merge(a,b,by="set_id")
cat("\nmiR-130a rho: ssGSEA vs GSVA Pearson", round(cor(m$rho_ss,m$rho_gv),3),
    "; sign agreement", round(mean(sign(m$rho_ss)==sign(m$rho_gv)),3), "\n")
fwrite(m, file.path(RES,"metabolic_mir130a_rho_method_concordance.csv"))
cat("DONE 23\n")
