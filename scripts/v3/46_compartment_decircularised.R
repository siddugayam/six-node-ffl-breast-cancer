#!/usr/bin/env Rscript
## PART 2H, decisive control.  Rebuild the "anti-correlated target" sets for miR-29 and
## miR-130a-3p from CAF-ADJUSTED partial Spearman correlations instead of raw bulk
## correlations, so that the compositional axis cannot drive set membership.
suppressPackageStartupMessages({library(data.table)})
REV <- "/path/to/revision"; OUT <- file.path(REV,"results/v3")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")
mir <- readRDS(file.path(REV,"data/brca_mirna_expr_canonical.rds"))
ge  <- readRDS(file.path(REV,"data/brca_gene_expr.rds"))
stopifnot(sum(duplicated(rownames(ge)))==0)
tum <- function(x) x[substr(x,14,15)=="01"]
cs <- intersect(tum(colnames(mir)), tum(colnames(ge)))
D <- fread(file.path(OUT,"deconv_all_celltype_estimates.csv"))
caf <- D[method=="InstaPrism_Wu" & grepl("^CAFs$", feature)]
if (nrow(caf)==0) caf <- D[grepl("CAF", feature)][method==method[1]]
msg("CAF mediator used:", caf$method[1], "/", caf$feature[1], "| n samples:", nrow(caf))
epi <- D[method==caf$method[1] & feature=="Cancer Epithelial"]
msg("epithelial feature:", if(nrow(epi)) epi$feature[1] else "none")
cs <- intersect(cs, caf$sample)
msg("samples with miRNA + gene + CAF estimate:", length(cs))
cv <- setNames(caf$value, caf$sample)[cs]
ev <- if (nrow(epi)) setNames(epi$value, epi$sample)[cs] else NULL
mir29 <- colMeans(mir[grep("hsa-miR-29[abc]$", rownames(mir)), cs, drop=FALSE], na.rm=TRUE)
m130  <- mir["hsa-miR-130a", cs]
msg("miR-29 rows:", paste(grep("hsa-miR-29[abc]$", rownames(mir), value=TRUE), collapse=","))
G <- ge[, cs, drop=FALSE]
keep <- rowMeans(!is.na(G) & G > 0) >= 0.25
G <- G[keep, ]; msg("genes retained:", nrow(G))
stopifnot(mean(grepl("^[A-Za-z]", rownames(G))) > 0.95)   # rule 4 guard
rk <- function(m) t(apply(m, 1, rank, na.last="keep"))
Gr <- rk(G)
resid_of <- function(Y, Z) { Q <- qr(cbind(1, Z)); t(qr.resid(Q, t(Y))) }
Zc <- cbind(rank(cv), if(!is.null(ev)) rank(ev) else NULL)
Gres <- resid_of(Gr, Zc)
partial <- function(x) {
  xr <- rank(x); xres <- qr.resid(qr(cbind(1, Zc)), xr)
  as.numeric(cor(t(Gres), xres))
}
raw_rho <- function(x) as.numeric(cor(t(Gr), rank(x)))
res <- data.table(gene=rownames(G),
                  rho29_raw=raw_rho(mir29), rho29_adj=partial(mir29),
                  rho130_raw=raw_rho(m130), rho130_adj=partial(m130))
n <- length(cs); k <- ncol(Zc)
pv <- function(r, df) 2*pt(-abs(r*sqrt(df/(1-r^2))), df)
res[, p29_raw := pv(rho29_raw, n-2)][, p29_adj := pv(rho29_adj, n-2-k)]
res[, p130_raw := pv(rho130_raw, n-2)][, p130_adj := pv(rho130_adj, n-2-k)]
for (cc in c("p29_raw","p29_adj","p130_raw","p130_adj"))
  res[[sub("^p","fdr",cc)]] <- p.adjust(res[[cc]], "BH")
fwrite(res, file.path(OUT,"sc_verify_partial_correlations.csv"))
msg("raw vs CAF-adjusted rho, Spearman across genes: miR-29", round(cor(res$rho29_raw,res$rho29_adj,method="spearman"),3),
    "| miR-130a", round(cor(res$rho130_raw,res$rho130_adj,method="spearman"),3))
## restrict to database-validated targets, as the original sets did
a29 <- fread(file.path(OUT,"screens_mir29_target_set.csv")); a130 <- fread(file.path(OUT,"mir130a_target_set.csv"))
v29 <- a29[tier %in% c("STRONG_lowthroughput","WEAK_highthroughput"), gene]
v130 <- a130[tier %in% c("STRONG_lowthroughput","WEAK_highthroughput"), gene]
mk <- function(genes, rhoc, fdrc) res[gene %in% genes & get(fdrc) < 0.05 & get(rhoc) < 0, gene]
sets <- list(
  mir29_anticorr_RAW  = mk(v29,  "rho29_raw",  "fdr29_raw"),
  mir29_anticorr_ADJ  = mk(v29,  "rho29_adj",  "fdr29_adj"),
  mir130a_anticorr_RAW= mk(v130, "rho130_raw", "fdr130_raw"),
  mir130a_anticorr_ADJ= mk(v130, "rho130_adj", "fdr130_adj"))
for (nm in names(sets)) msg("set", nm, ":", length(sets[[nm]]), "genes")
saveRDS(sets, "/path/to/scratch/decirc_sets.rds")
for (nm in names(sets)) fwrite(data.table(gene=sets[[nm]]), paste0("/path/to/scratch/set_",nm,".csv"))
msg("overlap RAW n ADJ: miR-29", length(intersect(sets$mir29_anticorr_RAW, sets$mir29_anticorr_ADJ)),
    "| miR-130a", length(intersect(sets$mir130a_anticorr_RAW, sets$mir130a_anticorr_ADJ)))
