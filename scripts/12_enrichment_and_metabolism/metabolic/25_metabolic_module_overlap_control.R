## 25_metabolic_module_overlap_control.R ----------------------------------
## The FFL module scores correlate with metabolic scores at rho up to 0.87.
## Those sets SHARE GENES. Two controls that remove the shared-gene artefact:
##  (a) metabolism-free module scores (module genes minus every gene that
##      appears in any KEGG/Reactome/Hallmark/Human-GEM metabolic set)
##  (b) per-pair leave-shared-out: recompute the metabolic score without the
##      genes it shares with that module, for the strongest pairs.
suppressPackageStartupMessages({library(data.table); library(GSVA); library(BiocParallel)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
set.seed(1)
E <- readRDS(file.path(ROOT,"data/brca_gene_expr.rds")); sdv <- apply(E,1,sd); Ef <- E[sdv>0,]
ss <- readRDS(file.path(RES,"metabolic_ssgsea_scores.rds"))
sets <- readRDS(file.path(RES,"metabolic_genesets.rds"))
ph <- readRDS(file.path(ROOT,"data/brca_pheno.rds")); ph <- ph[match(colnames(ss), ph$sample),]
tum <- ph$sample[ph$sample_type=="Primary Tumor"]
inv <- fread(file.path(RES,"metabolic_geneset_inventory.csv"))
bp <- MulticoreParam(workers=12, RNGseed=1)
sp <- function(x,y){ct<-suppressWarnings(cor.test(x,y,method="spearman",exact=FALSE));c(rho=unname(ct$estimate),p=ct$p.value)}

mod_ids <- grep("^FFLCLASS\\|", names(sets), value=TRUE)
met_ids <- grep("^(KEGG|KEGG_EXTRA|REACTOME|HALLMARK|HUMANGEM|CURATED)\\|", names(sets), value=TRUE)
met_union <- unique(unlist(sets[grep("^(KEGG|REACTOME|HALLMARK|HUMANGEM)\\|", names(sets), value=TRUE)]))
cat("union of all metabolic-collection genes:", length(met_union), "\n")

## ---- (a) metabolism-free module scores ----------------------------------
free <- lapply(sets[mod_ids], function(g) setdiff(g, met_union))
cat("\nmodule sizes before/after removing metabolic genes:\n")
print(data.table(module=mod_ids, n_full=lengths(sets[mod_ids]), n_metabfree=lengths(free)))
free <- free[lengths(free) >= 10]
names(free) <- paste0(names(free), "|METABFREE")
ssf <- gsva(ssgseaParam(Ef, free, alpha=0.25, normalize=TRUE, minSize=5), BPPARAM=bp, verbose=FALSE)

pairsA <- rbindlist(lapply(met_ids, function(a) rbindlist(lapply(names(free), function(m) {
  m0 <- sub("\\|METABFREE$","",m)
  r_full <- sp(ss[a,tum], ss[m0,tum]); r_free <- sp(ss[a,tum], ssf[m,tum])
  ov <- length(intersect(sets[[a]], sets[[m0]]))
  data.table(metabolic_set=a, module=m0, n_shared_genes=ov,
             rho_full=r_full["rho"], p_full=r_full["p"],
             rho_metabfree=r_free["rho"], p_metabfree=r_free["p"],
             attenuation=r_full["rho"]-r_free["rho"])
}))))
pairsA[, `:=`(FDR_full=p.adjust(p_full,"BH"), FDR_metabfree=p.adjust(p_metabfree,"BH"))]
pairsA <- merge(pairsA, inv[,.(metabolic_set=set_id, collection)], by="metabolic_set", all.x=TRUE)
setorder(pairsA, -rho_full)
fwrite(pairsA, file.path(RES,"metabolic_module_overlap_control.csv"))
cat("\n-- (a) top 15 module-metabolic pairs, full vs metabolism-free module score --\n")
print(pairsA[1:15, .(metabolic_set=substr(metabolic_set,1,52), module=sub("FFLCLASS\\|","",module),
      shared=n_shared_genes, rho_full=round(rho_full,3), rho_free=round(rho_metabfree,3),
      FDR_free=signif(FDR_metabfree,3))])
cat("\nmedian |rho| full:", round(median(abs(pairsA$rho_full)),3),
    " metabolism-free:", round(median(abs(pairsA$rho_metabfree)),3),
    "; Pearson(rho_full, rho_free) =", round(cor(pairsA$rho_full, pairsA$rho_metabfree),3), "\n")

## ---- (b) per-pair leave-shared-out on the metabolic side -----------------
top <- pairsA[order(-abs(rho_full))][1:60]
lso <- Map(function(a,m) setdiff(sets[[a]], sets[[m]]), top$metabolic_set, top$module)
names(lso) <- paste0(top$metabolic_set, "@@", top$module)
lso <- lso[lengths(lso) >= 5]
ssl <- gsva(ssgseaParam(Ef, lso, alpha=0.25, normalize=TRUE, minSize=5), BPPARAM=bp, verbose=FALSE)
pairsB <- rbindlist(lapply(rownames(ssl), function(k) {
  a <- sub("@@.*","",k); m <- sub(".*@@","",k)
  r <- sp(ssl[k,tum], ss[m,tum]); rf <- sp(ss[a,tum], ss[m,tum])
  data.table(metabolic_set=a, module=m, n_genes_after_removal=length(lso[[k]]),
             n_shared_removed=length(sets[[a]])-length(lso[[k]]),
             rho_full=rf["rho"], rho_leave_shared_out=r["rho"], p_lso=r["p"])
}))
pairsB[, FDR_lso := p.adjust(p_lso,"BH")]
pairsB <- pairsB[order(-abs(rho_full))]
fwrite(pairsB, file.path(RES,"metabolic_module_leave_shared_out.csv"))
cat("\n-- (b) leave-shared-out on the metabolic side (top 15) --\n")
print(pairsB[1:15, .(metabolic_set=substr(metabolic_set,1,50), module=sub("FFLCLASS\\|","",module),
      removed=n_shared_removed, rho_full=round(rho_full,3), rho_lso=round(rho_leave_shared_out,3),
      FDR=signif(FDR_lso,3))])
cat("DONE 25\n")
