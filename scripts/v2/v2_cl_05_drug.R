#!/usr/bin/env Rscript
## ============================================================================
## STROMA-FREE TEST -- part E (EXPLORATORY): does a collagen / miR-29-target
## module score predict drug sensitivity in breast cancer cell lines?
## GDSC2 (local GDSC2_dose_response.xlsx) primary; PRISM secondary screen 2nd.
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(readxl)})
options(warn=1); set.seed(20260908)
REV <- "/path/to/revision"
DM  <- file.path(REV,"data/depmap24q4"); OUT <- file.path(REV,"results/v2")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

W <- readRDS(file.path(DM,"cl_workspace.rds")); E<-W$E
mod <- fread(file.path(DM,"Model.csv"))
ed  <- fread(file.path(REV,"data/canonical_edges.tsv"))

## ---- module scores ---------------------------------------------------------
ECM <- intersect(c("COL1A1","COL3A1","COL1A2","COL5A1","COL5A2","COL6A3","COL11A1",
                   "FN1","SPARC","LOX"), colnames(E))
m29 <- unique(ed[edge_type=="miRNA_target" & source %in% c("hsa-miR-29a","hsa-miR-29b",
                                                           "hsa-miR-29c"), target])
m29 <- intersect(m29, colnames(E))
msg("ECM module genes:", length(ECM), "| miR-29 target module genes:", length(m29))
ids <- W$breast_ca
sc  <- data.table(ModelID=ids,
                  collagen_score = rowMeans(scale(E[ids,ECM,drop=FALSE])),
                  mir29_target_score = rowMeans(scale(E[ids,m29,drop=FALSE])))
sc[, COL1A1 := E[ids,"COL1A1"]][, COL3A1 := E[ids,"COL3A1"]]
fwrite(sc, file.path(OUT,"celllines_module_scores.csv"))
msg("collagen vs miR-29-target module score spearman:",
    round(cor(sc$collagen_score, sc$mir29_target_score, method="spearman"),3))

## ---- GDSC2 -----------------------------------------------------------------
g <- as.data.table(read_excel("/path/to/home/Desktop/DD/R_GPR/data/depmap/GDSC2_dose_response.xlsx"))
msg("GDSC2 rows:", nrow(g), "| lines:", uniqueN(g$SANGER_MODEL_ID), "| drugs:", uniqueN(g$DRUG_ID))
map <- mod[SangerModelID!="", .(SANGER_MODEL_ID=SangerModelID, ModelID)]
g <- merge(g, map, by="SANGER_MODEL_ID")
gb <- g[ModelID %in% ids]
msg("GDSC2 breast cancer lines matched to expression:", uniqueN(gb$ModelID),
    "| drug-line measurements:", nrow(gb), "| drugs:", uniqueN(gb$DRUG_ID))
gb <- merge(gb, sc, by="ModelID")

runsc <- function(dt, scorecol, valcol, minn=15){
  res <- dt[, {
    x <- get(scorecol); y <- get(valcol)
    ok <- is.finite(x)&is.finite(y)
    if(sum(ok)<minn) .(n=sum(ok), rho=NA_real_, p=NA_real_) else {
      ct <- suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE))
      .(n=sum(ok), rho=unname(ct$estimate), p=ct$p.value) }
  }, by=.(DRUG_ID, DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME)]
  res <- res[is.finite(rho)]
  res[, q := p.adjust(p,"BH")][order(p)]
}
r_col <- runsc(gb,"collagen_score","LN_IC50");   r_col[, score := "collagen"]
r_m29 <- runsc(gb,"mir29_target_score","LN_IC50"); r_m29[, score := "miR29_target"]
r_c1  <- runsc(gb,"COL1A1","LN_IC50");           r_c1[,  score := "COL1A1"]
allr <- rbind(r_col,r_m29,r_c1)
fwrite(allr, file.path(OUT,"celllines_gdsc2_module_drug_associations.csv"))
msg("GDSC2 drugs tested per score:", nrow(r_col))
for(s in unique(allr$score)) msg(" ", s, ": FDR<0.05 =", sum(allr[score==s]$q<0.05),
                                 "| FDR<0.10 =", sum(allr[score==s]$q<0.10))
print(allr[q<0.10][order(score,p)][,.(score,DRUG_NAME,PUTATIVE_TARGET,PATHWAY_NAME,n,rho,p,q)])

## pathway-level: is any GDSC pathway class enriched among collagen-associated drugs?
pw <- r_col[, .(n_drugs=.N, median_rho=median(rho),
                p=if(.N>=5) wilcox.test(rho)$p.value else NA_real_), by=PATHWAY_NAME]
pw[, q := p.adjust(p,"BH")][order(p)]
fwrite(pw, file.path(OUT,"celllines_gdsc2_pathway_level.csv"))
print(head(pw,12))

## ---- PRISM secondary (fewer breast lines; reported for completeness) -------
pr <- fread("/path/to/home/Desktop/DD/R_GPR/data/depmap/prism_secondary_dose_response.csv",
            select=c("depmap_id","auc","name","moa","screen_id","broad_id"))
pr <- pr[depmap_id %in% ids & is.finite(auc)]
msg("PRISM secondary: breast lines matched:", uniqueN(pr$depmap_id),
    "| compounds:", uniqueN(pr$broad_id))
pr <- merge(pr, sc[,.(depmap_id=ModelID, collagen_score, mir29_target_score)], by="depmap_id")
pr2 <- pr[, .(auc=median(auc)), by=.(depmap_id, broad_id, name, moa, collagen_score, mir29_target_score)]
prr <- pr2[, { ok <- is.finite(collagen_score)&is.finite(auc)
  if(sum(ok)<15) .(n=sum(ok),rho=NA_real_,p=NA_real_) else {
    ct<-suppressWarnings(cor.test(collagen_score[ok],auc[ok],method="spearman",exact=FALSE))
    .(n=sum(ok),rho=unname(ct$estimate),p=ct$p.value)}}, by=.(broad_id,name,moa)]
prr <- prr[is.finite(rho)][, q := p.adjust(p,"BH")][order(p)]
fwrite(prr, file.path(OUT,"celllines_prism_collagen_drug_associations.csv"))
msg("PRISM compounds tested:", nrow(prr), "| FDR<0.05:", sum(prr$q<0.05),
    "| FDR<0.10:", sum(prr$q<0.10))
print(head(prr,10))
msg("done")
