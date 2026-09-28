#!/usr/bin/env Rscript
## PART 1C) Does a miR-29 target score, a miR-130a target score or a collagen score
## predict sensitivity to any compound class in BREAST cell lines?
## Two independent resources: PRISM Repurposing 19Q4 secondary screen (AUC) and
## GDSC2 (LN_IC50).  Spearman association per compound, BH-corrected within resource.
suppressPackageStartupMessages({library(data.table)})
set.seed(20260909)
REV <- "/path/to/revision"
OUT <- file.path(REV,"results/v3"); CA <- file.path(REV,"cache/v7")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

S <- fread(file.path(OUT,"screens_celline_module_scores.csv"))
scores <- c("COLLAGEN","MIR29_TARGET","MIR130A_TARGET","MIR29_STRONG","MIR130A_STRONG")
msg("breast lines with module scores:", nrow(S))

## ---------------- PRISM ----------------
P <- fread(file.path(CA,"prism_secondary_dose_response.csv"),
           select=c("broad_id","depmap_id","ccle_name","screen_id","auc","ic50","name","moa","target","phase"))
msg("PRISM secondary rows:", nrow(P), "| compounds:", uniqueN(P$broad_id), "| lines:", uniqueN(P$depmap_id))
Pb <- P[depmap_id %in% S$ModelID]
msg("PRISM rows in breast lines with a module score:", nrow(Pb),
    "| breast lines:", uniqueN(Pb$depmap_id), "| compounds:", uniqueN(Pb$broad_id))
Pb <- Pb[is.finite(auc)]
## one AUC per compound x line (median over screen_ids / replicate entries)
Pb <- Pb[, .(auc=median(auc), name=name[1], moa=moa[1], target=target[1], phase=phase[1]),
         by=.(broad_id, depmap_id)]
cnt <- Pb[, .N, by=broad_id]
keepc <- cnt[N>=15, broad_id]
msg("compounds profiled in >=15 breast lines:", length(keepc), "of", nrow(cnt))
Pb <- Pb[broad_id %in% keepc]
W <- dcast(Pb, broad_id ~ depmap_id, value.var="auc")
meta <- unique(Pb[, .(broad_id, name, moa, target, phase)])

run_assoc <- function(mat_dt, id_col, sc, S, resource, value_name) {
  ids <- setdiff(names(mat_dt), id_col)
  ids <- intersect(ids, S$ModelID)
  x <- setNames(S[[sc]], S$ModelID)[ids]
  M <- as.matrix(mat_dt[, ..ids])
  res <- rbindlist(lapply(seq_len(nrow(M)), function(i){
    y <- M[i,]; ok <- is.finite(y) & is.finite(x)
    if (sum(ok) < 15) return(NULL)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
    data.table(resource=resource, score=sc, id=mat_dt[[id_col]][i], n_lines=sum(ok),
               rho=unname(ct$estimate), p=ct$p.value)
  }))
  res
}
PR <- rbindlist(lapply(scores, function(sc) run_assoc(W, "broad_id", sc, S, "PRISM_secondary_AUC", "auc")))
PR[, fdr := p.adjust(p, "BH"), by=score]
PR <- merge(PR, meta, by.x="id", by.y="broad_id", all.x=TRUE)
setorder(PR, score, p)
fwrite(PR, file.path(OUT,"screens_prism_associations.csv"))
msg("PRISM associations tested:", nrow(PR), "| by score:"); print(PR[, .(n_compounds=.N, n_FDR05=sum(fdr<0.05), min_fdr=min(fdr)), by=score])
msg("PRISM: compounds surviving BH FDR<0.05 (negative rho = higher score -> lower AUC = MORE sensitive):")
print(PR[fdr<0.05][order(score, p)][, .(score, name, moa, phase, n_lines, rho, p, fdr)], nrows=60)

## ---------------- GDSC2 ----------------
suppressPackageStartupMessages(library(readxl))
gf <- file.path(CA,"GDSC2_fitted_dose_response.xlsx")
if (file.exists(gf) && file.size(gf) > 5e6) {
  G <- as.data.table(read_excel(gf))
  msg("GDSC2 rows:", nrow(G), "| columns:", paste(head(names(G),20), collapse=","))
  cmp <- fread(file.path(CA,"cmp_model_list.csv"))
  mod <- fread(file.path(REV,"data/depmap/Model_24Q4.csv"))
  ## map GDSC SANGER_MODEL_ID -> DepMap ModelID via Model_24Q4 SangerModelID
  map <- mod[SangerModelID != "" & !is.na(SangerModelID), .(SangerModelID, ModelID)]
  G[, SangerModelID := SANGER_MODEL_ID]
  G <- merge(G, map, by="SangerModelID")
  Gb <- G[ModelID %in% S$ModelID]
  msg("GDSC2 rows mapped into breast lines with a module score:", nrow(Gb),
      "| lines:", uniqueN(Gb$ModelID), "| drugs:", uniqueN(Gb$DRUG_ID))
  Gb <- Gb[is.finite(LN_IC50)]
  Gb <- Gb[, .(LN_IC50=median(LN_IC50), DRUG_NAME=DRUG_NAME[1], PUTATIVE_TARGET=PUTATIVE_TARGET[1],
               PATHWAY_NAME=PATHWAY_NAME[1]), by=.(DRUG_ID, ModelID)]
  cnt2 <- Gb[, .N, by=DRUG_ID]; keepd <- cnt2[N>=15, DRUG_ID]
  msg("GDSC2 drugs profiled in >=15 breast lines:", length(keepd), "of", nrow(cnt2))
  Gb <- Gb[DRUG_ID %in% keepd]
  WG <- dcast(Gb, DRUG_ID ~ ModelID, value.var="LN_IC50")
  metaG <- unique(Gb[, .(DRUG_ID, DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME)])
  GR <- rbindlist(lapply(scores, function(sc) run_assoc(WG, "DRUG_ID", sc, S, "GDSC2_LN_IC50", "ln_ic50")))
  GR[, fdr := p.adjust(p, "BH"), by=score]
  GR <- merge(GR, metaG, by.x="id", by.y="DRUG_ID", all.x=TRUE)
  setorder(GR, score, p)
  fwrite(GR, file.path(OUT,"screens_gdsc2_associations.csv"))
  msg("GDSC2 associations tested:", nrow(GR), "| by score:")
  print(GR[, .(n_drugs=.N, n_FDR05=sum(fdr<0.05), min_fdr=min(fdr)), by=score])
  msg("GDSC2: drugs surviving BH FDR<0.05 (negative rho = higher score -> lower IC50 = MORE sensitive):")
  print(GR[fdr<0.05][order(score, p)][, .(score, DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME, n_lines, rho, p, fdr)], nrows=80)
  ## pathway-level: is any MOA/pathway enriched among the strongest associations?
  GR[, dir := fifelse(rho<0, "sensitising","resistance")]
  pw <- GR[score=="COLLAGEN", .(n=.N, med_rho=median(rho), n_p05=sum(p<0.05)), by=PATHWAY_NAME][order(-n_p05)]
  msg("GDSC2 COLLAGEN score: nominal (p<0.05) hits by pathway"); print(head(pw, 20))
  fwrite(GR[, .(resource, score, id, DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME, n_lines, rho, p, fdr)],
         file.path(OUT,"screens_gdsc2_associations.csv"))
} else {
  msg("GDSC2 file absent or incomplete -- GDSC arm NOT run")
}
msg("DONE 22")
