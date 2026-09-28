## 29_metabolic_summary.R --------------------------------------------------
## Consolidation: composition-adjusted miR-130a metabolic associations,
## CCLE replication of the TCGA top hits, and the headline summary tables.
suppressPackageStartupMessages({library(data.table)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
set.seed(1)
ss  <- readRDS(file.path(RES,"metabolic_ssgsea_scores.rds"))
ph  <- readRDS(file.path(ROOT,"data/brca_pheno.rds")); ph <- ph[match(colnames(ss), ph$sample),]
MI  <- readRDS(file.path(ROOT,"data/brca_mirna_expr.rds"))
tum <- ph$sample[ph$sample_type=="Primary Tumor"]; smp <- intersect(tum, colnames(MI))
cors<- fread(file.path(RES,"metabolic_feature_correlations.csv"))
tvn <- fread(file.path(RES,"metabolic_tumour_vs_normal.csv"))
cox <- fread(file.path(RES,"metabolic_survival_cox.csv"))
ccle<- fread(file.path(RES,"metabolic_ccle_mir130a_correlations.csv"))
inv <- fread(file.path(RES,"metabolic_geneset_inventory.csv"))

## ---- A. composition-adjusted miR-130a associations ----------------------
ctx <- ss[c("CONTROL|CONTROL_EPITHELIAL","CONTROL|CONTROL_CAF_FULL",
            "CONTROL|CONTROL_IMMUNE_FULL","CONTROL|CONTROL_PROLIFERATION"), smp]
met <- grep("^(KEGG|KEGG_EXTRA|REACTOME|HALLMARK|HUMANGEM|CURATED)\\|", rownames(ss), value=TRUE)
x <- rank(MI["hsa-miR-130a-3p", smp])
Zc <- apply(ctx, 1, rank)   # n_samples x 4 covariates
rx <- resid(lm(x ~ Zc))
adj <- rbindlist(lapply(met, function(p) {
  ry <- resid(lm(rank(ss[p,smp]) ~ Zc))
  ct <- cor.test(rx, ry)
  raw <- suppressWarnings(cor.test(ss[p,smp], MI["hsa-miR-130a-3p",smp], method="spearman", exact=FALSE))
  data.table(set_id=p, rho_raw=unname(raw$estimate), p_raw=raw$p.value,
             r_partial=unname(ct$estimate), p_partial=ct$p.value)
}))
adj[, `:=`(FDR_raw=p.adjust(p_raw,"BH"), FDR_partial=p.adjust(p_partial,"BH"))]
adj <- merge(adj, inv[,.(set_id, collection, n_genes_in_universe)], by="set_id")
adj <- adj[order(-abs(r_partial))]
fwrite(adj, file.path(RES,"metabolic_mir130a_composition_adjusted.csv"))
cat("== miR-130a metabolic associations, adjusted for epithelial + CAF + immune + proliferation ==\n")
cat("significant raw:", adj[FDR_raw<0.05,.N], " partial:", adj[FDR_partial<0.05,.N], "of", nrow(adj), "\n")
print(adj[1:15, .(set_id=substr(set_id,1,58), rho_raw=round(rho_raw,3),
    r_partial=round(r_partial,3), FDR_partial=signif(FDR_partial,3))])
cat("\n-- most negative after adjustment --\n")
print(adj[order(r_partial)][1:10, .(set_id=substr(set_id,1,58), rho_raw=round(rho_raw,3),
    r_partial=round(r_partial,3), FDR_partial=signif(FDR_partial,3))])

## ---- B. does the TCGA miR-130a signal replicate in stroma-free CCLE? ----
tc <- cors[score=="ssGSEA" & feature=="hsa-miR-130a-3p", .(set_id, rho_TCGA=rho, FDR_TCGA=FDR)]
rep <- merge(tc, ccle[, .(set_id, rho_CCLE, p_CCLE, FDR_CCLE)], by="set_id")
rep <- merge(rep, adj[, .(set_id, r_partial, FDR_partial)], by="set_id", all.x=TRUE)
top <- rep[FDR_TCGA<0.05][order(-abs(rho_TCGA))][1:25]
top[, replicates_sign := sign(rho_CCLE)==sign(rho_TCGA)]
top[, nominal_CCLE := p_CCLE < 0.05]
fwrite(rep, file.path(RES,"metabolic_mir130a_tcga_vs_ccle.csv"))
cat("\n== top 25 TCGA miR-130a metabolic hits, checked in 50 stroma-free CCLE breast lines ==\n")
print(top[, .(set_id=substr(set_id,1,52), rho_TCGA=round(rho_TCGA,3), r_partial=round(r_partial,3),
   rho_CCLE=round(rho_CCLE,3), p_CCLE=signif(p_CCLE,2), same_sign=replicates_sign)])
cat("\nsame sign:", sum(top$replicates_sign), "/25 ; nominally significant in CCLE:",
    sum(top$nominal_CCLE), "/25\n")

## ---- C. the three claims, one row each ---------------------------------
g <- function(id, ep="OS") {
  a <- tvn[set_id==id & score=="ssGSEA"]
  b <- cors[score=="ssGSEA" & feature=="hsa-miR-130a-3p" & set_id==id]
  c_<- cox[score=="ssGSEA" & endpoint==ep & set_id==id]
  data.table(set_id=id, n=a$n_genes_in_universe, d_TvN=a$cohens_d, FDR_TvN=a$FDR_limma,
             rho_mir130a=b$rho, FDR_mir130a=b$FDR,
             endpoint=ep, HR=c_$HR_per_SD, p_cox=c_$p, FDR_cox=c_$FDR)
}
claim_ids <- c("HALLMARK|HALLMARK_GLYCOLYSIS", grep("hsa00010", rownames(ss), value=TRUE),
  "CURATED|CURATED_GLYCOLYSIS_CORE","CURATED|CURATED_ANTIFOLATE_DETERMINANTS",
  "CURATED|CURATED_ONE_CARBON_FOLATE", grep("hsa01523|hsa00670", rownames(ss), value=TRUE),
  grep("hsa04933", rownames(ss), value=TRUE))
claims <- rbindlist(lapply(claim_ids, function(i) rbindlist(lapply(c("OS","PFI","DSS"), function(e) g(i,e)))))
fwrite(claims, file.path(RES,"metabolic_three_claims_summary.csv"))
cat("\n== the three claims: score-level summary ==\n")
print(claims[endpoint=="OS", .(set_id=substr(set_id,1,50), n, d_TvN=round(d_TvN,3),
   FDR_TvN=signif(FDR_TvN,2), rho_130a=round(rho_mir130a,3), FDR_130a=signif(FDR_mir130a,2))])
cat("\nsurvival for the same sets (age+stage adjusted):\n")
print(claims[, .(set_id=substr(set_id,1,42), endpoint, HR=round(HR,3), p=signif(p_cox,2), FDR=signif(FDR_cox,2))])

## ---- D. headline ranked tables -----------------------------------------
t1 <- tvn[score=="ssGSEA"][order(-cohens_d)][, rank_up := .I]
fwrite(t1[, .(rank_up, set_id, collection, n=n_genes_in_universe, cohens_d, t, p_limma, FDR_limma,
              auc_tumour_higher, mean_tumour, mean_normal)],
       file.path(RES,"metabolic_tumour_vs_normal_RANKED.csv"))
cat("\n== TOP 20 metabolic programmes UP in breast tumour ==\n")
print(t1[1:20, .(rank_up, set_id=substr(set_id,1,60), d=round(cohens_d,2), FDR=signif(FDR_limma,2))])
cat("\n== TOP 20 metabolic programmes DOWN in breast tumour ==\n")
print(t1[order(cohens_d)][1:20, .(set_id=substr(set_id,1,60), d=round(cohens_d,2), FDR=signif(FDR_limma,2))])

pro <- cox[score=="ssGSEA" & FDR<0.05][order(endpoint, p)]
fwrite(pro, file.path(RES,"metabolic_survival_significant.csv"))
cat("\n== prognostic metabolic programmes (FDR<0.05, age+stage adjusted), counts ==\n")
print(cox[score=="ssGSEA", .(sig=sum(FDR<0.05), adverse=sum(FDR<0.05 & HR_per_SD>1),
      protective=sum(FDR<0.05 & HR_per_SD<1)), by=endpoint])
cat("\nsets significant for ALL THREE endpoints:\n")
allthree <- cox[score=="ssGSEA" & FDR<0.05, .N, by=set_id][N==3]
p3 <- cox[score=="ssGSEA" & set_id %in% allthree$set_id]
print(dcast(p3, set_id ~ endpoint, value.var="HR_per_SD")[order(OS)][
  , .(set_id=substr(set_id,1,58), OS=round(OS,2), PFI=round(PFI,2), DSS=round(DSS,2))])
cat("\nDONE 29\n")
