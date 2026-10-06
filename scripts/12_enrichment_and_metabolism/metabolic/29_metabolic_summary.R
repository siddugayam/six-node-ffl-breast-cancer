## 29_metabolic_summary.R --------------------------------------------------
## Consolidation: the three metabolic claims, one row each, and the headline
## summary tables.
suppressPackageStartupMessages({library(data.table)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
set.seed(1)
ss  <- readRDS(file.path(RES,"metabolic_ssgsea_scores.rds"))
tvn <- fread(file.path(RES,"metabolic_tumour_vs_normal.csv"))
cox <- fread(file.path(RES,"metabolic_survival_cox.csv"))

## ---- A. the three claims, one row each ---------------------------------
g <- function(id, ep="OS") {
  a <- tvn[set_id==id & score=="ssGSEA"]
  c_<- cox[score=="ssGSEA" & endpoint==ep & set_id==id]
  data.table(set_id=id, n=a$n_genes_in_universe, d_TvN=a$cohens_d, FDR_TvN=a$FDR_limma,
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
   FDR_TvN=signif(FDR_TvN,2))])
cat("\nsurvival for the same sets (age+stage adjusted):\n")
print(claims[, .(set_id=substr(set_id,1,42), endpoint, HR=round(HR,3), p=signif(p_cox,2), FDR=signif(FDR_cox,2))])

## ---- B. headline ranked tables -----------------------------------------
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
