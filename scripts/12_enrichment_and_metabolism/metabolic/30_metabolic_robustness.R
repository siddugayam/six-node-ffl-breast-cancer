## 30_metabolic_robustness.R ----------------------------------------------
## Three robustness checks on the headline metabolic result:
##  1. within-PAM50-subtype miR-130a metabolic correlations
##  2. the two cell lines the group will use, on the CCLE metabolic scores
##  3. is the miR-130a metabolic axis independent of miR-130a's own PFI effect?
suppressPackageStartupMessages({library(data.table); library(survival)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
set.seed(1)
ss <- readRDS(file.path(RES,"metabolic_ssgsea_scores.rds"))
ph <- readRDS(file.path(ROOT,"data/brca_pheno.rds")); ph <- ph[match(colnames(ss), ph$sample),]
MI <- readRDS(file.path(ROOT,"data/brca_mirna_expr.rds"))
cl <- fread(file.path(ROOT,"data/brca_clinicalMatrix.tsv"))
sv <- as.data.table(readRDS(file.path(ROOT,"data/brca_survival.rds")))
tum <- ph$sample[ph$sample_type=="Primary Tumor"]; smp <- intersect(tum, colnames(MI))
met <- grep("^(KEGG|KEGG_EXTRA|REACTOME|HALLMARK|HUMANGEM|CURATED)\\|", rownames(ss), value=TRUE)
sp <- function(x,y){ct<-suppressWarnings(cor.test(x,y,method="spearman",exact=FALSE));c(rho=unname(ct$estimate),p=ct$p.value)}

## ---- 1. within-subtype -------------------------------------------------
scol <- names(cl)[1]
cl[, samp := get(scol)]
pam <- setNames(cl$PAM50Call_RNAseq, cl$samp)[smp]
cat("PAM50 calls available for", sum(!is.na(pam) & pam!=""), "of", length(smp), "paired tumours\n")
print(table(pam, useNA="ifany"))
sub_res <- rbindlist(lapply(c("LumA","LumB","Basal","Her2"), function(g) {
  s <- smp[!is.na(pam) & pam==g]
  if (length(s) < 60) return(NULL)
  rbindlist(lapply(met, function(p) {
    r <- sp(ss[p,s], MI["hsa-miR-130a-3p",s])
    data.table(subtype=g, n=length(s), set_id=p, rho=r["rho"], p=r["p"])
  }))
}))
sub_res[, FDR := p.adjust(p,"BH"), by=subtype]
fwrite(sub_res, file.path(RES,"metabolic_mir130a_within_subtype.csv"))
key <- c("REACTOME|R-HSA-390918_PEROXISOMAL_LIPID_METABOLISM",
         "REACTOME|R-HSA-2022854_KERATAN_SULFATE_BIOSYNTHESIS",
         "HUMANGEM|KERATAN_SULFATE_BIOSYNTHESIS",
         "KEGG|hsa00601_GLYCOSPHINGOLIPID_BIOSYNTHESIS_LACTO_AND_NEOLACTO_SERIES",
         "HUMANGEM|BETA_OXIDATION_OF_EVEN_CHAIN_FATTY_ACIDS_MITOCHONDRIAL_",
         "REACTOME|R-HSA-8978868_FATTY_ACID_METABOLISM",
         "HALLMARK|HALLMARK_PEROXISOME","CURATED|CURATED_GLYCOLYSIS_CORE",
         "CURATED|CURATED_ANTIFOLATE_DETERMINANTS",
         "KEGG_EXTRA|hsa04933_AGE_RAGE_SIGNALING_PATHWAY_IN_DIABETIC_COMPLICATIONS")
key <- intersect(key, rownames(ss))
cat("\n== miR-130a vs key metabolic scores, all tumours and within PAM50 subtype ==\n")
allr <- rbindlist(lapply(key, function(p){r<-sp(ss[p,smp],MI["hsa-miR-130a-3p",smp])
        data.table(subtype="ALL", n=length(smp), set_id=p, rho=r["rho"], p=r["p"], FDR=NA_real_)}))
w <- dcast(rbind(allr, sub_res[set_id %in% key]), set_id ~ subtype, value.var="rho")
print(w[, lapply(.SD, function(z) if(is.numeric(z)) round(z,3) else substr(z,1,52))])
wp <- dcast(rbind(allr, sub_res[set_id %in% key]), set_id ~ subtype, value.var="p")
fwrite(w, file.path(RES,"metabolic_mir130a_within_subtype_key.csv"))
cat("\np-values:\n"); print(wp[, lapply(.SD, function(z) if(is.numeric(z)) signif(z,2) else substr(z,1,52))])

## ---- 2. the two planned cell lines -------------------------------------
cc <- readRDS(file.path(RES,"metabolic_ccle_ssgsea_scores.rds"))
cat("\n== CCLE metabolic scores: MDA-MB-231 (miR-130a 203.3) vs MCF7 (miR-130a 5.4) ==\n")
cat("(prediction from the TCGA axis: the miR-130a-HIGH line should score HIGHER on\n",
    " keratan-sulfate/glycan biosynthesis and LOWER on peroxisomal/fatty-acid oxidation)\n")
lines2 <- c("MDAMB231_BREAST","MCF7_BREAST")
stopifnot(all(lines2 %in% colnames(cc)))
tab2 <- rbindlist(lapply(intersect(key, rownames(cc)), function(p) {
  data.table(set_id=p, MDAMB231=cc[p,"MDAMB231_BREAST"], MCF7=cc[p,"MCF7_BREAST"],
    pct_MDAMB231=round(100*mean(cc[p,]<=cc[p,"MDAMB231_BREAST"]),1),
    pct_MCF7=round(100*mean(cc[p,]<=cc[p,"MCF7_BREAST"]),1),
    diff_231_minus_MCF7=cc[p,"MDAMB231_BREAST"]-cc[p,"MCF7_BREAST"])
}))
print(tab2[, .(set_id=substr(set_id,1,55), pct_MDAMB231, pct_MCF7,
               diff=round(diff_231_minus_MCF7,4))])
fwrite(tab2, file.path(RES,"metabolic_ccle_two_lines.csv"))

## ---- 3. miR-130a's own survival effect vs the metabolic axis ------------
sv2 <- sv[sample %in% smp]
sv2[, stage := fcase(grepl("Stage IV", ajcc_pathologic_tumor_stage),"IV",
                     grepl("Stage III",ajcc_pathologic_tumor_stage),"III",
                     grepl("Stage II", ajcc_pathologic_tumor_stage),"II",
                     grepl("Stage I($| |A|B|C)",ajcc_pathologic_tumor_stage),"I", default=NA_character_)]
sv2[, stage := factor(stage, levels=c("I","II","III","IV"))]
sv2 <- sv2[!is.na(stage) & !is.na(age_at_initial_pathologic_diagnosis)]
s <- sv2$sample
m130 <- as.numeric(scale(MI["hsa-miR-130a-3p", s]))
ks   <- as.numeric(scale(ss["HUMANGEM|KERATAN_SULFATE_BIOSYNTHESIS", s]))
fao  <- as.numeric(scale(ss["HUMANGEM|BETA_OXIDATION_OF_EVEN_CHAIN_FATTY_ACIDS_MITOCHONDRIAL_", s]))
cat("\n== Cox: miR-130a and the metabolic axis it tracks, same", nrow(sv2), "tumours ==\n")
for (ep in c("OS","PFI","DSS")) {
  ev <- sv2[[ep]]; ti <- sv2[[paste0(ep,".time")]]
  ok <- !is.na(ev)&!is.na(ti)&ti>0
  y <- Surv(ti[ok], ev[ok]); a <- sv2$age_at_initial_pathologic_diagnosis[ok]; g <- sv2$stage[ok]
  f1 <- coxph(y ~ m130[ok] + a + g)
  f2 <- coxph(y ~ ks[ok] + a + g)
  f3 <- coxph(y ~ m130[ok] + ks[ok] + fao[ok] + a + g)
  cf1 <- summary(f1)$coefficients[1,]; cf2 <- summary(f2)$coefficients[1,]
  cf3 <- summary(f3)$coefficients
  cat(sprintf("%s (events %d): miR-130a alone HR %.3f p %.3g | keratan-sulfate alone HR %.3f p %.3g | joint: miR-130a HR %.3f p %.3g, KS HR %.3f p %.3g, FAO HR %.3f p %.3g\n",
      ep, sum(ev[ok]==1), exp(cf1[1]), cf1[5], exp(cf2[1]), cf2[5],
      exp(cf3[1,1]), cf3[1,5], exp(cf3[2,1]), cf3[2,5], exp(cf3[3,1]), cf3[3,5]))
}
cat("DONE 30\n")
