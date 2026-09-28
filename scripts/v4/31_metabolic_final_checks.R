## 31_metabolic_final_checks.R --------------------------------------------
## Reconciliation and final headline numbers.
suppressPackageStartupMessages({library(data.table); library(survival)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
set.seed(1)
ss <- readRDS(file.path(RES,"metabolic_ssgsea_scores.rds"))
ph <- readRDS(file.path(ROOT,"data/brca_pheno.rds")); ph <- ph[match(colnames(ss), ph$sample),]
MI <- readRDS(file.path(ROOT,"data/brca_mirna_expr.rds"))
sv <- as.data.table(readRDS(file.path(ROOT,"data/brca_survival.rds")))
tum <- ph$sample[ph$sample_type=="Primary Tumor"]; smp <- intersect(tum, colnames(MI))

## 1. reconcile the previously reported miR-130a PFI HR (0.840 per SD, p=0.048)
s <- intersect(sv$sample, smp); sv2 <- sv[match(s, sample)]
x <- as.numeric(scale(MI["hsa-miR-130a-3p", s]))
out <- rbindlist(lapply(c("OS","PFI","DSS"), function(ep) {
  ev <- sv2[[ep]]; ti <- sv2[[paste0(ep,".time")]]; ok <- !is.na(ev)&!is.na(ti)&ti>0
  f <- coxph(Surv(ti[ok],ev[ok]) ~ x[ok]); cf <- summary(f)$coefficients[1,]
  data.table(endpoint=ep, model="unadjusted", n=sum(ok), events=sum(ev[ok]==1),
             HR=exp(cf[1]), p=cf[5])
}))
cat("== miR-130a-3p per-SD Cox, UNADJUSTED (reconciles the previously reported HR 0.840 / p 0.048) ==\n")
print(out[, .(endpoint, n, events, HR=round(HR,3), p=signif(p,3))])
fwrite(out, file.path(RES,"metabolic_mir130a_survival_reconciliation.csv"))

## 2. headline counts for parts B and C
tvn <- fread(file.path(RES,"metabolic_tumour_vs_normal.csv"))[score=="ssGSEA"]
cors<- fread(file.path(RES,"metabolic_feature_correlations.csv"))[score=="ssGSEA"]
cat("\n== headline counts ==\n")
cat("metabolic score sets:", nrow(tvn), " significant tumour-vs-normal (BH<0.05):",
    sum(tvn$FDR_limma<0.05), " up:", sum(tvn$FDR_limma<0.05 & tvn$cohens_d>0),
    " down:", sum(tvn$FDR_limma<0.05 & tvn$cohens_d<0), "\n")
cat("pathway x feature pairs tested:", nrow(cors), " significant (BH<0.05 across all pairs):",
    sum(cors$FDR<0.05), "\n")
for (f in c("hsa-miR-130a-3p","FAMILY_miR130_301_seed","hsa-miR-29a-3p","hsa-miR-29b-3p",
            "hsa-miR-29c-3p","FAMILY_miR29_3p","hsa-let-7b-5p","FAMILY_let7_5p")) {
  d <- cors[feature==f & !grepl("^FFLCLASS|^CONTROL", set_id)]
  cat(sprintf("%-24s n=%d  FDR<0.05: %d   strongest: %s rho=%.3f\n", f, nrow(d), sum(d$FDR<0.05),
      substr(d[which.max(abs(rho)), set_id],1,50), d[which.max(abs(rho)), rho]))
}

## 3. how much of the AGE-RAGE score is explained by the CAF score?
tumall <- tum
m <- lm(ss["KEGG_EXTRA|hsa04933_AGE_RAGE_SIGNALING_PATHWAY_IN_DIABETIC_COMPLICATIONS", tumall] ~
        ss["CONTROL|CONTROL_CAF_FULL", tumall])
cat("\nAGE-RAGE score ~ CAF score (tumours only): R^2 =", round(summary(m)$r.squared,3), "\n")
m2 <- lm(ss["KEGG_EXTRA|hsa04933_AGE_RAGE_SIGNALING_PATHWAY_IN_DIABETIC_COMPLICATIONS", tumall] ~
         ss["CONTROL|CONTROL_CAF_FULL", tumall] + ss["CONTROL|CONTROL_IMMUNE_FULL", tumall] +
         ss["CONTROL|CONTROL_EPITHELIAL", tumall] + ss["CONTROL|CONTROL_PROLIFERATION", tumall])
cat("AGE-RAGE ~ CAF + immune + epithelial + proliferation: R^2 =", round(summary(m2)$r.squared,3), "\n")

## 4. output inventory
f <- list.files(RES, pattern="^metabolic_", full.names=TRUE)
inv <- data.table(file=basename(f), bytes=file.size(f))
inv[grepl("\\.csv$", file), rows := sapply(file, function(z)
     tryCatch(nrow(fread(file.path(RES,z), select=1L, showProgress=FALSE)), error=function(e) NA_integer_))]
setorder(inv, file)
fwrite(inv, file.path(RES,"metabolic_OUTPUT_INVENTORY.csv"))
cat("\n== outputs written (", nrow(inv), "files ) ==\n"); print(inv)
cat("DONE 31\n")
