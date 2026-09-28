#!/usr/bin/env Rscript
## Independent re-derivation of PART 2J: within-CAF TF -> COL1A1, with and without
## adjustment for study and CAF-subtype (myCAF-like) composition.
suppressPackageStartupMessages({library(data.table)})
REV <- "/path/to/revision"; OUT <- file.path(REV,"results/v3")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")
D <- fread(file.path(OUT,"sc_hbca_donor_compartment_pseudobulk.csv"))
FR <- fread(file.path(OUT,"sc_hbca_donor_caf_subtype_fractions.csv"))
setnames(FR, 1, "donor")
msg("donor pseudobulk rows:", nrow(D), "| compartments:", paste(unique(D$compartment), collapse=","))
C <- D[compartment=="CAF" & n_cells>=50]
msg("CAF donors with >=50 cells:", nrow(C))
M <- merge(C, FR, by="donor")
msg("merged with subtype fractions:", nrow(M))
myc <- intersect(c("CAFs (LAMP5)","CAFs (COL11A1+)"), names(M))
msg("myCAF-like subsets used:", paste(myc, collapse=" + "))
M[, myCAF_frac := rowSums(.SD), .SDcols=myc]
lg <- function(v) log2(v+1)
y <- lg(M$COL1A1)
ct <- cor.test(M$myCAF_frac, y, method="spearman", exact=FALSE)
msg("myCAF-like fraction vs CAF COL1A1: rho =", round(ct$estimate,3), "p =", signif(ct$p.value,3))
st <- factor(M$study)
TFS <- c("NFKB1","RELA","SP1","ETS1","MYC","STAT3","HIF1A","TP53","MKL1","SRF","RUNX1","RUNX2",
         "TWIST1","SNAI2","ZEB1","JUN","FOS","EGR1","CEBPB","YAP1","WWTR1","SMAD3","SMAD4","TGFB1","TGFBR2")
TFS <- intersect(TFS, names(M))
R <- rbindlist(lapply(TFS, function(tf){
  x <- lg(M[[tf]])
  r0 <- cor.test(x, y, method="spearman", exact=FALSE)
  rx <- residuals(lm(rank(x) ~ st)); ry <- residuals(lm(rank(y) ~ st))
  r2 <- cor.test(rx, ry)
  rx1 <- residuals(lm(rank(x) ~ st + rank(M$myCAF_frac)))
  ry1 <- residuals(lm(rank(y) ~ st + rank(M$myCAF_frac)))
  r1 <- cor.test(rx1, ry1)
  data.table(TF=tf, n=nrow(M), rho_raw=unname(r0$estimate), p_raw=r0$p.value,
             rho_study=unname(r2$estimate), p_study=r2$p.value,
             rho_full=unname(r1$estimate), p_full=r1$p.value)
}))
R[, fdr_full := p.adjust(p_full, "BH")]
setorder(R, -rho_full)
print(R, digits=3, nrows=40)
fwrite(R, file.path(OUT,"sc_verify_within_caf_tf_col1a1.csv"))
dep <- fread(file.path(OUT,"sc_tf_collagen_within_caf_composition_adjusted.csv"))
m <- merge(R, dep, by="TF")
msg("max |d rho_raw| =", signif(max(abs(m$rho_raw - m$rho_raw.y)),3),
    "| max |d rho_study| =", signif(max(abs(m$rho_study - m$rho_study_adj)),3),
    "| max |d rho_full| =", signif(max(abs(m$rho_full - m$rho_study_and_myCAFfrac_adj)),3))
msg("TFs positive at BH<0.05 after full adjustment:", paste(R[rho_full>0 & fdr_full<0.05, TF], collapse=", "))
