## S8 (analyses/six_node_pattern): per-circuit clinical test on BHAT6 composite and MODEL6-architecture node sets.
## Repeats the comparison of every six-node circuit with its own three-node core (EA6: 1,989 circuits,
## scripts/10_survival_and_clinical/cox_models/13_circuit_level_3node_vs_6node.R) on the six-node objects of S0/S1.
##   circuits  the 2,498 BHAT6 composite node sets and the 71 MODEL6-architecture (FULL) node sets
##             (S1/obs listings).  For a node set with several role assignments the first one listed is
##             used.  Core = (TF1, miR1, G1); extension = (TF2, miR2, G2).  In BHAT6 the core's T1 -> G1
##             arc is not a compulsory edge (Bhat Fig. 1d), so the BHAT6 "core" is the role triple only.
##   data      TCGA-BRCA primary tumours with gene + miRNA expression and survival (as scripts 12 and 13):
##             data/brca_gene_expr.rds, brca_mirna_expr_canonical.rds, brca_pheno.rds,
##             brca_survival_clean.rds; endpoints OS and PFI
##   PRIMARY (the analysis plan's specification) GSVA scores (kcdf Gaussian, minSize 2, as script 12), z-scaled;
##             Cox models adjusted for age and AJCC stage (complete cases):
##             core model ~ core + age + stage; six-node model ~ six + age + stage;
##             nested model ~ core + extension + age + stage.
##             dC = C(six-node model) - C(core model) (Harrell); LRT core vs nested; BH across circuits.
##   SECONDARY (the stored EA6 method, script 13) mean-z scores, univariable Cox, the same dC and LRT.
## Output: s8_circuits.csv (one row per circuit x endpoint x method; summary rows only), s8_summary.csv.
suppressPackageStartupMessages({ library(data.table); library(survival); library(GSVA) })
REV <- "/path/to/revision"
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
S1 <- file.path(HERE, "..", "S1")
set.seed(1234)
names587 <- readLines(file.path(S1, "node_names.txt"))[1:587]
rd <- function(f, cols) { d <- fread(f, header = FALSE, col.names = cols); d[, lapply(.SD, function(i) names587[i + 1])] }
BH <- rd(file.path(S1, "obs/nolegacy_bhat6_composite_instances.tsv"), c("miR1", "miR2", "TF1", "TF2", "G1", "G2"))
MF <- rd(file.path(S1, "obs/nolegacy_model6_full_instances.tsv"), c("TF1", "TF2", "miR1", "miR2", "G1", "G2"))
first_per_set <- function(D, fam) {
  D[, key := apply(.SD, 1, function(r) paste(sort(r), collapse = ";")), .SDcols = c("TF1", "TF2", "miR1", "miR2", "G1", "G2")]
  D <- D[!duplicated(key)]; D[, family := fam]; D
}
CIR <- rbind(first_per_set(BH, "BHAT6_composite"), first_per_set(MF, "MODEL6_architecture"), use.names = TRUE)
cat("node sets: ", paste(names(table(CIR$family)), table(CIR$family), collapse = "; "), "\n")

gex <- readRDS(file.path(REV, "data/brca_gene_expr.rds")); mex <- readRDS(file.path(REV, "data/brca_mirna_expr_canonical.rds"))
pheno <- readRDS(file.path(REV, "data/brca_pheno.rds")); surv <- as.data.table(readRDS(file.path(REV, "data/brca_survival_clean.rds")))
tumour <- pheno$sample[pheno$sample_type == "Primary Tumor"]
common <- Reduce(intersect, list(colnames(gex), colnames(mex), tumour, surv$sample_id))
X <- rbind(gex[, common, drop = FALSE], mex[, common, drop = FALSE])
v <- apply(X, 1, function(z) stats::sd(z, na.rm = TRUE)); X <- X[is.finite(v) & v > 0, , drop = FALSE]
cat("samples:", length(common), "; features:", nrow(X), "\n")
cl <- surv[sample_id %in% common]; setkey(cl, sample_id); cl <- cl[common]; stopifnot(identical(cl$sample_id, common))

have <- rownames(X)
CIR[, covered := TF1 %in% have & TF2 %in% have & miR1 %in% have & miR2 %in% have & G1 %in% have & G2 %in% have]
cat("circuits with all six members measured:", sum(CIR$covered), "of", nrow(CIR), "\n")
CIR <- CIR[covered == TRUE]
CIR[, core := paste(TF1, miR1, G1, sep = ";")][, ext := paste(TF2, miR2, G2, sep = ";")]
setsl <- unique(c(CIR$key, CIR$core, CIR$ext))
gsl <- setNames(lapply(setsl, function(s) strsplit(s, ";")[[1]]), setsl)
cat("distinct gene sets to score:", length(gsl), "\n")
t0 <- Sys.time()
GS <- gsva(gsvaParam(X, gsl, kcdf = "Gaussian", minSize = 2, maxSize = Inf), verbose = FALSE)
cat("GSVA done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min; ", paste(dim(GS), collapse = " x "), "\n")
Z <- t(scale(t(X))); MZ <- t(sapply(gsl, function(s) colMeans(Z[s, , drop = FALSE])))
sc <- function(M, s) as.numeric(scale(M[s, ]))

res <- list(); k <- 0L
for (ep in c("OS", "PFI")) {
  tt <- cl[[paste0(ep, ".time")]]; ee <- cl[[ep]]
  okb <- is.finite(tt) & tt > 0 & is.finite(ee)
  okc <- okb & !is.na(cl$age) & !is.na(cl$stage_group)
  cat(ep, ": univariable n =", sum(okb), "events", sum(ee[okb] == 1), "; adjusted n =", sum(okc), "events", sum(ee[okc] == 1), "\n")
  for (meth in c("GSVA_adjusted_age_stage", "meanZ_univariable")) {
    M <- if (meth == "GSVA_adjusted_age_stage") GS else MZ
    ok <- if (meth == "GSVA_adjusted_age_stage") okc else okb
    base <- data.frame(time = tt[ok], ev = ee[ok], age = cl$age[ok], stage = droplevels(cl$stage_group[ok]))
    for (i in seq_len(nrow(CIR))) {
      d <- base; d$core <- sc(M, CIR$core[i])[ok]; d$ext <- sc(M, CIR$ext[i])[ok]; d$six <- sc(M, CIR$key[i])[ok]
      if (meth == "GSVA_adjusted_age_stage") {
        f3 <- coxph(Surv(time, ev) ~ core + age + stage, data = d); f6 <- coxph(Surv(time, ev) ~ six + age + stage, data = d)
        fn <- coxph(Surv(time, ev) ~ core + ext + age + stage, data = d)
      } else {
        f3 <- coxph(Surv(time, ev) ~ core, data = d); f6 <- coxph(Surv(time, ev) ~ six, data = d)
        fn <- coxph(Surv(time, ev) ~ core + ext, data = d)
      }
      lrt <- anova(f3, fn)[2, "Pr(>|Chi|)"]
      c3 <- summary(f3)$concordance[1]; c6 <- summary(f6)$concordance[1]
      k <- k + 1L
      res[[k]] <- data.table(method = meth, endpoint = ep, family = CIR$family[i], members = CIR$key[i], core = CIR$core[i],
                             extension = CIR$ext[i], n = nrow(d), events = sum(d$ev == 1),
                             HR_six = unname(exp(coef(f6)["six"])), p_six = summary(f6)$coefficients["six", "Pr(>|z|)"],
                             C_core = c3, C_six = c6, deltaC = c6 - c3, LRT_p_extension_adds = lrt)
    }
  }
}
R <- rbindlist(res)
R[, q_BH := p.adjust(LRT_p_extension_adds, "BH"), by = .(method, endpoint, family)]
fwrite(R, file.path(HERE, "s8_circuits.csv"))
S <- R[, .(circuits = .N, median_deltaC = median(deltaC), pct_six_beats_core = 100 * mean(deltaC > 0),
           n_LRT_nominal_p05 = sum(LRT_p_extension_adds < 0.05), n_q_lt_005 = sum(q_BH < 0.05)), by = .(method, endpoint, family)]
fwrite(S, file.path(HERE, "s8_summary.csv")); print(S)
