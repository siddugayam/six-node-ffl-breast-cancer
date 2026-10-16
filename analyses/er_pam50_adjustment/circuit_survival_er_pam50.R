## circuit_survival_er_pam50.R -- the six-node circuit survival test of analyses/six_node_pattern/S8/s8_circuit_survival.R
## (GSVA scores of each six-node set, its three-node core and its extension; Cox models for overall survival adjusted for
## age and AJCC stage; likelihood-ratio test of core vs core + extension; BH within family), repeated with oestrogen-
## receptor status and PAM50 added (see SETTINGS.md). Data preparation and scoring are copied from s8_circuit_survival.R.
## Three fits per circuit: published (age + stage, published sample), same tumours (age + stage on the tumours with ER
## and PAM50), ER and PAM50 (age + stage + ER + PAM50 on those tumours). Overall survival only.
## Root of the repository: environment variable FFL_ROOT. No sample-level values are written.
## usage: FFL_ROOT=<repository checkout> Rscript circuit_survival_er_pam50.R
suppressPackageStartupMessages({ library(data.table); library(survival); library(GSVA) })
REV <- Sys.getenv("FFL_ROOT"); stopifnot(nzchar(REV))
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
S1 <- file.path(REV, "analyses", "six_node_pattern", "S1")
set.seed(1234)
names587 <- readLines(file.path(S1, "node_names.txt"))[1:587]
rd <- function(f, cols) { d <- fread(f, header = FALSE, col.names = cols); d[, lapply(.SD, function(i) names587[i + 1])] }
BH <- rd(file.path(S1, "obs/noexemplar_bhat6_composite_instances.tsv"), c("miR1", "miR2", "TF1", "TF2", "G1", "G2"))
MF <- rd(file.path(S1, "obs/noexemplar_model6_full_instances.tsv"), c("TF1", "TF2", "miR1", "miR2", "G1", "G2"))
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
cat("GSVA done in", sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "mins"))), "min; ", paste(dim(GS), collapse = " x "), "\n")
Z <- t(scale(t(X))); MZ <- t(sapply(gsl, function(s) colMeans(Z[s, , drop = FALSE])))
sc <- function(M, s) as.numeric(scale(M[s, ]))

## ---- ER and PAM50 ------------------------------------------------------------------------------------------------
clin <- fread(file.path(REV, "data/brca_clinicalMatrix.tsv"), sep = "\t", quote = "")
mm <- match(common, clin$sampleID); cat("tumours with a clinical-matrix row:", sum(!is.na(mm)), "of", length(common), "\n")
fac <- function(x, lev) factor(ifelse(x %in% lev, x, NA), levels = lev)
ER <- fac(clin$breast_carcinoma_estrogen_receptor_status[mm], c("Positive", "Negative"))
PAM50 <- fac(clin$PAM50Call_RNAseq[mm], c("LumA", "LumB", "Her2", "Basal", "Normal"))

tt <- cl[["OS.time"]]; ee <- cl[["OS"]]
okc <- is.finite(tt) & tt > 0 & is.finite(ee) & !is.na(cl$age) & !is.na(cl$stage_group)
okE <- okc & !is.na(ER) & !is.na(PAM50)
cat("OS: published sample n =", sum(okc), "events", sum(ee[okc] == 1), "; with ER and PAM50 n =", sum(okE), "events", sum(ee[okE] == 1), "\n")
fit_one <- function(d, cov) {
  f3 <- coxph(as.formula(paste("Surv(time, ev) ~ core +", cov)), data = d)
  f6 <- coxph(as.formula(paste("Surv(time, ev) ~ six +", cov)), data = d)
  fn <- coxph(as.formula(paste("Surv(time, ev) ~ core + ext +", cov)), data = d)
  s6 <- summary(f6)$coefficients["six", ]
  data.table(n = nrow(d), events = sum(d$ev == 1), HR_six = unname(exp(s6["coef"])),
             HR_lo = unname(exp(s6["coef"] - 1.96 * s6["se(coef)"])), HR_hi = unname(exp(s6["coef"] + 1.96 * s6["se(coef)"])),
             p_six = unname(s6["Pr(>|z|)"]), C_core = summary(f3)$concordance[1], C_six = summary(f6)$concordance[1],
             LRT_p_extension_adds = anova(f3, fn)[2, "Pr(>|Chi|)"])
}
FITS <- list(published = list(ok = okc, cov = "age + stage"), same_tumours = list(ok = okE, cov = "age + stage"),
             er_pam50 = list(ok = okE, cov = "age + stage + ER + PAM50"))
res <- list(); k <- 0L
for (fn in names(FITS)) {
  ok <- FITS[[fn]]$ok
  base <- data.frame(time = tt[ok], ev = ee[ok], age = cl$age[ok], stage = droplevels(cl$stage_group[ok]),
                     ER = droplevels(ER[ok]), PAM50 = droplevels(PAM50[ok]))
  for (i in seq_len(nrow(CIR))) {
    d <- base; d$core <- sc(GS, CIR$core[i])[ok]; d$ext <- sc(GS, CIR$ext[i])[ok]; d$six <- sc(GS, CIR$key[i])[ok]
    k <- k + 1L
    res[[k]] <- cbind(data.table(fit = fn, covariates = FITS[[fn]]$cov, family = CIR$family[i], members = CIR$key[i],
                                 core = CIR$core[i], extension = CIR$ext[i]), fit_one(d, FITS[[fn]]$cov))
  }
  cat(fn, "done\n")
}
R <- rbindlist(res)
R[, q_BH := p.adjust(LRT_p_extension_adds, "BH"), by = .(fit, family)]

## control: the published fit reproduces the stored S8 table (overall survival, GSVA, age + stage)
P <- fread(file.path(REV, "analyses/six_node_pattern/S8/six_node_pattern_survival_all_circuits.csv"))
P <- P[method == "GSVA_adjusted_age_stage" & endpoint == "OS"]
Q <- merge(R[fit == "published", .(family, members, n, events, HR_six, p_six, LRT_p_extension_adds, q_BH)],
           P[, .(family, members, n_p = n, events_p = events, HR_p = HR_six, p_six_p = p_six, LRT_p = LRT_p_extension_adds, q_p = q_BH)],
           by = c("family", "members"))
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
ctrl <- data.table(rows_compared = nrow(Q), rows_published = nrow(P), n_equal = all(Q$n == Q$n_p & Q$events == Q$events_p),
                   max_rel_HR = rel(Q$HR_six, Q$HR_p), max_rel_p_six = rel(Q$p_six, Q$p_six_p),
                   max_rel_LRT = rel(Q$LRT_p_extension_adds, Q$LRT_p), max_rel_q = rel(Q$q_BH, Q$q_p))
ctrl[, reproduced := rows_compared == rows_published & n_equal & max_rel_HR < 1e-6 & max_rel_p_six < 1e-6 & max_rel_LRT < 1e-6 & max_rel_q < 1e-6]
fwrite(ctrl, file.path(HERE, "control_reproduction.csv")); print(ctrl); stopifnot(ctrl$reproduced)

## plain names for the two circuit families
R[, family := c(BHAT6_composite = "six-node composite pattern", MODEL6_architecture = "six-node model architecture")[family]]
stopifnot(!anyNA(R$family))
fwrite(R, file.path(HERE, "circuits_er_pam50.csv"))
SIX <- R[fit == "published" & family == "six-node composite pattern" & q_BH < 0.05, members]
cat("published patterns with q < 0.05:", length(SIX), "\n")
W <- dcast(R[members %in% SIX & family == "six-node composite pattern"], members ~ fit,
           value.var = c("n", "events", "HR_six", "HR_lo", "HR_hi", "p_six", "LRT_p_extension_adds", "q_BH"))
fwrite(W, file.path(HERE, "sixteen_patterns_before_after.csv"))
S <- R[, .(circuits = .N, n = n[1], events = events[1], n_LRT_p05 = sum(LRT_p_extension_adds < 0.05), n_q_lt_005 = sum(q_BH < 0.05),
           n_HR_p05 = sum(p_six < 0.05), sixteen_LRT_p05 = sum(LRT_p_extension_adds[members %in% SIX] < 0.05),
           sixteen_q_lt_005 = sum(q_BH[members %in% SIX] < 0.05), sixteen_HR_p05 = sum(p_six[members %in% SIX] < 0.05),
           median_deltaC = median(C_six - C_core)), by = .(fit, family)]
fwrite(S, file.path(HERE, "circuit_summary.csv")); print(S)
