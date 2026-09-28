options(width=250)
setwd("/path/to/revision")
OUT <- "results/v2"
E <- read.csv(file.path(OUT,"newcohorts_per_edge.csv"), stringsAsFactors=FALSE)
S <- read.csv(file.path(OUT,"newcohorts_summary.csv"), stringsAsFactors=FALSE)
M <- read.csv(file.path(OUT,"newcohorts_meta.csv"), stringsAsFactors=FALSE)

cat("\n================ COHORT INVENTORY ================\n")
inv <- unique(E[, c("cohort","accession","n")])
inv <- inv[!duplicated(inv$cohort), ]
print(inv, row.names=FALSE)

cat("\n================ (a) HUB DIRECTION CONCORDANCE ================\n")
a <- S[grepl("^a[123]", S$analysis), c("cohort","analysis","n_tested","n_concordant","pct","binom_p")]
print(a, row.names=FALSE, digits=3)

cat("\n================ (b) COL1A1 ~ COL3A1 ================\n")
b <- E[E$edge=="COL1A1~COL3A1", c("cohort","n","rho_raw","p_raw","rho_cafA","p_cafA","rho_cafB","pct_retained_cafA")]
print(b, row.names=FALSE, digits=3)

cat("\n================ (c) TF -> COL1A1 ================\n")
for(tf in c("ETS1","NFKB1","RELA","SP1")){
  cc <- E[E$edge==paste0(tf,"->COL1A1"), c("cohort","n","rho_raw","p_raw","rho_cafA","p_cafA","rho_cafB")]
  cc$inverts <- is.finite(cc$rho_raw) & is.finite(cc$rho_cafA) & cc$rho_raw > 0 & cc$rho_cafA < 0
  cat("\n--", tf, "-> COL1A1 --\n"); print(cc, row.names=FALSE, digits=3)
}

cat("\n================ (d) SURVIVAL: modules ================\n")
d <- S[grepl("^d_", S$analysis), c("cohort","analysis","n","pct","binom_p","extra")]
names(d)[names(d)=="pct"] <- "HR_per_SD_univariate"
names(d)[names(d)=="binom_p"] <- "p_univariate"
print(d, row.names=FALSE, digits=3)

cat("\n================ (d) SURVIVAL: hubs passing FDR<0.05 per cohort ================\n")
h <- E[E$class=="cox_hub" & is.finite(E$q_uni) & E$q_uni < 0.05,
       c("cohort","target","source","n","rho_raw","p_raw","q_uni","rho_cafA","q_adj")]
names(h) <- c("cohort","endpoint","hub","n","HR_uni","p_uni","q_uni","HR_adj","q_adj")
print(h[order(h$cohort, h$q_uni), ], row.names=FALSE, digits=3)

cat("\n================ (e) RANDOM-EFFECTS META-ANALYSIS (tumour cohorts) ================\n")
m <- M[M$row_type=="RE_pooled" & M$pool=="tumour_only" & !grepl("^COXMETA", M$edge),
       c("edge","adjustment","k","n","rho","ci_lo_rho","ci_hi_rho","pooled_z","p_pooled","I2","tau2","p_Q")]
print(m, row.names=FALSE, digits=3)

cat("\n================ (e) META-ANALYSIS of Cox HRs, FDR<0.05 ================\n")
cx <- M[M$row_type=="RE_pooled" & grepl("^COXMETA", M$edge), ]
cx <- cx[order(cx$p_pooled), c("edge","k","n","rho","ci_lo_rho","ci_hi_rho","I2","p_pooled","q_pooled")]
names(cx)[names(cx)=="rho"] <- "HR_per_SD"
print(cx[is.na(cx$q_pooled) | cx$q_pooled < 0.05, ], row.names=FALSE, digits=3)

cat("\n================ normal(GTEx) vs pooled tumour ================\n")
nv <- S[grepl("^e_normal_vs_tumour", S$analysis), c("analysis","n","n_tested","pct","binom_p","extra")]
print(nv, row.names=FALSE, digits=3)
