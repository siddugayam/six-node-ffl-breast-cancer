#!/usr/bin/env Rscript
# Is the TCGA "no hub is prognostic" null a power artefact?
# Downsample METABRIC to TCGA's event count and re-run the identical hub Cox screen.
suppressMessages({library(survival)})
set.seed(42)
RES  <- "/path/to/revision/INBOX_2026-09-26/rerun/A3_power_downsample/"
DATA <- "/path/to/revision/data/"
mb <- readRDS(paste0(DATA,"metabric.rds")); M <- mb$M; cl <- as.data.frame(mb$cl)
rownames(cl) <- cl$PATIENT_ID; cl <- cl[colnames(M),]
cl$OS_ev <- as.numeric(cl$OS_STATUS=="1:DECEASED"); cl$OS_t <- as.numeric(cl$OS_MONTHS)
hubs <- read.csv("/path/to/revision/results/network_topology_hubs.csv", stringsAsFactors=FALSE)
hg <- intersect(hubs$name[hubs$is_hub & hubs$type!="miRNA"], rownames(M))
ok <- !is.na(cl$OS_ev) & !is.na(cl$OS_t) & cl$OS_t > 0
X <- t(scale(t(M[hg, ok]))); ev <- cl$OS_ev[ok]; tt <- cl$OS_t[ok]
cat("METABRIC full: n =", sum(ok), " events =", sum(ev), " hubs =", length(hg), "\n")
tcga_events <- 151
frac <- tcga_events/sum(ev)
NREP <- 200
nsig <- integer(NREP); nnom <- integer(NREP); nev <- integer(NREP)
for (r in 1:NREP) {
  idx <- sample(seq_len(ncol(X)), size=round(frac*ncol(X)))
  p <- sapply(hg, function(g) {
    m <- try(coxph(Surv(tt[idx], ev[idx]) ~ X[g, idx]), silent=TRUE)
    if (inherits(m,"try-error")) NA else summary(m)$coefficients[1,"Pr(>|z|)"] })
  q <- p.adjust(p, "BH")
  nsig[r] <- sum(q < 0.05, na.rm=TRUE); nnom[r] <- sum(p < 0.05, na.rm=TRUE); nev[r] <- sum(ev[idx])
}
out <- data.frame(
  metric=c("METABRIC_full_n","METABRIC_full_events","TCGA_OS_events_matched","downsampled_mean_n",
           "downsampled_mean_events","reps",
           "mean_hubs_BHq05_downsampled","median_hubs_BHq05_downsampled","pct_reps_with_zero_hubs_BHq05",
           "mean_hubs_nominal_p05_downsampled",
           "hubs_BHq05_METABRIC_full_OS","hubs_BHq05_TCGA_full_OS","hubs_nominal_p05_TCGA_full_OS"),
  value=c(sum(ok), sum(ev), tcga_events, round(frac*ncol(X)), mean(nev), NREP,
          mean(nsig), median(nsig), 100*mean(nsig==0), mean(nnom), 19, 0, 8))
write.csv(out, paste0(RES,"metabric_power_downsampling.csv"), row.names=FALSE)
print(out, row.names=FALSE)
## ---- added for INBOX_2026-09-26: per-subsample counts (the loop above is unchanged) ----
per <- data.frame(rep=1:NREP, n_subsample=round(frac*ncol(X)), events=nev,
                  hubs_BHq05=nsig, hubs_nominal_p05=nnom)
write.csv(per, paste0(RES,"metabric_power_downsampling_per_subsample.csv"), row.names=FALSE)
cat("hubs reaching BH q<0.05 per subsample: median", median(nsig), " IQR", quantile(nsig,.25), "-", quantile(nsig,.75),
    " mean", mean(nsig), " range", min(nsig), "-", max(nsig), "\n")
cat("fraction of subsamples with >=1 hub at q<0.05:", mean(nsig>=1), " (", sum(nsig>=1), "of", NREP, ")\n")
cat("events per subsample: mean", mean(nev), " median", median(nev), " range", min(nev), "-", max(nev), "\n")
cat("R", R.version.string, " survival", as.character(packageVersion("survival")), " RNGkind", paste(RNGkind(), collapse="/"), "\n")
