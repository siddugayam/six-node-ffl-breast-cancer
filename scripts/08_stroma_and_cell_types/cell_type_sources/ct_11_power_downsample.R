#!/usr/bin/env Rscript
# Is the TCGA "no hub is prognostic" null a power artefact?
# Downsample METABRIC to TCGA's event count and re-run the identical hub Cox screen.
suppressMessages({library(survival)})
set.seed(42)
RES  <- "/path/to/revision/results/multiomics/"
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
