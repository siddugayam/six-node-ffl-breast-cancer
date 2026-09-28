## ==========================================================================
## N4_normal_vs_tumour.R -- contrast the GTEx normal-breast correlations with
## the random-effects pooled tumour estimate for the same edge, and append the
## rows to newcohorts_summary.csv.
## ==========================================================================
setwd("/path/to/revision")
OUT <- "results/v2"
E <- read.csv(file.path(OUT,"newcohorts_per_edge.csv"), stringsAsFactors=FALSE)
M <- read.csv(file.path(OUT,"newcohorts_meta.csv"), stringsAsFactors=FALSE)
S <- read.csv(file.path(OUT,"newcohorts_summary.csv"), stringsAsFactors=FALSE)
fz <- function(r) 0.5*log((1+r)/(1-r))
EDGES <- c("COL1A1~COL3A1","ETS1->COL1A1","NFKB1->COL1A1","RELA->COL1A1","SP1->COL1A1")
new <- list()
for(ed in EDGES) for(aj in c("raw","cafA")){
  g <- E[E$edge==ed & E$cohort=="GTEx_breast", ]
  pm <- M[M$edge==ed & M$adjustment==aj & M$pool=="tumour_only" & M$row_type=="RE_pooled", ]
  if(!nrow(g) || !nrow(pm)) next
  rg <- if(aj=="raw") g$rho_raw else g$rho_cafA
  sg <- if(aj=="raw") g$se_raw  else g$se_cafA
  if(!is.finite(rg) || !is.finite(sg)) next
  zdiff <- (fz(rg) - pm$fisher_z) / sqrt(sg^2 + pm$se_z^2)
  p <- 2*pnorm(-abs(zdiff))
  cat(sprintf("%-14s %-4s  GTEx normal rho=%+.3f (n=%d)  vs pooled tumour rho=%+.3f (k=%d)  z=%+.2f p=%.3g\n",
              ed, aj, rg, g$n, pm$rho, pm$k, zdiff, p))
  new[[length(new)+1]] <- data.frame(
    cohort="GTEx_breast_vs_pooled_tumour", accession="GTEx v8 vs 7-cohort RE pool",
    platform="mixed", n=g$n,
    analysis=paste0("e_normal_vs_tumour_contrast_", gsub("[^A-Za-z0-9]","_",ed), "_", aj),
    n_tested=pm$k, n_concordant=NA, pct=rg, binom_p=p,
    extra=sprintf("GTEx_normal_rho=%+.4f; pooled_tumour_rho=%+.4f; z_diff=%+.3f", rg, pm$rho, zdiff),
    stringsAsFactors=FALSE)
}
if(length(new)){
  S2 <- rbind(S, do.call(rbind,new))
  write.csv(S2, file.path(OUT,"newcohorts_summary.csv"), row.names=FALSE)
  cat("\nappended", length(new), "contrast rows ->", nrow(S2), "summary rows\n")
}
