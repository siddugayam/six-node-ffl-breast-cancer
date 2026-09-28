## Pairwise subtype contrasts + within-subtype pooled estimates (TCGA-BRCA, METABRIC)
setwd("/path/to/revision")
OUT <- "results/v2"; options(width=200)
S <- read.csv(file.path(OUT,"subtype_stratified.csv"), stringsAsFactors=FALSE)
fz <- function(r) 0.5*log((1+r)/(1-r))
ztest <- function(r1,n1,r2,n2){
  if (any(!is.finite(c(r1,n1,r2,n2))) || n1<6 || n2<6) return(c(z=NA,p=NA))
  z <- (fz(r1)-fz(r2))/sqrt(1/(n1-3)+1/(n2-3)); c(z=z, p=2*pnorm(-abs(z)))
}
out <- list()
for (co in unique(S$cohort)){
  Sc <- S[S$cohort==co, ]
  sts <- setdiff(unique(Sc$stratum), "ALL")
  for (ax in unique(Sc$axis)){
    d <- Sc[Sc$axis==ax, ]
    if (!("Basal" %in% d$stratum)) next
    b <- d[d$stratum=="Basal", ]
    rest <- d[d$stratum %in% setdiff(sts,"Basal"), ]
    ## pooled rest (fixed effect on Fisher z)
    w <- rest$n-3; zb <- sum(w*fz(rest$rho))/sum(w); nr <- sum(rest$n)
    se_rest <- 1/sqrt(sum(w))
    zz <- (fz(b$rho)-zb)/sqrt(1/(b$n-3)+se_rest^2)
    ## within-subtype pooled over ALL subtypes
    wa <- d$n-3; za <- sum(wa*fz(d$rho))/sum(wa)
    allr <- Sc$rho[Sc$axis==ax & Sc$stratum=="ALL"][1]
    alln <- Sc$n[Sc$axis==ax & Sc$stratum=="ALL"][1]
    row <- data.frame(cohort=co, axis=ax,
      rho_marginal_ALL=allr, n_ALL=alln,
      rho_pooled_within_subtype=tanh(za),
      shrinkage_pct=100*(1 - abs(tanh(za))/max(abs(allr),1e-9)),
      rho_Basal=b$rho, n_Basal=b$n,
      rho_pooled_nonBasal=tanh(zb), n_nonBasal=nr,
      z_Basal_vs_rest=zz, p_Basal_vs_rest=2*pnorm(-abs(zz)), stringsAsFactors=FALSE)
    for (st in setdiff(sts,"Basal")){
      o <- d[d$stratum==st,]
      t <- if (nrow(o)) ztest(b$rho,b$n,o$rho,o$n) else c(z=NA,p=NA)
      row[[paste0("p_Basal_vs_",gsub("-","",st))]] <- unname(t["p"])
    }
    out[[length(out)+1]] <- row
  }
}
allnm <- unique(unlist(lapply(out, names)))
out <- lapply(out, function(x){ for (nm in setdiff(allnm, names(x))) x[[nm]] <- NA; x[allnm] })
CO <- do.call(rbind, out)
rownames(CO) <- NULL
CO$p_Basal_vs_rest_BH <- p.adjust(CO$p_Basal_vs_rest, method="BH")
write.csv(CO, file.path(OUT,"subtype_basal_contrast_v3.csv"), row.names=FALSE)
cat("=== BASAL vs REST, and marginal vs within-subtype-pooled ===\n")
print(CO[,c("cohort","axis","rho_marginal_ALL","rho_pooled_within_subtype","shrinkage_pct",
            "rho_Basal","n_Basal","rho_pooled_nonBasal","z_Basal_vs_rest","p_Basal_vs_rest","p_Basal_vs_rest_BH")], digits=3)
cat("\nrows:", nrow(CO), "\n")
