#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table)})
OUT <- "/path/to/revision/results/multiomics"
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")
L <- fread(file.path(OUT,"lincs_TF_signatures_collagen.csv"))
Lcd <- L[value_col=="CD-coefficient"]
BREAST <- c("MCF7","HCC1806","MDAMB231","BT20","HS578T","SKBR3")
msg("breast-derived LINCS lines present:", paste(intersect(BREAST, unique(Lcd$cell_line)), collapse=","))
B <- Lcd[cell_line %in% BREAST]
msg("breast-line TF perturbation signatures:", nrow(B))
print(B[, .N, by=.(tf, pert_type, cell_line)][order(tf,pert_type)])

res <- rbindlist(lapply(c("COL1A1","COL3A1"), function(g){
  rbindlist(lapply(split(B, by=c("tf","pert_type"), drop=TRUE), function(s){
    v <- s[[paste0(g,"_value")]]; p <- s[[paste0(g,"_pctile")]]
    ok <- is.finite(v); v <- v[ok]; p <- p[ok]
    if(!length(v)) return(NULL)
    data.table(tf=s$tf[1], perturbation=s$pert_type[1], cell_lines=paste(sort(unique(s$cell_line)),collapse=","),
      readout_gene=g, n_signatures=length(v), median_CD=median(v), mean_CD=mean(v),
      median_pctile=median(p), n_up=sum(v>0), n_down=sum(v<0),
      binom_p = if(length(v)>=3) binom.test(sum(v>0), length(v), 0.5)$p.value else NA_real_)
  }))
}))
res[, expected_if_TF_activates := fifelse(perturbation=="Overexpression","up","down")]
res[, observed := fifelse(median_CD>0,"up","down")]
res[, supports_activation := expected_if_TF_activates==observed]
setorder(res, readout_gene, tf, perturbation)
print(res)
fwrite(res, file.path(OUT,"lincs_TF_collagen_breastlines.csv"))
msg("wrote lincs_TF_collagen_breastlines.csv rows:", nrow(res))
