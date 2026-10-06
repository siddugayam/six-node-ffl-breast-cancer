# Adds Benjamini-Hochberg q values across the miRNAs tested to the pooled survival meta-analysis
# (results/v6/mirna_meta_pooled.csv) and to the forest-plot table.
suppressPackageStartupMessages(library(data.table))
setwd("/path/to/revision"); OUT<-"results/v6"
M<-fread(file.path(OUT,"mirna_meta_pooled.csv"))
M[row_type=="RE_pooled", q_BH := p.adjust(p,"BH"), by=analysis]
fwrite(M, file.path(OUT,"mirna_meta_pooled.csv"))
F<-fread(file.path(OUT,"mirna_meta_forest_table.csv"))
F <- merge(F, M[row_type=="RE_pooled" & analysis=="primary_endpoint_univariate", .(miRNA,q_BH)],
           by="miRNA", all.x=TRUE)
fwrite(F, file.path(OUT,"mirna_meta_forest_table.csv"))
cat("=== primary univariate pooled, with BH q across the 11 miRNAs ===\n")
print(as.data.frame(M[row_type=="RE_pooled"&analysis=="primary_endpoint_univariate"][order(p),
  .(miRNA,k,n,nevent,HR=round(HR,3),lo=round(lo,3),hi=round(hi,3),p=signif(p,3),
    q_BH=signif(q_BH,3),I2=round(I2,1),tau2=signif(tau2,3),p_Q=signif(p_Q,3))]))
cat("\n=== no-TCGA pooled, with BH q ===\n")
print(as.data.frame(M[row_type=="RE_pooled"&analysis=="primary_endpoint_univariate_noTCGA"][order(p),
  .(miRNA,k,n,nevent,HR=round(HR,3),lo=round(lo,3),hi=round(hi,3),p=signif(p,3),
    q_BH=signif(q_BH,3),I2=round(I2,1))]))
cat("\n=== per-cohort rows for the two headline miRNAs ===\n")
C<-fread(file.path(OUT,"mirna_meta_percohort_cox.csv"))
print(as.data.frame(C[model=="univariate"&primary==1&miRNA %in% c("miR-29a","miR-195"),
  .(miRNA,cohort,accession,endpoint,probe,n,nevent,HR=round(HR,3),lo=round(lo,3),
    hi=round(hi,3),p=signif(p,3))][order(miRNA,cohort)]))
