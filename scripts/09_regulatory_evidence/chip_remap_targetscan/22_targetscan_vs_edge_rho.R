## Does TargetScan site quality predict which network miRNA-target edges are
## sign-concordant (rho < 0) in TCGA-BRCA?  All inputs recomputed in this run.
suppressWarnings(suppressMessages({library(data.table)}))
BASE <- "/path/to/revision"
rho <- fread(file.path(BASE,"results/v2/v2d_edge_rho_recomputed.csv"))
ts  <- fread(file.path(BASE,"results/v2/targetscan_sites.csv"))
ts  <- ts[record_type=="network_edge"]
setnames(ts,c("network_miRNA","gene"),c("source","ts_gene"))
cat("network_edge rows in targetscan_sites.csv:",nrow(ts),"\n")
E <- fread(file.path(BASE,"data/canonical_edges.tsv"))[edge_type=="miRNA_target"]
ts[, target := E$target ]                      # same order, same file, 1:1
stopifnot(nrow(ts)==nrow(E), all(ts$source==E$source))
d <- merge(rho[edge_type=="miRNA_target"], ts, by=c("source","target"), all.x=TRUE)
cat("merged edges:",nrow(d)," with rho:",sum(!is.na(d$rho)),"\n")
d <- d[!is.na(rho)]
d[, concordant := as.integer(rho < 0)]
d[, has_site   := as.integer(n_sites > 0)]
d[, has_cons   := as.integer(n_conserved_sites > 0)]
d[, has_8mer   := as.integer(n_8mer > 0)]
d[, ctx := fifelse(is.na(total_context_pp), 0, total_context_pp)]
cat("\n=== n =",nrow(d)," measurable miRNA_target edges ===\n")
cat("overall concordance (rho<0):",round(mean(d$concordant),4),"\n")
cat("mean rho:",round(mean(d$rho),5),"\n\n")

## 1. concordance by presence/absence of a predicted site
for (v in c("has_site","has_cons","has_8mer")) {
  t <- table(d[[v]], d$concordant)
  pr <- prop.table(t,1)[,2]
  bt <- prop.test(t[,2], rowSums(t))
  cat(sprintf("%-9s : no-site %d edges %.4f concordant | site %d edges %.4f concordant | diff %.4f | prop.test p=%.3g\n",
      v, sum(d[[v]]==0), pr[1], sum(d[[v]]==1), pr[2], pr[2]-pr[1], bt$p.value))
}
cat("\n")
## 2. mean rho by site class
for (v in c("has_site","has_cons","has_8mer")) {
  a <- d[get(v)==0, rho]; b <- d[get(v)==1, rho]
  w <- wilcox.test(a,b)
  cat(sprintf("%-9s : mean rho no=%+.5f (n=%d) vs yes=%+.5f (n=%d) ; Wilcoxon p=%.3g\n",
      v, mean(a), length(a), mean(b), length(b), w$p.value))
}
## 3. is context++ score associated with rho / concordance?
cat("\n--- context++ (total, more negative = better predicted repression) ---\n")
sub <- d[!is.na(total_context_pp)]
cat("edges with a context++ score:",nrow(sub),"\n")
ct <- cor.test(sub$total_context_pp, sub$rho, method="spearman")
cat(sprintf("Spearman(total context++, rho) = %+.4f  p=%.3g  n=%d\n", ct$estimate, ct$p.value, nrow(sub)))
cat("  (positive rho-vs-score correlation = better sites give MORE NEGATIVE rho)\n")
gl <- glm(concordant ~ total_context_pp, data=sub, family=binomial)
s <- summary(gl)$coefficients
cat(sprintf("logistic concordant ~ context++ : beta=%+.4f SE=%.4f z=%+.2f p=%.3g\n",
    s[2,1],s[2,2],s[2,3],s[2,4]))
## quintiles of context++ among edges that have a site
sub[, q := cut(total_context_pp, breaks=quantile(total_context_pp, probs=seq(0,1,0.2)),
               include.lowest=TRUE, labels=paste0("Q",1:5))]
qt <- sub[, .(n=.N, min_ctx=round(min(total_context_pp),3), max_ctx=round(max(total_context_pp),3),
              concordance=round(mean(concordant),4), mean_rho=round(mean(rho),5)), by=q][order(q)]
cat("\nQuintiles of total context++ (Q1 = strongest predicted sites):\n"); print(qt)
tq <- prop.test(c(sum(sub[q=="Q1",concordant]), sum(sub[q=="Q5",concordant])),
                c(nrow(sub[q=="Q1"]), nrow(sub[q=="Q5"])))
cat(sprintf("Q1 vs Q5 concordance: %.4f vs %.4f, prop.test p=%.3g\n",
    mean(sub[q=="Q1",concordant]), mean(sub[q=="Q5",concordant]), tq$p.value))
## 4. site-count dose response
dr <- d[, .(n=.N, concordance=round(mean(concordant),4), mean_rho=round(mean(rho),5)),
        by=.(n_sites=pmin(n_sites,4))][order(n_sites)]
cat("\nConcordance by number of predicted sites (capped at 4+):\n"); print(dr)
dc <- d[, .(n=.N, concordance=round(mean(concordant),4), mean_rho=round(mean(rho),5)),
        by=.(n_conserved=pmin(n_conserved_sites,3))][order(n_conserved)]
cat("\nConcordance by number of CONSERVED sites (capped at 3+):\n"); print(dc)
## 5. best single site type
d[, best_type := fifelse(n_8mer>0,"8mer", fifelse(n_7mer_m8>0,"7mer-m8",
                  fifelse(n_7mer_A1>0,"7mer-A1","none/6mer-only")))]
bt <- d[, .(n=.N, concordance=round(mean(concordant),4), mean_rho=round(mean(rho),5)),
        by=best_type][order(-concordance)]
cat("\nConcordance by best site type present:\n"); print(bt)
fwrite(d[,.(source,target,ts_gene,best_arm,rho,p,concordant,n_sites,n_conserved_sites,
            n_8mer,n_7mer_m8,n_7mer_A1,total_context_pp,cum_weighted_context_pp,
            aggregate_PCT,best_site_context_pp,best_site_context_pp_percentile,best_type)],
       file.path(BASE,"results/v2/targetscan_concordance_edges.csv"))
out <- rbind(
 data.table(metric="n_edges_measurable", value=nrow(d)),
 data.table(metric="overall_concordance", value=round(mean(d$concordant),4)),
 data.table(metric="concordance_no_site", value=round(mean(d[has_site==0,concordant]),4)),
 data.table(metric="concordance_with_site", value=round(mean(d[has_site==1,concordant]),4)),
 data.table(metric="concordance_with_conserved_site", value=round(mean(d[has_cons==1,concordant]),4)),
 data.table(metric="concordance_with_8mer", value=round(mean(d[has_8mer==1,concordant]),4)),
 data.table(metric="spearman_ctx_vs_rho", value=round(unname(ct$estimate),4)),
 data.table(metric="spearman_ctx_vs_rho_p", value=signif(ct$p.value,3)),
 data.table(metric="logistic_beta_ctx", value=round(s[2,1],4)),
 data.table(metric="logistic_p_ctx", value=signif(s[2,4],3)))
fwrite(out, file.path(BASE,"results/v2/targetscan_concordance_stats_v2d.csv"))
cat("\nWROTE targetscan_concordance_edges.csv and targetscan_concordance_stats_v2d.csv\n")
