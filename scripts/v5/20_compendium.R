suppressMessages({library(data.table)})
REV <- "/path/to/revision"; R5 <- file.path(REV,"results/v5")
lg <- function(...) cat(...,"\n")
P     <- fread(file.path(R5,"node_prioritisation_full.csv"))
top   <- P[!is.na(priority)][, head(.SD, 10), by=type]
edges <- fread(file.path(REV,"data/canonical_edges.tsv"))
cores <- fread(file.path(REV,"results/v2/ffl_cores_coherence_corrected.csv"))
ec    <- unique(fread(file.path(REV,"results/edge_correlation.csv"))[, .(source,target,rho,p)])
ev    <- fread(file.path(REV,"data/edge_evidence_tier.tsv"))
setnames(ev, grep("tier", names(ev), value=TRUE)[1], "tier")

pos <- rbindlist(lapply(top$name, function(n){
  cr <- cores[regulator==n | intermediate==n | target==n]
  data.table(name=n, ffl_total=nrow(cr),
    as_regulator=cores[regulator==n,.N], as_intermediate=cores[intermediate==n,.N],
    as_target=cores[target==n,.N],
    classes=paste(sort(unique(cr$class)), collapse="; "),
    coherent_frac=if (nrow(cr[coherence!="unresolved"])) round(mean(substr(cr[coherence!="unresolved"]$coherence,1,1)=="C"),3) else NA_real_,
    out_degree=edges[source==n,.N], in_degree=edges[target==n,.N])
}))

part <- rbindlist(lapply(top$name, function(n){
  e <- edges[source==n | target==n][, .(source,target,edge_type)]
  e[, partner := fifelse(source==n, target, source)]
  e[, role    := fifelse(source==n, "regulates", "regulated by")]
  e <- merge(e, ec, by=c("source","target"), all.x=TRUE)
  e <- merge(e, ev[, .(source,target,tier)], by=c("source","target"), all.x=TRUE)
  e <- unique(e)[order(-abs(rho))][1:min(8,.N)]
  e[, node := n][, .(node, partner, role, edge_type, tier, rho=round(rho,3), p=signif(p,3))]
}), fill=TRUE)
fwrite(part, file.path(R5,"node_top_partners.csv"))

C <- merge(top, pos, by="name"); setorder(C, type, -priority)
keep <- intersect(c("name","type","priority","rank_within_type","ffl_total","as_regulator",
  "as_intermediate","as_target","classes","coherent_frac","out_degree","in_degree",
  "betweenness","frac_strong","n_edges","logFC","de_fdr","frac_same_dir","prot_rho",
  "meth_delta","mut_freq","surv_minp","chronos"), names(C))
fwrite(C[, ..keep], file.path(R5,"node_compendium_table.csv"))
lg("wrote node_compendium_table.csv:", nrow(C), "nodes,", length(keep), "fields")
lg("wrote node_top_partners.csv:", nrow(part), "rows")
for (ty in c("TF","Gene","miRNA")) {
  lg(sprintf("\n===== %s =====", ty))
  print(C[type==ty, .(name, ffl=ffl_total, reg=as_regulator, tgt=as_target, coh=coherent_frac,
        logFC=round(logFC,2), meth=round(meth_delta,3), surv=signif(surv_minp,2),
        chr=round(chronos,2), strong=round(frac_strong,2))])
}
