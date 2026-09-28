## ---------------------------------------------------------------------------
## v7 / BioNet  step 4: (a) edge-resampling stability of module membership
##                      (b) is the overlap with the prioritised 30 anything more
##                          than shared dependence on the same DE p-values?
## ---------------------------------------------------------------------------
suppressMessages({library(BioNet); library(igraph); library(data.table); library(graph)})
set.seed(20260912)
ROOT <- "/path/to/revision"
OUT  <- file.path(ROOT,"results/v7"); P <- function(...) file.path(OUT, paste0("bionet_", ...))
gp <- readRDS(P("graph_and_pvals.rds")); nod <- gp$nod; G <- gp$G; gi <- gp$gi
fits <- readRDS(P("bum_fits.rds")); tauf <- function(fb,f) BioNet:::fdrThreshold(f,fb)
pri30 <- fread(file.path(ROOT,"results/v5/tables/Table4_prioritised_30.csv"))$name
FDR <- fread(P("module_summary.csv"))[run=="primary"]$fdr

sc <- setNames(ifelse(nod$is_mir,
        (fits$mir_all$a-1)*(log(nod$pval)-log(tauf(fits$mir_all,FDR))),
        (fits$gene_gw$a-1)*(log(nod$pval)-log(tauf(fits$gene_gw,FDR)))), nod$name)
mb0 <- graph::nodes(runFastHeinz(G, sc[graph::nodes(G)]))
k <- length(mb0); cat("primary module n =", k, " at FDR", format(FDR,digits=3), "\n")

## ---- (a) edge-resampling stability (90% of edges, 200 replicates) --------
B <- 200; keepfrac <- 0.90
cnt <- setNames(integer(length(sc)), names(sc)); sizes <- integer(B)
el <- as_edgelist(gi)
for(b in 1:B){
  idx <- sample.int(nrow(el), floor(keepfrac*nrow(el)))
  gb  <- graph_from_edgelist(el[idx,,drop=FALSE], directed=FALSE)
  gb  <- add_vertices(gb, length(setdiff(V(gi)$name, V(gb)$name)),
                      name=setdiff(V(gi)$name, V(gb)$name))
  Gb  <- rmSelfLoops(igraph::as_graphnel(simplify(gb)))
  mbb <- runFastHeinz(Gb, sc[graph::nodes(Gb)])
  mbb <- if(is.null(mbb)) character(0) else graph::nodes(mbb)
  sizes[b] <- length(mbb); cnt[mbb] <- cnt[mbb] + 1L
}
stab <- data.table(name=names(cnt), type=nod$type[match(names(cnt),nod$name)],
                   pval=nod$pval[match(names(cnt),nod$name)],
                   selection_freq=cnt/B, in_primary=names(cnt) %in% mb0,
                   in_prioritised30=names(cnt) %in% pri30)
setorder(stab, -selection_freq)
fwrite(stab, P("edge_resample_stability.csv"))
cat(sprintf("edge-resample module sizes: median %d (IQR %d-%d)\n", median(sizes),
            quantile(sizes,.25), quantile(sizes,.75)))
cat("selection frequency of the paper's key nodes:\n")
print(stab[name %in% c("COL1A1","COL3A1","FN1","SPARC","POSTN","THBS2","TGFBI","EZH2",
                       "ETS1","NFKB1","RELA","SP1","hsa-miR-29a","hsa-miR-101"),
           .(name,type,pval,selection_freq,in_primary)])

## ---- (b) is the pri30 overlap more than shared use of the same p-values? --
ord   <- names(sort(sc, decreasing=TRUE))
topk  <- ord[1:k]                                   # DE ranking alone, no network
cat(sprintf("\ntop-%d by DE score alone: overlap with pri30 = %d ; BioNet module = %d ; Jaccard(module,topk) = %.3f\n",
            k, length(intersect(topk,pri30)), length(intersect(mb0,pri30)),
            length(intersect(mb0,topk))/length(union(mb0,topk))))

## p-value-rank-stratified null: draw k nodes matched to the module's DE-rank deciles
rk  <- rank(-sc); dec <- cut(rk, breaks=quantile(rk, seq(0,1,0.1)), include.lowest=TRUE, labels=FALSE)
tab <- table(dec[match(mb0, names(sc))])
nullov <- replicate(10000, {
  s <- unlist(lapply(names(tab), function(d) sample(names(sc)[dec==as.integer(d)], tab[[d]])))
  length(intersect(s, pri30))})
obs <- length(intersect(mb0, pri30))
cat(sprintf("DE-rank-matched null: observed overlap %d, null mean %.2f (sd %.2f), p = %.4f\n",
            obs, mean(nullov), sd(nullov), (sum(nullov>=obs)+1)/(length(nullov)+1)))
## unconditional null for comparison
nullov2 <- replicate(10000, length(intersect(sample(names(sc), k), pri30)))
cat(sprintf("unconditional null   : observed overlap %d, null mean %.2f, p = %.4f\n",
            obs, mean(nullov2), (sum(nullov2>=obs)+1)/(length(nullov2)+1)))
fwrite(data.table(null=c(rep("DE_rank_matched",length(nullov)), rep("unconditional",length(nullov2))),
                  overlap=c(nullov,nullov2)), P("pri30_overlap_nulls.csv"))

## ---- (c) where do the collagen nodes sit? --------------------------------
key <- c("COL1A1","COL3A1","FN1","SPARC","POSTN","THBS2","TGFBI","ASPN","SPP1","FBN1",
         "hsa-miR-29a","hsa-miR-101","ETS1","NFKB1","RELA","SP1")
kd <- data.table(name=key, type=nod$type[match(key,nod$name)], pval=nod$pval[match(key,nod$name)],
                 score_at_primary=sc[key], rank_by_score=rank(-sc)[key],
                 in_primary_module=key %in% mb0)
setorder(kd, rank_by_score); fwrite(kd, P("key_stroma_nodes_status.csv")); print(kd)
