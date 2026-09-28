#!/usr/bin/env Rscript
# v7 / 10 : compile the headline numbers quoted in the write-up.
suppressMessages({library(igraph); library(data.table); library(mclust)})
setwd("/path/to/revision")
g <- readRDS("results/v7/cm2_graph_undirected.rds"); nm <- V(g)$name; N <- length(nm)
P <- readRDS("results/v7/cm2_all_partitions.rds"); s <- fread("results/v7/cm2_run_summary.csv")
nodes <- fread("data/canonical_nodes.tsv"); ntype <- setNames(nodes$type, nodes$name)
prior30 <- fread("results/v5/tables/Table4_prioritised_30.csv")$name
EMT <- intersect(unique(fread("results/v7/cm2_genesets_msigdb.csv")[
  gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", gene_symbol]), nm)
ECM <- intersect(c("COL1A1","COL3A1","COL7A1","COL18A1","FN1","SPARC","POSTN","LOX","MMP2",
  "MMP14","THBS1","THBS2","TGFB2","TGFBI","SERPINE1","ACTA2","TAGLN","VIM","CDH11","FBN1",
  "CTGF","PDGFRB","ITGA5","ITGB1"), nm)
cat("igraph", as.character(packageVersion("igraph")), "| mclust",
    as.character(packageVersion("mclust")), "| MCL pkg",
    as.character(packageVersion("MCL")), "| R", R.version.string, "\n\n")

## cross-implementation agreement for the same algorithm
same <- data.table(
  algorithm = c("MCL I=2.2", "MCL I=2.5", "Louvain", "Fast-greedy", "Infomap",
                "Label propagation", "Leading eigenvector"),
  igraph_run = c("MCL_I2.2","MCL_I2.5","Louvain_r1.00","FastGreedy","Infomap_t10",
                 "LabelProp_best","LeadingEigen"),
  cm2_run = c("cm2_mcl_I22","cm2_mcl_I25","cm2_louvain_cm2","cm2_fastgreedy_cm2",
              "cm2_infomap_cm2_t10","cm2_labelprop_cm2","cm2_leadeig_cm2"))
same[, ARI := mapply(function(a,b) round(adjustedRandIndex(P[[a]], P[[b]]), 3), igraph_run, cm2_run)]
same[, k_igraph := sapply(igraph_run, function(r) s[run==r, k_clusters])]
same[, k_cm2    := sapply(cm2_run,    function(r) s[run==r, k_clusters])]
fwrite(same, "results/v7/cm2_cross_implementation_ARI.csv")
cat("== same algorithm, two implementations ==\n"); print(same)

## cross-algorithm ARI among the informative modularity-based partitions
core <- c("MCL_I2.5","cm2_mcl_I25","Louvain_r1.00","cm2_louvain_cm2","Walktrap_s4",
          "Infomap_t10","LabelProp_best","FastGreedy","cm2_fastgreedy_cm2","cm2_glay_cm2",
          "LeadingEigen","Leiden_mod_r1.0","EdgeBetweenness")
A <- outer(core, core, Vectorize(function(a,b) adjustedRandIndex(P[[a]], P[[b]])))
dimnames(A) <- list(core, core)
off <- A[upper.tri(A)]
cat(sprintf("\ncross-algorithm ARI over %d runs: median %.3f, IQR %.3f-%.3f, max %.3f\n",
            length(core), median(off), quantile(off, .25), quantile(off, .75), max(off)))
## excluding same-algorithm pairs
pairsame <- outer(sub("^cm2_","",core), sub("^cm2_","",core), "==")
off2 <- A[upper.tri(A) & !pairsame]
cat(sprintf("excluding same-algorithm pairs: median %.3f, max %.3f, n pairs %d\n",
            median(off2), max(off2), length(off2)))

## leading modules of the two reference partitions
qc <- function(m) { mm <- ecount(g); Ad <- as_adjacency_matrix(g, sparse=TRUE); dt <- degree(g)
  sapply(sort(unique(m)), function(cl){ i <- which(m==cl); sum(Ad[i,i])/2/mm - (sum(dt[i])/(2*mm))^2 }) }
for (k in c("Louvain_r1.00","cm2_louvain_cm2","Leiden_mod_r1.0")) {
  m <- P[[k]]; q <- qc(m); tb <- table(m)
  d <- data.table(cluster=as.integer(names(tb)), size=as.integer(tb),
                  q=round(q[match(as.integer(names(tb)), sort(unique(m)))],4))
  d[, `:=`(nG=sapply(cluster,function(c)sum(ntype[names(m)[m==c]]=="Gene")),
           nTF=sapply(cluster,function(c)sum(ntype[names(m)[m==c]]=="TF")),
           nmiR=sapply(cluster,function(c)sum(ntype[names(m)[m==c]]=="miRNA")),
           ECM=sapply(cluster,function(c)length(intersect(names(m)[m==c],ECM))),
           EMT=sapply(cluster,function(c)length(intersect(names(m)[m==c],EMT))),
           p30=sapply(cluster,function(c)length(intersect(names(m)[m==c],prior30))),
           hubs=sapply(cluster,function(c){mem<-names(m)[m==c]
             paste(names(sort(degree(g,mem),decreasing=TRUE))[1:6],collapse=", ")}),
           hasCOL=sapply(cluster,function(c)"COL1A1"%in%names(m)[m==c]))]
  setorder(d,-size)
  cat(sprintf("\n== leading modules, %s (Q=%.3f) ==\n", k, modularity(g,m)))
  print(d, nrows=30)
  fwrite(d, sprintf("results/v7/cm2_modules_%s.csv", k))
}
