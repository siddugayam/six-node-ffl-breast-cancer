#!/usr/bin/env Rscript
# v7 / 11 : one summary figure for the convergent-validity test.
suppressMessages({library(igraph); library(data.table); library(mclust)})
setwd("/path/to/revision")
dir.create("figures/v7", showWarnings = FALSE, recursive = TRUE)
g <- readRDS("results/v7/cm2_graph_undirected.rds"); nm <- V(g)$name
P <- readRDS("results/v7/cm2_all_partitions.rds"); s <- fread("results/v7/cm2_run_summary.csv")
tab <- fread("results/v7/cm2_TABLE_method_battery.csv")

core <- c(MCL="MCL_I2.5", `MCL (cM2)`="cm2_mcl_I25", Louvain="Louvain_r1.00",
          `Louvain (cM2)`="cm2_louvain_cm2", Walktrap="Walktrap_s4", Infomap="Infomap_t10",
          `LabelProp`="LabelProp_best", FastGreedy="FastGreedy",
          `FastGreedy (cM2)`="cm2_fastgreedy_cm2", `GLay (cM2)`="cm2_glay_cm2",
          LeadEigen="LeadingEigen", Leiden="Leiden_mod_r1.0", `Girvan-Newman`="EdgeBetweenness")
A <- outer(core, core, Vectorize(function(a,b) adjustedRandIndex(P[[a]], P[[b]])))
dimnames(A) <- list(names(core), names(core))

draw <- function() {
par(mfrow = c(2,2), mar = c(5,5,4,2), cex.lab = 1.05)

## A ---- modularity vs k
s[, alg := sub("_.*","", run)]
cols <- setNames(rainbow(length(unique(s$alg)), v=.8), sort(unique(s$alg)))
plot(s$k_clusters, s$modularity, log="x", pch=ifelse(s$degenerate,1,19),
     col=cols[s$alg], xlab="number of clusters (log)", ylab="modularity Q",
     main="A  61 runs of the battery on one network")
abline(h=0.3, lty=2, col="grey40"); text(2, 0.31, "Q = 0.3", pos=4, cex=.8, col="grey30")
legend("topright", legend=names(cols), col=cols, pch=19, cex=.65, ncol=2, bty="n")
mtext("open = degenerate (k<3, giant >60%, or all singletons)", side=3, line=0.1, cex=.65)

## B ---- ARI heatmap
par(mar=c(9,9,4,2))
image(seq_len(nrow(A)), seq_len(ncol(A)), A[, ncol(A):1], col=hcl.colors(50,"Blues",rev=TRUE),
      zlim=c(0,1), axes=FALSE, xlab="", ylab="", main="B  adjusted Rand index between methods")
axis(1, at=seq_len(nrow(A)), labels=rownames(A), las=2, cex.axis=.7)
axis(2, at=seq_len(ncol(A)), labels=rev(colnames(A)), las=2, cex.axis=.7)
for (i in seq_len(nrow(A))) for (j in seq_len(ncol(A)))
  text(i, ncol(A)-j+1, sprintf("%.2f", A[i,j]), cex=.5,
       col=ifelse(A[i,j] > .5, "white", "grey20"))

## C ---- consensus co-clustering of the stroma set
co <- fread("results/v7/cm2_stroma_coclustering.csv"); rn <- co$node
C <- as.matrix(co[,-1]); rownames(C) <- rn
ord <- hclust(as.dist(1-C))$order; C <- C[ord, ord]
par(mar=c(9,9,4,2))
image(seq_len(nrow(C)), seq_len(ncol(C)), C[, ncol(C):1], col=hcl.colors(50,"Reds",rev=TRUE),
      zlim=c(0,1), axes=FALSE, xlab="", ylab="",
      main="C  co-clustering of ECM genes + miR-29 family")
axis(1, at=seq_len(nrow(C)), labels=rownames(C), las=2, cex.axis=.55)
axis(2, at=seq_len(ncol(C)), labels=rev(colnames(C)), las=2, cex.axis=.55)

## D ---- COL1A1 module
par(mar=c(5,14,4,2))
t2 <- tab[order(COL1A1_module_size)]
bp <- barplot(t2$COL1A1_module_size, horiz=TRUE, names.arg=t2$method, las=1, cex.names=.6,
              col=ifelse(t2$degenerate, "grey80", ifelse(t2$miR29a_with_COL1A1,"#c0392b","#2980b9")),
              xlab="size of the module containing COL1A1 (nodes)",
              main="D  no method isolates a stroma module")
abline(v=587, lty=3); text(587, max(bp), "whole network", srt=90, pos=2, cex=.6)
legend("bottomright", c("COL1A1 + COL3A1 + miR-29a together","COL1A1 alone/without miR-29a",
                        "degenerate partition"),
       fill=c("#c0392b","#2980b9","grey80"), cex=.6, bty="n")
}
pdf("figures/v7/fig_v7_clustermaker2_battery.pdf", width = 12, height = 10); draw(); invisible(dev.off())
png("figures/v7/fig_v7_clustermaker2_battery.png", width = 2400, height = 2000, res = 190); draw(); invisible(dev.off())
cat("wrote figures/v7/fig_v7_clustermaker2_battery.pdf\n")
