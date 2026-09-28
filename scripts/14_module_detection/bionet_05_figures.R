## v7 / BioNet step 5: supplementary figure (print, light surface)
suppressMessages({library(BioNet); library(igraph); library(data.table); library(graph)})
ROOT <- "/path/to/revision"
OUT <- file.path(ROOT,"results/v7"); FIG <- file.path(ROOT,"figures/v7")
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
P <- function(...) file.path(OUT, paste0("bionet_", ...))

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUT<-"#8a8983"
BLUE<-"#2a78d6"; RED<-"#e34948"; ORNG<-"#eb6834"; AQUA<-"#1baf7a"; GRIDC<-"#e6e5e1"

gp <- readRDS(P("graph_and_pvals.rds")); nod <- gp$nod; gi <- gp$gi
fits <- readRDS(P("bum_fits.rds")); fb <- fits$gene_gw
scan <- fread(P("fdr_scan.csv")); S <- fread(P("module_summary.csv"))
mem  <- fread(P("module_primary_members.csv")); mm <- fread(P("fdr_scan_with_members.csv"))
enr  <- fread(P("primary_module_signature_enrichment.csv"))
pri30 <- fread(file.path(ROOT,"results/v5/tables/Table4_prioritised_30.csv"))$name
PFDR <- S[run=="primary"]$fdr

hdr <- function(a,b){ title(main=a, adj=0, col.main=INK, font.main=2, cex.main=1.05, line=1.7)
                      mtext(b, adj=0, line=0.5, cex=0.72, col=INK2) }

pdf(file.path(FIG,"bionet_supplementary.pdf"), width=11, height=9.5)
par(bg=SURF, mfrow=c(2,2), mar=c(4.8,5.0,3.4,1.6), col.axis=INK2, col.lab=INK2,
    fg=MUT, family="sans", las=1, cex.axis=0.85)

## -- A: BUM fit ------------------------------------------------------------
h <- hist(fb$pvalues, breaks=50, plot=FALSE)
plot(h, freq=FALSE, border=NA, col="#dfe9f7", main="", xlab="DE p-value (TCGA-BRCA, 20,250 genes)",
     ylab="density", ylim=c(0,6), xaxs="i")
xs <- seq(1e-4,1,length.out=4000)
lines(xs, fb$lambda+(1-fb$lambda)*fb$a*xs^(fb$a-1), col=BLUE, lwd=2)
pi0 <- fb$lambda+(1-fb$lambda)*fb$a
abline(h=pi0, col=ORNG, lwd=2, lty=2)
text(0.30, 3.1, "fitted beta-uniform mixture", col=BLUE, cex=0.8, pos=4)
text(0.30, pi0+0.45, sprintf("null component, %.1f%%", 100*pi0), col=ORNG, cex=0.8, pos=4)
text(0.055, 5.6, "first bin clipped (density 37)", col=INK2, cex=0.7, pos=4)
hdr("A  Beta-uniform mixture fit",
    sprintf("lambda = %.3f, a = %.3f; at most %.0f%% of the transcriptome is null", fb$lambda, fb$a, 100*pi0))

## -- B: FDR vs module size -------------------------------------------------
plot(scan$fdr, scan$module_n, log="x", type="n", xlab="FDR used for node scoring",
     ylab="nodes in maximum-scoring module", ylim=c(0,660), xlim=c(1e-104,1), panel.first={
       abline(h=seq(0,600,100), col=GRIDC); abline(v=10^seq(-100,0,20), col=GRIDC)})
lastT <- function(col) min(scan$fdr[scan[[col]]])   # smallest FDR at which the node is still in
for(k in list(c("has_COL3A1","COL3A1"), c("has_COL1A1","COL1A1"), c("has_miR29a","miR-29a"))){
  f <- lastT(k[1]); abline(v=f, col=RED, lty=3)
}
lines(scan$fdr, scan$module_n, col=BLUE, lwd=2)
points(scan$fdr, scan$module_n, col=BLUE, pch=19, cex=0.5)
abline(h=587, col=MUT, lty=3); text(1e-101, 587, "whole network (587)", pos=4, cex=0.72, col=INK2)
abline(h=58, col=AQUA, lty=2)
points(PFDR, 58, pch=21, bg=AQUA, col=SURF, cex=1.7, lwd=2)
text(PFDR, 58, "primary\nFDR 1.6e-54\n58 nodes", pos=3, cex=0.72, col=INK, offset=0.7)
points(1e-3, 444, pch=21, bg=ORNG, col=SURF, cex=1.6, lwd=2)
text(1e-19, 500, "BioNet vignette default\nFDR 1e-3: 444 nodes (76%)", pos=2, cex=0.72, col=INK)
text(lastT("has_COL1A1"), 95, "COL1A1 lost below", srt=90, adj=c(0,-0.45), cex=0.67, col=RED)
text(lastT("has_COL3A1"), 95, "COL3A1 lost below", srt=90, adj=c(0,-0.45), cex=0.67, col=RED)
text(lastT("has_miR29a"), 95, "miR-29a lost below", srt=90, adj=c(0,-0.45), cex=0.67, col=RED)
hdr("B  Module size is set by the chosen FDR, not by the data",
    "no FDR gives both an interpretable size and the collagen axis")

## -- C: signature enrichment of the primary module -------------------------
pick <- c("FFL_3node","FFL_higher_order_only","MS_exemplar_module","NABA_CORE_MATRISOME",
          "NABA_COLLAGENS","FARMER_STROMAL","HALLMARK_EMT","ESTIMATE_STROMAL","TRIULZI_ECM",
          "HALLMARK_TGF_BETA","CHANG_CSR_LATE_UP","FINAK_SDPP")
e <- enr[signature %in% pick][match(pick, signature)]
e <- e[!is.na(signature)]; e[, lf := log2(pmax(fold, 0.06))]
setorder(e, lf); y <- 1:nrow(e)
par(mar=c(4.8,12.5,3.4,1.6))
plot(NA, xlim=c(-5.4,2.9), ylim=c(0.4,nrow(e)+0.6), yaxt="n", xlab="log2 fold enrichment in the 58-node module",
     ylab="", panel.first={abline(v=0, col=MUT); abline(h=y, col=GRIDC)})
segments(0, y, e$lf, y, col=ifelse(e$fold>=1, RED, BLUE), lwd=2)
points(e$lf, y, pch=19, cex=1.3, col=ifelse(e$fold>=1, RED, BLUE))
axis(2, at=y, labels=gsub("_"," ",e$signature), cex.axis=0.72, tick=FALSE)
text(e$lf, y, sprintf(" %d/%d", e$n_overlap, e$n_ref), pos=ifelse(e$fold>=1,4,2), cex=0.66, col=INK2)
text(-5.35, nrow(e)+0.45, "depleted", col=BLUE, cex=0.72, pos=4)
text(2.85, nrow(e)+0.45, "enriched", col=RED, cex=0.72, pos=2)
hdr("C  What the module is made of",
    "universe = 364 gene/TF network nodes; only FFL 3-node reaches q < 0.05")

## -- D: fate of the stroma axis -------------------------------------------
par(mar=c(4.8,8.2,3.4,1.6))
kn <- c("COL1A1","COL3A1","FN1","SPARC","POSTN","THBS2","hsa-miR-29a","hsa-miR-101","EZH2")
inmat <- sapply(kn, function(g) sapply(strsplit(mm$members,";"), function(v) g %in% v))
plot(NA, xlim=c(-104,2), ylim=c(0.5,length(kn)+0.7), yaxt="n", xlab="log10 FDR", ylab="",
     panel.first=abline(h=1:length(kn), col=GRIDC))
axis(2, at=1:length(kn), labels=rev(kn), cex.axis=0.8, tick=FALSE)
for(i in seq_along(kn)){ yy <- length(kn)-i+1
  points(log10(mm$fdr), rep(yy, nrow(mm)), pch=19, cex=1.1,
         col=ifelse(inmat[,i], BLUE, "#ebeae7")) }
abline(v=log10(PFDR), col=AQUA, lwd=2, lty=2)
text(log10(PFDR), length(kn)+0.6, "primary FDR", col=AQUA, cex=0.72, pos=2)
hdr("D  Membership of the paper's key nodes across FDR",
    "filled = present in the maximum-scoring module at that FDR")

## -- page 2: the module ----------------------------------------------------
par(mfrow=c(1,1), mar=c(0.5,0.5,3.6,0.5))
sub <- induced_subgraph(gi, mem$name)
V(sub)$lfc <- mem$logFC[match(V(sub)$name, mem$name)]
V(sub)$pri <- V(sub)$name %in% pri30
V(sub)$tp  <- mem$type[match(V(sub)$name, mem$name)]
set.seed(11); L <- layout_with_fr(sub, niter=4000)
plot(sub, layout=L, rescale=TRUE, vertex.size=ifelse(V(sub)$pri, 6.5, 4.5),
     vertex.color=ifelse(V(sub)$lfc>0, RED, BLUE),
     vertex.frame.color=ifelse(V(sub)$pri, INK, SURF), vertex.frame.width=ifelse(V(sub)$pri,2,1),
     vertex.shape=ifelse(V(sub)$tp=="miRNA","square","circle"),
     vertex.label=V(sub)$name, vertex.label.cex=0.62, vertex.label.color=INK,
     vertex.label.dist=0.95, vertex.label.degree=-pi/2, vertex.label.family="sans",
     edge.color="#dcdbd7", edge.width=0.8, asp=0.72)
hdr("E  The 58-node maximum-scoring module (BioNet, FDR = 1.6e-54)",
    "red = up in tumour, blue = down; square = miRNA; black ring = one of the 30 prioritised nodes. No collagen node is present.")
dev.off()
cat("figure written:", file.path(FIG,"bionet_supplementary.pdf"), "\n")
