## ---------------------------------------------------------------------------
## v7 / BioNet  (Beisser et al. 2010, Bioconductor 1.68.0)
## Step 1: build undirected graph, attach DE p-values, fit beta-uniform mixture,
##         scan FDR -> module size.
## ---------------------------------------------------------------------------
suppressMessages({library(BioNet); library(igraph); library(data.table); library(graph)})
set.seed(1)
ROOT <- "/path/to/revision"
OUT  <- file.path(ROOT,"results/v7"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
FIG  <- file.path(ROOT,"figures/v7"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
P <- function(...) file.path(OUT, paste0("bionet_", ...))

## ---- 1. network -----------------------------------------------------------
E <- fread(file.path(ROOT,"data/canonical_edges.tsv"))
nodes <- unique(rbind(E[,.(name=source,type=source_type)], E[,.(name=target,type=target_type)]), by="name")
setorder(nodes, name)
gi <- simplify(graph_from_data_frame(E[source!=target,.(source,target)], directed=FALSE,
                                     vertices=nodes), remove.multiple=TRUE, remove.loops=TRUE)
G  <- igraph::as_graphnel(gi)                 # BioNet works on graphNEL, undirected
G  <- rmSelfLoops(G)
cat("graph: V=", numNodes(G), " E=", numEdges(G), " components=", components(gi)$no, "\n")

## ---- 2. p-values ----------------------------------------------------------
DEg <- fread(file.path(ROOT,"results/v2/v2_DE_genes.csv"))
DEm <- fread(file.path(ROOT,"results/v2/v2_DE_mirnas.csv"))
pg_all <- setNames(DEg$P.Value, DEg$feature)       # 20250 genome-wide, as specified
pm_all <- setNames(DEm$P.Value, DEm$feature)       # 495 miRNAs, same limma pipeline

nod <- copy(nodes); nod[, is_mir := type=="miRNA"]
nod[, pval := ifelse(is_mir, pm_all[name], pg_all[name])]
nod[, has_p := !is.na(pval)]
cat("nodes with p:", sum(nod$has_p), "/", nrow(nod),
    " (missing:", paste(nod[has_p==FALSE]$name, collapse=","), ")\n")
nod[has_p==FALSE, pval := 1]     # neutral / maximally non-significant, keeps topology intact

## ---- 3. beta-uniform mixture fits ----------------------------------------
fitsafe <- function(p, nm){
  fb <- fitBumModel(p, plot=FALSE)
  cat(sprintf("BUM %-22s n=%5d  lambda=%.4f  a=%.4f  negLL=%.2f\n", nm, length(p), fb$lambda, fb$a, fb$negLL))
  fb
}
fb_gene_gw  <- fitsafe(pg_all, "genes genome-wide")                       # PRIMARY (as pre-specified)
fb_gene_net <- fitsafe(pg_all[intersect(names(pg_all), nod[is_mir==FALSE]$name)], "genes network-only")
fb_mir_all  <- fitsafe(pm_all, "miRNAs all")
fb_mir_net  <- fitsafe(pm_all[intersect(names(pm_all), nod[is_mir==TRUE]$name)], "miRNAs network-only")
fb_pooled   <- fitsafe(c(pg_all, pm_all), "genes+miRNAs pooled")

saveRDS(list(gene_gw=fb_gene_gw, gene_net=fb_gene_net, mir_all=fb_mir_all,
             mir_net=fb_mir_net, pooled=fb_pooled), P("bum_fits.rds"))

## goodness of fit: BioNet's own plot + a KS-style check on the uniform part
pdf(file.path(FIG,"bionet_bum_fits.pdf"), width=9, height=6)
par(mfrow=c(2,3))
for(nm in c("gene_gw","gene_net","mir_all","mir_net","pooled")){
  fb <- list(gene_gw=fb_gene_gw, gene_net=fb_gene_net, mir_all=fb_mir_all,
             mir_net=fb_mir_net, pooled=fb_pooled)[[nm]]
  hist(fb, main=nm)
}
dev.off()

## ---- 4. tau (score threshold) as a function of FDR ------------------------
## Beisser et al.: tau = ( (lambda - fdr*lambda) / (fdr - fdr*lambda) )^(1/(a-1)) ... via BioNet:::fdrThreshold
tauf <- function(fb, fdr) BioNet:::fdrThreshold(fdr, fb)

## ---- 5. FDR scan ----------------------------------------------------------
fdrs <- c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15,
          1e-18,1e-20,1e-25,1e-30,1e-40,1e-50,1e-60,1e-80,1e-100)

scoreAll <- function(fdr, fbG, fbM){
  tg <- tauf(fbG, fdr); tm <- tauf(fbM, fdr)
  s <- numeric(nrow(nod))
  ig <- !nod$is_mir; im <- nod$is_mir
  s[ig] <- (fbG$a - 1) * (log(nod$pval[ig]) - log(tg))
  s[im] <- (fbM$a - 1) * (log(nod$pval[im]) - log(tm))
  setNames(s, nod$name)
}

scan_res <- rbindlist(lapply(fdrs, function(f){
  sc <- scoreAll(f, fb_gene_gw, fb_mir_all)
  sc <- sc[nodes(G)]
  mod <- runFastHeinz(G, sc)
  mn  <- if(is.null(mod)) character(0) else nodes(mod)
  data.table(fdr=f,
             tau_gene=tauf(fb_gene_gw,f), tau_mir=tauf(fb_mir_all,f),
             n_pos=sum(sc>0), n_pos_gene=sum(sc[nod[match(names(sc),nod$name)]$is_mir==FALSE]>0),
             n_pos_mir=sum(sc[nod[match(names(sc),nod$name)]$is_mir==TRUE]>0),
             module_n=length(mn), module_score=sum(sc[mn]),
             n_mir=sum(mn %in% nod[is_mir==TRUE]$name),
             n_neg_in_mod=sum(sc[mn]<0),
             has_COL1A1="COL1A1" %in% mn, has_COL3A1="COL3A1" %in% mn,
             has_FN1="FN1" %in% mn, has_EZH2="EZH2" %in% mn,
             has_miR29a="hsa-miR-29a" %in% mn, has_miR101="hsa-miR-101" %in% mn,
             members=paste(sort(mn), collapse=";"))
}))
fwrite(scan_res[, !"members"], P("fdr_scan.csv"))
fwrite(scan_res, P("fdr_scan_with_members.csv"))
print(scan_res[, .(fdr, tau_gene, tau_mir, n_pos, module_n, n_mir, n_neg_in_mod, module_score,
                   has_COL1A1, has_COL3A1, has_EZH2, has_miR29a, has_miR101)])
saveRDS(list(nod=nod, G=G, gi=gi), P("graph_and_pvals.rds"))
