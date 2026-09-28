## ---------------------------------------------------------------------------
## v7 / BioNet  step 3: sensitivity of the maximum-scoring module to every
## analyst choice (BUM fitting set, miRNA handling, missing p-values, size rule)
## ---------------------------------------------------------------------------
suppressMessages({library(BioNet); library(igraph); library(data.table); library(graph)})
set.seed(1)
ROOT <- "/path/to/revision"
OUT  <- file.path(ROOT,"results/v7"); P <- function(...) file.path(OUT, paste0("bionet_", ...))
gp <- readRDS(P("graph_and_pvals.rds")); nod <- gp$nod; G <- gp$G; gi <- gp$gi
fits <- readRDS(P("bum_fits.rds"))
tauf <- function(fb, fdr) BioNet:::fdrThreshold(fdr, fb)
pri30 <- fread(file.path(ROOT,"results/v5/tables/Table4_prioritised_30.csv"))$name
sets  <- readRDS(file.path(ROOT,"cache/v5/farmer_signature_sets.rds"))
ecm_core <- intersect(nod$name, unique(c(sets$NABA_COLLAGENS, sets$NABA_CORE_MATRISOME)))

mkScores <- function(fdr, fbG, fbM, dat){
  tg <- tauf(fbG,fdr); tm <- tauf(fbM,fdr)
  setNames(ifelse(dat$is_mir, (fbM$a-1)*(log(dat$pval)-log(tm)),
                              (fbG$a-1)*(log(dat$pval)-log(tg))), dat$name)
}
runMod <- function(fdr, fbG, fbM, graph, dat){
  sc <- mkScores(fdr, fbG, fbM, dat)[graph::nodes(graph)]
  m  <- runFastHeinz(graph, sc)
  if(is.null(m)) character(0) else graph::nodes(m)
}
hyp <- function(mb, ref, uni){
  x<-length(intersect(mb,ref)); m<-length(intersect(ref,uni)); k<-length(intersect(mb,uni)); N<-length(uni)
  c(n_ov=x, exp=k*m/N, fold=if(x==0) 0 else (x/k)/(m/N), p=phyper(x-1,m,N-m,k,lower.tail=FALSE))
}
describe <- function(lab, mb, uni){
  h30 <- hyp(mb,pri30,uni); he <- hyp(mb,ecm_core,uni)
  data.table(scenario=lab, module_n=length(mb), frac_net=length(mb)/length(uni),
    n_miRNA=sum(mb %in% nod[is_mir==TRUE]$name),
    ov_pri30=h30["n_ov"], exp_pri30=round(h30["exp"],2), fold_pri30=round(h30["fold"],2), p_pri30=h30["p"],
    n_ecm=he["n_ov"], fold_ecm=round(he["fold"],2), p_ecm=he["p"],
    COL1A1="COL1A1" %in% mb, COL3A1="COL3A1" %in% mb, FN1="FN1" %in% mb,
    SPARC="SPARC" %in% mb, POSTN="POSTN" %in% mb, THBS2="THBS2" %in% mb,
    EZH2="EZH2" %in% mb, miR29a="hsa-miR-29a" %in% mb, miR101="hsa-miR-101" %in% mb,
    members=paste(sort(mb), collapse=";"))
}

## reference: the primary run's FDR
prim_fdr <- fread(P("module_summary.csv"))[run=="primary"]$fdr
cat("primary FDR from step 2:", format(prim_fdr, digits=3), "\n")
uni_all <- nod$name
R <- list()

## ---- S0 reference ---------------------------------------------------------
R[["S0 primary (genome-wide gene BUM + all-miRNA BUM, p=1 for 11 unmeasured)"]] <-
  runMod(prim_fdr, fits$gene_gw, fits$mir_all, G, nod)

## ---- S1 BUM fitted on network nodes only ---------------------------------
R[["S1 BUM fitted on the 587 network nodes only"]] <-
  runMod(prim_fdr, fits$gene_net, fits$mir_net, G, nod)

## ---- S2 pooled BUM for genes and miRNAs ----------------------------------
R[["S2 single pooled BUM for genes+miRNAs"]] <-
  runMod(prim_fdr, fits$pooled, fits$pooled, G, nod)

## ---- S3 gene/TF-only subnetwork (miRNAs deleted) -------------------------
gi_g <- induced_subgraph(gi, V(gi)$name[!(V(gi)$name %in% nod[is_mir==TRUE]$name)])
gi_g <- induced_subgraph(gi_g, V(gi_g)$name[components(gi_g)$membership==which.max(components(gi_g)$csize)])
Gg   <- rmSelfLoops(igraph::as_graphnel(gi_g)); nodg <- nod[name %in% graph::nodes(Gg)]
cat("gene/TF-only subnetwork: V=",numNodes(Gg)," E=",numEdges(Gg),"\n")
R[["S3 gene/TF-only subnetwork (miRNAs removed)"]] <- runMod(prim_fdr, fits$gene_gw, fits$mir_all, Gg, nodg)

## ---- S4 drop the 11 nodes with no measurement ----------------------------
keep <- nod[has_p==TRUE]$name
gi_k <- induced_subgraph(gi, keep)
gi_k <- induced_subgraph(gi_k, V(gi_k)$name[components(gi_k)$membership==which.max(components(gi_k)$csize)])
Gk   <- rmSelfLoops(igraph::as_graphnel(gi_k)); nodk <- nod[name %in% graph::nodes(Gk)]
R[["S4 unmeasured nodes deleted rather than scored at p=1"]] <- runMod(prim_fdr, fits$gene_gw, fits$mir_all, Gk, nodk)

## ---- S5 alternative size rules -------------------------------------------
msz <- function(f, fbG=fits$gene_gw, fbM=fits$mir_all) length(runMod(f, fbG, fbM, G, nod))
bisect <- function(tmax, lo=-110, hi=-2){
  while((hi-lo)>0.05){ mid<-(lo+hi)/2; if(msz(10^mid)<=tmax) lo<-mid else hi<-mid }
  10^lo
}
f30  <- bisect(30);  f100 <- bisect(100)
cat(sprintf("size rule <=30 nodes -> FDR %.3g ; <=100 nodes -> FDR %.3g\n", f30, f100))
R[[sprintf("S5a size rule <=30 nodes (FDR %.2g)", f30)]]   <- runMod(f30,  fits$gene_gw, fits$mir_all, G, nod)
R[[sprintf("S5b size rule <=100 nodes (FDR %.2g)", f100)]] <- runMod(f100, fits$gene_gw, fits$mir_all, G, nod)

out <- rbindlist(lapply(names(R), function(k) describe(k, R[[k]],
        if(grepl("^S3", k)) nodg$name else if(grepl("^S4", k)) nodk$name else uni_all)))
fwrite(out, P("sensitivity.csv"))
print(out[, .(scenario, module_n, n_miRNA, ov_pri30, exp_pri30, fold_pri30, p_pri30,
              n_ecm, fold_ecm, p_ecm, COL1A1, COL3A1, EZH2, miR29a, miR101)])

## ---- S6: how much of the module is just "the k most significant nodes"? --
sc0 <- mkScores(prim_fdr, fits$gene_gw, fits$mir_all, nod)
mb0 <- R[[1]]; k <- length(mb0)
topk <- names(sort(sc0, decreasing=TRUE))[1:k]
cat(sprintf("\nJaccard(primary module, top-%d nodes by node score) = %.3f ; |intersect| = %d/%d\n",
            k, length(intersect(mb0,topk))/length(union(mb0,topk)), length(intersect(mb0,topk)), k))
fwrite(data.table(topk_by_score=topk, in_module=topk %in% mb0), P("primary_vs_topk_by_score.csv"))
