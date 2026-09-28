## ---------------------------------------------------------------------------
## v7 / BioNet  step 2: primary + sensitivity maximum-scoring modules
## ---------------------------------------------------------------------------
suppressMessages({library(BioNet); library(igraph); library(data.table); library(graph)})
set.seed(1)
ROOT <- "/path/to/revision"
OUT  <- file.path(ROOT,"results/v7"); FIG <- file.path(ROOT,"figures/v7")
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
P <- function(...) file.path(OUT, paste0("bionet_", ...))

gp   <- readRDS(P("graph_and_pvals.rds")); nod <- gp$nod; G <- gp$G; gi <- gp$gi
fits <- readRDS(P("bum_fits.rds"))
tauf <- function(fb, fdr) BioNet:::fdrThreshold(fdr, fb)

DEg <- fread(file.path(ROOT,"results/v2/v2_DE_genes.csv"))
DEm <- fread(file.path(ROOT,"results/v2/v2_DE_mirnas.csv"))
lfc <- c(setNames(DEg$logFC, DEg$feature), setNames(DEm$logFC, DEm$feature))
qv  <- c(setNames(DEg$adj.P.Val, DEg$feature), setNames(DEm$adj.P.Val, DEm$feature))

pri30 <- fread(file.path(ROOT,"results/v5/tables/Table4_prioritised_30.csv"))$name
stopifnot(length(pri30)==30, all(pri30 %in% nod$name))

## ---- scoring helper -------------------------------------------------------
## Beisser/Dittrich score:  s(p) = (a-1) * ( log p - log tau(fdr) )
mkScores <- function(fdr, fbG, fbM, dat=nod){
  tg <- tauf(fbG, fdr); tm <- tauf(fbM, fdr)
  s <- ifelse(dat$is_mir, (fbM$a-1)*(log(dat$pval)-log(tm)), (fbG$a-1)*(log(dat$pval)-log(tg)))
  setNames(s, dat$name)
}
runMod <- function(fdr, fbG, fbM, graph=G, dat=nod){
  sc <- mkScores(fdr, fbG, fbM, dat)
  sc <- sc[graph::nodes(graph)]
  m  <- runFastHeinz(graph, sc)
  list(fdr=fdr, scores=sc, module=m, members=if(is.null(m)) character(0) else graph::nodes(m))
}

## ---- pre-specified size rule: module <= 10% of network nodes (<= 58) ------
target_max <- floor(0.10*numNodes(G))
msize <- function(f) length(runMod(f, fits$gene_gw, fits$mir_all)$members)
## bracket from the coarse scan (1e-50 -> 67 nodes, 1e-60 -> 45 nodes), then bisect on log10
lo <- -60; hi <- -50; trace <- list()
sz_lo <- msize(10^lo); sz_hi <- msize(10^hi)
trace[[length(trace)+1]] <- data.table(log10fdr=c(lo,hi), module_n=c(sz_lo,sz_hi))
stopifnot(sz_lo <= target_max, sz_hi > target_max)
while((hi-lo) > 0.05){
  mid <- (lo+hi)/2; s_mid <- msize(10^mid)
  trace[[length(trace)+1]] <- data.table(log10fdr=mid, module_n=s_mid)
  if(s_mid <= target_max) lo <- mid else hi <- mid
}
primary_fdr <- 10^lo
cat(sprintf("size rule: module <= %d nodes (10%% of %d); primary FDR = %.3g -> %d nodes\n",
            target_max, numNodes(G), primary_fdr, msize(primary_fdr)))
fwrite(rbindlist(trace)[order(-log10fdr)], P("fdr_bisection_trace.csv"))

## ---- the three reported runs ---------------------------------------------
runs <- list(
  primary     = list(fdr=primary_fdr, lab=sprintf("PRIMARY  FDR=%.2g (size rule <=%d nodes)", primary_fdr, target_max)),
  sens_vign   = list(fdr=1e-3,        lab="SENSITIVITY  FDR=1e-3 (BioNet vignette default)"),
  sens_mid    = list(fdr=1e-25,       lab="SENSITIVITY  FDR=1e-25 (intermediate)")
)

hyper <- function(members, ref, universe){
  x <- length(intersect(members, ref)); m <- length(intersect(ref, universe))
  k <- length(intersect(members, universe)); n <- length(universe)-m
  data.table(n_module=k, n_ref=m, n_overlap=x, expected=k*m/length(universe),
             fold=(x/k)/(m/length(universe)),
             p_hyper=phyper(x-1, m, n, k, lower.tail=FALSE),
             p_depletion=phyper(x, m, n, k, lower.tail=TRUE))
}

universe <- nod$name
sets <- readRDS(file.path(ROOT,"cache/v5/farmer_signature_sets.rds"))
ecm_core <- sort(intersect(nod$name, unique(c(sets$NABA_COLLAGENS, sets$NABA_CORE_MATRISOME))))
cat("core-matrisome/collagen nodes present in network (n=",length(ecm_core),"): ",
    paste(ecm_core, collapse=", "), "\n", sep="")

res_members <- list(); summ <- list()
for(nm in names(runs)){
  r  <- runMod(runs[[nm]]$fdr, fits$gene_gw, fits$mir_all)
  mb <- r$members
  dt <- data.table(name=mb)[nod, on="name", nomatch=0][
        , .(name, type, pval, score=r$scores[name], logFC=lfc[name], adj_P=qv[name],
            in_prioritised30 = name %in% pri30,
            in_NABA_collagen = name %in% sets$NABA_COLLAGENS,
            in_core_matrisome= name %in% sets$NABA_CORE_MATRISOME,
            in_FARMER_STROMAL= name %in% sets$FARMER_STROMAL,
            in_HALLMARK_EMT  = name %in% sets$HALLMARK_EMT)]
  setorder(dt, -score)
  dt[, score_rank := 1:.N]
  fwrite(dt, P("module_", nm, "_members.csv"))
  res_members[[nm]] <- dt

  h30 <- hyper(mb, pri30, universe)
  hec <- hyper(mb, ecm_core, universe)
  summ[[nm]] <- data.table(run=nm, label=runs[[nm]]$lab, fdr=runs[[nm]]$fdr,
      module_n=length(mb), frac_of_network=length(mb)/length(universe),
      n_gene=sum(dt$type=="Gene"), n_TF=sum(dt$type=="TF"), n_miRNA=sum(dt$type=="miRNA"),
      module_score=sum(r$scores[mb]), n_negative_members=sum(r$scores[mb]<0),
      overlap_pri30=h30$n_overlap, expected_pri30=h30$expected, fold_pri30=h30$fold, p_pri30=h30$p_hyper,
      n_ecm=hec$n_overlap, expected_ecm=hec$expected, fold_ecm=hec$fold, p_ecm=hec$p_hyper,
      COL1A1=("COL1A1" %in% mb), COL3A1=("COL3A1" %in% mb), FN1=("FN1" %in% mb),
      SPARC=("SPARC" %in% mb), POSTN=("POSTN" %in% mb), THBS2=("THBS2" %in% mb),
      EZH2=("EZH2" %in% mb), miR29a=("hsa-miR-29a" %in% mb), miR101=("hsa-miR-101" %in% mb),
      ETS1=("ETS1" %in% mb), NFKB1=("NFKB1" %in% mb), RELA=("RELA" %in% mb), SP1=("SP1" %in% mb),
      pri30_members=paste(sort(intersect(mb,pri30)),collapse=";"),
      members=paste(sort(mb), collapse=";"))
  if(!is.null(r$module)) saveRDS(r$module, P("module_", nm, ".rds"))
}
S <- rbindlist(summ); fwrite(S, P("module_summary.csv"))
print(S[, .(run, fdr, module_n, frac_of_network, n_gene, n_TF, n_miRNA,
            overlap_pri30, expected_pri30, fold_pri30, p_pri30, n_ecm, fold_ecm, p_ecm)])

## ---- full signature enrichment for the primary module --------------------
prim <- res_members$primary$name
gene_universe <- nod[is_mir==FALSE]$name          # gene sets are protein-coding only
prim_g <- intersect(prim, gene_universe)
enr <- rbindlist(lapply(names(sets), function(s){
  cbind(data.table(signature=s), hyper(prim_g, sets[[s]], gene_universe))}))
enr[, q_BH := p.adjust(p_hyper, "BH")]; setorder(enr, p_hyper)
fwrite(enr, P("primary_module_signature_enrichment.csv"))
cat("\n--- signature enrichment, PRIMARY module (universe = 364 gene/TF network nodes) ---\n")
print(enr[1:15, .(signature, n_module, n_ref, n_overlap, expected, fold, p_hyper, q_BH)])

## same for the vignette-default module, for contrast
vig_g <- intersect(res_members$sens_vign$name, gene_universe)
enr2 <- rbindlist(lapply(names(sets), function(s){
  cbind(data.table(signature=s), hyper(vig_g, sets[[s]], gene_universe))}))
enr2[, q_BH := p.adjust(p_hyper, "BH")]; setorder(enr2, p_hyper)
fwrite(enr2, P("sens_vign_module_signature_enrichment.csv"))
cat("\n--- signature enrichment, FDR=1e-3 module ---\n")
print(enr2[1:8, .(signature, n_module, n_ref, n_overlap, expected, fold, p_hyper, q_BH)])

cat("\n--- PRIMARY module members, top 40 by node score ---\n")
print(res_members$primary[1:min(40,.N), .(score_rank,name,type,logFC,pval,score,in_prioritised30,in_core_matrisome)])
