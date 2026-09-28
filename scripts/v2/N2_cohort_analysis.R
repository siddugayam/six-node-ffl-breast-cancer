## ==========================================================================
## N2_cohort_analysis.R  --  replication of the TCGA-BRCA findings in
## independent transcriptomic cohorts.
##   (a) direction concordance of the TCGA-significant hubs
##   (b) COL1A1~COL3A1 unadjusted and CAF-adjusted
##   (c) ETS1/NFKB1/RELA/SP1 -> COL1A1 unadjusted and CAF-adjusted
##   (d) Cox models: 50 protein-coding hubs (FDR), collagen module, 4-TF module
## Outputs newcohorts_summary.csv / newcohorts_per_edge.csv (+ inputs for meta)
## ==========================================================================
suppressPackageStartupMessages({
  library(data.table); library(matrixStats); library(survival); library(limma)
})
setwd("/path/to/revision")
OUT <- "results/v2"; CA <- "cache/newcohorts"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
set.seed(20260908)

cohorts <- readRDS(file.path(CA,"cohorts_base.rds"))
cat("cohorts loaded:", paste(names(cohorts), collapse=", "), "\n\n")

## ------------------------------------------------------------ gene sets ---
hubs_tab <- read.csv("results/network_topology_hubs.csv", stringsAsFactors=FALSE)
HUBS <- hubs_tab$name[hubs_tab$is_hub & hubs_tab$type != "miRNA"]
stopifnot(length(HUBS)==50, !any(grepl("^hsa-", HUBS)))
cat("protein-coding hubs:", length(HUBS), "\n")

cafA_genes <- readLines(file.path(CA,"cafA_signature_genes.txt"))
cafB_genes <- c("DCN","LUM","FAP","THY1")
nodes <- read.delim("data/canonical_nodes.tsv", stringsAsFactors=FALSE)$name
stopifnot(!any(grepl("^COL", c(cafA_genes,cafB_genes))),
          !any(c(cafA_genes,cafB_genes) %in% nodes))
cat("CAF-A signature genes (ESTIMATE stromal, collagens+network nodes removed):",
    length(cafA_genes), "\n")
TFS <- c("ETS1","NFKB1","RELA","SP1")
COLMOD <- c("COL1A1","COL3A1")

## ------------------------------------------------------------- helpers ----
zrow <- function(M){ s <- rowSds(M, na.rm=TRUE); s[s==0] <- NA; (M - rowMeans(M, na.rm=TRUE))/s }
score_of <- function(X, genes){
  g <- intersect(genes, rownames(X)); if(length(g) < 2) return(list(v=NULL, n=length(g)))
  Z <- zrow(X[g,,drop=FALSE]); Z <- Z[rowSums(is.na(Z))==0,,drop=FALSE]
  list(v=colMeans(Z), n=nrow(Z))
}
rk <- function(v) { r <- rank(v, na.last="keep"); (r-mean(r,na.rm=TRUE))/sd(r,na.rm=TRUE) }
## Spearman zero-order and partial; returns rho, n, p, and Fisher-z SE
pct_ret <- function(r0, r1){  ## meaningless when the raw rho is ~0
  if(!is.finite(r0) || !is.finite(r1) || abs(r0) < 0.05) return(NA_real_)
  100*r1/r0
}
sp <- function(x, y, z=NULL){
  ok <- is.finite(x) & is.finite(y) & (if(is.null(z)) TRUE else is.finite(z))
  x <- x[ok]; y <- y[ok]; n <- sum(ok)
  if(n < 10) return(c(rho=NA, n=n, p=NA, se=NA))
  a <- rk(x); b <- rk(y)
  if(is.null(z)){ r <- cor(a,b); k <- 0 } else {
    c0 <- rk(z[ok])
    rxy <- cor(a,b); rxz <- cor(a,c0); ryz <- cor(b,c0)
    r <- (rxy - rxz*ryz)/sqrt((1-rxz^2)*(1-ryz^2)); k <- 1
  }
  df <- n - 2 - k
  tt <- r*sqrt(df/(1-r^2))
  c(rho=r, n=n, p=2*pt(-abs(tt), df), se=1/sqrt(n-3-k))
}
binom_pct <- function(nsucc, ntot){
  if(ntot==0) return(c(pct=NA, p=NA))
  c(pct=100*nsucc/ntot, p=binom.test(nsucc, ntot, 0.5)$p.value)
}

SUM  <- list()   ## cohort-level summary rows
EDGE <- list()   ## per-edge rows (b + c + hub cox)
addS <- function(...) SUM[[length(SUM)+1]] <<- data.frame(..., stringsAsFactors=FALSE)
addE <- function(...) EDGE[[length(EDGE)+1]] <<- data.frame(..., stringsAsFactors=FALSE)

## ================= TCGA reference quantities (computed in THIS run) ========
Xt <- cohorts$TCGA_BRCA$X
hub_t <- intersect(HUBS, rownames(Xt))
cat("hubs measurable in TCGA:", length(hub_t), "\n")
Rt <- t(apply(Xt[hub_t,,drop=FALSE], 1, rank))
Zt <- t(scale(t(Rt)))
Ct <- (Zt %*% t(Zt)) / (ncol(Zt)-1)
pairs_idx <- which(upper.tri(Ct), arr.ind=TRUE)
tcga_pairs <- data.frame(g1=hub_t[pairs_idx[,1]], g2=hub_t[pairs_idx[,2]],
                         rho=Ct[upper.tri(Ct)], stringsAsFactors=FALSE)
nT <- ncol(Xt)
tt <- tcga_pairs$rho*sqrt((nT-2)/(1-tcga_pairs$rho^2))
tcga_pairs$p <- 2*pt(-abs(tt), nT-2)
tcga_pairs$q <- p.adjust(tcga_pairs$p, "BH")
tcga_pairs$sig <- tcga_pairs$q < 0.05
cat("TCGA hub-hub pairs:", nrow(tcga_pairs), " significant (BH q<0.05):",
    sum(tcga_pairs$sig), "\n")

## TCGA tumour-vs-normal logFC for the hubs, recomputed here with limma
Nt <- cohorts$TCGA_BRCA$normals
Xall <- cbind(Xt, Nt)
grp <- factor(c(rep("Tumor", ncol(Xt)), rep("Normal", ncol(Nt))), levels=c("Normal","Tumor"))
keepv <- rowSds(Xall) > 0
Xall <- Xall[keepv,]
stopifnot(sum(duplicated(rownames(Xall)))==0)
fit <- eBayes(lmFit(Xall, model.matrix(~grp)))
tt2 <- topTable(fit, coef="grpTumor", number=Inf, sort.by="none")
## SILENT-FAILURE GUARD: topTable must return gene symbols, not integer indices
stopifnot(!all(grepl("^[0-9]+$", rownames(tt2))))
stopifnot(mean(rownames(tt2) %in% rownames(Xall)) == 1)
tcga_de <- data.frame(gene=rownames(tt2), logFC=tt2$logFC, q=tt2$adj.P.Val,
                      stringsAsFactors=FALSE)
cat("TCGA DE recomputed: n_tumor =", ncol(Xt), " n_normal =", ncol(Nt),
    " genes =", nrow(tcga_de), " sig(q<0.05,|lfc|>1) =",
    sum(tcga_de$q<0.05 & abs(tcga_de$logFC)>1), "\n")
hub_de <- tcga_de[match(hub_t, tcga_de$gene), ]
hub_de_sig <- hub_de[!is.na(hub_de$q) & hub_de$q < 0.05, ]
cat("hubs DE-significant in TCGA (q<0.05):", nrow(hub_de_sig), "\n")

## TCGA hub Cox (reference for Cox-direction concordance)
tcga_cox <- NULL
{
  ph <- cohorts$TCGA_BRCA$pheno
  ok <- is.finite(ph$OS_time) & ph$OS_time > 0 & !is.na(ph$OS_event)
  res <- lapply(hub_t, function(g){
    v <- scale(rank(Xt[g, ok]))[,1]
    f <- try(coxph(Surv(ph$OS_time[ok], ph$OS_event[ok]) ~ v), silent=TRUE)
    if(inherits(f,"try-error")) return(c(NA,NA))
    s <- summary(f)$coefficients; c(s[1,1], s[1,5])
  })
  tcga_cox <- data.frame(gene=hub_t, beta=sapply(res,`[`,1), p=sapply(res,`[`,2))
  tcga_cox$q <- p.adjust(tcga_cox$p, "BH")
  cat("TCGA hub Cox OS: n =", sum(ok), " events =", sum(ph$OS_event[ok]),
      " hubs q<0.05:", sum(tcga_cox$q<0.05, na.rm=TRUE), "\n\n")
}
saveRDS(list(pairs=tcga_pairs, de=tcga_de, cox=tcga_cox, hubs=HUBS),
        file.path(CA,"tcga_reference.rds"))

## ================================= per-cohort analysis =====================
for(nm in names(cohorts)){
  co <- cohorts[[nm]]; X <- co$X; ph <- co$pheno
  cat("=========================================================\n")
  cat("### ", nm, " | ", co$accession, " | ", co$platform, " | n =", ncol(X), "\n")
  gpres <- rownames(X)

  ## ---- CAF scores -------------------------------------------------------
  sa <- score_of(X, cafA_genes); sb <- score_of(X, cafB_genes)
  cafA <- sa$v; cafB <- sb$v
  cat("CAF-A genes used:", sa$n, " CAF-B genes used:", sb$n, "\n")
  cA_col1 <- if(!is.null(cafA) && "COL1A1" %in% gpres) sp(cafA, X["COL1A1",])["rho"] else NA
  cA_col3 <- if(!is.null(cafA) && "COL3A1" %in% gpres) sp(cafA, X["COL3A1",])["rho"] else NA
  cAB     <- if(!is.null(cafA) && !is.null(cafB)) sp(cafA, cafB)["rho"] else NA
  cat(sprintf("CAF-A vs COL1A1 rho=%.3f  vs COL3A1 rho=%.3f ; cor(CAF-A,CAF-B)=%.3f\n",
              cA_col1, cA_col3, cAB))

  ## ---- (a1) hub co-expression sign concordance --------------------------
  hp <- intersect(hub_t, gpres)
  pr <- tcga_pairs[tcga_pairs$g1 %in% hp & tcga_pairs$g2 %in% hp, ]
  prs <- pr[pr$sig, ]
  rho_ext <- if(nrow(prs)) mapply(function(a,b) sp(X[a,], X[b,])["rho"], prs$g1, prs$g2) else numeric(0)
  conc <- sum(sign(rho_ext) == sign(prs$rho), na.rm=TRUE)
  bp <- binom_pct(conc, sum(is.finite(rho_ext)))
  rho_all <- if(nrow(pr)) mapply(function(a,b) sp(X[a,], X[b,])["rho"], pr$g1, pr$g2) else numeric(0)
  sprho <- if(length(rho_all)>3) cor(rho_all, pr$rho, method="spearman", use="complete.obs") else NA
  cat(sprintf("(a1) hub co-expression: %d hubs measurable, %d TCGA-sig pairs tested, "
              , length(hp), sum(is.finite(rho_ext))))
  cat(sprintf("%.1f%% same sign, binom p=%.3g, rho(effect sizes)=%.3f\n",
              bp["pct"], bp["p"], sprho))
  addS(cohort=nm, accession=co$accession, platform=co$platform, n=ncol(X),
       analysis="a1_hub_coexpression_sign_concordance",
       n_tested=sum(is.finite(rho_ext)), n_concordant=conc,
       pct=unname(bp["pct"]), binom_p=unname(bp["p"]),
       extra=sprintf("spearman_effectsize_rho=%.3f; hubs_measurable=%d", sprho, length(hp)))

  ## ---- (b) COL1A1 ~ COL3A1 ---------------------------------------------
  if(all(COLMOD %in% gpres)){
    r0 <- sp(X["COL1A1",], X["COL3A1",])
    rA <- if(!is.null(cafA)) sp(X["COL1A1",], X["COL3A1",], cafA) else c(rho=NA,n=NA,p=NA,se=NA)
    rB <- if(!is.null(cafB)) sp(X["COL1A1",], X["COL3A1",], cafB) else c(rho=NA,n=NA,p=NA,se=NA)
    cat(sprintf("(b) COL1A1~COL3A1  raw rho=%.3f (p=%.3g, n=%d) | CAF-A adj=%.3f (p=%.3g) | CAF-B adj=%.3f\n",
                r0["rho"], r0["p"], r0["n"], rA["rho"], rA["p"], rB["rho"]))
    addE(cohort=nm, accession=co$accession, n=unname(r0["n"]), edge="COL1A1~COL3A1",
         source="COL1A1", target="COL3A1", class="gene_gene",
         rho_raw=unname(r0["rho"]), p_raw=unname(r0["p"]), se_raw=unname(r0["se"]),
         rho_cafA=unname(rA["rho"]), p_cafA=unname(rA["p"]), se_cafA=unname(rA["se"]),
         rho_cafB=unname(rB["rho"]), p_cafB=unname(rB["p"]), se_cafB=unname(rB["se"]),
         pct_retained_cafA=pct_ret(unname(r0["rho"]), unname(rA["rho"])),
         q_uni=NA_real_, q_adj=NA_real_, value_type="spearman_rho")
  } else cat("(b) SKIPPED: COL1A1/COL3A1 not both measurable\n")

  ## ---- (c) TF -> COL1A1 -------------------------------------------------
  if("COL1A1" %in% gpres){
    for(tf in TFS){
      if(!(tf %in% gpres)){ cat("(c)", tf, "not measurable\n"); next }
      r0 <- sp(X[tf,], X["COL1A1",])
      rA <- if(!is.null(cafA)) sp(X[tf,], X["COL1A1",], cafA) else c(rho=NA,n=NA,p=NA,se=NA)
      rB <- if(!is.null(cafB)) sp(X[tf,], X["COL1A1",], cafB) else c(rho=NA,n=NA,p=NA,se=NA)
      inv <- is.finite(r0["rho"]) && is.finite(rA["rho"]) && sign(r0["rho"]) != sign(rA["rho"])
      cat(sprintf("(c) %-6s -> COL1A1  raw=%+.3f (p=%.3g) | CAF-A=%+.3f (p=%.3g) | CAF-B=%+.3f | sign_inverts=%s\n",
                  tf, r0["rho"], r0["p"], rA["rho"], rA["p"], rB["rho"], inv))
      addE(cohort=nm, accession=co$accession, n=unname(r0["n"]),
           edge=paste0(tf,"->COL1A1"), source=tf, target="COL1A1", class="TF_target",
           rho_raw=unname(r0["rho"]), p_raw=unname(r0["p"]), se_raw=unname(r0["se"]),
           rho_cafA=unname(rA["rho"]), p_cafA=unname(rA["p"]), se_cafA=unname(rA["se"]),
           rho_cafB=unname(rB["rho"]), p_cafB=unname(rB["p"]), se_cafB=unname(rB["se"]),
           pct_retained_cafA=pct_ret(unname(r0["rho"]), unname(rA["rho"])),
           q_uni=NA_real_, q_adj=NA_real_, value_type="spearman_rho")
    }
  }
  ## also the miR-29 negative-control analogue is not available (no miRNA assay)

  ## ---- (d) survival -----------------------------------------------------
  has_surv <- !is.null(ph) && all(c("OS_time","OS_event") %in% names(ph)) &&
              sum(is.finite(ph$OS_time) & !is.na(ph$OS_event)) > 30
  endpoints <- list()
  if(!is.null(ph)){
    for(ep in c("OS","RFS","DFS","DMFS","MFS","DRFS")){
      tcol <- paste0(ep,"_time"); ecol <- paste0(ep,"_event")
      if(all(c(tcol,ecol) %in% names(ph))){
        okk <- is.finite(ph[[tcol]]) & ph[[tcol]] > 0 & !is.na(ph[[ecol]])
        if(sum(okk) > 30 && sum(ph[[ecol]][okk]) >= 10)
          endpoints[[ep]] <- list(t=ph[[tcol]], e=ph[[ecol]], ok=okk)
      }
    }
  }
  if(!length(endpoints)) cat("(d) no usable survival endpoint\n")
  for(ep in names(endpoints)){
    E <- endpoints[[ep]]; okk <- E$ok
    ## available covariates
    cvs <- c()
    for(cv in c("age","grade","stage","ER","node","size")){
      if(cv %in% names(ph)){
        v <- ph[[cv]][okk]
        if(mean(!is.na(v)) > 0.7 && length(unique(na.omit(v))) > 1) cvs <- c(cvs, cv)
      }
    }
    dd <- data.frame(time=E$t[okk], ev=E$e[okk])
    for(cv in cvs) dd[[cv]] <- ph[[cv]][okk]
    cmpl <- complete.cases(dd)
    cat(sprintf("(d) %s: n=%d events=%d ; covariates used: %s (adjusted n=%d)\n",
                ep, nrow(dd), sum(dd$ev), if(length(cvs)) paste(cvs,collapse="+") else "NONE",
                sum(cmpl)))
    fml_adj <- if(length(cvs)) paste("~ v +", paste(cvs, collapse=" + ")) else "~ v"
    ## 50 hubs
    hp2 <- intersect(hub_t, gpres)
    resl <- lapply(hp2, function(g){
      v <- scale(rank(X[g, okk]))[,1]
      d1 <- cbind(dd, v=v)
      f <- try(coxph(as.formula(paste("Surv(time,ev)", fml_adj)), data=d1[cmpl,]), silent=TRUE)
      fu <- try(coxph(Surv(time,ev) ~ v, data=d1), silent=TRUE)
      out <- c(NA,NA,NA,NA,NA,NA)
      if(!inherits(fu,"try-error")){ s <- summary(fu)$coefficients; out[1:3] <- c(s["v",1], s["v",5], s["v",3]) }
      if(!inherits(f,"try-error")){ s <- summary(f)$coefficients; out[4:6] <- c(s["v",1], s["v",5], s["v",3]) }
      out
    })
    hb <- data.frame(gene=hp2, beta_uni=sapply(resl,`[`,1), p_uni=sapply(resl,`[`,2),
                     se_uni=sapply(resl,`[`,3),
                     beta_adj=sapply(resl,`[`,4), p_adj=sapply(resl,`[`,5),
                     se_adj=sapply(resl,`[`,6),
                     stringsAsFactors=FALSE)
    hb$q_uni <- p.adjust(hb$p_uni,"BH"); hb$q_adj <- p.adjust(hb$p_adj,"BH")
    cat(sprintf("     hubs tested=%d  FDR<0.05 univariate=%d  FDR<0.05 adjusted=%d\n",
                nrow(hb), sum(hb$q_uni<0.05,na.rm=TRUE), sum(hb$q_adj<0.05,na.rm=TRUE)))
    for(i in seq_len(nrow(hb)))
      addE(cohort=nm, accession=co$accession, n=sum(cmpl),
           edge=paste0("COX_",ep,"_",hb$gene[i]), source=hb$gene[i], target=ep,
           class="cox_hub",
           rho_raw=exp(hb$beta_uni[i]), p_raw=hb$p_uni[i], se_raw=hb$se_uni[i],
           rho_cafA=exp(hb$beta_adj[i]), p_cafA=hb$p_adj[i], se_cafA=hb$se_adj[i],
           rho_cafB=NA, p_cafB=NA, se_cafB=NA,
           pct_retained_cafA=NA, q_uni=hb$q_uni[i], q_adj=hb$q_adj[i],
           value_type="hazard_ratio_per_SD (raw=univariate, cafA cols=covariate-adjusted)")
    ## (a2) Cox direction concordance vs TCGA (only for non-TCGA cohorts, OS)
    if(nm != "TCGA_BRCA" && !is.null(tcga_cox)){
      ## TCGA OS gives 0 hubs at FDR<0.05, so the reference set is the
      ## nominally significant one (stated explicitly in `extra`).
      ref <- tcga_cox[!is.na(tcga_cox$p) & tcga_cox$p < 0.05, ]
      mm <- match(ref$gene, hb$gene)
      keep2 <- !is.na(mm) & is.finite(hb$beta_uni[mm])
      nc <- sum(sign(hb$beta_uni[mm][keep2]) == sign(ref$beta[keep2]))
      bp2 <- binom_pct(nc, sum(keep2))
      mm2 <- match(tcga_cox$gene, hb$gene)
      ok3 <- !is.na(mm2) & is.finite(hb$beta_uni[mm2]) & is.finite(tcga_cox$beta)
      sr <- if(sum(ok3) > 5) cor(hb$beta_uni[mm2][ok3], tcga_cox$beta[ok3], method="spearman") else NA
      cat(sprintf("     (a2) Cox-direction concordance vs TCGA-OS nominally-sig hubs (n=%d): %.1f%% p=%.3g ; rho(beta,all %d hubs)=%.3f\n",
                  sum(keep2), bp2["pct"], bp2["p"], sum(ok3), sr))
      addS(cohort=nm, accession=co$accession, platform=co$platform, n=sum(okk),
           analysis=paste0("a2_hub_cox_direction_vs_TCGA_", ep),
           n_tested=sum(keep2), n_concordant=nc, pct=unname(bp2["pct"]),
           binom_p=unname(bp2["p"]),
           extra=sprintf("reference = %d TCGA hubs with OS nominal p<0.05 (0 survive FDR); spearman_beta_all_hubs=%.3f (n=%d)",
                         nrow(ref), sr, sum(ok3)))
      ## better powered: sign agreement over ALL measurable hubs
      nc_all <- sum(sign(hb$beta_uni[mm2][ok3]) == sign(tcga_cox$beta[ok3]))
      bp4 <- binom_pct(nc_all, sum(ok3))
      cat(sprintf("     (a2b) Cox-beta sign agreement over ALL %d measurable hubs: %.1f%% p=%.3g\n",
                  sum(ok3), bp4["pct"], bp4["p"]))
      addS(cohort=nm, accession=co$accession, platform=co$platform, n=sum(okk),
           analysis=paste0("a2b_hub_cox_sign_all_hubs_", ep),
           n_tested=sum(ok3), n_concordant=nc_all, pct=unname(bp4["pct"]),
           binom_p=unname(bp4["p"]),
           extra=sprintf("spearman of Cox log-HR vs TCGA OS = %.3f", sr))
    }
    ## module scores
    for(mdn in c("collagen_module","fourTF_module")){
      gs <- if(mdn=="collagen_module") COLMOD else TFS
      s <- score_of(X, gs); if(is.null(s$v)) next
      v <- scale(rank(s$v[okk]))[,1]
      d1 <- cbind(dd, v=v)
      fu <- try(coxph(Surv(time,ev) ~ v, data=d1), silent=TRUE)
      fa <- try(coxph(as.formula(paste("Surv(time,ev)", fml_adj)), data=d1[cmpl,]), silent=TRUE)
      gu <- if(inherits(fu,"try-error")) c(NA,NA,NA) else summary(fu)$coefficients["v",c(1,5,3)]
      ga <- if(inherits(fa,"try-error")) c(NA,NA,NA) else summary(fa)$coefficients["v",c(1,5,3)]
      cat(sprintf("     %s (%d genes) %s: HR/SD uni=%.3f p=%.3g | adj=%.3f p=%.3g\n",
                  mdn, s$n, ep, exp(gu[1]), gu[2], exp(ga[1]), ga[2]))
      addS(cohort=nm, accession=co$accession, platform=co$platform, n=sum(cmpl),
           analysis=paste0("d_", mdn, "_cox_", ep),
           n_tested=s$n, n_concordant=NA, pct=exp(gu[1]), binom_p=gu[2],
           extra=sprintf("HR_adj=%.3f p_adj=%.3g covars=%s", exp(ga[1]), ga[2],
                         if(length(cvs)) paste(cvs,collapse="+") else "none"))
      addE(cohort=nm, accession=co$accession, n=sum(cmpl),
           edge=paste0("COX_",ep,"_",mdn), source=mdn, target=ep, class="cox_module",
           rho_raw=exp(gu[1]), p_raw=gu[2], se_raw=gu[3],
           rho_cafA=exp(ga[1]), p_cafA=ga[2], se_cafA=ga[3],
           rho_cafB=NA, p_cafB=NA, se_cafB=NA, pct_retained_cafA=NA,
           q_uni=NA_real_, q_adj=NA_real_,
           value_type="hazard_ratio_per_SD (raw=univariate, cafA cols=covariate-adjusted)")
    }
  }
  ## ---- (a3) cross-study hub DE direction, non-tumour reference only -----
  if(nm == "GTEx_breast"){
    ## rank-percentile within each sample, then tumour(TCGA) vs normal(GTEx)
    pct_of <- function(M) apply(M, 2, function(v) rank(v)/length(v))
    gcom <- intersect(rownames(Xt), rownames(X))
    Pt <- pct_of(Xt[gcom,,drop=FALSE]); Pg <- pct_of(X[gcom,,drop=FALSE])
    d <- rowMedians(Pt) - rowMedians(Pg); names(d) <- gcom
    ref <- hub_de_sig[hub_de_sig$gene %in% gcom, ]
    nc3 <- sum(sign(d[ref$gene]) == sign(ref$logFC))
    bp3 <- binom_pct(nc3, nrow(ref))
    cat(sprintf("(a3) cross-study TCGA-tumour vs GTEx-normal hub direction: %d/%d = %.1f%% concordant, binom p=%.3g\n",
                nc3, nrow(ref), bp3["pct"], bp3["p"]))
    addS(cohort=nm, accession=co$accession, platform=co$platform, n=ncol(X),
         analysis="a3_hub_DE_direction_TCGAtumour_vs_GTExnormal_CROSS_STUDY",
         n_tested=nrow(ref), n_concordant=nc3, pct=unname(bp3["pct"]),
         binom_p=unname(bp3["p"]),
         extra="cross-study/cross-platform: within-sample rank percentiles; interpret with caution")
  }

  ## ---- sign-inversion replication flag for the two headline TF axes ------
  for(tf in c("ETS1","NFKB1")){
    ee <- Filter(function(z) z$cohort==nm && z$edge==paste0(tf,"->COL1A1"), EDGE)
    if(!length(ee)) next
    z <- ee[[1]]
    inv <- is.finite(z$rho_raw) && is.finite(z$rho_cafA) && z$rho_raw > 0 && z$rho_cafA < 0
    addS(cohort=nm, accession=co$accession, platform=co$platform, n=z$n,
         analysis=paste0("c_sign_inversion_", tf, "_to_COL1A1"),
         n_tested=1, n_concordant=as.integer(inv), pct=NA, binom_p=z$p_cafA,
         extra=sprintf("rho_raw=%+.4f -> rho_CAFA=%+.4f (inverts=%s); rho_CAFB=%+.4f",
                       z$rho_raw, z$rho_cafA, inv, z$rho_cafB))
  }

  ## record CAF diagnostics
  addS(cohort=nm, accession=co$accession, platform=co$platform, n=ncol(X),
       analysis="caf_score_diagnostics", n_tested=sa$n, n_concordant=sb$n,
       pct=unname(cA_col1), binom_p=NA,
       extra=sprintf("CAFA_vs_COL3A1=%.3f cor_CAFA_CAFB=%.3f", cA_col3, cAB))
  cat("\n")
}

SUMd <- do.call(rbind, SUM); EDGEd <- do.call(rbind, EDGE)
write.csv(SUMd, file.path(OUT,"newcohorts_summary.csv"), row.names=FALSE)
write.csv(EDGEd, file.path(OUT,"newcohorts_per_edge.csv"), row.names=FALSE)
cat("\nWROTE", file.path(OUT,"newcohorts_summary.csv"), nrow(SUMd), "rows\n")
cat("WROTE", file.path(OUT,"newcohorts_per_edge.csv"), nrow(EDGEd), "rows\n")
