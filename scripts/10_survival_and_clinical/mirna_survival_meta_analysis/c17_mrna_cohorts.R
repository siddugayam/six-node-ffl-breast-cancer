## ==========================================================================
## c17_mrna_cohorts.R -- Part B
## Seven additional breast-cancer mRNA cohorts with survival, none of which is
## already in the study: MAINZ, TRANSBIG, VDX, UPP, UNT, NKI (Bioconductor
## breastCancer* experiment packages) and GSE1456 (Stockholm) from GEO.
## Per-SD Cox for the 20 prioritised protein-coding nodes + RE meta-analysis,
## plus the CAF-adjustment sign test on the named TF->collagen axes.
## ==========================================================================
suppressPackageStartupMessages({library(Biobase); library(survival); library(metafor)
  library(data.table)})
setwd("/path/to/revision")
OUT<-"results/v6"
TF   <- c("E2F1","EZH2","GATA3","BRCA1","JUN","EGR2","ESR1","SREBF1","DNMT1","E2F3")
GENE <- c("CCND2","COL1A1","STAT5A","MYBL2","FN1","PDGFRB","MET","CXCL12","MMP14","PLAU")
NODES <- c(TF, GENE)
EXTRA <- c("COL3A1","ETS1","NFKB1","RELA","SP1","POSTN")
WANT  <- unique(c(NODES, EXTRA))
CAFG  <- readLines("results/v2/cafA_gene_list_v3.txt")

## ---------------------------------------------------------------- loaders --
COH <- list()
mk <- function(name, accession, platform, X, pd, endpoints, covars, note){
  COH[[name]] <<- list(name=name, accession=accession, platform=platform, X=X, pd=pd,
                       endpoints=endpoints, covars=covars, note=note)
}
collapse_probes <- function(E, sym){
  ok <- !is.na(sym) & sym!="" & !grepl("///|,", sym)
  E <- E[ok,,drop=FALSE]; sym <- sym[ok]
  mu <- rowMeans(E, na.rm=TRUE); o <- order(sym, -mu)
  E <- E[o,,drop=FALSE]; sym <- sym[o]
  keep <- !duplicated(sym); E <- E[keep,,drop=FALSE]; rownames(E) <- sym[keep]; E
}
PKG <- list(MAINZ=c("breastCancerMAINZ","mainz","GSE11121","Affymetrix HG-U133A (GPL96)"),
            TRANSBIG=c("breastCancerTRANSBIG","transbig","GSE7390","Affymetrix HG-U133A (GPL96)"),
            VDX=c("breastCancerVDX","vdx","GSE2034 + GSE5327","Affymetrix HG-U133A (GPL96)"),
            UPP=c("breastCancerUPP","upp","GSE3494","Affymetrix HG-U133A+B (GPL96/GPL97)"),
            UNT=c("breastCancerUNT","unt","GSE2990","Affymetrix HG-U133A+B (GPL96/GPL97)"),
            NKI=c("breastCancerNKI","nki","van de Vijver / NKI295","Agilent two-colour 25k (log ratio)"))
EP <- list(MAINZ=list(DMFS=c("t.dmfs","e.dmfs")),
           TRANSBIG=list(DMFS=c("t.dmfs","e.dmfs"),RFS=c("t.rfs","e.rfs"),OS=c("t.os","e.os")),
           VDX=list(DMFS=c("t.dmfs","e.dmfs")),
           UPP=list(RFS=c("t.rfs","e.rfs")),
           UNT=list(RFS=c("t.rfs","e.rfs"),DMFS=c("t.dmfs","e.dmfs")),
           NKI=list(DMFS=c("t.dmfs","e.dmfs"),RFS=c("t.rfs","e.rfs"),OS=c("t.os","e.os")))
for(k in names(PKG)){
  p<-PKG[[k]][1]; nm<-PKG[[k]][2]
  suppressPackageStartupMessages(library(p, character.only=TRUE))
  data(list=nm, package=p); e <- get(nm)
  fd <- fData(e)
  sym <- if("Gene.symbol" %in% colnames(fd)) fd$Gene.symbol else fd$HUGO.gene.symbol
  if(k=="NKI"){ s2 <- fd$NCBI.gene.symbol; sym[is.na(sym)|sym==""] <- s2[is.na(sym)|sym==""] }
  X <- collapse_probes(exprs(e), sym)
  mk(k, PKG[[k]][3], PKG[[k]][4], X, pData(e), EP[[k]],
     intersect(c("age","grade","node","er"), colnames(pData(e))),
     sprintf("Bioconductor %s", p))
  cat(sprintf("%-9s n=%4d genes=%5d nodes present=%2d/20\n", k, ncol(X), nrow(X),
              sum(NODES %in% rownames(X))))
}
## ---- GSE1456 (Stockholm) from the GEO series matrix, GPL96 ---------------
parse_sm <- function(path){
  L <- readLines(gzfile(path)); b<-grep("^!series_matrix_table_begin",L); en<-grep("^!series_matrix_table_end",L)
  dt <- fread(text=paste(L[(b+1):(en-1)],collapse="\n"), sep="\t", header=TRUE)
  X <- as.matrix(dt[,-1]); rownames(X) <- gsub('"','',dt[[1]])
  colnames(X) <- gsub('"','',colnames(X))
  ch <- grep("^!Sample_characteristics_ch1", L, value=TRUE)
  gsm <- gsub('"','',strsplit(grep("^!Sample_geo_accession",L,value=TRUE)[1],"\t")[[1]][-1])
  pd <- data.frame(gsm=gsm, stringsAsFactors=FALSE)
  for(line in ch){
    v <- gsub('"','',strsplit(line,"\t")[[1]][-1])
    key <- unique(sub(":.*","",v[grepl(":",v)]))
    if(length(key)!=1) next
    pd[[key]] <- ifelse(grepl(":",v), trimws(sub("^[^:]*:","",v)), NA)
  }
  list(X=X, pd=pd)
}
s <- parse_sm("cache/v6/cohorts/matrix/GSE1456-GPL96_series_matrix.txt.gz")
fdm <- fData(get("mainz"))
sym <- fdm$Gene.symbol[match(rownames(s$X), fdm$probe)]
Xs <- collapse_probes(s$X, sym)
pd <- s$pd
pd$t.rfs <- as.numeric(pd$SURV_RELAPSE)*365.25; pd$e.rfs <- as.numeric(pd$RELAPSE)
pd$t.os  <- as.numeric(pd$SURV_DEATH)*365.25;  pd$e.os  <- as.numeric(pd$DEATH)
pd$t.bcss<- as.numeric(pd$SURV_DEATH)*365.25;  pd$e.bcss<- as.numeric(pd$DEATH_BC)
pd$grade <- suppressWarnings(as.numeric(pd$ELSTON))
mk("STK","GSE1456 (Stockholm)","Affymetrix HG-U133A (GPL96)", Xs, pd,
   list(RFS=c("t.rfs","e.rfs"), OS=c("t.os","e.os"), BCSS=c("t.bcss","e.bcss")),
   "grade", "GEO series matrix, GPL96 arm")
cat(sprintf("%-9s n=%4d genes=%5d nodes present=%2d/20\n","STK",ncol(Xs),nrow(Xs),sum(NODES %in% rownames(Xs))))

PRIMARY <- c(MAINZ="DMFS", TRANSBIG="DMFS", VDX="DMFS", UPP="RFS", UNT="RFS",
             NKI="DMFS", STK="RFS")
## ---------------------------------------------------------------- Cox ------
rows<-list(); inv<-list(); INVL<-list()
for(cn in names(COH)){
  co<-COH[[cn]]; X<-co$X; pd<-co$pd
  stopifnot(ncol(X)==nrow(pd))
  cafg <- intersect(CAFG, rownames(X))
  cafs <- colMeans(t(scale(t(X[cafg,,drop=FALSE]))), na.rm=TRUE)
  colscore <- colMeans(t(scale(t(X[intersect(c("COL1A1","COL3A1"),rownames(X)),,drop=FALSE]))),na.rm=TRUE)
  for(nd in c(NODES,"COL3A1","collagen_module")){
    x <- if(nd=="collagen_module") colscore else if(nd %in% rownames(X)) as.numeric(X[nd,]) else NULL
    if(is.null(x)) next
    z <- as.numeric(scale(x))
    for(ep in names(co$endpoints)){
      tv<-co$endpoints[[ep]][1]; ev<-co$endpoints[[ep]][2]
      d <- data.frame(z=z, time=suppressWarnings(as.numeric(pd[[tv]])),
                      evn=suppressWarnings(as.numeric(pd[[ev]])), caf=cafs)
      for(cv in co$covars) d[[cv]] <- suppressWarnings(as.numeric(pd[[cv]]))
      for(model in c("univariate","adjusted","CAFadjusted")){
        vars <- switch(model, univariate="z", adjusted=c("z",co$covars), CAFadjusted=c("z","caf"))
        dd <- d[,c("z","time","evn", setdiff(vars,"z")),drop=FALSE]
        dd <- dd[complete.cases(dd)&is.finite(dd$time)&dd$time>0,,drop=FALSE]
        if(nrow(dd)<40 || sum(dd$evn)<10) next
        fit<-try(coxph(as.formula(paste0("Surv(time,evn) ~ ",paste(vars,collapse=" + "))),data=dd),silent=TRUE)
        if(inherits(fit,"try-error")) next
        sm<-summary(fit)$coefficients
        rows[[length(rows)+1]]<-data.frame(cohort=cn,accession=co$accession,platform=co$platform,
          node=nd,endpoint=ep,model=model,n=nrow(dd),nevent=sum(dd$evn),
          logHR=sm["z","coef"],se=sm["z","se(coef)"],HR=exp(sm["z","coef"]),
          lo=exp(sm["z","coef"]-1.96*sm["z","se(coef)"]),hi=exp(sm["z","coef"]+1.96*sm["z","se(coef)"]),
          p=sm["z","Pr(>|z|)"],primary=as.integer(ep==PRIMARY[[cn]]),stringsAsFactors=FALSE)
      }
    }
  }
  ## named-axis CAF-adjustment sign test (v2 protocol)
  ax <- list(c("ETS1","COL1A1"),c("NFKB1","COL1A1"),c("RELA","COL1A1"),c("SP1","COL1A1"),
             c("COL1A1","COL3A1"),c("EZH2","COL1A1"))
  for(a in ax){
    if(!all(a %in% rownames(X))) next
    u<-as.numeric(X[a[1],]); v<-as.numeric(X[a[2],])
    r0<-suppressWarnings(cor(u,v,method="spearman",use="pairwise"))
    ru<-residuals(lm(rank(u)~rank(cafs))); rv<-residuals(lm(rank(v)~rank(cafs)))
    r1<-cor(ru,rv)
    inv[[length(inv)+1]]<-data.frame(cohort=cn,accession=co$accession,n=ncol(X),
      axis=paste(a,collapse="~"),rho_raw=r0,rho_CAFadj=r1,
      inverts=as.integer(sign(r0)!=sign(r1) & r0>0),stringsAsFactors=FALSE)
  }
  inv2 <- data.frame(cohort=cn,accession=co$accession,platform=co$platform,n=ncol(X),
    genes=nrow(X), nodes_measurable=sum(NODES %in% rownames(X)),
    caf_genes=length(cafg), endpoints=paste(names(co$endpoints),collapse=";"),
    primary=PRIMARY[[cn]], note=co$note, stringsAsFactors=FALSE)
  INVL[[cn]]<-inv2
}
COX<-rbindlist(rows); fwrite(COX, file.path(OUT,"cohorts_mrna_percohort_cox.csv"))
INV<-rbindlist(INVL); fwrite(INV, file.path(OUT,"cohorts_mrna_inventory.csv"))
AX <- rbindlist(Filter(function(z) "axis" %in% names(z), inv))
fwrite(AX, file.path(OUT,"cohorts_mrna_named_axis_cafadj.csv"))
cat("\nCox fits:",nrow(COX),"cohorts:",length(unique(COX$cohort)),"\n")
print(as.data.frame(INV[,.(cohort,accession,n,genes,nodes_measurable,endpoints,primary)]))

## ------------------------------------------------------- meta-analysis -----
pool<-function(df,label){
  o<-list()
  for(nd in unique(df$node)){ s<-df[node==nd]; if(nrow(s)<3) next
    m<-try(rma(yi=s$logHR,sei=s$se,method="DL"),silent=TRUE); if(inherits(m,"try-error")) next
    w<-weights(m)
    for(i in seq_len(nrow(s))) o[[length(o)+1]]<-data.frame(analysis=label,node=nd,row_type="study",
      cohort=s$cohort[i],accession=s$accession[i],endpoint=s$endpoint[i],n=s$n[i],nevent=s$nevent[i],
      HR=s$HR[i],lo=s$lo[i],hi=s$hi[i],p=s$p[i],weight_pct=w[i],k=nrow(s),I2=NA,tau2=NA,Q=NA,p_Q=NA)
    o[[length(o)+1]]<-data.frame(analysis=label,node=nd,row_type="RE_pooled",
      cohort="POOLED (DerSimonian-Laird)",accession="",endpoint="primary per cohort",
      n=sum(s$n),nevent=sum(s$nevent),HR=exp(m$b[1]),lo=exp(m$ci.lb),hi=exp(m$ci.ub),p=m$pval,
      weight_pct=100,k=m$k,I2=m$I2,tau2=m$tau2,Q=m$QE,p_Q=m$QEp)
  }
  rbindlist(o)
}
M <- rbindlist(list(pool(COX[model=="univariate"&primary==1],"newcohorts_primary_univariate"),
                    pool(COX[model=="CAFadjusted"&primary==1],"newcohorts_primary_CAFadjusted")))
M[row_type=="RE_pooled", q := p.adjust(p,"BH"), by=analysis]
fwrite(M, file.path(OUT,"cohorts_mrna_meta.csv"))
cat("\n===== 7 NEW mRNA cohorts, pooled per-SD HR (primary endpoint, univariate) =====\n")
print(as.data.frame(M[analysis=="newcohorts_primary_univariate"&row_type=="RE_pooled"][order(p),
  .(node,k,n,nevent,HR=round(HR,3),lo=round(lo,3),hi=round(hi,3),p=signif(p,3),q=signif(q,3),I2=round(I2,1))]))
cat("\n===== same, after adjustment for CAF score =====\n")
print(as.data.frame(M[analysis=="newcohorts_primary_CAFadjusted"&row_type=="RE_pooled"][order(p),
  .(node,k,n,nevent,HR=round(HR,3),p=signif(p,3),q=signif(q,3),I2=round(I2,1))]))
cat("\n===== named-axis CAF adjustment =====\n")
print(as.data.frame(AX[,.(cohort,axis,rho_raw=round(rho_raw,3),rho_CAFadj=round(rho_CAFadj,3),inverts)]))
