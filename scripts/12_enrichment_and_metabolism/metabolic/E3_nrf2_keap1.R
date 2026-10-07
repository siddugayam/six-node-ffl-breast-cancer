#!/usr/bin/env Rscript
# =============================================================================
# E3_nrf2_keap1.R
# PART C: KEAP1-NFE2L2 enrichment in the 5- and 6-node networks and its link
# to therapy resistance. Tested four ways:
#   1. Are NFE2L2 / KEAP1 (or any canonical NRF2 target) actually IN the network?
#      What genes carry the reported Reactome KEAP1-NFE2L2 enrichment?
#   2. NFE2L2 target-signature GSVA score: tumour vs normal, and vs module scores,
#      benchmarked against 1,000 random same-size signatures (the same control as
#      for the NF-kB activity score).
#   3. NFE2L2 / KEAP1 expression and mutation in TCGA-BRCA (MC3).
#   4. BRCA vs LUAD / LUSC mutation frequency, since the KEAP1-NFE2L2
#      resistance literature is lung-derived.
# Output: results/multiomics/nrf2_keap1_brca.csv
# =============================================================================
suppressPackageStartupMessages({library(data.table)})
ROOT <- "/path/to/revision"
OUT  <- file.path(ROOT,"results/multiomics")
MC3  <- "/path/to/home/Desktop/DD/R_GPR/ML/external_cohorts/MC3/mc3_nonsilentGene.xena.gz"
PHE  <- "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/PANCAN/TCGA_phenotype_denseDataOnlyDownload.tsv.gz"
GMT  <- "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/PANCAN/c2.all.v2026.1.Hs.symbols.gmt"
con <- file(file.path(ROOT,"logs/E3_nrf2_keap1.log"), open="wt")
sink(con, split=TRUE); sink(con, type="message")
cat("=== E3_nrf2_keap1.R ", format(Sys.time()), " ===\n")
set.seed(42)
rows <- list()

## ==================== 1. NETWORK MEMBERSHIP =================================
cat("\n########## 1. IS THE NRF2 AXIS ACTUALLY IN THE NETWORK? ##########\n")
nodes <- fread(file.path(ROOT,"data/canonical_nodes.tsv"))
edges <- fread(file.path(ROOT,"data/canonical_edges.tsv"))
read_gmt <- function(f){ ll <- strsplit(readLines(f),"\t")
  setNames(lapply(ll, function(x) unique(x[-c(1,2)])), sapply(ll,`[`,1)) }
c2 <- read_gmt(GMT)
KEAP1_REACTOME <- c2[["REACTOME_KEAP1_NFE2L2_PATHWAY"]]
NRF2_CORE <- c("NQO1","GCLC","GCLM","TXNRD1","TXN","SLC7A11","GSR","G6PD","PGD",
               "ME1","AKR1B10","AKR1C1","AKR1C2","AKR1C3","SRXN1","PRDX1","FTL",
               "FTH1","HMOX1","OSGIN1","CBR1","GPX2","ABCC2","SQSTM1","GSTP1",
               "GSTM3","UGT1A6","EPHX1","CES1","TALDO1","IDH1")
netg <- nodes[type %in% c("TF","Gene")]$name
for (g in c("NFE2L2","KEAP1","CUL3","MAF","MAFG","MAFK")){
  cat(sprintf("  %-8s nodes=%d  edges=%d\n", g, sum(nodes$name==g),
              sum(edges$source==g | edges$target==g)))
}
ov <- intersect(netg, KEAP1_REACTOME)
cat("\n  Reactome KEAP1-NFE2L2 pathway size:", length(KEAP1_REACTOME), "genes\n")
cat("  Network genes overlapping it (n=", length(ov), "): ",
    paste(sort(ov), collapse=", "), "\n", sep="")
cat("  Of those, canonical NRF2 TARGET genes: ",
    paste(intersect(ov, NRF2_CORE), collapse=", "), "\n")
cat("  Of those, NFE2L2 or KEAP1 themselves: ",
    paste(intersect(ov, c("NFE2L2","KEAP1")), collapse=", "), " <- (empty means none)\n")
mod <- readRDS(file.path(ROOT,"data/ffl_module_sets.rds"))
cat("\n  Manuscript exemplar module (MS_MODULE, n=", length(mod$MS_MODULE), "): ",
    paste(mod$MS_MODULE, collapse=", "), "\n", sep="")
cat("  Its overlap with Reactome KEAP1-NFE2L2: ",
    paste(intersect(mod$MS_MODULE, KEAP1_REACTOME), collapse=", "), "\n")
cat("  Its overlap with canonical NRF2 targets: ",
    paste(intersect(mod$MS_MODULE, NRF2_CORE), collapse=", "),
    " <- (empty means none)\n")
rows$membership <- data.table(
  analysis="network_membership",
  item=c("NFE2L2_in_network","KEAP1_in_network","reactome_KEAP1_NFE2L2_size",
         "network_genes_in_reactome_set","canonical_NRF2_targets_in_network",
         "MS_module_genes_in_reactome_set"),
  value=c(sum(nodes$name=="NFE2L2"), sum(nodes$name=="KEAP1"), length(KEAP1_REACTOME),
          length(ov), length(intersect(netg,NRF2_CORE)),
          length(intersect(mod$MS_MODULE,KEAP1_REACTOME))),
  detail=c("","","", paste(sort(ov),collapse="|"),
           paste(intersect(netg,NRF2_CORE),collapse="|"),
           paste(intersect(mod$MS_MODULE,KEAP1_REACTOME),collapse="|")))

## ==================== 2. NRF2 SIGNATURE SCORE ===============================
cat("\n\n########## 2. NFE2L2 TARGET SIGNATURE SCORE IN TCGA-BRCA ##########\n")
S  <- readRDS(file.path(OUT,"metabolic_gsva_scores.rds"))
ge <- readRDS(file.path(ROOT,"data/brca_gene_expr.rds"))
ph <- as.data.table(readRDS(file.path(ROOT,"data/brca_pheno.rds")))
ph <- ph[sample %in% colnames(ge)]; ge <- ge[, ph$sample]
is_tum <- ph$sample_type=="Primary Tumor"
tum <- ph[sample_type=="Primary Tumor"]$sample
nrf2_sets <- grep("NFE2L2|NRF2|KEAP1|OXIDATIVE_STRESS|REACTIVE_OXYGEN", rownames(S), value=TRUE)
cat("NRF2-related score rows available:", paste(nrf2_sets, collapse=", "), "\n\n")
tvn <- rbindlist(lapply(nrf2_sets, function(p){
  w <- wilcox.test(S[p,is_tum], S[p,!is_tum], exact=FALSE)
  data.table(pathway=p, n_tumour=sum(is_tum), n_normal=sum(!is_tum),
             mean_tumour=mean(S[p,is_tum]), mean_normal=mean(S[p,!is_tum]),
             delta=mean(S[p,is_tum])-mean(S[p,!is_tum]), p=w$p.value)
}))
tvn[, FDR := p.adjust(p,"BH")]
cat("--- tumour vs normal ---\n")
print(tvn[order(p)][, .(pathway, delta=round(delta,4), p=signif(p,3), FDR=signif(FDR,3),
       direction=fifelse(FDR<0.05 & delta>0,"UP_in_tumour",
                 fifelse(FDR<0.05 & delta<0,"DOWN_in_tumour","ns")))])
rows$tvn <- tvn[, .(analysis="NRF2_score_tumour_vs_normal", item=pathway,
                    value=delta, detail=paste0("p=",signif(p,3)," FDR=",signif(FDR,3)))]

cat("\n--- NFE2L2 and KEAP1 mRNA, tumour vs normal (limma table) ---\n")
de <- fread(file.path(ROOT,"results/BRCA_DEX_genes.csv"))
print(de[feature %in% c("NFE2L2","KEAP1","CUL3","NQO1","GCLC","GCLM","TXNRD1","SLC7A11","SQSTM1")][
      order(adj.P.Val)][, .(feature, logFC=round(logFC,3), P.Value=signif(P.Value,3),
                            adj.P.Val=signif(adj.P.Val,3))])

## --- relation to module scores, WITH a random-signature null ---------------
cat("\n\n--- NRF2 signature score vs FFL MODULE SCORES (tumours) ---\n")
ms <- readRDS(file.path(ROOT,"data/tcga_module_scores.rds"))
msg <- ms$gsva; common <- intersect(colnames(msg), tum)
cat("shared tumours:", length(common), "\n")
mc <- rbindlist(lapply(rownames(msg), function(m){
  rbindlist(lapply(nrf2_sets, function(p){
    ct <- suppressWarnings(cor.test(as.numeric(msg[m,common]), as.numeric(S[p,common]),
                                    method="spearman", exact=FALSE))
    data.table(module=m, pathway=p, n=length(common), rho=unname(ct$estimate), p=ct$p.value)
  }))
}))
mc[, FDR := p.adjust(p,"BH")]
print(mc[module %in% c("MS_exemplar_module","FFL_3node","FFL_higher_order_only","FFL_6node_all")][
      order(pathway,module)][, .(module, pathway, n, rho=round(rho,3), FDR=signif(FDR,3))], nrows=60)
rows$modcor <- mc[, .(analysis="NRF2_score_vs_module_score", item=paste0(module,"::",pathway),
                      value=rho, detail=paste0("n=",n," FDR=",signif(FDR,3)))]

cat("\n--- RANDOM-SIGNATURE CONTROL for the NRF2-target score ---\n")
cat("Any coherent gene signature correlates with any other in bulk tumours.\n")
cat("Benchmarking the curated NRF2 target score against 1,000 random signatures\n")
cat("of the SAME SIZE drawn from the same expression matrix.\n")
sets_tab <- fread(file.path(OUT,"metabolic_pathway_genesets.csv"))
nrf2_genes <- strsplit(sets_tab[pathway=="NFE2L2_CORE_TARGETS_curated_here"]$genes, "\\|")[[1]]
k <- length(nrf2_genes); cat("NRF2 target set size:", k, "\n")
zscore <- function(M) (M - rowMeans(M)) / pmax(matrixStats::rowSds(M), 1e-9)
gz <- zscore(ge[, common, drop=FALSE])
obs_score <- colMeans(gz[intersect(nrf2_genes, rownames(gz)), , drop=FALSE])
pool <- rownames(gz)[matrixStats::rowSds(ge[,common,drop=FALSE]) > 0]
targets <- c("MS_exemplar_module","FFL_higher_order_only","FFL_3node")
nullres <- rbindlist(lapply(targets, function(m){
  obs <- suppressWarnings(cor(obs_score, as.numeric(msg[m,common]), method="spearman"))
  nullrho <- replicate(1000, {
    rs <- colMeans(gz[sample(pool, k), , drop=FALSE])
    suppressWarnings(cor(rs, as.numeric(msg[m,common]), method="spearman"))
  })
  data.table(module=m, observed_rho=obs, null_mean=mean(nullrho), null_sd=sd(nullrho),
             null_q025=quantile(nullrho,.025), null_q975=quantile(nullrho,.975),
             emp_p_two_sided=(sum(abs(nullrho)>=abs(obs))+1)/1001)
}))
print(nullres[, .(module, observed_rho=round(observed_rho,3), null_mean=round(null_mean,3),
                  null_sd=round(null_sd,3), null_95CI=paste0("[",round(null_q025,3),", ",
                  round(null_q975,3),"]"), emp_p=signif(emp_p_two_sided,3))])
rows$null <- nullres[, .(analysis="NRF2_score_vs_module_RANDOM_SIGNATURE_NULL", item=module,
                         value=observed_rho,
                         detail=paste0("null_mean=",round(null_mean,3),
                                       " null_sd=",round(null_sd,3),
                                       " emp_p=",signif(emp_p_two_sided,3)))]

## ==================== 3+4. MUTATION: BRCA vs LUNG ===========================
cat("\n\n########## 3+4. KEAP1 / NFE2L2 MUTATION: BRCA vs LUNG (MC3) ##########\n")
mc3 <- fread(cmd=paste("zcat", shQuote(MC3)), header=TRUE)
setnames(mc3, 1, "gene")
cat("MC3 matrix:", nrow(mc3), "genes x", ncol(mc3)-1, "samples\n")
phe <- fread(cmd=paste("zcat", shQuote(PHE)), header=TRUE)
setnames(phe, 1, "sample")
# NB: use EXACT names. grep("sample_type") matches "sample_type_id" first, which
# silently yields zero Primary Tumor rows.
dc <- "_primary_disease"; sc <- "sample_type"
stopifnot(all(c(dc,sc) %in% names(phe)))
cat("phenotype cols used:", dc, "|", sc, "\n")
cat("sample_type counts:\n"); print(head(sort(table(phe[[sc]]), decreasing=TRUE), 4))
phe <- phe[get(sc)=="Primary Tumor"]
lookup <- list(BRCA="breast invasive carcinoma",
               LUAD="lung adenocarcinoma",
               LUSC="lung squamous cell carcinoma",
               HNSC="head & neck squamous cell carcinoma",
               ESCA="esophageal carcinoma",
               UCEC="uterine corpus endometrioid carcinoma",
               OV  ="ovarian serous cystadenocarcinoma",
               LIHC="liver hepatocellular carcinoma")
mc3_samples <- setdiff(names(mc3), "gene")
genes_of_interest <- c("KEAP1","NFE2L2","CUL3","TP53","PIK3CA")
G <- mc3[gene %in% genes_of_interest]
mut <- rbindlist(lapply(names(lookup), function(ct){
  s <- intersect(phe[tolower(get(dc))==lookup[[ct]]]$sample, mc3_samples)
  if (!length(s)) return(NULL)
  rbindlist(lapply(genes_of_interest, function(g){
    v <- as.numeric(G[gene==g, ..s])
    data.table(cohort=ct, gene=g, n_samples=length(s), n_mutated=sum(v==1, na.rm=TRUE),
               freq=mean(v==1, na.rm=TRUE))
  }))
}))
cat("\n--- non-silent mutation frequency by cohort ---\n")
print(dcast(mut, cohort + n_samples ~ gene, value.var="freq")[
      , lapply(.SD, function(x) if(is.numeric(x)) round(x,4) else x)])
cat("\n--- KEAP1 / NFE2L2 / CUL3 detail ---\n")
print(mut[gene %in% c("KEAP1","NFE2L2","CUL3")][order(gene,-freq)][
      , .(gene, cohort, n_samples, n_mutated, freq=round(freq,4))])
# formal BRCA vs lung test
for (g in c("KEAP1","NFE2L2")){
  b <- mut[cohort=="BRCA" & gene==g]
  for (lg in c("LUAD","LUSC")){
    l <- mut[cohort==lg & gene==g]
    ft <- fisher.test(matrix(c(b$n_mutated, b$n_samples-b$n_mutated,
                               l$n_mutated, l$n_samples-l$n_mutated), nrow=2))
    cat(sprintf("  %s  BRCA %d/%d (%.2f%%) vs %s %d/%d (%.2f%%)  Fisher p=%.3g OR=%.3f\n",
        g, b$n_mutated, b$n_samples, 100*b$freq, lg, l$n_mutated, l$n_samples,
        100*l$freq, ft$p.value, unname(ft$estimate)))
    rows[[paste0("fisher_",g,"_",lg)]] <- data.table(
      analysis="KEAP1_NFE2L2_mutation_BRCA_vs_lung", item=paste0(g,"_BRCA_vs_",lg),
      value=unname(ft$estimate),
      detail=sprintf("BRCA %d/%d vs %s %d/%d Fisher p=%.3g",
                     b$n_mutated,b$n_samples,lg,l$n_mutated,l$n_samples,ft$p.value))
  }
}
# any NRF2-pathway alteration (KEAP1 or NFE2L2 or CUL3) per cohort
cat("\n--- ANY of KEAP1/NFE2L2/CUL3 mutated, per cohort ---\n")
anyalt <- rbindlist(lapply(names(lookup), function(ct){
  s <- intersect(phe[tolower(get(dc))==lookup[[ct]]]$sample, mc3_samples)
  if (!length(s)) return(NULL)
  m <- as.matrix(G[gene %in% c("KEAP1","NFE2L2","CUL3"), ..s])
  data.table(cohort=ct, n_samples=length(s), n_any=sum(colSums(m==1, na.rm=TRUE)>0),
             freq_any=mean(colSums(m==1, na.rm=TRUE)>0))
}))
print(anyalt[order(-freq_any)][, .(cohort, n_samples, n_any, freq_any=round(freq_any,4))])
rows$mut <- mut[, .(analysis="mutation_frequency_MC3", item=paste0(gene,"::",cohort),
                    value=freq, detail=paste0(n_mutated,"/",n_samples))]
rows$anyalt <- anyalt[, .(analysis="any_KEAP1_NFE2L2_CUL3_mutation", item=cohort,
                          value=freq_any, detail=paste0(n_any,"/",n_samples))]

## ==================== write =================================================
res <- rbindlist(rows, use.names=TRUE, fill=TRUE)
res[, evidence_class := "TCGA-BRCA expression (GSVA inference) + MC3 somatic mutation"]
fwrite(res, file.path(OUT,"nrf2_keap1_brca.csv"))
cat("\n\nWROTE", file.path(OUT,"nrf2_keap1_brca.csv"), " rows:", nrow(res), "\n")
fwrite(mut, file.path(OUT,"nrf2_keap1_mutation_by_cohort.csv"))
cat("WROTE nrf2_keap1_mutation_by_cohort.csv  rows:", nrow(mut), "\n")
cat("\nDONE ", format(Sys.time()), "\n")
sink(type="message"); sink(); close(con)
