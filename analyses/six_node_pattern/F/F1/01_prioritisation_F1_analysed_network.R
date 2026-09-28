# ==============================================================================
# Pre-specified composite prioritisation of network nodes.
#
#
# This script answers that. The criteria below are fixed BEFORE inspecting the ranking, every
# component is rank-normalised to [0,1] so no single scale dominates, missing components are
# handled by renormalising over the components a node actually has, and the FULL ranked table is
# written out so the selection can be audited rather than taken on trust.
#
# Seven evidence domains, weighted equally:
#   T  topology    - FFL core participation, degree, betweenness
#   E  evidence    - fraction of the node's edges validated by low-throughput assay
#   D  expression  - |log2FC| and significance, tumour vs normal (TCGA-BRCA)
#   R  replication - direction concordance across three independent GEO cohorts
#   M  multi-omic  - protein/mRNA concordance, promoter methylation, copy number, mutation
#   C  clinical    - strength of survival association (OS/PFI), FDR-aware
#   F  functional  - CRISPR essentiality in 53 breast cancer cell lines
# ==============================================================================
suppressMessages({library(data.table)})
REV <- "/path/to/revision"; R5 <- "/path/to/revision/analyses/six_node_pattern/F/F1/out_analysed_network"  # F1 sandbox: outputs here only
dir.create(R5, showWarnings=FALSE)
lg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")
rn <- function(x) { r <- data.table::frank(x, na.last="keep", ties.method="average")
                    (r-1)/(max(r,na.rm=TRUE)-1) }          # rank-normalise to [0,1]

nodes <- fread(file.path(REV,"data/canonical_nodes.tsv"))
edges <- fread(file.path(REV,"data/canonical_edges.tsv"))
edges <- edges[edge_type != "miRNA_miRNA"]            # F1 sandbox: analysed network, 6,829 edges
stopifnot(nrow(edges) == 6829)
lg("nodes", nrow(nodes), "edges", nrow(edges))
P <- nodes[, .(name, type)]

# ---- T: topology -------------------------------------------------------------
deg <- rbind(edges[, .(name=source)], edges[, .(name=target)])[, .(degree=.N), by=name]
P <- merge(P, deg, by="name", all.x=TRUE); P[is.na(degree), degree := 0]
top <- fread(file.path(REV,"results/network_topology_hubs.csv"))
nc <- intersect(c("name","node"), names(top))[1]; setnames(top, nc, "name")
bc <- grep("betw", names(top), value=TRUE, ignore.case=TRUE)[1]
if (!is.na(bc)) P <- merge(P, top[, .(name, betweenness=get(bc))], by="name", all.x=TRUE)
# FFL participation in ANY role (regulator, intermediate or target), recomputed here so that
# target genes are counted too - the deposited participation file counts regulators only.
cores <- fread(file.path(REV,"results/v2/ffl_cores_coherence_corrected.csv"))
allroles <- rbind(cores[, .(name=regulator)], cores[, .(name=intermediate)], cores[, .(name=target)])
fps <- allroles[, .(ffl_cores=.N), by=name]
P <- merge(P, fps, by="name", all.x=TRUE)
P[is.na(ffl_cores), ffl_cores := 0]
# also record the role breakdown for the manuscript table
roleb <- merge(merge(cores[, .(as_regulator=.N), by=.(name=regulator)],
                     cores[, .(as_intermediate=.N), by=.(name=intermediate)], by="name", all=TRUE),
               cores[, .(as_target=.N), by=.(name=target)], by="name", all=TRUE)
P <- merge(P, roleb, by="name", all.x=TRUE)
for (cc in c("as_regulator","as_intermediate","as_target")) P[is.na(get(cc)), (cc) := 0]
P[, T_topology := rowMeans(cbind(rn(degree), rn(betweenness), rn(ffl_cores)), na.rm=TRUE)]

# ---- E: edge evidence quality ------------------------------------------------
ev <- fread(file.path(REV,"data/edge_evidence_tier.tsv"))
tcol <- grep("tier", names(ev), value=TRUE)[1]
evn <- rbind(ev[, .(name=source, tier=get(tcol))], ev[, .(name=target, tier=get(tcol))])
evs <- evn[, .(n_edges=.N, frac_strong=mean(tier=="strong", na.rm=TRUE)), by=name]
P <- merge(P, evs, by="name", all.x=TRUE)
P[, E_evidence := rn(frac_strong)]

# ---- D: differential expression ---------------------------------------------
de <- fread(file.path(REV,"results/BRCA_DEX_ALL_nodes.csv"))
P <- merge(P, de[, .(name=Gene, logFC, de_fdr=adj.P.Val)], by="name", all.x=TRUE)
P[, D_expression := rowMeans(cbind(rn(abs(logFC)), rn(-log10(pmax(de_fdr,1e-300)))), na.rm=TRUE)]

# ---- R: independent replication ---------------------------------------------
ext <- fread(file.path(REV,"results/external_GEO_DE_vs_TCGA.csv"))
rep <- ext[!is.na(logFC_ext) & !is.na(logFC_tcga),
           .(n_coh=.N, frac_same_dir=mean(sign(logFC_ext)==sign(logFC_tcga))), by=.(name=feature)]
P <- merge(P, rep, by="name", all.x=TRUE)
P[, R_replication := rn(frac_same_dir)]

# ---- M: multi-omic support ---------------------------------------------------
mm <- function(p, key, val, new) {
  f <- file.path(REV,p); if (!file.exists(f)) return(NULL)
  d <- fread(f); if (!all(c(key,val) %in% names(d))) return(NULL)
  setnames(d[, c(key,val), with=FALSE], c("name", new))
}
pr <- mm("results/multiomics/cptac_mRNA_protein_concordance.csv","gene","rho","prot_rho")
if (!is.null(pr)) P <- merge(P, pr, by="name", all.x=TRUE)
me <- mm("results/multiomics/brca_methylation_nodes.csv","node","delta_beta","meth_delta")
if (!is.null(me)) P <- merge(P, me, by="name", all.x=TRUE)
mu <- fread(file.path(REV,"results/multiomics/brca_mutation_frequency.csv"))
mk <- names(mu)[1]; mv <- grep("freq", names(mu), value=TRUE)[1]
if (!is.na(mv)) { setnames(mu, c(mk,mv), c("name","mut_freq")); P <- merge(P, mu[, .(name, mut_freq)], by="name", all.x=TRUE) }
P[, M_multiomic := rowMeans(cbind(rn(abs(prot_rho)), rn(abs(meth_delta)), rn(mut_freq)), na.rm=TRUE)]

# ---- C: clinical -------------------------------------------------------------
sv <- fread(file.path(REV,"results/survival_cox_hubs.csv"))
sk <- intersect(c("feature","gene","name"), names(sv))[1]; setnames(sv, sk, "name")
pc <- grep("^p$|p_val|pvalue", names(sv), value=TRUE)[1]
ec <- grep("endpoint|outcome", names(sv), value=TRUE)[1]
svb <- if (!is.na(ec)) sv[get(ec) %in% c("OS","PFI")] else sv
svm <- svb[, .(surv_minp = min(get(pc), na.rm=TRUE)), by=name]
P <- merge(P, svm, by="name", all.x=TRUE)
P[, C_clinical := rn(-log10(pmax(surv_minp, 1e-300)))]

# ---- F: functional -----------------------------------------------------------
dp <- fread(file.path(REV,"results/multiomics/depmap_essentiality_nodes.csv"))
dk <- names(dp)[1]; dv <- grep("mean_chronos_breast", names(dp), value=TRUE)[1]
if (!is.na(dv)) { setnames(dp, c(dk,dv), c("name","chronos")); P <- merge(P, dp[, .(name, chronos)], by="name", all.x=TRUE) }
P[, F_functional := rn(-chronos)]                    # more negative Chronos = more essential

# ---- composite ---------------------------------------------------------------
dom <- c("T_topology","E_evidence","D_expression","R_replication","M_multiomic","C_clinical","F_functional")
P[, n_domains := rowSums(!is.na(.SD)), .SDcols=dom]
P[, priority  := rowMeans(.SD, na.rm=TRUE), .SDcols=dom]
P[n_domains < 3, priority := NA]                      # require >=3 domains to be rankable
setorder(P, -priority)
P[, rank_overall := seq_len(.N)]
P[, rank_within_type := seq_len(.N), by=type]
fwrite(P, file.path(R5,"node_prioritisation_full.csv"))
lg("wrote node_prioritisation_full.csv:", nrow(P), "nodes;", sum(!is.na(P$priority)), "rankable")

show <- c("name","type","priority","n_domains","ffl_cores","as_regulator","as_target","degree",
          "logFC","de_fdr","frac_strong","surv_minp","meth_delta")
for (ty in c("TF","Gene","miRNA")) {
  lg(sprintf("\n================ TOP 10 %s ================", ty))
  print(P[type==ty & !is.na(priority)][1:10, ..show][, lapply(.SD, function(x) if(is.numeric(x)) signif(x,3) else x)])
}
fwrite(P[!is.na(priority)][, head(.SD,10), by=type], file.path(R5,"top10_per_class.csv"))
lg("\nwrote top10_per_class.csv")
lg("domain coverage:"); print(P[, lapply(.SD, function(x) sum(!is.na(x))), .SDcols=dom])
