#!/usr/bin/env Rscript
## C) Assemble perturbation evidence for TF -> COL1A1 / COL3A1
suppressPackageStartupMessages({library(data.table)})
REV <- "/path/to/revision"
OUT <- file.path(REV,"results/multiomics")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

## ---------- (i) LINCS L1000 ----------
L <- fread(file.path(OUT,"lincs_TF_signatures_collagen.csv"))
msg("LINCS signature rows:", nrow(L))
Lcd <- L[value_col=="CD-coefficient"]                     # individual signatures
Lmean <- L[grepl("mean_sig", local_id)]                   # CRISPR-KO consensus
msg("individual CD signatures:", nrow(Lcd), "| consensus mean sigs:", nrow(Lmean))

stats <- rbindlist(lapply(c("COL1A1","COL3A1"), function(g){
  rbindlist(lapply(split(Lcd, by=c("tf","pert_type"), drop=TRUE), function(s){
    v <- s[[paste0(g,"_value")]]; p <- s[[paste0(g,"_pctile")]]
    v <- v[is.finite(v)]; p <- p[is.finite(p)]
    if(length(v) < 3) return(NULL)
    wt <- wilcox.test(v, mu=0)
    bt <- binom.test(sum(v>0), length(v), 0.5)
    data.table(resource="LINCS L1000 (SigCom LINCS, 2021 CD-coefficient signatures)",
      tf=s$tf[1], perturbation=s$pert_type[1], readout_gene=g,
      n_signatures=length(v), n_cell_lines=uniqueN(s$cell_line),
      mean_CD=mean(v), median_CD=median(v), median_pctile=median(p),
      frac_up=mean(v>0), wilcoxon_p_vs0=wt$p.value, binom_p_vs_half=bt$p.value)
  }))
}))
## expected direction under the manuscript model (TF activates collagen):
stats[, expected_direction := fifelse(perturbation=="Overexpression","up (positive CD)","down (negative CD)")]
stats[, observed_direction := fifelse(median_CD>0,"up (positive CD)","down (negative CD)")]
stats[, supports_activation := expected_direction==observed_direction]
stats[, wilcoxon_q := p.adjust(wilcoxon_p_vs0, method="BH")]
setorder(stats, readout_gene, tf, perturbation)
print(stats[, .(tf,perturbation,readout_gene,n_signatures,n_cell_lines,
                median_CD=signif(median_CD,3), median_pctile=round(median_pctile,3),
                frac_up=round(frac_up,3), p=signif(wilcoxon_p_vs0,3), q=signif(wilcoxon_q,3),
                supports_activation)])
fwrite(stats, file.path(OUT,"lincs_TF_collagen_stats.csv"))
msg("wrote lincs_TF_collagen_stats.csv rows:", nrow(stats))

## consensus mean signature table
cons <- Lmean[, .(resource="LINCS L1000 CRISPR-KO consensus signature (SigCom LINCS)",
                  tf, perturbation=pert_type, signature=local_id, n_genes,
                  COL1A1_value, COL1A1_pctile, COL3A1_value, COL3A1_pctile)]
print(cons)

## ---------- (iii) ChIP-seq / curated / GEO libraries ----------
M <- fread(file.path(OUT,"chipseq_TF_collagen_membership.csv"))
msum <- M[, .(n_sets=.N, n_sets_containing_gene=sum(present),
              sets_containing_gene=paste(set_name[present], collapse=" ; ")),
          by=.(resource, tf, readout_gene)]
setorder(msum, readout_gene, resource, tf)
fwrite(msum, file.path(OUT,"chipseq_TF_collagen_summary.csv"))
msg("wrote chipseq_TF_collagen_summary.csv rows:", nrow(msum))

## ---------- TRRUST v2 raw (mode of regulation + PMIDs) ----------
tr <- fread(file.path(REV,"data/db/trrust_human.tsv"), header=FALSE,
            col.names=c("TF","target","mode","pmids"))
trsub <- tr[TF %in% c("NFKB1","RELA","SP1","ETS1") & target %like% "^COL"]
print(trsub)
fwrite(trsub, file.path(OUT,"trrust_TF_collagen_records.csv"))
msg("wrote trrust_TF_collagen_records.csv rows:", nrow(trsub))

## ---------- unified evidence table ----------
rows <- list()
add <- function(...) rows[[length(rows)+1]] <<- data.table(...)

for(i in seq_len(nrow(stats))){
  s <- stats[i]
  add(resource      = s$resource,
      evidence_type = "perturbation transcriptomics (L1000 CD coefficients)",
      tf            = s$tf,
      perturbation  = s$perturbation,
      readout_gene  = s$readout_gene,
      n_units       = s$n_signatures,
      effect        = signif(s$median_CD,4),
      effect_units  = "median CD coefficient across signatures",
      percentile    = round(s$median_pctile,3),
      p_value       = signif(s$wilcoxon_p_vs0,3),
      q_value       = signif(s$wilcoxon_q,3),
      direction     = s$observed_direction,
      supports_TF_activates_collagen = s$supports_activation,
      detail        = paste0(s$n_signatures," signatures, ",s$n_cell_lines,
                             " cell lines; frac up = ", round(s$frac_up,3)))
}
for(i in seq_len(nrow(cons))){
  cc <- cons[i]
  for(g in c("COL1A1","COL3A1")){
    v <- cc[[paste0(g,"_value")]]; pc <- cc[[paste0(g,"_pctile")]]
    add(resource="LINCS L1000 CRISPR-KO consensus (SigCom LINCS mean signature)",
        evidence_type="perturbation transcriptomics (consensus signature)",
        tf=cc$tf, perturbation="CRISPR Knockout", readout_gene=g, n_units=1L,
        effect=signif(v,4), effect_units="consensus signature value",
        percentile=round(pc,3), p_value=NA_real_, q_value=NA_real_,
        direction=ifelse(v>0,"up (positive)","down (negative)"),
        supports_TF_activates_collagen = (v<0),
        detail=cc$signature)
  }
}
for(i in seq_len(nrow(msum))){
  s <- msum[i]
  add(resource=s$resource, evidence_type=
        fifelse(s$resource %in% c("ChEA_2022","ENCODE_TF_ChIP-seq_2015",
                                  "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"),
                "TF ChIP-seq target set",
        fifelse(s$resource=="TRRUST_Transcription_Factors_2019","curated literature TF-target",
        fifelse(s$resource=="ARCHS4_TFs_Coexp","bulk RNA-seq co-expression (ARCHS4)",
        fifelse(s$resource=="TF_Perturbations_Followed_by_Expression","GEO TF perturbation signature",
        fifelse(s$resource=="LINCS_L1000_CRISPR_KO_Consensus_Sigs","LINCS CRISPR-KO consensus up/down set",
                "text co-occurrence"))))),
      tf=s$tf, perturbation=NA_character_, readout_gene=s$readout_gene,
      n_units=s$n_sets, effect=s$n_sets_containing_gene,
      effect_units="number of gene sets containing the readout gene",
      percentile=NA_real_, p_value=NA_real_, q_value=NA_real_,
      direction=NA_character_,
      supports_TF_activates_collagen=NA,
      detail=ifelse(s$sets_containing_gene=="", "gene absent from all sets", s$sets_containing_gene))
}
for(i in seq_len(nrow(trsub))){
  s <- trsub[i]
  add(resource="TRRUST v2 (raw, human)", evidence_type="curated literature TF-target",
      tf=s$TF, perturbation=NA_character_, readout_gene=s$target, n_units=1L,
      effect=NA_real_, effect_units=paste("mode of regulation:", s$mode),
      percentile=NA_real_, p_value=NA_real_, q_value=NA_real_,
      direction=s$mode,
      supports_TF_activates_collagen = (s$mode=="Activation"),
      detail=paste("PMIDs:", s$pmids))
}
ev <- rbindlist(rows, fill=TRUE)
fwrite(ev, file.path(OUT,"perturbation_TF_collagen.csv"))
msg("wrote perturbation_TF_collagen.csv rows:", nrow(ev))
msg("DONE")
