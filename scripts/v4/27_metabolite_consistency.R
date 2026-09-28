## 27_metabolite_consistency.R --------------------------------------------
## Part F: metabolite-level evidence. Re-derived independently from the cached
## Metabolomics Workbench ST000054 / AN000092 mwTab already on disk
## (cache/metabolomics/an92.mwtab.txt) -- no re-download. Then asked whether the
## MEASURED metabolite changes agree with the TRANSCRIPTIONAL metabolic scores.
suppressPackageStartupMessages({library(data.table)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
set.seed(1)
mw <- readLines(file.path(ROOT,"cache/metabolomics/an92.mwtab.txt"), warn=FALSE)
i0 <- grep("^MS_METABOLITE_DATA_START$", mw); i1 <- grep("^MS_METABOLITE_DATA_END$", mw)
stopifnot(length(i0)==1, length(i1)==1)
blk <- mw[(i0+1):(i1-1)]
hdr <- strsplit(blk[1], "\t")[[1]]; fac <- strsplit(blk[2], "\t")[[1]]
stopifnot(hdr[1]=="Samples", fac[1]=="Factors")
samples <- hdr[-1]; factors <- fac[-1]
grp_raw <- trimws(sub(".*treatment:([^|]*)\\|.*", "\\1", factors))
pool <- grepl("pool aliquot:yes", factors)
cat("ST000054 / AN000092 samples in data block:", length(samples), "\n")
print(table(grp_raw, pool))
dat <- fread(text=paste(blk[-c(1,2)], collapse="\n"), header=FALSE, sep="\t")
mets <- dat[[1]]; X <- as.matrix(dat[, -1]); rownames(X) <- mets; colnames(X) <- samples
storage.mode(X) <- "numeric"
cat("rows in data block:", nrow(X), "\n")
## drop GC-MS internal standards (C08..C30 IST) and the phthalate contaminant
is_ist <- grepl("^C[0-9]+ IST$", rownames(X)); is_cont <- grepl("dioctylphtalate", rownames(X))
cat("internal standards dropped:", sum(is_ist), " contaminant rows dropped:", sum(is_cont), "\n")
X <- X[!is_ist & !is_cont, , drop=FALSE]
cat("biological metabolite features analysed:", nrow(X), "\n")

grp <- ifelse(grp_raw %in% c("Basal","HER2","LumA","LumB"), "InvasiveTumour",
       ifelse(grp_raw == "DCIS", "DCIS",
       ifelse(grp_raw == "Reduction", "Normal", NA)))
use <- !pool & !is.na(grp)
cat("\nnon-pool biological samples used:\n"); print(table(grp[use]))

test_grp <- function(g1, g2) {
  a <- X[, use & grp==g1, drop=FALSE]; b <- X[, use & grp==g2, drop=FALSE]
  rbindlist(lapply(rownames(X), function(m) {
    x <- a[m,]; y <- b[m,]
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    wt <- suppressWarnings(wilcox.test(x, y))
    data.table(metabolite=m, contrast=paste(g1,"vs",g2), n1=length(x), n2=length(y),
      median1=median(x), median2=median(y),
      log2FC=log2(median(x)/median(y)), p=wt$p.value)
  }))
}
## GLOBAL-SCALING CHECK: GC-TOF-MS peak heights are not total-signal normalised.
cat("\n-- global scaling check on raw peak heights --\n")
raw_lfc <- test_grp("InvasiveTumour","Normal")
cat("metabolites with log2FC > 0 (raw):", sum(raw_lfc$log2FC>0, na.rm=TRUE), "of", nrow(raw_lfc),
    "; median log2FC =", round(median(raw_lfc$log2FC, na.rm=TRUE),3), "\n")
cat("median per-sample total peak height, tumour vs normal:",
    round(median(colSums(X[, use & grp=="InvasiveTumour"], na.rm=TRUE))), "vs",
    round(median(colSums(X[, use & grp=="Normal"], na.rm=TRUE))), "\n")
X_raw <- X
## median-normalise each sample (standard for GC-MS peak heights) and repeat
X <- sweep(X, 2, apply(X, 2, median, na.rm=TRUE), "/")
norm_lfc <- test_grp("InvasiveTumour","Normal")
cat("after per-sample median normalisation: log2FC>0 in", sum(norm_lfc$log2FC>0, na.rm=TRUE),
    "of", nrow(norm_lfc), "; median log2FC =", round(median(norm_lfc$log2FC, na.rm=TRUE),3), "\n")
X <- X_raw   ## primary analysis stays on the raw peak heights, as published
mres_norm <- norm_lfc; mres_norm[, FDR := p.adjust(p,"BH")]
fwrite(mres_norm, file.path(RES,"metabolic_metabolite_ST000054_mediannorm.csv"))
cat("median-normalised: significant at FDR<0.05:", mres_norm[FDR<0.05, .N], "->",
    paste(mres_norm[FDR<0.05][order(p), metabolite], collapse=", "), "\n\n")

mres <- rbind(test_grp("InvasiveTumour","Normal"), test_grp("DCIS","Normal"))
mres[, FDR := p.adjust(p,"BH"), by=contrast]
setorder(mres, contrast, p)
fwrite(mres, file.path(RES,"metabolic_metabolite_ST000054.csv"))
cat("\n== ST000054 InvasiveTumour vs Normal: significant at FDR<0.05 ==\n")
print(mres[contrast=="InvasiveTumour vs Normal" & FDR<0.05,
  .(metabolite, n1, n2, log2FC=round(log2FC,3), p=signif(p,3), FDR=signif(FDR,3))])
cat("\nDCIS vs Normal significant at FDR<0.05:",
    mres[contrast=="DCIS vs Normal" & FDR<0.05, .N], "\n")

## cross-check against the earlier pass (results/multiomics/breast_metabolomics.csv)
old <- fread(file.path(ROOT,"results/multiomics/breast_metabolomics.csv"))
o <- old[accession=="ST000054" & contrast=="InvasiveTumour vs Normal", .(metabolite, log2FC_old=log2FC, p_old=p)]
m <- merge(mres[contrast=="InvasiveTumour vs Normal"], o, by="metabolite")
cat("\ncross-check vs earlier pass: n matched", nrow(m),
    " Pearson(log2FC) =", round(cor(m$log2FC, m$log2FC_old, use="complete.obs"),4),
    " max |dp| =", signif(max(abs(m$p - m$p_old), na.rm=TRUE),3), "\n")

## ---- what was measured? ------------------------------------------------
cat("\n== all", nrow(X), "measured metabolites ==\n"); print(sort(rownames(X)))
glycolytic <- c("glucose","glucose-6-phosphate","fructose-6-phosphate","fructose-1,6",
  "dihydroxyacetone","glyceraldehyde-3","3-phosphoglycerate","2-phosphoglycerate",
  "phosphoenolpyruvate","pyruv","lactate","lactic acid")
hit <- unlist(lapply(glycolytic, function(p) grep(p, rownames(X), value=TRUE, ignore.case=TRUE)))
cat("\nglycolytic intermediates present in the platform:",
    if (length(hit)) paste(unique(hit), collapse=", ") else "NONE", "\n")
folate <- grep("folate|folic|thymid|dUMP|dTMP|methionine|homocyst|serine|glycine",
               rownames(X), value=TRUE, ignore.case=TRUE)
cat("one-carbon / folate-cycle metabolites present:",
    if (length(folate)) paste(folate, collapse=", ") else "NONE", "\n")

## ---- consistency with the transcriptional scores ------------------------
## Curated metabolite -> pathway-score mapping, only for metabolites whose
## pathway assignment is unambiguous. Everything else is left unmapped and said so.
tvn <- fread(file.path(RES,"metabolic_tumour_vs_normal.csv"))[score=="ssGSEA"]
map <- rbind(
 ## --- unambiguous single-pathway assignments ---------------------------
 data.table(metabolite=c("citric acid","malic acid","alpha-ketoglutarate"),
            score_set="CURATED|CURATED_TCA_CYCLE", assignment="unambiguous"),
 data.table(metabolite="glutamic acid 2TMS",
            score_set="REACTOME|R-HSA-8964539_GLUTAMATE_AND_GLUTAMINE_METABOLISM", assignment="unambiguous"),
 data.table(metabolite=c("aspartate minor","N-acetyl-L-aspartic acid 1"),
            score_set="KEGG|hsa00250_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM", assignment="unambiguous"),
 data.table(metabolite=c("serine","serine minor","threonine","threonine minor"),
            score_set="KEGG|hsa00260_GLYCINE_SERINE_AND_THREONINE_METABOLISM", assignment="unambiguous"),
 data.table(metabolite=c("isoleucine","valine"),
            score_set="KEGG|hsa00280_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION", assignment="unambiguous"),
 data.table(metabolite=c("proline","trans-4-hydroxy-L-proline minor1","creatine major dehydrated TMS3"),
            score_set="KEGG|hsa00330_ARGININE_AND_PROLINE_METABOLISM", assignment="unambiguous"),
 data.table(metabolite="lysine", score_set="KEGG|hsa00310_LYSINE_DEGRADATION", assignment="unambiguous"),
 data.table(metabolite="oxoproline", score_set="KEGG|hsa00480_GLUTATHIONE_METABOLISM", assignment="unambiguous"),
 data.table(metabolite="urea", score_set="KEGG|hsa00220_ARGININE_BIOSYNTHESIS", assignment="unambiguous"),
 data.table(metabolite="beta-hydroxybutyric acid",
            score_set="REACTOME|R-HSA-74182_KETONE_BODY_METABOLISM", assignment="unambiguous"),
 data.table(metabolite="succinate semialdehyde meox1",
            score_set="KEGG|hsa00650_BUTANOATE_METABOLISM", assignment="unambiguous"),
 data.table(metabolite="oxalic acid",
            score_set="KEGG|hsa00630_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM", assignment="unambiguous"),
 data.table(metabolite="sorbitol", score_set="KEGG|hsa00051_FRUCTOSE_AND_MANNOSE_METABOLISM", assignment="unambiguous"),
 data.table(metabolite="ribitol", score_set="KEGG|hsa00040_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS", assignment="unambiguous"),
 data.table(metabolite="cholesterol", score_set="REACTOME|R-HSA-191273_CHOLESTEROL_BIOSYNTHESIS", assignment="unambiguous"),
 data.table(metabolite=c("glycerol","1-monopalmitin","1-monostearin"),
            score_set="KEGG|hsa00561_GLYCEROLIPID_METABOLISM", assignment="unambiguous"),
 data.table(metabolite=c("pinitol"), score_set="KEGG|hsa00562_INOSITOL_PHOSPHATE_METABOLISM", assignment="unambiguous"),
 ## --- free fatty acids: pool size is ambiguous between synthesis and
 ##     degradation, so both arms are scored and the pair is FLAGGED --------
 data.table(metabolite=c("palmitic acid","stearic acid","oleic acid","linoleic acid",
                         "myristic acid","palmitoleic acid","lauric acid","behenic acid",
                         "arachidiic acid","heptadecanoic acid","pentadecanoic acid",
                         "caprylic acid","pelargonic acid"),
            score_set="KEGG|hsa00061_FATTY_ACID_BIOSYNTHESIS", assignment="ambiguous_FA_pool"),
 data.table(metabolite=c("palmitic acid","stearic acid","oleic acid","linoleic acid",
                         "myristic acid","palmitoleic acid","lauric acid","behenic acid",
                         "arachidiic acid","heptadecanoic acid","pentadecanoic acid",
                         "caprylic acid","pelargonic acid"),
            score_set="KEGG|hsa00071_FATTY_ACID_DEGRADATION", assignment="ambiguous_FA_pool"))
## deliberately NOT mapped (exogenous / xenobiotic / plant-derived / no clear human
## pathway on this platform): alpha tocopherol, benzoic acid, salicylic acid,
## chlorogenic acid, shikimic acid, maleimide, aminomalonic acid,
## 2-hydroxybutanoic acid, octadecanol, phosphate.
mm <- mres[contrast=="InvasiveTumour vs Normal"]
## match metabolite names case/whitespace-insensitively
mm[, key := tolower(trimws(metabolite))]; map[, key := tolower(trimws(metabolite))]
cons <- merge(map[, .(key, metabolite_mapped=metabolite, score_set, assignment)], mm[, .(key, log2FC, p, FDR)], by="key")
cons <- merge(cons, tvn[, .(score_set=set_id, score_d=cohens_d, score_FDR=FDR_limma,
                            score_available=TRUE)], by="score_set", all.x=TRUE)
cons[, metabolite_direction := fifelse(log2FC>0,"UP in tumour","DOWN in tumour")]
cons[, transcript_direction := fifelse(score_d>0,"UP in tumour","DOWN in tumour")]
cons[, agree := metabolite_direction == transcript_direction]
cons[, metabolite_significant := FDR < 0.05]
setorder(cons, -metabolite_significant, FDR)
fwrite(cons, file.path(RES,"metabolic_metabolite_vs_transcript_consistency.csv"))
cat("\n== metabolite vs transcriptional-score direction ==\n")
print(cons[assignment=="unambiguous", .(metabolite=metabolite_mapped, score=substr(score_set,1,46),
   met_log2FC=round(log2FC,2), met_FDR=signif(FDR,2), met_sig=metabolite_significant,
   score_d=round(score_d,2), score_FDR=signif(score_FDR,2), agree)])
pri <- cons[assignment=="unambiguous" & !is.na(agree)]
cat("\n-- PRIMARY test, unambiguous assignments only --\n")
cat("direction agreement:", sum(pri$agree), "/", nrow(pri), "\n")
bt <- binom.test(sum(pri$agree), nrow(pri), 0.5)
cat("binomial vs 50%: p =", signif(bt$p.value,3), "\n")
prs <- pri[metabolite_significant==TRUE]
cat("restricted to metabolites significant at FDR<0.05:", sum(prs$agree), "/", nrow(prs), "\n")
amb <- cons[assignment=="ambiguous_FA_pool" & !is.na(agree)]
cat("\n-- free-fatty-acid pool (ambiguous synthesis/degradation), reported separately --\n")
cat("vs FATTY_ACID_BIOSYNTHESIS score: agree", amb[grepl("hsa00061",score_set), sum(agree)], "/", amb[grepl("hsa00061",score_set), .N], "\n")
cat("vs FATTY_ACID_DEGRADATION score:  agree", amb[grepl("hsa00071",score_set), sum(agree)], "/", amb[grepl("hsa00071",score_set), .N], "\n")
cat("\nmetabolite features measured:", nrow(mm),
    "| mapped unambiguously:", uniqueN(cons[assignment=="unambiguous", key]),
    "| mapped as ambiguous FA pool:", uniqueN(cons[assignment=="ambiguous_FA_pool", key]),
    "| deliberately unmapped:", nrow(mm) - uniqueN(cons$key), "\n")
## repeat the concordance on the median-normalised metabolite fold-changes
mn <- mres_norm[, .(key=tolower(trimws(metabolite)), log2FC_norm=log2FC, FDR_norm=FDR)]
cons2 <- merge(cons, mn, by="key")
cons2[, agree_norm := (log2FC_norm>0) == (score_d>0)]
p2 <- cons2[assignment=="unambiguous" & !is.na(agree_norm)]
cat("\n-- concordance after per-sample median normalisation of the metabolite data --\n")
cat("direction agreement:", sum(p2$agree_norm), "/", nrow(p2),
    "  binomial p =", signif(binom.test(sum(p2$agree_norm), nrow(p2), 0.5)$p.value,3), "\n")
fwrite(cons2, file.path(RES,"metabolic_metabolite_vs_transcript_consistency.csv"))
cat("DONE 27\n")
