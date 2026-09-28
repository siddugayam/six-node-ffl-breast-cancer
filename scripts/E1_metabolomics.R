#!/usr/bin/env Rscript
# =============================================================================
# E1_metabolomics.R  --  MEASURED metabolomics (not expression-inferred)
#
# Search performed (all endpoints verified live):
#  * Metabolomics Workbench REST. Exhaustive cross-study scan via the metstat
#    context (/rest/metstat/;;;Human;;;;) = 1,493 human study-source rows,
#    plus a full catalogue pull (/rest/study/study_id/ST/summary, 4,546 studies).
#    RESULT: ST000054 is the ONLY human breast-TISSUE tumour-vs-normal study with
#    a retrievable quantitative data table. Everything else with "breast" is a
#    cell line, blood/plasma/serum, breast milk, or adipose lipidomics.
#  * MetaboLights (EBI search + WS). Best hit MTBLS552 (72 IDC tissue samples,
#    malignant/benign/normal) -- BUT its MAF sample columns are EMPTY and only raw
#    .CDF files are deposited, so NO quantitative table exists. Recorded as a
#    failed retrieval, NOT silently dropped.
#
# DATASETS ACTUALLY ANALYSED
#  ST000054 / AN000092  breast TISSUE, GC-TOF-MS, 52 metabolites
#                       80 invasive tumour + 25 DCIS vs 20 reduction-mammoplasty
#  ST000355 / AN000580+AN000581  PLASMA, GC-MS + LC-MS, 135 BC vs 76 control
#  ST000356 / AN000582+AN000583  SERUM,  GC-MS + LC-MS, 103 BC vs 31 control
#
# Output: results/multiomics/breast_metabolomics.csv
# =============================================================================
suppressPackageStartupMessages({library(data.table)})

ROOT  <- "/path/to/revision"
CACHE <- file.path(ROOT,"cache/metabolomics")
OUT   <- file.path(ROOT,"results/multiomics")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
con <- file(file.path(ROOT,"logs/E1_metabolomics.log"), open="wt")
sink(con, split=TRUE); sink(con, type="message")
cat("=== E1_metabolomics.R ", format(Sys.time()), " ===\n")

## ---------------------------------------------------------------------------
## CURATED metabolite -> pathway class map (exact names as they appear in the
## deposited tables; curated by hand from the printed panels, not regex-guessed)
## ---------------------------------------------------------------------------
GLYCOLYSIS <- c("Lactate","Lactic acid","Pyruvate","pyruvic acid","Glucose",
                "Glucopyranose","D-Fructose","Fructose","Glyceraldehyde",
                "3-Phosphoglyceric acid, TMS","Glucose 6-phosphate",
                "Glyceraldehyde 3-phosphate","Dihydroxyacetone phosphate",
                "Glyceric acid, TMS","Glyceric acid")
ONE_CARBON <- c("Glycine","Serine","serine","serine minor","sarcosine","Sarcosine",
                "L-methionine","Methionine","Betaine","Choline","choline",
                "Cystathionine","Glycerophosphocholine","phospholic choline",
                "Homocysteic acid","Homoserine","L-Homoserine","Phosphoserine",
                "N,N-Dimethylglycine","N-formyl-glycine","N-Formyl-L-methionine",
                "Methionine sulfoxide")
# folate-dependent nucleotide synthesis == the "antifolate resistance" claim
ANTIFOLATE_NUC <- c("dUMP","Dihydrouracil","Dihydrothymine","uracil","Uracil",
                    "thymine","Ribothymidine","5-Methyldeoxycytidine","Cytidine",
                    "cytosine","Tetrahydroneopterin",
                    "5'-Phosphoribosyl-N-formylglycinamide","5-Phosphoribosylamine")
TCA <- c("Citric acid","Malate","Malic acid","Succinate","Succinic acid",
         "a-ketoglutarate","alpha-ketoglutarate","-Ketoglutaric acid",
         "alpha-Ketoglutaric acid","Oxaloacetate","2-Butenedioic acid")
PPP <- c("6-Phosphogluconic acid","D-Gluconic acid","Gluconate","Ribitol",
         "D-Ribofuranose","Arabitol","Threitol","Erythrose")
REDOX <- c("Glutathione","Cysteinylglycine","Dehydroascorbic acid","Hypotaurine",
           "Taurine","5-Oxoproline","5-oxoproline","oxoproline",
           "2-pyrrolidone-5-carboxylic acid","Pyrrolidonecarboxylic acid")

# case-insensitive exact matching: the deposited tables are inconsistent about
# capitalisation ("Malic acid" in ST000356 vs "malic acid" in ST000054).
classify <- function(m){
  m <- tolower(trimws(m))
  if (m %in% tolower(GLYCOLYSIS))     return("glycolysis")
  if (m %in% tolower(ONE_CARBON))     return("one_carbon_folate")
  if (m %in% tolower(ANTIFOLATE_NUC)) return("folate_dependent_nucleotide")
  if (m %in% tolower(TCA))            return("TCA_OXPHOS")
  if (m %in% tolower(PPP))            return("pentose_phosphate")
  if (m %in% tolower(REDOX))          return("redox_glutathione")
  "other"
}

## ---------------------------------------------------------------------------
## generic Wilcoxon tester on a long table (sample, group, metabolite, abundance)
## ---------------------------------------------------------------------------
test_pair <- function(dt, g1, g2){
  d <- dt[group %in% c(g1,g2)]
  res <- d[, {
    x <- abundance[group==g1]; y <- abundance[group==g2]
    x <- x[is.finite(x)];      y <- y[is.finite(y)]
    if (length(x)>=5 && length(y)>=5 && sd(c(x,y))>0) {
      wt <- suppressWarnings(wilcox.test(x, y, exact=FALSE))
      mx <- median(x); my <- median(y)
      .(n_group1=length(x), n_group2=length(y),
        median_group1=mx, median_group2=my,
        log2FC = if (is.finite(mx) && is.finite(my) && mx>0 && my>0) log2(mx/my) else NA_real_,
        W=unname(wt$statistic), p=wt$p.value)
    } else .(n_group1=length(x), n_group2=length(y), median_group1=NA_real_,
             median_group2=NA_real_, log2FC=NA_real_, W=NA_real_, p=NA_real_)
  }, by=.(metabolite)]
  res[, contrast := paste0(g1," vs ",g2)]
  res[, FDR := p.adjust(p,"BH")]
  res[]
}

read_mw <- function(analysis_id){
  f <- file.path(CACHE, paste0("dt_",analysis_id,".txt"))
  if (analysis_id=="AN000092") f <- file.path(CACHE,"dt2.txt")
  stopifnot(file.exists(f))
  d <- fread(f, sep="\t", header=TRUE, na.strings=c("","NA"))
  setnames(d, 1:2, c("sample","class")); d
}

all_res <- list()

## ============================ ST000054  BREAST TISSUE ======================
cat("\n########## ST000054 / AN000092 : BREAST TISSUE (GC-TOF-MS) ##########\n")
d <- read_mw("AN000092")
d[, pool := grepl("pool aliquot:yes", class)]
d[, treatment := trimws(sub(".*treatment:([^|]*)\\|.*","\\1", class))]
cat("non-pool sample counts:\n"); print(d[pool==FALSE, .N, by=treatment][order(-N)])
d2 <- d[pool==FALSE & treatment %in% c("Basal","LumA","LumB","HER2","DCIS","Reduction")]
d2[, group := fifelse(treatment=="Reduction","Normal",
             fifelse(treatment=="DCIS","DCIS","InvasiveTumour"))]
mets <- setdiff(names(d), c("sample","class","pool","treatment"))
L <- melt(d2[, c("sample","group",mets), with=FALSE], id.vars=c("sample","group"),
          variable.name="metabolite", value.name="abundance")
L[, abundance := as.numeric(abundance)]; L[, metabolite := as.character(metabolite)]
r <- rbind(test_pair(L,"InvasiveTumour","Normal"), test_pair(L,"DCIS","Normal"))
r[, `:=`(repository="Metabolomics Workbench", accession="ST000054",
         analysis_id="AN000092", matrix_type="breast tissue",
         platform="GC-TOF-MS (peak height)")]
all_res$ST000054 <- r
cat("\nn metabolites tested:", uniqueN(r$metabolite), "\n")

## ============================ ST000355  PLASMA =============================
cat("\n########## ST000355 : PLASMA (breast cancer vs control) ##########\n")
for (aid in c("AN000580","AN000581")) {
  d <- read_mw(aid)
  d[, group := fifelse(grepl("Diagnosis:Control|Diagnosis:control", class), "Control", "BreastCancer")]
  cat(aid, "groups:\n"); print(d[, .N, by=group])
  mets <- setdiff(names(d), c("sample","class","group"))
  L <- melt(d[, c("sample","group",mets), with=FALSE], id.vars=c("sample","group"),
            variable.name="metabolite", value.name="abundance")
  L[, abundance := as.numeric(abundance)]; L[, metabolite := as.character(metabolite)]
  r <- test_pair(L,"BreastCancer","Control")
  r[, `:=`(repository="Metabolomics Workbench", accession="ST000355", analysis_id=aid,
           matrix_type="plasma (systemic, NOT tumour tissue)",
           platform=ifelse(aid=="AN000580","GC-MS (peak intensity)","LC-MS (peak intensity)"))]
  all_res[[aid]] <- r
}

## ============================ ST000356  SERUM ==============================
cat("\n########## ST000356 : SERUM (breast cancer vs control) ##########\n")
for (aid in c("AN000582","AN000583")) {
  d <- read_mw(aid)
  d[, group := fifelse(grepl("Diagnosis:Control|Diagnosis:control", class), "Control", "BreastCancer")]
  cat(aid, "groups:\n"); print(d[, .N, by=group])
  mets <- setdiff(names(d), c("sample","class","group"))
  L <- melt(d[, c("sample","group",mets), with=FALSE], id.vars=c("sample","group"),
            variable.name="metabolite", value.name="abundance")
  L[, abundance := as.numeric(abundance)]; L[, metabolite := as.character(metabolite)]
  r <- test_pair(L,"BreastCancer","Control")
  r[, `:=`(repository="Metabolomics Workbench", accession="ST000356", analysis_id=aid,
           matrix_type="serum (systemic, NOT tumour tissue)",
           platform=ifelse(aid=="AN000582","GC-MS (peak intensity)","LC-MS (peak intensity)"))]
  all_res[[aid]] <- r
}

## ============================ combine ======================================
A <- rbindlist(all_res, fill=TRUE)
A[, pathway_class := sapply(metabolite, classify)]
A[, direction := fifelse(is.na(p),NA_character_,
      fifelse(p<0.05 & !is.na(log2FC) & log2FC>0,"UP_in_group1",
      fifelse(p<0.05 & !is.na(log2FC) & log2FC<0,"DOWN_in_group1","ns")))]
A[, evidence_class := "MEASURED metabolite abundance"]
setcolorder(A, c("repository","accession","analysis_id","matrix_type","platform",
                 "contrast","metabolite","pathway_class","n_group1","n_group2",
                 "median_group1","median_group2","log2FC","W","p","FDR",
                 "direction","evidence_class"))
setorder(A, accession, analysis_id, contrast, p)
fwrite(A, file.path(OUT,"breast_metabolomics.csv"))
cat("\nWROTE", file.path(OUT,"breast_metabolomics.csv"), " rows:", nrow(A), "\n")

## ============================ report =======================================
show <- function(cls, title){
  cat("\n=========== ", title, " ===========\n")
  s <- A[pathway_class==cls & !is.na(p)]
  if (!nrow(s)) { cat("  NO metabolites of this class measured in any retrieved dataset.\n"); return(invisible()) }
  print(s[order(p), .(accession, matrix_type=substr(matrix_type,1,14), contrast,
                      metabolite, n1=n_group1, n2=n_group2,
                      log2FC=round(log2FC,3), p=signif(p,3), FDR=signif(FDR,3), direction)],
        nrows=100)
}
show("glycolysis","GLYCOLYSIS  (manuscript claim i: let-7b-HK2 modulates glycolysis)")
show("one_carbon_folate","ONE-CARBON / FOLATE CYCLE (claim ii)")
show("folate_dependent_nucleotide","FOLATE-DEPENDENT NUCLEOTIDE SYNTHESIS (antifolate claim ii)")
show("TCA_OXPHOS","TCA / OXPHOS")
show("pentose_phosphate","PENTOSE PHOSPHATE")
show("redox_glutathione","REDOX / GLUTATHIONE (relevant to NRF2 claim iii)")

cat("\n=========== CLASS COVERAGE PER DATASET ===========\n")
print(dcast(A[!duplicated(paste(accession,analysis_id,metabolite))],
            accession + analysis_id + matrix_type ~ pathway_class,
            value.var="metabolite", fun.aggregate=length))

cat("\n=========== TISSUE-ONLY SUMMARY (ST000054) ===========\n")
t54 <- A[accession=="ST000054" & contrast=="InvasiveTumour vs Normal" & !is.na(p)]
cat("metabolites tested:", nrow(t54), "; FDR<0.05:", sum(t54$FDR<0.05), "\n")
print(t54[FDR<0.05][order(p), .(metabolite, pathway_class, log2FC=round(log2FC,3),
                               p=signif(p,3), FDR=signif(FDR,3))])
cat("\nGLYCOLYTIC INTERMEDIATES IN THE ONLY BREAST-TISSUE DATASET:",
    nrow(A[accession=="ST000054" & pathway_class=="glycolysis"]), "\n")

cat("\n########## SAMPLE OVERLAP WITH TRANSCRIPTOMICS ##########\n")
cat("NONE. ST000054, ST000355 and ST000356 have no matched RNA-seq/array data for the\n")
cat("same subjects, and none of their sample IDs are TCGA barcodes. The requested\n")
cat("metabolite-vs-network-gene-expression correlation is therefore NOT computable.\n")
cat("No gene-expression proxy was substituted for it.\n")
cat("\nDONE ", format(Sys.time()), "\n")
sink(type="message"); sink(); close(con)
