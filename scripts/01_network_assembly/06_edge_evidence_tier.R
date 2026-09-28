#!/usr/bin/env Rscript
# E) Evidence tier for every miRNA_target edge in the canonical network.
#    The manuscript pooled miRWalk (PREDICTED) with miRTarBase (VALIDATED) and treated the
#    union as one homogeneous edge set. Here every miRNA-target pair is checked against
#    multiMiR's *validated* tables (miRTarBase / TarBase / miRecords) and tiered.
#
#  IMPORTANT SEMANTICS (verified against the retrieved records, not assumed):
#    * TarBase support_type "negative"/"positive" encodes the DIRECTION of regulation
#      (miRNA down- vs up-regulates the target); it does NOT mean the interaction was
#      refuted. e.g. hsa-miR-34a-5p -> AXL, Luciferase + Western blot, support_type
#      "negative". These rows are therefore kept as supporting evidence.
#    * Only miRTarBase "Non-Functional MTI" flags an interaction reported as NOT functional.
#    * Assay tokens are classified individually after splitting on "//" and ";".
#         strong = low-throughput causal assays (luciferase/reporter, Western/immunoblot,
#                  qRT-PCR, Northern blot, ELISA, IHC/IF/ICC, in-situ hybridisation)
#         weak   = high-throughput binding/expression assays (HITS-CLIP, PAR-CLIP, CLASH,
#                  qCLASH, microarray, RNA-Seq, degradome, pSILAC, proteomics, AGO-IP,
#                  biotin pull-down, chimeric fragments) or miRTarBase "Functional MTI (Weak)"
suppressPackageStartupMessages({library(multiMiR); library(data.table)})
BASE  <- "/path/to/revision"
CACHE <- file.path(BASE,"cache","multimir"); dir.create(CACHE, showWarnings=FALSE, recursive=TRUE)
source(file.path(BASE,"scripts","01_network_assembly","00_mirna_canon.R"))
LOG <- file.path(BASE,"logs","06_edge_evidence_tier.log")
say <- function(...) { m <- paste0(format(Sys.time(),"%H:%M:%S")," | ",paste0(...,collapse=""))
                       cat(m,"\n"); cat(m,"\n",file=LOG,append=TRUE); flush.console() }
pr <- function(x){ print(x); sink(LOG,append=TRUE); print(x); sink() }
cat("", file=LOG)

edges <- fread(file.path(BASE,"data","canonical_edges.tsv"), sep="\t", header=TRUE)
mt <- edges[edge_type=="miRNA_target", .(source, target, edge_type)]
say("miRNA_target edges in canonical network = ", nrow(mt),
    " ; distinct miRNAs=", uniqueN(mt$source), " ; distinct targets=", uniqueN(mt$target))
tg <- sort(unique(mt$target))

## ---- fetch validated interactions, batched + cached ----
B <- 20; batches <- split(tg, ceiling(seq_along(tg)/B))
say("multiMiR validated query: ", length(batches), " batches of <= ", B, " genes")
all <- list()
for (i in seq_along(batches)) {
  f <- file.path(CACHE, sprintf("validated_batch_%03d.rds", i))
  if (file.exists(f)) d <- readRDS(f) else {
    ok <- FALSE; att <- 0
    while (!ok && att < 4) { att <- att + 1
      d <- tryCatch({ r <- get_multimir(org="hsa", target=batches[[i]], table="validated",
                                        summary=FALSE, legacy.out=FALSE); as.data.table(r@data) },
                    error=function(e){ say("  batch ",i," attempt ",att," ERROR: ",conditionMessage(e)); NULL })
      if (!is.null(d)) ok <- TRUE else Sys.sleep(20) }
    if (!ok) stop("multiMiR failed for batch ", i)
    saveRDS(d, f); Sys.sleep(1) }
  all[[i]] <- d
}
V <- rbindlist(all, fill=TRUE)
say("total validated rows retrieved = ", nrow(V))
say("  rows per database:"); pr(V[, .N, by=database][order(-N)])
saveRDS(V, file.path(CACHE,"validated_all.rds"))

## ---- normalise ----
V[, mir_canon := canon_mirna(mature_mirna_id)]
V[, gene := toupper(target_symbol)]
V[, st_l := tolower(trimws(support_type))]
say("support_type by database:"); pr(V[, .N, by=.(database, support_type)][order(database,-N)])

## refuted only when miRTarBase explicitly says Non-Functional
V[, refuted := grepl("non-functional", st_l)]
say("rows flagged refuted (miRTarBase Non-Functional MTI) = ", sum(V$refuted))
Vp <- V[refuted == FALSE]
say("supporting rows kept = ", nrow(Vp))

STRONG_RE <- paste0("^(.*luciferase.*|.*reporter assay.*|luciferasereporterassay|",
                    "western ?blot.*|westernblotting|immunoblot.*|qrt-?pcr|qrtpcr|q-?pcr|",
                    "rt-?pcr|northern ?blot|elisa|immunohistochemistry.*|immunofluorescence|",
                    "immunocytochemistry|in situ hybridi[sz]ation)$")
WEAK_RE   <- paste0("^(.*clip.*|.*clash.*|microarray.*|biotin-.*|rna-?seq|srna-?seq|rpf-?seq|",
                    "rip-?seq|sequencing|next generation sequencing.*|degradome.*|psilac|",
                    "proteomics|chimeric fragments|ago-ip|immunoprecipitaion.*|",
                    "immunoprecipitation.*|chip.*|trap|flow.*|facs|other|)$")
classify <- function(expstr) {
  tk <- tolower(trimws(unlist(strsplit(expstr, "//|;"))))
  list(strong = any(grepl(STRONG_RE, tk)), weak = any(grepl(WEAK_RE, tk) & nzchar(tk)))
}
utok <- unique(Vp$experiment)
cls  <- rbindlist(lapply(utok, function(x){ z <- classify(x); data.table(experiment=x, s=z$strong, w=z$weak) }))
Vp <- merge(Vp, cls, by="experiment", all.x=TRUE)
Vp[, strong_hit := s | st_l == "functional mti"]
Vp[, weak_hit   := w | grepl("functional mti \\(weak\\)", st_l)]
say("supporting rows with a STRONG assay token = ", sum(Vp$strong_hit),
    " ; with only weak/high-throughput = ", sum(!Vp$strong_hit))
## unclassified assay strings (sanity check)
unc <- Vp[strong_hit==FALSE & weak_hit==FALSE, .N, by=experiment][order(-N)]
say("supporting rows whose assay string matched NEITHER list = ", sum(unc$N),
    " over ", nrow(unc), " distinct strings"); pr(head(unc, 15))

agg <- Vp[, .(databases     = paste(sort(unique(database)), collapse=";"),
              n_pubmed      = uniqueN(unlist(strsplit(paste(pubmed_id, collapse=","), ","))),
              any_strong    = any(strong_hit),
              support_types = paste(sort(unique(trimws(unlist(strsplit(paste(experiment, collapse="//"), "//|;"))))), collapse=";")),
          by=.(source=mir_canon, target=gene)]
agg[, support_types := gsub("^;", "", substr(support_types, 1, 400))]
say("distinct supporting (canonical miRNA, target) pairs in multiMiR = ", nrow(agg))

## ---- join onto the canonical miRNA_target edges ----
mt[, target := toupper(target)]
res <- merge(mt, agg, by=c("source","target"), all.x=TRUE)
res[, validated := !is.na(databases)]
res[, tier := fifelse(validated & any_strong, "strong", fifelse(validated, "weak", "predicted_only"))]
res[is.na(databases), databases := ""]; res[is.na(support_types), support_types := ""]
out <- res[, .(source, target, edge_type, validated, databases, support_types, tier)]
setorder(out, -validated, source, target)
fp <- file.path(BASE,"data","edge_evidence_tier.tsv")
fwrite(out, fp, sep="\t", quote=TRUE)
say("WROTE ", fp, " rows=", nrow(out))
fwrite(res, file.path(BASE,"results","miRNA_target_evidence_extended.tsv"), sep="\t", quote=TRUE)

say("=== HEADLINE NUMBERS ===")
say("miRNA_target edges total = ", nrow(out))
say("experimentally validated = ", sum(out$validated), " (", sprintf("%.1f",100*mean(out$validated)), "%)")
say("--- tier breakdown ---")
pr(out[, .(n=.N, pct=round(100*.N/nrow(out),1)), by=tier][order(-n)])
say("--- database combination among validated edges ---")
pr(out[validated==TRUE, .N, by=databases][order(-N)])
say("--- per-database coverage (exact token match) ---")
for (db in c("mirtarbase","tarbase","mirecords")) {
  hit <- sapply(strsplit(out$databases, ";"), function(x) db %in% x)
  say("   ", db, " : ", sum(hit), " edges") }
say("edges validated ONLY by weak/high-throughput evidence = ", sum(out$tier=="weak"))
say("edges with NO experimental validation (prediction-only)  = ", sum(out$tier=="predicted_only"))
## how many distinct miRNAs / genes have at least one validated edge
say("miRNAs with >=1 validated target edge = ", uniqueN(out[validated==TRUE, source]), " / ", uniqueN(out$source))
say("targets with >=1 validated miRNA edge = ", uniqueN(out[validated==TRUE, target]), " / ", uniqueN(out$target))
fwrite(out[, .N, by=tier], file.path(BASE,"results","miRNA_target_tier_summary.tsv"), sep="\t", quote=FALSE)
say("DONE")
