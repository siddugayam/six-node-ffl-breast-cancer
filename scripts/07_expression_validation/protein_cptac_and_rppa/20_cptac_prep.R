#!/usr/bin/env Rscript
# 20_cptac_prep.R -- harmonise CPTAC-BRCA proteome / phosphoproteome / RNA / miRNA
# Output: data-side RDS bundle in results/multiomics/cptac_bundle.rds (kept inside multiomics/)
suppressMessages({library(data.table); library(org.Hs.eg.db); library(AnnotationDbi)})

REV  <- "/path/to/revision"
OUT  <- file.path(REV, "results", "multiomics")
RAW  <- "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(REV, "logs", "cptac_prep.log")
con <- file(LOG, open = "wt")
say <- function(...) { m <- paste0(format(Sys.time(),"%H:%M:%S")," | ",paste0(...,collapse="")); cat(m,"\n"); cat(m,"\n",file=con); flush(con) }
say("=== 20_cptac_prep.R START ===")

f_prot   <- file.path(RAW,"Proteome_BCM_GENCODE_v34_harmonized_v1/Proteome_BCM_GENCODE_v34_harmonized_v1/BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt")
f_phos   <- file.path(RAW,"Phosphoproteome_BCM_GENCODE_v34_harmonized_v1/Phosphoproteome_BCM_GENCODE_v34_harmonized_v1/BRCA_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt")
f_rna    <- file.path(RAW,"RNA_BCM_v1/RNA_BCM_v1/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt")
f_mir    <- file.path(RAW,"miRNA_BCM_v1/BRCA_miRNAseq_mature_miRNA_RPM_log2_Tumor.txt")
f_meta   <- file.path(RAW,"Clinical_meta_data_v1/Clinical_meta_data_v1/BRCA_meta.txt")
for (f in c(f_prot,f_phos,f_rna,f_mir,f_meta)) if(!file.exists(f)) stop("MISSING: ", f)

rd <- function(f) { d <- fread(f, sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
                    rn <- as.character(d[[1]]); d <- d[,-1,drop=FALSE]
                    M <- as.matrix(d); rownames(M) <- rn; storage.mode(M) <- "double"; M }
P  <- rd(f_prot);  say("proteome  ", nrow(P), " x ", ncol(P))
PH <- rd(f_phos);  say("phospho   ", nrow(PH)," x ", ncol(PH))
R  <- rd(f_rna);   say("RNA       ", nrow(R), " x ", ncol(R))
MI <- rd(f_mir);   say("miRNA     ", nrow(MI)," x ", ncol(MI))

## ---------------- sample harmonisation ----------------
say("--- sample ID harmonisation (raw IDs are CPTAC case IDs like 11BR047) ---")
say("prot ids ex: ", paste(head(colnames(P),4), collapse=","))
say("rna  ids ex: ", paste(head(colnames(R),4), collapse=","))
say("mir  ids ex: ", paste(head(colnames(MI),4),collapse=","))
say("phos ids ex: ", paste(head(colnames(PH),4),collapse=","))
say("non-alnum in any id? ", any(grepl("[^A-Za-z0-9]", c(colnames(P),colnames(R),colnames(MI),colnames(PH)))))
s_all  <- sort(Reduce(intersect, list(colnames(R), colnames(P), colnames(MI))))
s_rp   <- sort(intersect(colnames(R), colnames(P)))
s_pph  <- sort(intersect(colnames(P), colnames(PH)))
s_allp <- sort(Reduce(intersect, list(colnames(R), colnames(P), colnames(MI), colnames(PH))))
say("n RNA=",ncol(R)," protein=",ncol(P)," phospho=",ncol(PH)," miRNA=",ncol(MI))
say("RNA+protein             n = ", length(s_rp))
say("protein+phospho         n = ", length(s_pph))
say("RNA+protein+miRNA       n = ", length(s_all))
say("RNA+protein+miRNA+phos  n = ", length(s_allp))
say("in protein but not RNA: ", paste(setdiff(colnames(P), colnames(R)), collapse=","))
say("in RNA but not protein: ", paste(setdiff(colnames(R), colnames(P)), collapse=","))
say("in protein but not miRNA: ", paste(setdiff(colnames(P), colnames(MI)), collapse=","))

## ---------------- ENSG -> SYMBOL ----------------
ens_base <- function(x) sub("\\..*$","",x)
netnodes <- read.delim(file.path(REV,"data","canonical_nodes.tsv"), stringsAsFactors=FALSE)
gsym <- netnodes$name[netnodes$type %in% c("TF","Gene")]
# alias fallback for symbols org.Hs.eg.db does not carry under SYMBOL
alias_fix <- c(CTGF="CCN2", CYR61="CCN1", MKL1="MRTFA")
ask <- unique(c(gsym, unname(alias_fix)))
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=ask, keytype="SYMBOL", columns="ENSEMBL"))
map <- map[!is.na(map$ENSEMBL),]
# map back the aliases to the network name
back <- setNames(names(alias_fix), unname(alias_fix))
map$NETNAME <- ifelse(map$SYMBOL %in% names(back), back[map$SYMBOL], map$SYMBOL)
say("network gene/TF nodes: ", length(gsym), " ; with >=1 ENSEMBL id: ", length(unique(map$NETNAME)))
say("still unmapped: ", paste(setdiff(gsym, map$NETNAME), collapse=","))

collapse_to_symbol <- function(M, map, label) {
  eb <- ens_base(rownames(M))
  sub <- map[map$ENSEMBL %in% eb, ]
  res <- list(); nmulti <- 0
  for (g in unique(sub$NETNAME)) {
    ids <- unique(sub$ENSEMBL[sub$NETNAME==g]); rr <- which(eb %in% ids)
    if (!length(rr)) next
    if (length(rr)>1) { nmulti <- nmulti+1; v <- colMeans(M[rr,,drop=FALSE], na.rm=TRUE) } else v <- M[rr,]
    v[is.nan(v)] <- NA_real_
    res[[g]] <- v
  }
  X <- do.call(rbind, res); rownames(X) <- names(res)
  say(label, ": ", nrow(X), " network genes recovered (", nmulti, " collapsed from >1 ENSG row)")
  X
}
Pg <- collapse_to_symbol(P, map, "PROTEIN")
Rg <- collapse_to_symbol(R, map, "mRNA")

## ---------------- miRNA canonicalisation ----------------
source(file.path(REV,"scripts","01_network_assembly","00_mirna_canon.R"))
cm <- canon_mirna(rownames(MI))
say("miRNA rows ", nrow(MI), " -> ", length(unique(cm)), " canonical stems")
idx <- split(seq_len(nrow(MI)), cm)
MIc <- do.call(rbind, lapply(idx, function(i) if(length(i)==1) MI[i,] else colMeans(MI[i,,drop=FALSE], na.rm=TRUE)))
rownames(MIc) <- names(idx)
mirnodes <- netnodes$name[netnodes$type=="miRNA"]
say("network miRNA nodes ", length(mirnodes), " ; present in CPTAC canonical matrix ", length(intersect(mirnodes, rownames(MIc))))
say("network miRNAs ABSENT from CPTAC: ", paste(setdiff(mirnodes, rownames(MIc)), collapse=","))

meta <- fread(f_meta, sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
say("meta rows ", nrow(meta), " cols ", ncol(meta))

saveRDS(list(P=P, PH=PH, R=R, MI=MI, Pg=Pg, Rg=Rg, MIc=MIc, meta=meta,
             s_all=s_all, s_rp=s_rp, s_pph=s_pph, s_allp=s_allp, map=map),
        file.path(OUT,"cptac_bundle.rds"))
say("WROTE ", file.path(OUT,"cptac_bundle.rds"))

# small human-readable sample-overlap table
ov <- data.frame(assay=c("RNAseq","Proteome","Phosphoproteome","miRNAseq",
                         "RNA+protein","protein+phospho","RNA+protein+miRNA","RNA+protein+miRNA+phospho"),
                 n_samples=c(ncol(R),ncol(P),ncol(PH),ncol(MI),
                             length(s_rp),length(s_pph),length(s_all),length(s_allp)))
write.csv(ov, file.path(OUT,"cptac_sample_overlap.csv"), row.names=FALSE)
say("WROTE cptac_sample_overlap.csv rows=", nrow(ov))
say("=== DONE ===" ); close(con)
