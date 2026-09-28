## ============================================================================
## TASK C: CPTAC PAN-CANCER PROTEOME
## Replicate the miR-29 -> collagen PROTEIN-level result in every CPTAC tumour
## type that has both a mature-miRNA matrix and a harmonised proteome.
## Also compare mRNA-level vs protein-level correlation in the same samples.
## Everything computed in this run.
## ============================================================================
setwd("/path/to/revision")
OUT <- "results/v2"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

ROOT <- "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data"
PDIR <- file.path(ROOT, "Proteome_BCM_GENCODE_v34_harmonized_v1",
                        "Proteome_BCM_GENCODE_v34_harmonized_v1")
MDIR <- file.path(ROOT, "miRNA_BCM_v1")
RDIR <- file.path(ROOT, "RNA_BCM_v1")

cat("================ PART C : CPTAC PAN-CANCER ================\n")

## ------------------------------------------------- inventory the cohorts ---
mfiles <- list.files(MDIR, pattern="_miRNAseq_mature_miRNA_RPM_log2_Tumor.txt$", full.names=TRUE)
pfiles <- list.files(PDIR, pattern="_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt$", full.names=TRUE)
## matched bulk RNA: use the gene-level RSEM coding matrices only
rfiles <- list.files(RDIR, recursive=TRUE,
                     pattern="_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor\\.txt$",
                     full.names=TRUE)
rfiles <- rfiles[!grepl("MACOSX", rfiles)]
mcoh <- sub("_miRNAseq.*","", basename(mfiles))
pcoh <- sub("_proteomics.*","", basename(pfiles))
cat("miRNA cohorts   :", paste(sort(mcoh), collapse=", "), "\n")
cat("proteome cohorts:", paste(sort(pcoh), collapse=", "), "\n")
cat("RNA files found :", length(rfiles), "\n")
if(length(rfiles)) cat("  ", paste(basename(rfiles), collapse="\n   "), "\n")

## the two naming mismatches: ccRCC/CCRCC and HGSC/OV are the same cohorts
alias <- c(ccRCC="CCRCC", HGSC="OV")
pmatch_name <- function(x) if (x %in% names(alias)) alias[[x]] else x

## ------------------------------------------------------- gene annotation ---
ann <- read.delim(file.path(PDIR,"README","Gene_annotation_and_representable_isoform_mapping_table.txt"),
                  stringsAsFactors=FALSE, check.names=FALSE)
cat("\nannotation rows:", nrow(ann), "\n")
ens2sym <- setNames(ann$gene_name, ann$gene)          # versioned ENSG -> symbol
cat("ENSG(versioned) -> symbol entries:", length(ens2sym), "\n")

## -------------------------------------------------------------- helpers ----
read_mat <- function(f){
  d <- read.delim(f, stringsAsFactors=FALSE, check.names=FALSE, row.names=1)
  as.matrix(d)
}
to_symbol <- function(X){
  sym <- ens2sym[rownames(X)]
  ok  <- !is.na(sym) & sym != ""
  X <- X[ok,,drop=FALSE]; sym <- sym[ok]
  ## collapse duplicated symbols by the row with the most non-missing values
  if (any(duplicated(sym))){
    nobs <- rowSums(is.finite(X))
    ord  <- order(sym, -nobs)
    X <- X[ord,,drop=FALSE]; sym <- sym[ord]
    keep <- !duplicated(sym)
    X <- X[keep,,drop=FALSE]; sym <- sym[keep]
  }
  rownames(X) <- sym
  stopifnot(!any(duplicated(rownames(X))))
  X
}
sp <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 15) return(c(rho=NA, p=NA, n=sum(ok)))
  if (sd(x[ok])==0 || sd(y[ok])==0) return(c(rho=NA, p=NA, n=sum(ok)))
  tt <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
  c(rho=unname(tt$estimate), p=tt$p.value, n=sum(ok))
}
## Steiger test for two DEPENDENT correlations sharing variable x
## (r_xy1 vs r_xy2 with y1,y2 the mRNA and protein measurement of the same gene)
steiger <- function(r12, r13, r23, n){
  if (any(!is.finite(c(r12,r13,r23))) || n < 10) return(c(t=NA,p=NA))
  rm2 <- (r12^2 + r13^2)/2
  f   <- (1 - r23)/(2*(1 - rm2)); if (f > 1) f <- 1
  h   <- (1 - f*rm2)/(1 - rm2)
  z1 <- atanh(r12); z2 <- atanh(r13)
  z  <- (z1 - z2) * sqrt((n-3)/(2*(1-r23)*h))
  c(t=z, p=2*pnorm(-abs(z)))
}

ARMS <- list("hsa-miR-29a"=c("hsa-miR-29a-3p","hsa-miR-29a-5p"),
             "hsa-miR-29b"=c("hsa-miR-29b-1-5p","hsa-miR-29b-2-5p","hsa-miR-29b-3p"),
             "hsa-miR-29c"=c("hsa-miR-29c-3p","hsa-miR-29c-5p"))
COLS <- c("COL1A1","COL3A1","COL1A2","COL5A1","COL5A2","COL6A3","COL11A1","FN1","SPARC","LOX")
KEY  <- c("COL1A1","COL3A1")

## ------------------------------------------------------------- main loop ---
res <- list(); inv <- list()
for (ct in sort(mcoh)){
  pn <- pmatch_name(ct)
  pf <- pfiles[pcoh == pn]
  if (!length(pf)){ cat("\n[", ct, "] NO PROTEOME -> skipped\n");
    inv[[length(inv)+1]] <- data.frame(cohort=ct, proteome_name=pn, has_proteome=FALSE,
      n_mirna=NA, n_proteome=NA, n_overlap=NA, has_rna=FALSE, n_rna_overlap=NA); next }

  Mi <- read_mat(mfiles[mcoh==ct])
  Pr <- to_symbol(read_mat(pf[1]))
  ov <- intersect(colnames(Mi), colnames(Pr))
  cat(sprintf("\n[%s -> %s] miRNA %d x %d | proteome %d x %d | OVERLAP n=%d\n",
              ct, pn, nrow(Mi), ncol(Mi), nrow(Pr), ncol(Pr), length(ov)))

  ## optional matched RNA for the same samples
  rf <- rfiles[basename(rfiles) == paste0(ct, "_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt")]
  if (!length(rf))
    rf <- rfiles[basename(rfiles) == paste0(pn, "_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt")]
  Rn <- NULL; nrov <- NA
  if (length(rf)){
    Rn <- try(to_symbol(read_mat(rf[1])), silent=TRUE)
    if (inherits(Rn,"try-error")) Rn <- NULL else nrov <- length(intersect(ov, colnames(Rn)))
  }
  inv[[length(inv)+1]] <- data.frame(cohort=ct, proteome_name=pn, has_proteome=TRUE,
    n_mirna=ncol(Mi), n_proteome=ncol(Pr), n_overlap=length(ov),
    has_rna=!is.null(Rn), n_rna_overlap=nrov)

  if (length(ov) < 20){ cat("   overlap < 20 -> no statistics\n"); next }

  Mi <- Mi[, ov, drop=FALSE]; Pr <- Pr[, ov, drop=FALSE]
  Rs <- if (!is.null(Rn) && length(intersect(ov, colnames(Rn)))>=20) Rn[, ov[ov %in% colnames(Rn)], drop=FALSE] else NULL

  ## build the two miRNA mappings
  for (mp in c("arm_averaged","dominant_arm")){
    for (lab in names(ARMS)){
      a <- intersect(ARMS[[lab]], rownames(Mi))
      if (!length(a)) next
      mv <- if (mp=="arm_averaged") colMeans(Mi[a,,drop=FALSE], na.rm=TRUE)
            else Mi[a[which.max(rowMeans(Mi[a,,drop=FALSE], na.rm=TRUE))], ]
      for (g in COLS){
        if (!(g %in% rownames(Pr))) next
        pv <- Pr[g, ]
        e_prot <- sp(mv, pv)
        ## same samples, mRNA layer
        e_rna <- c(rho=NA,p=NA,n=NA); st <- c(t=NA,p=NA)
        if (!is.null(Rs) && g %in% rownames(Rs)){
          s2 <- colnames(Rs)
          e_rna <- sp(mv[s2], Rs[g, ])
          rmp <- sp(Rs[g,], Pr[g, s2])["rho"]     # mRNA-protein concordance of the gene
          st <- steiger(unname(e_rna["rho"]), unname(sp(mv[s2], Pr[g,s2])["rho"]),
                        unname(rmp), min(e_rna["n"], length(s2)))
        }
        res[[length(res)+1]] <- data.frame(cohort=ct, proteome_cohort=pn,
          mir=lab, gene=g, mir_mapping=mp, arms_used=paste(a, collapse="+"),
          n_protein=unname(e_prot["n"]), rho_protein=unname(e_prot["rho"]), p_protein=unname(e_prot["p"]),
          n_mrna=unname(e_rna["n"]), rho_mrna=unname(e_rna["rho"]), p_mrna=unname(e_rna["p"]),
          steiger_z=unname(st["t"]), steiger_p=unname(st["p"]),
          stringsAsFactors=FALSE)
      }
    }
  }
  d <- do.call(rbind, res)
  d <- d[d$cohort==ct & d$mir_mapping=="arm_averaged" & d$gene %in% KEY,]
  for (i in seq_len(nrow(d)))
    cat(sprintf("   %-12s -> %-7s protein rho=%+.3f (p=%.2e, n=%d)\n",
                d$mir[i], d$gene[i], d$rho_protein[i], d$p_protein[i], d$n_protein[i]))
}
INV <- do.call(rbind, inv); rownames(INV) <- NULL
RES <- do.call(rbind, res); rownames(RES) <- NULL

cat("\n================ COHORT INVENTORY ================\n"); print(INV)

## ------------------------------------------------------------- summaries ---
cat("\n================ miR-29 -> COL1A1 / COL3A1 AT PROTEIN LEVEL ================\n")
key <- RES[RES$gene %in% KEY & RES$mir_mapping=="arm_averaged" & is.finite(RES$rho_protein),]
for (g in KEY) for (m in names(ARMS)){
  d <- key[key$gene==g & key$mir==m,]
  if(!nrow(d)) next
  d <- d[order(d$rho_protein),]
  cat(sprintf("\n-- %s -> %s (arm-averaged) --\n", m, g))
  for (i in seq_len(nrow(d)))
    cat(sprintf("   %-8s n=%3d rho=%+.3f p=%.2e %s\n", d$cohort[i], d$n_protein[i],
                d$rho_protein[i], d$p_protein[i],
                ifelse(d$p_protein[i]<0.05, ifelse(d$rho_protein[i]<0,"SIG-neg","SIG-POS"), "")))
  k <- nrow(d); neg <- sum(d$rho_protein<0)
  bt <- binom.test(neg, k, 0.5, alternative="greater")
  z <- atanh(d$rho_protein); tt <- t.test(z)
  cat(sprintf("   => %d/%d cohorts negative; sign-test p=%.4g; mean Fisher-z=%.3f (t=%.2f, p=%.4g); median rho=%+.3f\n",
              neg, k, bt$p.value, mean(z), tt$statistic, tt$p.value, median(d$rho_protein)))
}

cat("\n================ PROTEIN vs mRNA IN THE SAME SAMPLES ================\n")
cmp <- RES[RES$mir_mapping=="arm_averaged" & is.finite(RES$rho_protein) & is.finite(RES$rho_mrna),]
if (nrow(cmp)){
  cmpk <- cmp[cmp$gene %in% KEY,]
  cat("testable pairs (all ECM genes):", nrow(cmp), "  key collagens:", nrow(cmpk), "\n")
  st <- sum(abs(cmp$rho_protein) > abs(cmp$rho_mrna))
  cat(sprintf("|rho| larger at PROTEIN in %d/%d pairs, binomial p=%.4g\n",
              st, nrow(cmp), binom.test(st, nrow(cmp), 0.5)$p.value))
  stk <- sum(abs(cmpk$rho_protein) > abs(cmpk$rho_mrna))
  cat(sprintf("   restricted to COL1A1/COL3A1: %d/%d, binomial p=%.4g\n",
              stk, nrow(cmpk), binom.test(stk, nrow(cmpk), 0.5)$p.value))
  cat("\nper-cohort (COL1A1/COL3A1, arm-averaged):\n")
  print(format(cmpk[,c("cohort","mir","gene","n_mrna","rho_mrna","n_protein","rho_protein","steiger_p")],
               digits=3), row.names=FALSE)
} else {
  cat("NO matched RNA matrices were found under", RDIR, "\n")
  cat("-> the mRNA-vs-protein comparison could not be run pan-cancer;\n")
  cat("   protein-level correlations are still reported above.\n")
}

write.csv(RES, file.path(OUT,"cptac_pancancer_mir29.csv"), row.names=FALSE)
write.csv(INV, file.path(OUT,"cptac_pancancer_inventory.csv"), row.names=FALSE)
cat("\nWROTE", file.path(OUT,"cptac_pancancer_mir29.csv"), nrow(RES), "rows\n")
cat("WROTE", file.path(OUT,"cptac_pancancer_inventory.csv"), nrow(INV), "rows\n")
cat("DONE PART C\n")
