## v7 / BioNet step 7: one consolidated headline table
suppressMessages(library(data.table))
ROOT <- "/path/to/revision"; OUT <- file.path(ROOT,"results/v7")
P <- function(...) file.path(OUT, paste0("bionet_", ...))
s <- fread(P("fdr_scan.csv")); S <- fread(P("module_summary.csv")); fb <- readRDS(P("bum_fits.rds"))$gene_gw
first <- function(col) { r <- s[get(col)==TRUE][which.min(fdr)]; sprintf("%s (module n=%d, %.0f%% of network)",
                          format(r$fdr, digits=2), r$module_n, 100*r$module_n/587) }
rows <- rbind(
 data.table(quantity="package / version", value="BioNet 1.68.0 (Bioconductor 3.21), RBGL 1.84.0, graph 1.86.0, igraph 2.2.2, R 4.5.1"),
 data.table(quantity="graph", value="587 nodes, 5,636 undirected edges (from 6,859 directed signed edges), 1 connected component, no self-loops"),
 data.table(quantity="p-values", value="genes/TFs 363/364 from results/v2/v2_DE_genes.csv; miRNAs 213/223 from results/v2/v2_DE_mirnas.csv; 11 unmeasured nodes scored at p=1"),
 data.table(quantity="BUM fit (genes, genome-wide n=20,250)",
            value=sprintf("lambda=%.4f, a=%.4f; upper bound on null fraction = %.1f%%", fb$lambda, fb$a, 100*(fb$lambda+(1-fb$lambda)*fb$a))),
 data.table(quantity="modules returned by runFastHeinz", value="1 (the method returns a single maximum-scoring connected subnetwork, not a partition)"),
 data.table(quantity="primary FDR (pre-specified rule: largest FDR giving <=10% of nodes)",
            value=sprintf("%.2e -> %d nodes (24 genes, 19 TFs, 15 miRNAs), module score %.0f", S[run=="primary"]$fdr, S[run=="primary"]$module_n, S[run=="primary"]$module_score)),
 data.table(quantity="sensitivity FDR 1e-3 (BioNet vignette default)", value=sprintf("%d nodes = %.0f%% of the network", S[run=="sens_vign"]$module_n, 100*S[run=="sens_vign"]$frac_of_network)),
 data.table(quantity="sensitivity FDR 1e-25", value=sprintf("%d nodes = %.0f%% of the network", S[run=="sens_mid"]$module_n, 100*S[run=="sens_mid"]$frac_of_network)),
 data.table(quantity="smallest FDR at which COL1A1 is still a member", value=first("has_COL1A1")),
 data.table(quantity="smallest FDR at which COL3A1 is still a member", value=first("has_COL3A1")),
 data.table(quantity="smallest FDR at which hsa-miR-101 is still a member", value=first("has_miR101")),
 data.table(quantity="smallest FDR at which hsa-miR-29a is still a member", value=first("has_miR29a")),
 data.table(quantity="overlap with the 30 prioritised nodes (primary)",
            value=sprintf("13/30 in a 58-node module; expected 2.96; 4.39-fold; hypergeometric p = 7.5e-07 (universe 587)")),
 data.table(quantity="same overlap under a DE-rank-matched null", value="observed 13, null mean 10.97 (sd 1.18), p = 0.083 - not significant once shared dependence on the same DE p-values is controlled"),
 data.table(quantity="module vs DE ranking alone", value="Jaccard 0.785 with the top-58 nodes by node score; the top-58 list alone already overlaps the prioritised 30 by 12"),
 data.table(quantity="collagen / matrisome content of the primary module", value="NABA collagens 0/4; NABA core matrisome 1/17 (FN1 only, expected 2.0); Farmer stroma-related 0/13; HALLMARK EMT 4/47 (expected 5.6)"),
 data.table(quantity="paper's module layers", value="FFL 3-node members 34/204 enriched 1.41x, p=7.7e-4, q=0.027; higher-order-only members 7/104, 0.57x, depletion p=0.038"),
 data.table(quantity="edge-resampling stability (200 x 90% of edges)", value="module size median 54 (IQR 52-57); selection frequency EZH2 1.00, FN1 1.00, NFKB1 0.65, SP1 0.59, RELA 0.02, COL1A1 0.00, COL3A1 0.00, ETS1 0.00, POSTN 0.00, SPARC 0.00, THBS2 0.00, TGFBI 0.00, miR-29a 0.00, miR-101 0.00"),
 data.table(quantity="compartment skew (Wu breast atlas pseudobulk)", value="module median log2(CAF/carcinoma) = 1.49 vs 0.73 for the network background, Wilcoxon p = 0.021; but 25 of the 26 CAF-skewed members are DOWN in tumour (only FN1 is up)")
)
fwrite(rows, P("HEADLINE_SUMMARY.csv")); print(rows)
