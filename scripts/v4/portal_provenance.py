"""Provenance record: exactly what was retrieved from each portal, how, and every caveat."""
import csv, os, time
RES = '/path/to/revision/results/v4/'
NOW = time.strftime('%Y-%m-%d', time.gmtime())
F = ["portal", "version_or_module", "base_url", "how_retrieved", "what_was_retrieved",
     "cohorts_and_n", "output_file", "caveats"]
R = [
 dict(portal="KM Plotter", version_or_module="breast cancer miRNA module (mirPower)",
   base_url="https://kmplot.com/analysis/index.php?p=service&cancer=breast_mirna",
   how_retrieved="multipart POST replay of the exact field set the site's own form submits "
     "(step=2). The field set was captured by driving the real form once in headless Chromium "
     "(scripts/v4/portal_km_probe.py) and verified to reproduce the browser result digit for "
     "digit (hsa-miR-130a/METABRIC/auto-cutoff: HR 0.81, 95% CI 0.65-0.99, p=0.0431 both ways). "
     "Per-arm n and event counts were read from the portal's own 'export plot data as text' file.",
   what_was_retrieved="136 analyses: 8 miRNAs x 4 datasets x {auto best cut-off, median split}, "
     "plus 13 subtype/nodal restrictions for hsa-miR-130a and hsa-miR-29a",
   cohorts_and_n="METABRIC n=1262 (400 OS events); TCGA n=1062 (148); GSE40267 n=85 (58); "
     "GSE19783 n=93 (24)",
   output_file="portals_kmplot_mirna_breast.csv",
   caveats="This module offers OVERALL SURVIVAL ONLY - there is no RFS or DMFS option for miRNA "
     "(the surv dropdown has a single entry, 'OS'). Its TCGA dataset is the same TCGA-BRCA cohort "
     "we analyse ourselves, so only METABRIC, GSE40267 and GSE19783 are independent of us. "
     "The 'auto select best cut-off' p-values are optimised over all cut-points; KM Plotter's own "
     "FDR column is reported alongside and for hsa-miR-130a in METABRIC it is 'over 50%'. "
     "Restrictions that returned no result did so because the portal reported an insufficient "
     "sample size - the exact message is stored in the 'error' column."),
 dict(portal="KM Plotter", version_or_module="breast cancer mRNA (gene chip) module",
   base_url="https://kmplot.com/analysis/index.php?p=service&cancer=breast",
   how_retrieved="same multipart POST replay; field set captured from a real browser submission "
     "(scripts/v4/portal_km_probe_mrna.py). probe_set_option=jetset (JetSet best probe set), "
     "array_quality=2 (exclude biased arrays), redundant_publication and propcheck on - i.e. the "
     "site's own recommended defaults.",
   what_was_retrieved="144 analyses: 16 genes x {RFS, OS, DMFS} x {median split, auto cut-off}, "
     "plus 7 ER/PAM50 restrictions on RFS for 8 genes",
   cohorts_and_n="pooled public microarray cohorts; RFS up to n=4929 (1513 events), "
     "OS n=1879 (438), DMFS n=2765 (649). Genes whose JetSet probe is HG-U133 Plus 2.0-only "
     "(PTEN, SMAD4, XIAP, SP1, ETS1) fall to n=2032 RFS / n=943 OS.",
   output_file="portals_kmplot_mrna_breast.csv",
   caveats="The probe actually used is recorded per row. RELA resolves to probe 201783_s_at, "
     "which the portal labels with the legacy alias NFKB3. Sample size differs between genes "
     "because it depends on which arrays carry the JetSet probe - do not compare HRs across "
     "genes without noting n."),
 dict(portal="UALCAN", version_or_module="TCGA-BRCA expression / promoter methylation; "
     "CPTAC breast proteome",
   base_url="https://ualcan.path.uab.edu/cgi-bin/",
   how_retrieved="direct GET of TCGAExResultNew2.pl, TCGA-methyl-Result.pl and CPTAC-Result.pl. "
     "UALCAN embeds the complete box-plot statistics (whisker low, Q1, median, Q3, whisker high) "
     "and the group sizes in the Highcharts configuration inside the returned HTML, so every "
     "number here is the portal's own value, parsed, not read off a picture.",
   what_was_retrieved="16 genes x 11 groupings (sample type, stage, race, gender, age, major "
     "subclasses, subclasses with TNBC types, menopause status, histology, nodal metastasis "
     "status, TP53 status) for expression; 10 groupings for promoter methylation; 17 for CPTAC "
     "protein. All pairwise t-test p-values captured separately.",
   cohorts_and_n="expression 114 normal vs 1097 primary tumour; methylation 97 vs 793; "
     "CPTAC breast 18 normal vs 125 primary tumour",
   output_file="portals_ualcan_expression.csv, portals_ualcan_methylation.csv, "
     "portals_ualcan_cptac_protein.csv, portals_ualcan_stats.csv",
   caveats="(1) UALCAN's CPTAC values are Z-scores centred on the whole cancer-type cohort, which "
     "is 125 tumours and 18 normals - so the TUMOUR median is ~0 by construction for every gene "
     "and the comparison is really 'where do the 18 normals sit', not a fold change. Do not quote "
     "a CPTAC 'fold change' from UALCAN. (2) TGFBR2 and CCND2 are not quantified in the CPTAC "
     "breast proteome; those queries return no chart, which is a genuine absence, not a fetch "
     "failure. (3) UALCAN's survival module renders its Kaplan-Meier curves as cairo SVGs in "
     "which the text is glyph paths, so the survival p-value is not machine-readable and was NOT "
     "extracted - survival evidence here comes from KM Plotter, bc-GenExMiner and GEPIA2."),
 dict(portal="bc-GenExMiner", version_or_module="v5.3 (site build dated 7 September 2026); "
     "exhaustive prognostic (mode=2), targeted differential expression (mode=8), "
     "targeted correlation (mode=4)",
   base_url="https://bcgenex.ico.unicancer.fr/BC-GEM/GEM-requete.php",
   how_retrieved="two-step POST (query form, then the site's 'Start analysis' validation step), "
     "replaying the body captured from a real browser submission "
     "(scripts/v4/portal_bcgem_probe.py). Results are read from the portal's own HTML tables; "
     "the correlation matrices are the portal's own downloadable CSV exports.",
   what_was_retrieved="prognosis: 16 genes x 27 ER/PR/HER2 population strata x {DMFS, OS, DFS} "
     "univariate Cox (p, HR, 95% CI, n, events, direction of good prognosis) = 81 rows per gene. "
     "differential expression: healthy vs tumour-adjacent vs tumour on the GTEx&TCGA RNA-seq "
     "source, with Welch p and Dunnett-Tukey-Kramer pairwise p. correlation: pairwise Pearson r "
     "over all 16 genes, all patients and stratified by ER / PAM50 / TNBC status.",
   cohorts_and_n="prognosis, all-ER/all-PR/all-HER2 stratum: DFS n up to 9080, DMFS up to 6495, "
     "OS up to 5451. correlation: n up to 10,275 patients per pair. "
     "differential expression source = GTEx & TCGA RNA-seq.",
   output_file="portals_bcgenexminer_prognosis.csv, portals_bcgenexminer_diffexpr.csv, "
     "portals_bcgenexminer_correlation.csv",
   caveats="(1) TLS: this host serves an INCOMPLETE certificate chain (openssl reports 'Verify "
     "return code: 21, unable to verify the first certificate'). The leaf certificate is a valid "
     "Sectigo OV certificate covering *.ico.unicancer.fr, so the site is authentic; it simply "
     "does not send the intermediate. curl/requests therefore need verification disabled for this "
     "host, and that is what the scripts do - deliberately, for this host only. (2) The Kaplan-"
     "Meier figures are cairo SVGs with glyph-path text, but the HR/p/CI/n are in HTML tables, "
     "which is what was harvested. (3) The differential-expression module's 'nature of the "
     "tissue' split is only available on the GTEx&TCGA RNA-seq source (Gene2=G2RT); on the "
     "microarray sources it returns 'less than 2 subtypes with more than 11 patients'. "
     "(4) p-values are printed as '< 0.0001' or '> 0.10', not as numbers."),
 dict(portal="GEPIA2", version_or_module="Survival Analysis and pan-cancer expression bar plot",
   base_url="http://gepia2.cancer-pku.cn/assets/PHP4/survival_zf.php , "
     "http://gepia2.cancer-pku.cn/assets/PHP2/barplot.php",
   how_retrieved="direct POST to the two endpoints GEPIA2's own front end calls (found by "
     "capturing the SPA's network traffic). Survival returns a PDF whose text layer carries "
     "Logrank p, HR(high), p(HR), n(high), n(low); those were extracted with pypdf, not read off "
     "an image. Expression medians come from the bar-plot response's own data arrays.",
   what_was_retrieved="16 genes x {OS, DFS} x {median split, 75/25 quartile split} on TCGA-BRCA; "
     "tumour and normal median TPM for all 31 GEPIA2 cancer types",
   cohorts_and_n="TCGA-BRCA, n(high)=n(low)=535 at median split (534 for XIAP/EZH2), "
     "268/268 at the quartile split",
   output_file="portals_gepia2_survival.csv, portals_gepia2_expression.csv",
   caveats="(1) GEPIA2's 'normal' is TCGA normals PLUS GTEx healthy breast, which is largely "
     "adipose and stroma. That makes GEPIA2's tumour-vs-normal direction differ from UALCAN's "
     "and ours for stromal/secreted genes; see the FLAG column of portals_CROSSCHECK_genes.csv. "
     "(2) GEPIA2 survival is TCGA-BRCA only (~1070 patients, few events), so it is much less "
     "powered than KM Plotter's pooled microarray cohorts and returns null for almost everything "
     "- absence of significance here is a power statement, not a contradiction. (3) The plot does "
     "not print a 95% CI, so only HR(high) and the two p-values are recorded."),
 dict(portal="TIMER2.0", version_or_module="published deconvolution matrix "
     "(infiltration_estimation_for_tcga.csv)",
   base_url="http://timer.cistrome.org/infiltration_estimation_for_tcga.csv.gz",
   how_retrieved="downloaded TIMER2.0's own precomputed estimate file (8.7 MB gz, 11,070 TCGA "
     "samples x 119 estimates from TIMER, CIBERSORT, CIBERSORT-ABS, QUANTISEQ, MCPCOUNTER, XCELL "
     "and EPIC) and computed Spearman correlations against our own TCGA-BRCA miRNA and gene "
     "matrices in R.",
   what_was_retrieved="529 feature x cell-type Spearman correlations, focused on the three "
     "cancer-associated-fibroblast estimators plus endothelial, uncharacterised and the six "
     "TIMER immune fractions",
   cohorts_and_n="1065 TCGA-BRCA primary tumours with paired miRNA; 1093 with gene expression",
   output_file="portals_timer2_infiltration_correlations.csv",
   caveats="The interactive TIMER site now redirects to TIMER3 (https://compbio.cn/timer3/), a "
     "Shiny application; its Gene_Correlation / Immune_Gene modules could not be driven reliably "
     "headlessly, so the numbers here are computed by us FROM TIMER2.0's published estimates "
     "rather than read off TIMER's own plot. That is stated explicitly because it is a different "
     "provenance class from the other rows. TIMER2.0's web module also offers a tumour-purity "
     "partial correlation, which is NOT replicated here - these are raw Spearman correlations."),
]
with open(RES + 'portals_PROVENANCE.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=F); w.writeheader(); w.writerows(R)
print("wrote portals_PROVENANCE.csv", len(R), "rows")

# inventory of every portals_* file actually produced
inv = []
for f in sorted(os.listdir(RES)):
    if f.startswith('portals_') and f.endswith('.csv'):
        n = sum(1 for _ in open(RES + f)) - 1
        inv.append(dict(file=f, rows=n, bytes=os.path.getsize(RES + f), retrieved_on=NOW))
with open(RES + 'portals_OUTPUT_INVENTORY.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=["file", "rows", "bytes", "retrieved_on"])
    w.writeheader(); w.writerows(inv)
for i in inv:
    print(f"{i['file']:52s} {i['rows']:7d} rows")
