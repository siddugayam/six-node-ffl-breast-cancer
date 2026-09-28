# -*- coding: utf-8 -*-
"""Curated literature notes and one-line verdicts for the 30 prioritised nodes.
Every [PMID:x] marker is verified against a live NCBI esummary record by 43_verify_pmids.py
before the compendium is written; an unverifiable marker aborts the build.
Claims are written to distinguish established consensus from contested findings."""

LIT = {

# ------------------------------------------------------------------ TFs
"E2F1": (
 "E2F1 is the transcriptional effector of the RB pathway and the switch that commits cells to "
 "S phase; its de-repression following RB1 loss or cyclin D-CDK4/6 hyperactivation is one of the "
 "best-established routes to unscheduled proliferation in epithelial cancers [PMID:31053804]. "
 "In breast cancer specifically, co-overexpression of E2F1 and cyclin A is associated with "
 "shorter disease-free survival in node-negative disease [PMID:17201163]. What is established is "
 "that E2F1 activity is elevated in proliferative breast tumours; what is not settled is whether "
 "E2F1 is oncogenic or tumour-suppressive in any given context, since it also drives apoptotic "
 "targets — and in our data it shows no prognostic effect at the mRNA level (OS HR 1.06, "
 "p = 0.44) despite a large fold change."),

"EZH2": (
 "EZH2 is the catalytic subunit of PRC2 and deposits H3K27me3. Its over-expression in breast "
 "cancer, and its association with aggressive, poorly differentiated and metastatic disease, "
 "was established two decades ago and has been repeatedly reproduced [PMID:14500907; "
 "PMID:33491744]. The specific regulatory relationship this study recovers — repression of EZH2 "
 "by miR-101, with genomic loss of the miR-101 loci driving EZH2 over-expression — was described "
 "in prostate and other cancers [PMID:19008416] and has since been reported in breast models "
 "[PMID:28034643]. EZH2 is pharmacologically tractable: tazemetostat is an approved EZH2 "
 "inhibitor, although its approvals are in follicular lymphoma and INI1-negative epithelioid "
 "sarcoma, not breast cancer [PMID:33035457]."),

"GATA3": (
 "GATA3 specifies and maintains the luminal mammary epithelial fate; loss of GATA3 in mouse "
 "models causes loss of differentiation and promotes dissemination [PMID:17129787; PMID:18242514]. "
 "It is one of the three genes most frequently mutated in breast cancer in TCGA, alongside PIK3CA "
 "and TP53 [PMID:23000897], which our data reproduce exactly (non-silent mutation frequency "
 "12.3 %, by far the highest of the 30 nodes). GATA3 is therefore an established lineage factor "
 "whose importance is not in dispute; what its mutations do functionally — loss of function, "
 "gain, or neomorphic activity depending on mutation class — remains actively debated."),

"BRCA1": (
 "BRCA1 is the canonical hereditary breast cancer gene and functions in homologous-recombination "
 "repair and replication-fork protection [PMID:22193408]. Its clinical importance is "
 "therapeutically actionable rather than merely prognostic: olaparib improves progression-free "
 "survival over chemotherapy in germline BRCA-mutated metastatic breast cancer [PMID:28578601]. "
 "Note that BRCA1 enters this ranking on evidence weight, not on network topology (2 FFL cores), "
 "and its somatic non-silent mutation rate in TCGA-BRCA is only 2.8 % — the hereditary "
 "contribution is germline and therefore invisible to the somatic analysis reported here."),

"JUN": (
 "JUN is a component of the AP-1 dimeric transcription factor, whose role in tumorigenesis is "
 "established but explicitly context-dependent — AP-1 can promote or suppress transformation "
 "depending on dimer composition and cell state [PMID:14668816]. In breast cancer the most "
 "clinically salient recent finding is that CDK4/6 inhibition reprograms the enhancer landscape "
 "by stimulating AP-1 transcriptional activity, making AP-1 a resistance node rather than a "
 "simple oncogene [PMID:33997789]. The strong down-regulation we observe (log2FC -1.58) is "
 "therefore not straightforwardly interpretable as loss of an oncogenic driver."),

"EGR2": (
 "EGR2 is genuinely under-characterised in breast cancer, and this is a strength of an unbiased "
 "ranking rather than a weakness: PubMed indexes 1,597 EGR2 records but only 8 that are also "
 "indexed under Breast Neoplasms (0.5 %; counts retrieved 2026-09-09). Most of what is "
 "established comes from other systems — EGR2 is the Krox20 myelination factor and a negative "
 "regulator of T-cell activation [PMID:15834410]. The small cancer literature reports "
 "tumour-suppressive functions: EGR2 induces apoptosis by direct transactivation of BNIP3L and "
 "BAK [PMID:12687019] and acts within the PTEN signalling pathway [PMID:11494141]; in breast "
 "cells an EGR2/CITED1 complex regulates ERBB2 [PMID:17938205], and EGR2 is a downstream target "
 "of SFRP1 signalling [PMID:28687085]. Its strong down-regulation here (log2FC -2.43) is "
 "consistent with, but does not establish, that role. Its perfect validated-edge fraction (1.00) "
 "rests on a single miRNA-target edge and should be read as an absence of predicted-only noise, "
 "not as depth of evidence."),

"ESR1": (
 "ESR1 encodes oestrogen receptor alpha, the single most consequential molecular axis in breast "
 "cancer and the target of all endocrine therapy. Activating mutations in the ligand-binding "
 "domain are an established mechanism of acquired endocrine resistance in metastatic disease "
 "[PMID:24185512; PMID:26122181]. ESR1 amplification was reported to be frequent in breast "
 "cancer [PMID:17417639], but that claim has been contested and is not consensus — consistent "
 "with our data, in which ESR1 shows copy gain in 15.9 % of tumours and loss in 30.8 %. ESR1 is "
 "not differentially expressed at the cohort level here (log2FC +0.37, FDR 0.26), which is "
 "expected: ER status stratifies the disease rather than shifting uniformly."),

"SREBF1": (
 "SREBF1 encodes SREBP-1, the master transcriptional regulator of fatty-acid and cholesterol "
 "synthesis. De novo lipogenesis is an established feature of breast and other carcinomas, and "
 "the downstream enzyme FASN is a long-standing candidate target [PMID:17882277]. Direct "
 "evidence that SREBP-1 supports breast cancer growth exists but is model-based rather than "
 "clinical [PMID:26148231]. SREBF1 participates in no 3-node FFL core in our network and enters "
 "the prioritised set purely on evidence weight; it is the one clear link from this regulatory "
 "census to tumour lipid metabolism."),

"DNMT1": (
 "DNMT1 is the maintenance DNA methyltransferase that propagates CpG methylation through "
 "replication [PMID:22641018]. It is over-expressed in breast cancer and has been proposed as a "
 "drug target particularly in triple-negative disease [PMID:32461152; PMID:35307730]. Its "
 "relevance here is internal and mechanistic rather than clinical: DNMT1 is the most "
 "CRISPR-essential transcription factor in the prioritised set (essential in 34 of 53 breast "
 "lines, mean Chronos -0.55), and the network contains several miRNA nodes carrying substantial "
 "promoter hypermethylation, so a candidate silencing enzyme and its plausible substrates are "
 "both present. The causal link between them is not tested by this study."),

"E2F3": (
 "E2F3 is a second activating E2F family member that overlaps E2F1 functionally. The breast "
 "cancer literature is thin but consistent: E2F3 has been reported to drive "
 "epithelial-mesenchymal transition, invasion and metastasis [PMID:34365840] and to support DNA "
 "damage repair, stem-like properties and therapy resistance [PMID:37499929]. Both are "
 "single-group findings in cell-line models and should be treated as suggestive rather than "
 "established. Its value in this analysis is corroborative — it reproduces the E2F1 signal "
 "independently, from a partly non-overlapping set of edges."),

# ------------------------------------------------------------------ Genes
"CCND2": (
 "Loss of cyclin D2 expression in the majority of breast cancers, and its association with "
 "promoter hypermethylation, is an old and well-replicated observation [PMID:11289162], "
 "subsequently proposed as a biomarker and drug target across lung and breast tumours "
 "[PMID:30308939]. Our data reproduce this exactly and independently: promoter delta-beta "
 "+0.128 (FDR 1.4e-31), methylation-expression rho = -0.29 (FDR 6.7e-16), and log2FC -1.40. "
 "This is the clearest instance in the prioritised set where the observed dysregulation has a "
 "documented, measured epigenetic mechanism rather than an inferred one."),

"COL1A1": (
 "COL1A1 encodes the alpha-1 chain of fibrillar type I collagen. That collagen deposition and "
 "crosslinking increase tissue stiffness, enhance integrin signalling and promote breast tumour "
 "progression is established [PMID:19931152], and COL1A1 itself has been reported to promote "
 "breast cancer metastasis [PMID:29906404]. The critical qualification, which this study "
 "quantifies, is compartmental: COL1A1 is overwhelmingly a fibroblast product [PMID:31980749], "
 "and we measure it as 146-fold enriched in cancer-associated fibroblasts over carcinoma cells. "
 "Bulk-tissue associations involving COL1A1 are therefore dominated by stromal content, and its "
 "CRISPR gene effect in carcinoma cells is essentially zero (essential in 0 of 53 lines)."),

"STAT5A": (
 "STAT5A is the prolactin/JAK2 effector controlling luminal differentiation and milk-protein "
 "gene expression. In breast cancer the established clinical finding is counter-intuitive and "
 "well replicated: nuclear, tyrosine-phosphorylated STAT5 marks favourable prognosis, and its "
 "loss predicts poor outcome and antioestrogen failure [PMID:15169792; PMID:21576635]. The "
 "apparent paradox — a growth-factor-pathway effector that is prognostically protective — is "
 "explicitly discussed in the field [PMID:23161573]. Our data agree in direction (log2FC -1.54; "
 "TCGA OS HR 0.82, p = 0.018; METABRIC OS HR 0.81, p = 1.9e-12)."),

"MYBL2": (
 "MYBL2 (B-Myb) is a MuvB/DREAM-complex transcription factor that drives late cell-cycle (G2/M) "
 "transcription and is a constituent of commercial proliferation signatures. Its amplification "
 "and over-expression in breast cancer, and the association with high grade and poor outcome, "
 "are established [PMID:32853735]. It carries the largest tumour-versus-normal fold change of "
 "any node in this study (log2FC +3.72, FDR 4.7e-89) and is among the most CRISPR-essential "
 "(essential in 20 of 53 breast lines). It is not, however, prognostic in our TCGA analysis "
 "(OS p = 0.56), a useful reminder that magnitude of dysregulation and clinical informativeness "
 "are different properties."),

"FN1": (
 "FN1 encodes fibronectin, a central organiser of the provisional matrix and a ligand for "
 "integrin adhesion. Its role in the tumour microenvironment is established [PMID:32266654], "
 "but the literature is explicit that fibronectin is context-dependent rather than simply "
 "pro-tumorigenic [PMID:31861892]. In our data FN1 behaves much like COL1A1: strongly "
 "up-regulated (log2FC +2.68), prognostic in METABRIC (OS HR 1.12, p = 3.0e-4), and completely "
 "non-essential in carcinoma cells (0 of 53 lines; mean Chronos +0.08) — the profile of a "
 "stromal product rather than a cell-autonomous dependency. Its one strong-tier repression edge "
 "from a prioritised miRNA, miR-101 -| FN1, is the strongest such edge in the set and "
 "strengthens under CAF adjustment (rho -0.296 to -0.338)."),

"PDGFRB": (
 "PDGFRB marks pericytes and fibroblasts rather than carcinoma cells, and stromal PDGFR-beta "
 "expression is an established adverse prognostic marker in breast cancer [PMID:19498003]; it is "
 "a standard marker in current CAF taxonomies [PMID:31980749]. Its presence in an unbiased "
 "prioritisation is itself a result: it demonstrates that the curated regulatory network carries "
 "a substantial non-epithelial component. Consistent with that, PDGFRB is the most heavily "
 "promoter-hypermethylated protein-coding node in the set (delta-beta +0.173, FDR 1.3e-44) while "
 "being only weakly down-regulated (log2FC -0.44), and 12 of its 13 sign-testable edges are "
 "sign-concordant (92.3 %) — the highest sign concordance of any of the 30 nodes."),

"MET": (
 "MET encodes the hepatocyte growth factor receptor; MET signalling drives scattering, invasion "
 "and therapy resistance, and transgenic Met is sufficient to induce mammary tumours with basal "
 "features in mice [PMID:19617568]. MET is clinically druggable, with an extensive inhibitor and "
 "antibody pipeline, although the approved indications are in lung and renal cancer rather than "
 "breast [PMID:32532746]. Our data are discordant with a simple oncogene reading: MET mRNA is "
 "strongly down-regulated in tumour versus normal (log2FC -2.48) and is not prognostic — a "
 "reminder that MET pathology in breast cancer is driven by activation and amplification in "
 "subsets, not by cohort-level transcript abundance."),

"CXCL12": (
 "The CXCL12-CXCR4 axis directing breast cancer metastatic homing is one of the foundational "
 "results in the chemokine-metastasis literature [PMID:11242036], and CXCL12 secretion by "
 "carcinoma-associated fibroblasts promoting tumour growth and angiogenesis is equally "
 "established [PMID:15882617]. Because the ligand is predominantly stromal and the receptor "
 "epithelial, bulk CXCL12 measurements report largely on stromal content — which is why our "
 "CXCL12 signal (log2FC -2.24; METABRIC OS HR 0.85, p = 1.3e-8) should be read as a "
 "microenvironmental readout, not as a carcinoma-cell property."),

"MMP14": (
 "MMP14 (MT1-MMP) is the principal pericellular collagenase and the physiological activator of "
 "pro-MMP2; the demonstration that tumour-cell traffic through three-dimensional collagen is "
 "controlled by MT1-MMP is a foundational and well-replicated result [PMID:15557125]. In breast "
 "cancer, MT1-MMP delivery to invadopodia is required for ERBB2-induced invasion [PMID:26899482]. "
 "MMP14 is therefore a direct effector of invasion rather than a correlate of it, and it is the "
 "most prognostic ECM node in our TCGA analysis (PFI HR 1.26, p = 0.012; DSS p = 0.0076)."),

"PLAU": (
 "PLAU encodes urokinase plasminogen activator. uPA and its inhibitor PAI-1 are among the very "
 "few molecular markers in breast cancer with level-of-evidence-1 clinical validation: a pooled "
 "analysis of 8,377 patients established independent prognostic value [PMID:11792750], and a "
 "prospective randomised adjuvant trial used uPA/PAI-1 to select node-negative patients for "
 "chemotherapy [PMID:11416112]. This is an external anchor for the prioritisation — a node whose "
 "clinical value was demonstrated prospectively, entirely independently of this analysis."),

# ------------------------------------------------------------------ miRNAs
"hsa-miR-21": (
 "miR-21 is the canonical oncomiR and was among the miRNAs first reported as deregulated in "
 "human breast cancer [PMID:16103053]; meta-analysis supports high miR-21 as a marker of "
 "unfavourable survival [PMID:26349663]. Its literature is by far the largest of any miRNA here "
 "(10,744 PubMed records, 539 also indexed under Breast Neoplasms). It functions in this study "
 "as a positive control: any credible breast cancer miRNA prioritisation must rank it first, and "
 "this one does (composite rank 1 of 219 miRNAs; log2FC +2.04, FDR 9.7e-124)."),

"hsa-miR-195": (
 "miR-195 is a miR-15/16-family member; its down-regulation in breast cancer and "
 "tumour-suppressive behaviour in breast models are reported and reasonably well replicated "
 "[PMID:21350001], including a role in tamoxifen resistance through the MIR497HG locus "
 "[PMID:36815359]. The distinctive contribution here is mechanistic rather than functional: "
 "miR-195 carries the largest promoter hypermethylation of any prioritised node "
 "(delta-beta +0.212, FDR 1.1e-45; across all 575 nodes with methylation data it is second only "
 "to its own genomic cluster partner hsa-miR-497 at +0.218), with a strong inverse "
 "methylation-expression correlation (rho = -0.50) and copy loss in 60.8 % of tumours — two "
 "independent silencing mechanisms converging on one node."),

"hsa-miR-204": (
 "miR-204 is comparatively under-studied in breast cancer (1,261 PubMed records overall, 44 also "
 "indexed under Breast Neoplasms). The breast-specific literature is small and recent; a "
 "miR-204/COX5A axis has been implicated in invasion and chemotherapy resistance in "
 "oestrogen-receptor-positive disease [PMID:32758616]. Its best-characterised biology lies "
 "elsewhere, in the retina and the TRPM3 locus. It reaches the top of our ranking on measured "
 "evidence rather than prior attention: log2FC -1.58 at FDR 7.5e-94, promoter hypomethylation "
 "delta-beta -0.174, and TCGA OS HR 0.83 (p = 0.045), with only a single 3-node FFL core."),

"hsa-miR-383": (
 "hsa-miR-383 is the most genuinely under-characterised node in the prioritised set, and we "
 "state that explicitly: PubMed indexes 273 miR-383 records in total and only 9 that are also "
 "indexed under Breast Neoplasms (3.3 %; counts retrieved 2026-09-09). Most of its literature is "
 "in glioma, medulloblastoma [PMID:23227829], hepatocellular and ovarian cancer, where it "
 "behaves as a tumour suppressor and represses the glycolytic enzyme LDHA [PMID:28043152]. The "
 "single substantial breast report restored miR-383-5p in MDA-MB-231 cells in combination with "
 "paclitaxel [PMID:34761351]. Notably, our network independently recovers the LDHA relationship "
 "(miR-383 -| LDHA, strong tier, rho = -0.160, FDR 5.1e-7) — external corroboration of an edge "
 "that no breast cancer study had reported. An under-studied node surfacing from an unbiased "
 "ranking is the strongest kind of candidate this analysis can produce."),

"hsa-miR-124": (
 "miR-124 is an established tumour-suppressive miRNA that is silenced by promoter "
 "hypermethylation in many carcinomas; in breast models it represses STAT3 [PMID:30896795] and "
 "its loss de-represses a MALAT1/CDK4/E2F1 axis [PMID:26918449]. Our data support the silencing "
 "mechanism strongly (delta-beta +0.182 across 47 probes, FDR 1.7e-48) but not the expression "
 "consequence — miR-124 is not differentially expressed in TCGA-BRCA (log2FC +0.13, FDR 0.10). "
 "Despite that it carries the best survival signal of any miRNA in the study (PFI HR 1.20, "
 "p = 3.0e-4), which is a direct argument against selecting nodes on fold change."),

"hsa-miR-155": (
 "miR-155 is an inflammation-associated oncomiR; its function in breast cancer through "
 "repression of SOCS1 is established [PMID:20354188], and systematic review supports an adverse "
 "prognostic association [PMID:32823863]. Its role in this study is evidential rather than "
 "biological: 19 of its 35 miRNA-target edges (54.3 %) carry low-throughput experimental "
 "validation — the highest validated fraction of any miRNA in the network with at least ten "
 "target edges — so it is the miRNA whose network neighbourhood can be trusted most. Its very "
 "strong reciprocal wiring with STAT1, SPI1 and IRF1 is, however, exactly the kind of mutual "
 "TF-miRNA pair that our own reciprocity audit flags as partly a database-curation artefact."),

"hsa-miR-429": (
 "miR-429 belongs to the miR-200 family, whose repression of the E-cadherin repressors ZEB1 and "
 "ZEB2 to enforce the epithelial state is one of the best-established regulatory circuits in "
 "cancer biology [PMID:18376396; PMID:18381893]. Its role in breast cancer specifically is "
 "contested: it has been reported both to inhibit migration and invasion [PMID:25405387] and to "
 "act as an oncogene in triple-negative disease by degrading DLC1 [PMID:37737712]. Our data "
 "recover the canonical miR-429 -| ZEB2 edge with the expected sign (rho = -0.301, FDR 1.0e-22), "
 "but that correlation collapses when tumour CAF content is partialled out (rho_adj = -0.097), "
 "so in bulk tissue it is substantially a compartment effect."),

"hsa-miR-141": (
 "miR-141 is the second miR-200-family member in the prioritised set and shares the ZEB1/ZEB2 "
 "targeting mechanism [PMID:18381893]; circulating miR-200c/miR-141 have been evaluated as "
 "outcome markers in breast cancer with mixed results [PMID:25885099]. It functions here as "
 "internal replication of the miR-429 result and carries the same caveat: miR-141 -| ZEB2 is "
 "rho = -0.395 raw but only -0.095 after CAF adjustment. It also supplies a strong-tier "
 "repression edge onto another prioritised node (miR-141 -| STAT5A, rho = -0.268, FDR 5.3e-18) "
 "which, unlike its ZEB2 edge, largely survives CAF adjustment (rho_adj = -0.182)."),

"hsa-miR-34a": (
 "miR-34a is the principal p53 effector miRNA [PMID:17554337] and the subject of the first miRNA "
 "mimic to enter clinical trial: MRX34, a liposomal miR-34a mimic, reached phase 1 in solid "
 "tumours but the trial was halted for immune-related adverse events [PMID:27917453; "
 "PMID:32238921]. That history is the reason to include it — it is the set's concrete precedent "
 "for, and cautionary tale about, therapeutic miRNA restoration. It is not differentially "
 "expressed here (log2FC +0.09, FDR 0.36) yet is prognostic (OS HR 1.23, p = 0.018), and 13 of "
 "its 32 target edges are experimentally validated."),

"hsa-miR-101": (
 "miR-101 represses EZH2, and genomic loss of the miR-101 loci leading to EZH2 over-expression "
 "is an established mechanism first shown in prostate cancer [PMID:19008416] and subsequently "
 "reported in breast models [PMID:28034643]; our copy-number data are consistent with it, with "
 "loss at hsa-mir-101-1 in 29.2 % of TCGA-BRCA tumours against gain in 13.4 %. This axis is the strongest single validated result "
 "in the present study: it is one of only two miRNA-target relationships that survive every "
 "control we applied, and it strengthens rather than collapses under CAF adjustment "
 "(rho = -0.190 to rho_adj = -0.230), the opposite of the behaviour expected from a "
 "stromal-content artefact — as does its second strong-tier edge onto a prioritised node, "
 "miR-101 -| FN1 (rho = -0.296 to -0.338). miR-101 is also the only miRNA in the set that is "
 "significant on both survival endpoints (OS HR 0.80, p = 0.0047; PFI HR 0.78, p = 0.0028), "
 "although neither survives genome-wide FDR correction."),
}

VERDICT = {
"E2F1":   "The one node on which topology and evidence agree — 9th by FFL participation and 1st on the seven-domain composite — which is what makes the disagreement everywhere else interpretable.",
"EZH2":   "Anchors the epigenetic programme, carries the study's second fully validated miRNA axis (miR-101 -| EZH2), and is directly druggable with an approved catalytic inhibitor (tazemetostat), albeit outside breast cancer.",
"GATA3":  "Shows that the evidence composite recovers established lineage biology (12.3 % mutation rate, mRNA-protein rho 0.86) that FFL participation alone would have discarded at 3 cores.",
"BRCA1":  "A sanity check on the prioritisation: any ranking of breast cancer regulators that omitted BRCA1 would be suspect, and it enters on evidence with only 2 FFL cores.",
"JUN":    "The strongest down-regulated TF hub (45 cores, every one as regulator), and evidence that high FFL participation does not imply up-regulation in tumour.",
"EGR2":   "An under-studied node surfaced by unbiased ranking rather than prior attention — 0 FFL cores and 8 breast cancer papers, yet the largest down-regulation of any transcription factor in the set and a METABRIC survival association at q = 3.2e-4.",
"ESR1":   "Demonstrates why differential expression is a poor prioritisation criterion: not differentially expressed at all (FDR 0.26) yet the dominant clinical axis in the disease, and the most promoter-hypomethylated TF here.",
"SREBF1": "The network's only link to lipid metabolism, and a second example of a node with zero FFL participation that the evidence weighting promotes.",
"DNMT1":  "Supplies a candidate mechanism for the promoter hypermethylation seen at the network's own miRNA nodes, and is the most CRISPR-essential TF in the set (34 of 53 lines).",
"E2F3":   "Independent internal corroboration of the proliferation programme, and the node whose edges retain the most signal after CAF adjustment (median |rho_adj| 0.296).",

"CCND2":  "The second most FFL-embedded gene in the network after VEGFA (136 cores, every one occupied as target) and the clearest case of a co-regulated target whose loss has a documented epigenetic mechanism.",
"COL1A1": "The anchor of the ECM programme and the node around which the study's principal methodological finding turns — 146-fold CAF-enriched, prognostic in METABRIC only after CAF adjustment, and essential in 0 of 53 carcinoma lines.",
"STAT5A": "A second demonstration that topology and evidence measure different things (0 FFL cores, 3rd on the composite), and the strongest METABRIC survival signal of any node here.",
"MYBL2":  "The proliferation programme's strongest single marker (log2FC +3.72) and proof that fold-change magnitude and prognostic value are dissociable.",
"FN1":    "Extends the ECM programme beyond collagen and shows the compartment caveat is general, not specific to COL1A1.",
"PDGFRB": "Included deliberately as a stromal marker; its arrival in an unbiased top ten is itself evidence that the curated network carries a substantial non-epithelial component.",
"MET":    "The most extensively drugged gene in the set (118 DGIdb interactions, 18 approved antineoplastics, second only to ESR1 across all 30 nodes) and a caution that transcript-level down-regulation does not exclude pathway activation.",
"CXCL12": "Connects the network to metastatic tropism through a foundational, well-replicated mechanism, while illustrating the ligand-stromal / receptor-epithelial compartment split.",
"MMP14":  "The functional endpoint of the ECM programme — a direct effector of invasion rather than a correlate — and the most prognostic ECM node in TCGA.",
"PLAU":   "An external anchor: the one node in the set whose clinical value was established prospectively in a randomised trial, entirely independently of this analysis.",

"hsa-miR-21":  "The positive control that validates the ranking procedure itself: the strongest miRNA differential expression in the network and rank 1 of 219 on the composite.",
"hsa-miR-195": "The cleanest case in the study of an epigenetically silenced miRNA — the largest promoter hypermethylation of any prioritised node, plus copy loss in 61 % of tumours.",
"hsa-miR-204": "A strongly dysregulated, clinically associated miRNA with almost no FFL participation — exactly the node a topology-only ranking discards.",
"hsa-miR-383": "The clearest novel candidate the study produces: entirely coherent FFL cores, 9 breast cancer papers, and independent recovery of its published LDHA target.",
"hsa-miR-124": "Dissociates prognostic value from differential expression — the best miRNA survival signal in the study with no significant fold change — and carries a documented silencing mechanism.",
"hsa-miR-155": "The miRNA whose network edges are most trustworthy (54 % validated), and simultaneously the clearest illustration of the reciprocal TF-miRNA curation artefact.",
"hsa-miR-429": "The most FFL-embedded miRNA in the set (31 cores) and the EMT axis — but its canonical ZEB2 edge is shown here to be largely a stromal-compartment effect.",
"hsa-miR-141": "Internal replication of the miR-429 EMT result, carrying the same CAF caveat, plus a strong-tier repression edge onto STAT5A that survives it.",
"hsa-miR-34a": "The set's precedent for therapeutic miRNA restoration, including the negative clinical experience with MRX34 that any translational claim must acknowledge.",
"hsa-miR-101": "Together with miR-29a -| collagen, the study's positive result — it strengthens rather than collapses under stromal adjustment, and is the only miRNA in the set significant on both OS and PFI.",
}
