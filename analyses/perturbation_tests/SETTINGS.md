# analyses/perturbation_tests — fixed settings (written 2026-09-27 22:40, before any dataset was searched or opened)

The brief's rules and decision rules are fixed as sent. This file records how the rules are put into practice.
Anything decided here is not changed after a result is seen. If a setting has to change, both versions are run and
both are reported, with the change logged at the end of this file.

## Common to P1–P3

- **Search.** NCBI GEO (db=gds, entry type GSE) is searched with the analysis plan's query terms through E-utilities. The date,
  query and hit count are logged per part. KnockTF 2.0 (P1a, P3b), the ENCODE portal (P1a) and the L1000 landmark
  list (P1a) are searched as named.
  - Every hit is listed in DATASETS.tsv, included or excluded with a reason.
  - Series without sample-level data in GEO (for example SuperSeries shells or non-expression data) are excluded as
    "not an expression dataset".
- **Dataset unit.** A dataset is one contrast: the perturbed arm(s) against the matched control in one cell type.
  - A series with several cell types gives one dataset per cell type that has its own matched control.
  - Several siRNAs, shRNAs or guides against the same gene in the same cell type are pooled into one perturbed arm.
  - A series with several time points or doses: the latest time point and the highest dose that have a matched control.
  - A series with a stimulus (for example TGF-β): the unstimulated contrast if present. Otherwise the stimulated one,
    flagged.
- **Replicates.** At least 2 samples in each arm of the contrast.
- **Raw data**, in this order of preference:
  1. RNA-seq: NCBI-generated raw counts (GEO "RNA-seq raw counts", GRCh38.p13), which NCBI computes from the SRA
     reads.
  2. RNA-seq: the submitter's raw count matrix.
  3. Arrays: probe-level raw data (Affymetrix CEL; Illumina non-normalised; Agilent raw).
  4. Otherwise processed values (series matrix or normalised tables), flagged "processed" in DATASETS.tsv and in
     every result row.
- **Differential expression.**
  - RNA-seq: edgeR filterByExpr (defaults), TMM normalisation, limma-voom, design `~ group`.
  - Arrays: RMA for Affymetrix CEL, neqc for Illumina, normexp + quantile for Agilent. For processed values, log2 if
    not already log scale (max > 50). Then limma, design `~ group`.
  - Probe to gene: the probe with the highest mean expression across all samples.
  - log2FC is perturbed minus control. SE = coefficient / moderated t. A gene removed by filterByExpr has no value
    (NOT FOUND in that dataset).
- **Expression calls.**
  - RNA-seq "expressed" = kept by filterByExpr.
  - Arrays "expressed" = mean log2 expression in the controls at or above the median of all genes.
- **Pooling.** metafor `rma(yi, sei, method = "REML")`, with 95 % CI from the default z-based interval. k = 1 is
  pooled the same way, which returns the dataset's own estimate and SE. Fold change = 2^log2FC.
- **Cell class.**
  - Fibroblast lineage: fibroblasts of any tissue, CAFs, myofibroblasts, hepatic or pancreatic stellate cells, and
    mesenchymal stromal/stem cells.
  - Carcinoma: cell lines derived from epithelial cancers.
  - Other: everything else, including normal epithelium, endothelium, immune cells, leukaemia, lymphoma, sarcoma,
    glioma, melanoma, neuroblastoma, pluripotent cells and HEK293.
  - The class is assigned from the GEO sample characteristics and the cell-line identity, before any expression value
    is read.
- **Seeds.** 20250908 for every random step (bootstrap, shuffles).

## P1a

- **TF knockdown verified:** TF log2FC ≤ −0.74 in the recomputed contrast. A CRISPR knockout without an mRNA drop is
  included only if the series text reports loss of the protein, and is flagged.
- **COL1A1 expressed in the controls:**
  - RNA-seq: mean TPM ≥ 10 in the controls. Uses NCBI's TPM table when it exists, otherwise TPM computed from counts
    and NCBI gene lengths.
  - Arrays: the COL1A1 gene (best probe) at or above the median of all genes in the controls.
- **Positive control.** The TF's TRRUST v2 targets annotated "Activation" are compared with expressed genes that are
  not TRRUST targets of the TF in any mode.
  - Measure: difference of the median log2FC, with a two-sided Wilcoxon rank-sum test.
  - TRRUST v2 file: the project's copy, recorded by md5.
- **Decision.**
  - Applied to the fibroblast-lineage REML pool as written in the analysis plan (0.80–1.25 band on the pooled fold change and
    its 95 % CI).
  - Other cell types are pooled and reported without a decision.
  - NFKB1, RELA and SP1 are reported with the same rule, marked "descriptive".
- **LINCS L1000:** used only if COL1A1 or COL3A1 is a landmark gene in the GSE92742 gene table (`pr_is_lm`).

## P1b

- **Gain:** a miR-29a, -29b or -29c mimic, or overexpression. **Loss:** an inhibitor, antagomir, sponge or knockout.
- **Repression effect:** −log2FC for gain, +log2FC for loss.
- **Collagen effect per dataset.**
  - The pair COL1A1/COL3A1 is analysed as one pseudo-gene: the per-sample mean of the two genes' normalised log2
    values (voom log-CPM for RNA-seq, with the mean of the two genes' voom weights).
  - It is added as a row before lmFit/eBayes, so its SE reflects the correlation between the two genes.
  - COL1A1 and COL3A1 are also reported singly.
- **Meta-regression.** `rma(yi, sei, mods = ~ fibroblast, method = "REML")` over fibroblast-lineage and carcinoma
  datasets only. The test is the z-test of the fibroblast coefficient.
- **Target shift.**
  - Site definitions: TargetScan 8.0 files in the project cache (`cache/direct/ts/`), human rows (Gene Tax ID 9606),
    miR-29-3p family.
  - "Site" = 8mer or 7mer-m8 in Conserved_Site_Context_Scores. "No site" = no miR-29 entry of any type in either the
    conserved or the non-conserved file.
  - Shift = median sign-aligned effect of the expressed site genes minus that of the expressed no-site genes.
  - SE from 1,000 bootstrap resamples of genes within each class.
- **Baseline expression:** controls' mean TPM (RNA-seq) or percentile among genes (arrays).
- **Positive control:** miR-101 gain datasets (same inclusion rules), with the EZH2 log2FC pooled by REML.

## P2

- **Clusters.**
  - Membership comes from the project's 10 kb co-transcription layer (`data/layer_miRNA_miRNA.tsv`, `cluster_id`), the
    definition behind the 70.1 %.
  - Where a named cluster has no `cluster_id` in the network, its miRBase 10 kb cluster is used, flagged.
- **Perturbed members:**
  - the members whose mimic or inhibitor was applied;
  - the members encoded by the deleted or overexpressed locus.
- **Classes.**
  - Network edges: the analysed network's miRNA→target edges (`data/canonical_edges.tsv`).
  - TargetScan: conserved 8mer/7mer-m8 sites of the members' families.
  - Six-node class: G1 of every BHAT6 composite instance (analyses/six_node_pattern `S1/obs`) whose M1 and M2 are both
    perturbed members.
- **Primary statistic:** difference of medians, targets of ≥ 2 members minus targets of exactly 1 member. SE comes from
  1,000 bootstrap resamples within each class, and the Wilcoxon p is reported.
- **Site adjustment:**
  - model: `effect ~ class + total_sites + log10(UTR length)`;
  - total sites and UTR length come from TargetScan 8.0 (human 3′ UTR of the representative transcript);
  - the class coefficient is pooled by REML.
- **Mouse datasets:** mouse genes mapped to human by one-to-one orthologues (MGI HOM_MouseHumanSequence).

## P3

- **P3a.**
  - Tiers come from `data/edge_evidence_tier.tsv` (`tier`: strong, weak, predicted).
  - Background: expressed genes with no edge from the perturbed miRNA and no TargetScan site (any type) for its family.
  - Shift = median sign-aligned effect of the tier's targets minus the background median, with SE from 1,000
    bootstrap resamples.
  - If more than 40 datasets qualify, the 40 with the most samples are kept (ties broken by GSE number) and the rest
    are listed.
- **P3b.**
  - Datasets: KnockTF 2.0 human datasets for network TFs, meeting P1a's criteria except the COL1A1 requirement.
  - Values: KnockTF's own differential expression, flagged "processed", unless the raw data are processed as in P1a.
  - The annotated mode is TRRUST's Activation or Repression. Edges with no mode enter the |log2FC| test only.

## P4

- **Replication settings:**
  - the 1,626 calls of `results/v3/seqreg_ext_occlusion_allsites.csv`;
  - 5 dinucleotide-preserving shuffles per call, with the shuffle function and seed of `scripts/09_regulatory_evidence/sequence_models/41_seqreg_ext_occlusion_null.py`;
  - the three point-mutation controls of `scripts/09_regulatory_evidence/sequence_models/45_seqreg_ext_splicedonor_control.py`.
- **Model and tracks.** Borzoi replicates 0–3 (borzoi-pytorch port), in bf16, batch size 1. Tracks: fibroblast RNA-seq
  and fibroblast CAGE, listed by index in P4/p4_tracks.tsv before any prediction.
- **Score.** Predicted coverage summed over the gene's exon bins, after undoing Borzoi's squashed scale as in the
  official scoring code.
- **Shared GPU.** No GPU work starts until the authors confirm it.

## Changes after this file was written

1. **P4 decision rule revised, 2026-09-27 22:44, before any P4 output existed (P4 launched 22:45:47).** The GPU run had not started, and no
   Borzoi prediction had been made.
   - **Why.** Criterion 1 as first sent required the most influential NF-κB, SP and ETS1/ETS2/GABPA calls at COL1A1
     to rank outside the top 10 % (rank > 62 of 618). Enformer itself fails that: its best SP call ranks 33 of 618
     (`results/v3/seqreg_ext_occlusion_calibrated.csv` line 3). This was reported to the authors, who replaced the rule.
   - **Revised rule, fixed from now.** REPLICATED if all three hold, NOT REPLICATED otherwise:
     1. At COL1A1, the most influential call overall is not an NF-κB, SP or ETS-family call, and the most
        influential call of each of these three families changes the prediction by at most one third as much as the
        most influential call overall.
     2. At COL3A1, the most influential call overlaps the exon-1 donor (+184 to +210).
     3. Destroying the donor alone changes the prediction more than substituting the GGAA core (+191 to +194) while
        sparing the donor.
   - **How it is applied.**
     - "Most influential" = the most negative mean change (the Enformer convention of
       `scripts/09_regulatory_evidence/sequence_models/44_seqreg_ext_synthesis.py`).
     - Families are those of that script, lines 14–16: NF-κB = NFKB1, NFKB2, RELA, RELB, REL; SP = SP1, SP2, SP3;
       ETS = ETS1, ETS2, ELK1, ELK4, GABPA.
     - "Overlaps the donor": the call's span [rel, rel + length − 1] intersects +184..+210.
     - Criterion 3 compares donor GT→AA with GGAA→CCTT, the first-listed variant of each in `45_`. The GT→CC and
       GGAA→TTTT variants are reported alongside.
     - The primary score is the fibroblast RNA-seq exon sum. The CAGE score is reported alongside.
     - The rule is applied to the deposited Enformer values first, as a sanity check.
2. **Clarification of the dataset-unit rules, 2026-09-27 22:56, before any expression value was read.** Some series
   are time courses of a stimulus rather than of the knockdown, for example GSE122918: WI-38 cells with siRNA,
   sampled 0, 72 and 144 h after RAS induction. In such series the stimulus rule takes precedence: the contrast at the
   unstimulated time point (0 h) is primary. The latest stimulated time point with a matched control is reported as a
   secondary result, flagged "stimulated", and does not enter the decision. The "latest time point" rule applies to
   time courses of the perturbation itself.
   Samples shared by two series (GSE144397 re-uses GSE122918's siRNA samples GSM3488178–GSM3488219) count once,
   under the series that first deposited them.
3. **Two further clarifications, 2026-09-27 22:57, before any expression value was read.**
   - All control samples of the same cell type and condition (for example two different control siRNAs) are pooled
     into the control arm, as several siRNAs against the TF are pooled into the perturbed arm. An untransfected
     parental control is accepted and flagged.
   - The COL1A1 expression criterion may be checked on the NCBI TPM table or on the series matrix before the raw data
     are fetched. A dataset that fails it is excluded without further processing. Every dataset that passes is
     analysed from raw data as above.
4. **Implementation detail, 2026-09-27 23:20, before any expression value was read.** RMA is run over the
   arrays of each contrast (perturbed plus control arrays only). Secondary contrasts, such as a later stimulated time
   point, are listed in P1a/p1a_secondary.tsv.

## P5 (written 2026-09-27 23:35, before any P5 expression value was read)

- **Atlases.**
  - FANTOM5 miRNA atlas (de Rie et al. 2017): the viewer's human data package, with mature-miRNA CPM as provided
    (human.srna.cpm.txt).
  - McCall et al. 2017, "Toward the human cellular microRNAome": the Bioconductor data package microRNAome. Values
    are as provided; if the package gives counts only, CPM = count / total miRNA counts of the sample × 10^6, stated
    as such.
- **Samples.**
  - Whole-cell samples only. FANTOM5 samples marked "(nuclear fraction)" or "(cytoplasmic fraction)" are excluded.
  - Fibroblasts: every sample whose cell type or description names a fibroblast.
  - Mammary fibroblasts: "Fibroblast - Mammary".
  - Mammary epithelial cells: "Mammary Epithelial Cell" (FANTOM5), or the package's mammary epithelial cell type.
  - If an atlas has no mammary epithelial cells, its other epithelial cells are used and flagged.
- **Values.** hsa-miR-29a-3p, -29b-3p, -29c-3p and -101-3p. The group value is the mean over its samples. Ratio =
  fibroblast mean / mammary-epithelial mean.
- **Decision per miRNA.**
  - FIBROBLAST-ENRICHED if the ratio is ≥ 2 in every atlas that has both groups.
  - EPITHELIAL-ENRICHED if it is ≤ 0.5 in every such atlas.
  - Otherwise NO CONSISTENT DIFFERENCE.
5. **Descriptive addition requested by the authors after the primary P1a result, 2026-09-27 23:45.** For
   GSE122918 at 0 h, each ETS1 siRNA is analysed separately (2 arrays each) against the same 4 controls. RMA runs over
   those 6 arrays. The results are labelled "descriptive, added after the primary result" and change neither the
   pooled estimate nor the decision.
6. **TargetScan implementation detail, 2026-09-27 23:48, before any target shift was computed.** Site classes
   come from TargetScan 8.0's "Summary Counts, all predictions" file, which counts sites per gene and miRNA family,
   human (Species ID 9606).
   - "Site": conserved 8mer + conserved 7mer-m8 ≥ 1 for the family.
   - "No site": total conserved + non-conserved sites = 0 for the family.
   These are the same site classes as in the P1b rule, counted by TargetScan itself. The project's cached context++
   files would need site-type codes decoded, so they are not used.
7. **P1b sensitivity analysis and a bug fix, 2026-09-28 00:07, after interim P1b values had been seen.**
   - **The analysis.** In several carcinoma RNA-seq datasets, filterByExpr drops COL3A1 (and in some COL1A1), so the
     COL1A1/COL3A1 pseudo-gene is undefined there and the dataset leaves the primary meta-regression.
     - The analysis as specified stays primary, and it alone carries the decision.
     - A second version keeps COL1A1 and COL3A1 whenever they have any counts ("force-keep"). It is reported alongside,
       labelled sensitivity.
   - **The bug.** For count data, _de.R wrote the "expressed" flag unnamed, so it was NA for every gene. This emptied
     the target-shift and positive-control gene sets of the RNA-seq datasets. It is fixed, and every RNA-seq contrast
     is re-run. The fix changes no rule.
   - **A second fix of the same kind.** _de2.R had dropped every probe with a missing value in any array, which removed
     COL3A1 and EZH2 from GSE65704. Rows with at least 2 values are now kept, and limma handles the missing values. No
     rule changes.

## P3 (written 2026-09-28 00:12, before any P3 value was computed)

- **P3b edges.**
  - The analysed network has 770 TF→target edges (`data/canonical_edges.tsv`, edge_type TF_target), all present in
    the project's TRRUST v2 file. The "TRRUST" and "all edges" sets are therefore identical, and "hTFtarget-only" is
    empty: NOT FOUND.
  - The annotated mode is TRRUST's (Activation / Repression). Unknown-mode edges enter the |log2FC| test only.
- **P3b datasets.**
  - KnockTF 2.0 human datasets whose TF is one of the 157 network TFs with TF→target edges.
  - Included if KnockTF's dataset page gives ≥ 2 control and ≥ 2 treated samples, and the TF's own Log2FC in KnockTF's
    differential expression is ≤ −0.74. KnockTF's values are processed and flagged as such.
- **P3b measures.**
  - Expressed genes: Mean_Control ≥ the dataset's median.
  - Difference of medians of |Log2FC|: the TF's network targets minus expressed genes that are not TRRUST targets of
    the TF in any mode. Two-sided Wilcoxon test; SE from 1,000 bootstrap resamples (seed 20250908).
  - Sign agreement: an Activation target is expected to fall on knockdown, a Repression target to rise; Log2FC = 0
    counts as disagreement.
  - Pooling: REML for the |log2FC| difference. Pooled agreement = all agreements / all signed targets, with a
    two-sided binomial test against 0.5.
  - Decision: SUPPORTED if the pooled difference is > 0 with a 95 % CI excluding 0, AND the pooled agreement is > 50 %
    with binomial p < 0.05. Otherwise NOT SUPPORTED.
- **P3a.**
  - Datasets: human gain or loss datasets for any network miRNA, found with the P1b search pattern applied to each
    network miRNA, plus the P1b and P2 datasets. If more than 40 qualify, the 40 with the most samples are kept.
  - The perturbed mature miRNA maps to its network node by removing the arm suffix (hsa-miR-29a-3p → hsa-miR-29a).
  - Tiers come from `data/edge_evidence_tier.tsv`.
  - Background: expressed genes with no network edge from that miRNA and no TargetScan site of any type for its family.
  - Shift = median sign-aligned effect of the tier's targets minus the background median, with bootstrap SE.
  - Per tier: REML pool; SUPPORTED if the pooled shift is > 0 with a 95 % CI excluding 0, otherwise NOT SUPPORTED.
  - Predicted minus strong: the per-dataset difference (SE = sqrt of the summed variances), pooled by REML.
8. **P2 implementation details, 2026-09-28 00:19, before any P2 class comparison was computed.**
   - **Member seeds.** The member's arms whose TargetScan family is conserved (Family Conservation ≥ 1 in
     miR_Family_Info). If no arm has a conserved family, all arms are used.
   - **3′ UTR length.** The longest human (9606) 3′ UTR per gene symbol in TargetScan 8.0 UTR_Sequences, gaps removed.
   - **Total sites.** Conserved plus non-conserved sites, summed over the members' seeds.
   - **Which datasets carry the decision.** The five primary clusters' human datasets. Human secondary-cluster datasets
     are reported separately.
   - **Unspecified mimics.** A "miR-29 mimic/ASO" whose members are not named (GSE159688, GSE100735, GSE53143) does not
     show that ≥ 2 members were perturbed, and is excluded.
9. **GSE140424 (SP1, adipocytes, other class), 2026-09-28 00:20.** oligo RMA of its HTA 2.0 CEL files
   crashed twice ("free(): invalid pointer"). The processed series-matrix values are used instead, flagged
   "processed (raw RMA failed)".
10. **P3b data completeness, 2026-09-28 00:25.** KnockTF's bulk file (knocktf_v2_main_human.txt) holds
    only the DataSet_01/02 series and 4 DataSet_03 series. The tables of the other 157 network-TF datasets
    (DataSet_03/04) are fetched from KnockTF's per-dataset endpoint, the table shown on each dataset page. A first pass
    on the bulk file alone (25 datasets with an estimate) is superseded by the complete run, and both are recorded.
11. **Bootstrap seeding, 2026-09-28 00:30.** Every bootstrap now uses its own generator, seeded from 20250908
    and the dataset and test (_seed.py), so an SE no longer depends on which datasets were processed before it. Interim
    values sent earlier from runs that shared one generator across datasets differ in the third decimal of a few SEs.
    No rule changes.
12. **Second reading of the P1a and P1b hits, 2026-09-28 00:55.** The first-pass screen only ranks hits for reading.
    Hits it had set aside (P1a 149, P1b 17) were then read: sample lists for every hit whose samples name one of the
    TFs or miRNAs, series titles for the rest. No rule changes. Every hit is in DATASETS.tsv with its reason.
    - New P1a candidates: GSE171677 (ETS1, RELA; primary CD4+ T cells), GSE247495 (NFKB1, RELA; Akata, with LMP1
      induced, so flagged "stimulated"), GSE247497 (NFKB1, RELA; GM12878). All are "other" cell types.
    - The two "check" entries were resolved: GSE27869 has one array per siRNA and is excluded. GSE50588's SP1 and RELA
      arms are matched to the non-silencing arrays of the same transfection plates.
    - These passed through the same COL1A1 screen and knockdown check as before (P1a/p1a_col1a1_screen_reread.tsv).
    - Two implementation fixes, neither changing a rule:
      - A sample missing from NCBI's counts is left out and noted, provided each arm keeps ≥ 2 samples.
      - The COL1A1 screen falls back to overlap-based probe mapping when a platform's annotation has no symbol
        column, and states why COL1A1 is missing when it is.
13. **Second reading of the P2 hits, 2026-09-28 01:00–01:45.** This was after the first-pass P2 values had been seen.
    The first pass had read only the hits the screen ranked highest. It had also excluded GSE137048 on its submitter
    table (1,880 genes), without checking NCBI's counts, which cover a valid contrast. Every hit was then read (P2/p2_reasons.tsv).
    - **Human primary**, added: GSE137048, GSE193482, GSE99050 and GSE99134.
    - **Human secondary**, added: GSE51217, GSE193479, GSE262478 and GSE20158. GSE262478 is labelled Mus musculus in
      GEO, but its samples are human MDA-MB-436 cells.
    - **Mouse**: the qualifying series in P2/p2_mouse_reread.py.
    - **Both versions are reported.** The first-pass files are kept as P2/*_first_pass.*.
    - **Mouse details, fixed before any of the added mouse contrasts was computed:**
      - One contrast per series and tissue or cell type. Where a series also has partial deletions (GSE63813), only
        the whole-cluster knockout is used.
      - A perturbation spanning several clusters takes the union of their members.
      - Knockout arms that differ only in an outcome (GSE57112, normo- and hypoglycaemic) are pooled.
      - Stimulated-only designs are flagged, and so are processed tables.
      - Series whose table columns cannot be matched to samples are not used: GSE250213 and GSE235186.
      - Series whose sample annotations conflict are not used: GSE95579.
14. **P3a selection procedure, 2026-09-28 01:10.** This was after the interim P3a values (P1b and P2 datasets only) had
    been seen, and before any of the 38 added contrasts was analysed. The selection procedure is in
    P3a/p3a_triage_note.txt. The 40-with-most-samples rule is unchanged, and the interim values are kept as
    P3/*_interim_P1b_P2_only.*.
15. **Annotation and parsing fixes, 2026-09-28 01:30–01:40.** None changes a rule.
    - HuGene 2.0 ST transcript clusters (GPL16686) take their symbols from hugene20sttranscriptcluster.db, because the
      platform table has no symbol column.
    - HuEx 1.0 ST arrays are summarised to transcript clusters and annotated with GPL5175.
    - GSE99050's CEL files are HuGene 2.0 ST. They are processed with the standard design and annotated with GPL16686,
      because its GEO platform GPL21970 is a Brainarray remapping.
    - The P4 score files hold each value as a NumPy repr, "np.float64(x)". p4_summarise.py parses the number back at
      full precision.
16. **P5 sample groups corrected, 2026-09-28 10:15, after the P5 result had been delivered.** A wording check of the
    placed P5 sentences against the sample tables found two kinds of sample in the fibroblast groups that are not
    whole-cell fibroblasts. Both versions are reported.
    - **FANTOM5.** The whole-cell rule above names "(nuclear fraction)" and "(cytoplasmic fraction)", but FANTOM5 also
      marks fractions "(cytosolic)" and "(nuclear)". These 18 samples are now excluded too, 24 in all. Two of them,
      "Fibroblast - Dermal (cytosolic)" and "(nuclear)", were among the 30 fibroblasts, which are now 28. The mammary
      groups are unchanged.
    - **microRNAome.** The CellType "iPSC_fibroblast" names a fibroblast, but the package classes its 32 samples as
      "Stem": they are induced pluripotent stem cells. They are no longer counted as fibroblasts, which are now 96 of
      128. Their mean miR-29a-3p is 175 CPM (P5/p5_other_cell_types.tsv line 800), against 1.314e+04 in the 96
      fibroblasts (P5/p5_summary.txt line 4).
    - **Effect.** miR-29c-3p changes from EPITHELIAL-ENRICHED (ratios 0.273 and 0.491) to NO CONSISTENT DIFFERENCE
      (0.293 and 0.614). miR-29a-3p, miR-29b-3p and miR-101-3p stay NO CONSISTENT DIFFERENCE. The mammary-fibroblast
      ratios do not change.
    - **Files.** P5/p5_atlas.py writes the corrected version (p5_summary.txt, p5_groups.tsv,
      p5_other_cell_types.tsv). With --as-delivered it reproduces the first version byte for byte
      (P5/*_as_delivered.*).
