# analyses/bhat_pattern_analysis: the twelve Bhat FFL networks, analysed with the paper's settings (report only)

**The request.** Run each of the paper's analyses on the twelve FFL networks of `DEPOSIT_ffl_networks_2026-09-28.zip`
(md5 182a1629…) and report the result beside the paper's value. The networks are three classes (miRNA-FFL, TF-FFL,
composite FFL) at three to six nodes, following Bhat et al. (2024) on the census graph. Nothing goes into the
manuscript.

**When settings were fixed.** Each section below was written before its analysis was run, and carries the time it was
written. A change made after seeing a result is logged under "Changes", with both versions reported. Analyses added
after the first results are marked "added".

**The networks.** The instances, node sets and edges come from `make_ffl_networks.py` (data/ffl_networks), imported
read-only. Every script checks that the rebuilt network files are byte-identical to the deposited ones (by md5).

**What a decision line judges.** Each analysis ends with one line: supports, contradicts, or leaves unchanged the
paper's six-node conclusion it bears on. The conclusion is quoted from the 27 Sep text (Results and Discussion), so the
current text may word it differently.

## T. Topology (written 2026-09-28 16:59)

- **T1, composition.** For each network: instances, distinct node sets, nodes by type, and edges by type and by source
  layer. These are taken from the deposited files.
- **T2, degree and hubs.** Each network is analysed as undirected, with the NetworkAnalyzer definitions of
  analyses/six_node_pattern_networks `topology_mcc.py`:
  - neighbours, and degree counted as edges;
  - betweenness within each connected component, normalised;
  - closeness and clustering;
  - cytoHubba MCC, the sum over the maximal cliques containing a node of (|C| − 1)!.
  - Reported: the degree distribution (minimum, quartiles, maximum, mean, and the share of nodes with degree 1), and
    the top 10 nodes by MCC, by neighbours and by betweenness.
  - Beside it, the same measures on the analysed network (the interaction table without the 30 miRNA–miRNA edges,
    587 nodes, 6,829 edges). Also each top-10 list's overlap with the analysed network's top 10.
- **T3, overlap between networks.** For every pair of networks, the Jaccard index of their node sets and of their edge
  sets. An undirected link is one unordered pair; a directed arc is an ordered pair.
- **T4, overlap with the typed cores.** The paper's 1,649 typed cores are `results/v5/tables/TableS2_ffl_cores.csv`,
  Table S4 in the paper. For each network:
  - the share of its nodes that belong to any typed core;
  - the share of its edges that some typed core uses (a core uses R→M, M→T and R→T, plus M→R for a composite).
- **T5, overlap with D1–D4 modules.** For each network, the share of distinct instance node sets that are D1–D4 modules
  of the same size on the census graph. The graph is the one of `scripts/03_ffl_census/03_ffl_census_nolegacymirna.py`
  `build_graph()`, including its TF↔miRNA contraction. The test is that script's `is_ffl` on the induced subgraph.
- **Paper values placed beside the results:**
  - Table 3: 5,833, 116,505, 1,506,263 and 19,401,840 modules at n = 3–6; typed cores 1,434 / 206 / 9.
  - analyses/six_node_pattern S0: no BHAT6 node set of any class passes D1–D4.
- **Decision line.** It judges the claim that the six-node module is the smallest size at which two TFs and two miRNAs
  act on two genes, and that it is a distinct object from the three-node cores. It reports whether the higher-order
  Bhat networks add nodes and edges beyond the cores, or re-use them.

## O. Over-representation under NULL-B, NULL-C and NULL-L (written 2026-09-28 17:58)

- **What is counted.** All twelve networks, as the authors asked (their correction of 28 Sep): the Bhat three-, four-,
  five- and six-node patterns in the miRNA-FFL, TF-FFL and composite classes. For each, both instances and distinct
  node sets are counted, with the definitions of `make_ffl_networks.py`.
- **Program.** A copy of the S1 null program (`S1/s1_null6.c` of analyses/six_node_pattern, md5 4395d3e8…).
  - Counters for the three-, four- and five-node patterns are added inside its six-node count. The six-node counts
    are its existing counters.
  - The randomisation, both random streams and the order of draws are unchanged. Counting draws no random numbers,
    so replicate r is the same graph as S1's replicate r.
  - S1's other columns (the three-node D1-D4 counts and MODEL6) are kept for the replicate checks.
  - CENSUS6c is not computed (census_reps = 0). It is not a Bhat count, and it draws no random numbers.
- **Inputs.** S1's inputs, copied unchanged: `graph_nolegacy_null.txt` (md5 527cba76…), `graph_nolegacy_labels.txt`
  (f327921c…) and `node_names.txt` (2154653c…). They hold the census graph without the legacy miRNA-miRNA arcs,
  uncontracted (9,226 arcs), with one layer label per arc.
- **Settings, as S1.** 1,000 replicates per null, seed 20250908, 100 swap attempts per rewired item. Each replicate
  restarts from the observed graph.
  - Primary nulls: B, C and L.
  - Also LT, LG and LM, as in S1: NULL-L with only the TRRUST, STRING or 10 kb layer rewired. They show which layer
    carries any NULL-L excess.
- **Checks, before any statistic.**
  1. The input graph equals the census graph of `make_ffl_networks.py`: the same arcs, each in its own layer.
  2. The observed graph reproduces the instance and node-set counts of all twelve networks, and the instance lists
     equal those of the network builder.
  3. In every replicate of every null, S1's columns (all except CENSUS6c and the timing) equal S1's stored runs,
     row for row.
  4. In every replicate, the new three-node composite count equals S1's `bhat3_comp` column.
  - If check 3 fails, the new counts are taken from a separate pass with the same seed and replicates, as the
    authors asked, and the failure is reported.
- **Statistics (S1's summariser).**
  - Reported: observed count, null mean and SD, observed / mean, and Z.
  - p_upper = (#{x ≥ o} + 1) / (R + 1) and p_lower = (#{x ≤ o} + 1) / (R + 1).
- **Decision rule (S1's, fixed).** Applied to the instances and to the node sets of each network.
  - Over-represented: p_upper < 0.05 under both NULL-C and NULL-L.
  - Explained by the three-node cores: p_upper < 0.05 under NULL-C only.
  - Depleted: p_lower < 0.05 under both.
  - Otherwise: none of the three.
- **Known limit, stated before running: NULL-L cannot test the three-node networks.**
  - NULL-L fixes every analysed-network arc. Of the TRRUST layer, it rewires only TF->TF, Gene->TF and Gene->Gene arcs.
  - The three-node Bhat counts use only fixed arcs, so they are the same in every NULL-L replicate.
  - For the three-node networks, the rule can therefore give at most "explained by the three-node cores".
- **Paper and earlier values placed beside the results.**
  - The paper: motif excess of 8-11 % under reciprocity-preserving nulls.
  - analyses/six_node_pattern S1: six-node composite 2.35x NULL-C (p 0.002) and 3.69x NULL-L; all three six-node classes
    over-represented.
- **Decision line.** It judges whether higher-order FFL patterns are over-represented beyond what the three-node cores
  explain, as the paper's six-node conclusions assume, or only as the cores are.

## S. Circuit survival: each circuit against its three-node core (written 2026-09-28 18:13)

- **Script.** A generalised copy of `S8/s8_circuit_survival.R` (analyses/six_node_pattern). Only its input changes: it
  takes the instance listings of each network (section O's observed listings, in the S1 program's order). Data,
  scores, models and correction are unchanged.
- **Networks.** The nine four-, five- and six-node networks.
  - The three-node networks are not tested, because there the circuit is its core (agreed with the authors).
- **Circuits.** The distinct node sets of each network.
  - For a node set with several role assignments, the first one listed is used, as in S8.
  - Core = (T1, M1, G1). Extension = the other roles: G2 at four nodes; M2 and G2 at five; T2, M2 and G2 at six.
  - As in S8, the Bhat core is a triple of roles: at four to six nodes T1 -> G1 is not a compulsory edge (T1 -> G2 is).
  - Only circuits whose members are all measured are tested.
- **The one-gene extension (adaptation, fixed now).** At four nodes the extension is a single gene, G2. GSVA cannot
  score a one-member set (minSize 2), so its score is the gene's z-scaled expression. That is also what mean-z gives.
- **Data (as S8).** TCGA-BRCA primary tumours with gene and miRNA expression and survival: 1,065 samples and 20,742
  features. Endpoints are OS and PFI.
- **Primary analysis (as S8).**
  - Scores: GSVA (kcdf Gaussian, minSize 2, maxSize Inf), z-scaled.
  - Cox models adjusted for age and AJCC stage (complete cases):
    - core model: ~ core + age + stage;
    - circuit model: ~ circuit + age + stage;
    - nested model: ~ core + extension + age + stage.
  - ΔC = C(circuit model) − C(core model), Harrell's C.
  - LRT of the core model against the nested model.
  - BH across the circuits of each network, within method × endpoint.
- **Secondary analysis (S8's copy of the EA6 method).** Mean-z scores, univariable Cox, and the same ΔC and LRT.
- **Parallel computation.**
  - GSVA scores are computed in chunks of gene sets, in parallel processes. A test on this machine showed that GSVA
    2.2.1 gives identical scores for a gene set whatever other sets are scored with it (maximum absolute
    difference 0).
  - The Cox models are fitted in parallel over circuits. Each fit is deterministic.
- **Check.** The six-node composite network must reproduce S8 row for row: 2,383 circuits, 16 at q < 0.05 with GSVA
  for OS.
- **Also reported, as Q8 (analyses/six_node_followups).** For each network: the circuits at q < 0.05 in the primary analysis,
  and the nodes they share.
- **Paper and earlier values placed beside the results.**
  - The paper: "none of the 1,989 six-node census circuits beats its core after correction" (EA6 method).
  - S8: six-node composite, 16 of 2,383 circuits at q < 0.05 (GSVA, OS), 2 (GSVA, PFI), 0 with mean-z.
  - Q8: 12 of those 16 contain HIF1A or STAT3.
- **Decision line.** It judges the paper's statement that no six-node circuit beats its three-node core after
  correction.

## E. Enrichment: GSEA against the tumour vs normal signature, the network-background test, and ORA (written 2026-09-28 18:17)

- **Sets.** For each of the twelve networks: its protein-coding members (TF and Gene), plus TF-only and gene-only
  subsets, as `scripts/12_enrichment_and_metabolism/gsea/02_build_genesets.R` builds its FFL-class sets.
  - Sets with fewer than 5 members are dropped. `ALL_NETWORK_NODES` (the network's protein-coding nodes) is kept as
    the reference.
  - Derived sets, following `scripts/12_enrichment_and_metabolism/gsea/07_ffl_class_gsea.R`:
    - for each four-, five- and six-node network, its members that are in no three-node network (`_not_3node`);
    - the union of the three three-node networks;
    - the union of the nine higher-order networks, and the part of it in no three-node network;
    - the network's protein-coding nodes that are in no Bhat network.
- **E1 GSEA (settings of `scripts/12_enrichment_and_metabolism/gsea/03_run_gsea.R` and `07`).**
  - Ranked list: `results/v4/rank_tumour_vs_normal.csv`, the limma moderated t of primary tumour against normal
    (20,250 genes), read as `read_rank()` does. Only this list is used; it is the one Results 3.8 reports.
  - fgseaMultilevel with nPermSimple 10,000 and eps 0. `set.seed(20260908)` and MulticoreParam(8, RNGseed
    20260908).
  - Network sets: minSize 10 and maxSize 500, the FFL-class setting. Derived sets: minSize 5 and maxSize 600.
  - padj = BH within each collection (the network sets; the derived sets), as in the paper.
  - **Machinery check.** The paper's own FFL-class collection (`cache/v4/genesets.rds`) is run through the same code.
    Its ES must equal `results/v4/gsea_all_results.csv`, and its NES must agree within 0.02 (fgsea's NES depends
    on random sampling).
- **E2 Network-background test** (`07`, lines 77-95).
  - For each set: the mean |t| of its genes, against 10,000 size-matched draws without replacement from the 363
    network protein-coding nodes that have a t. p = (1 + #{draw >= observed}) / 10,001.
  - Also the signed version, and the transcriptome-wide version (draws from all 20,250 genes).
  - Sets with fewer than 5 genes, or at least as large as the background, are not tested (as in `07`).
  - `set.seed(20260908)` before the draws.
- **E3 ORA (settings of `scripts/12_enrichment_and_metabolism/over_representation/09_enrichment_ora.R`).**
  - Input: the protein-coding members of each of the twelve networks, mapped to Entrez with org.Hs.eg.db.
  - Universes: A, the genes tested for differential expression plus the network's protein-coding nodes; B, the
    network's protein-coding nodes.
  - Tests: enrichGO for BP, MF and CC; enrichKEGG; ReactomePA enrichPathway.
  - minGSSize 10, maxGSSize 500, BH correction. Significance is p.adjust < 0.05, with pvalueCutoff = qvalueCutoff = 1 so
    the full table is returned. Low-support flag: Count < 3. Seed 1234.
  - **KEGG is a download.** enrichKEGG fetches the current KEGG pathway data from the KEGG REST service when it runs.
    This is listed in DOWNLOADS.tsv and may differ from the release the paper used. GO and Reactome come from the
    installed packages (org.Hs.eg.db 3.21.0, GO.db 3.21.0, reactome.db 1.92.0; clusterProfiler 4.16.0).
  - **Machinery check.** The paper's 3-Comp and 6-node motif sets are run through the same code. Their GO and Reactome
    rows must equal `results/enrichment_all_motifs_FULL.csv`. KEGG differences are reported as release drift.
- **Paper values placed beside the results.**
  - GSEA: NES -1.85 to -1.91 (FDR < 1.2e-8) for the paper's FFL classes.
  - Network-background p: 3e-4 to 1.2e-3 for the three-node classes; 0.994 to 0.9995 for the sets of nodes added
    only at higher order.
  - ORA: the paper's significant-term counts per motif set (`results/enrichment_low_support_summary.csv`).
- **Decision line.** It judges the paper's reading in Results 3.8 on two points:
  - FFL nodes are enriched among genes lower in tumour (negative NES);
  - against the network background, the three-node nodes are more extreme, while nodes added only at higher order
    are not.

## X. Sign concordance in TCGA and CPTAC against matched nulls (written 2026-09-28 18:21)

- **Edges tested.** For each network, its edges with a predicted sign: the `sign` column of the network's edge file,
  following the published rule (miRNA -> target -1; TF edges +1 or -1 from the TRRUST or TransmiR mode).
  - STRING and 10 kb links carry no sign, so they are excluded, as in the paper.
  - Analysed-network edges are tested with the paper's stored results.
  - Census-layer TRRUST arcs have no stored correlation. They are computed here with the same method and reported
    separately; they are never pooled with the analysed-network edges.
- **Consistency check.** For every analysed-network edge in the networks, the file's sign must equal `predicted_sign`
  in `results/edge_correlation.csv`. Any difference is reported, and the file's sign is used.
- **X1 TCGA (`scripts/07_expression_validation/edge_correlations/08_expression_validation.R`).**
  - Stored results: rho, predicted sign and concordance per edge from `results/edge_correlation.csv`. Null pairs
    from `results/edge_correlation_null.csv`: 10 matched random pairs per edge, of the same node types, matched on
    expression decile within platform and on degree; real edges excluded. Their `edge_row` indexes
    `data/canonical_edges.tsv`.
  - New TRRUST-layer arcs:
    - Spearman rho on the same paired primary tumours, computed as the z-scored-rank product of `08`, with p from
      the t approximation.
    - Null pairs made with `08`'s `get_cand` procedure (the same pool, deciles, degrees from
      `data/canonical_edges.tsv`, and 10 pairs per edge), after `set.seed(1234)`.
  - Strata (`08`'s, restricted to the network's edges): all sign-annotated edges; miRNA -> target (all, and by
    tier: strong, weak, predicted only); TF -> target Activation and Repression; TF -> miRNA Activation and
    Repression. Also TF -> target Activation and Repression of the TRRUST layer, reported separately.
  - Tests (`08`): concordance rate against the stratum's null concordance rate,
    `binom.test(k, n, p0, "greater")`, the binomial test against 0.5, and `08`'s stratified permutation
    (10,000 draws).
- **X2 CPTAC (`scripts/07_expression_validation/protein_cptac_and_rppa/21_cptac_protein_validation.R`).**
  - Stored per-edge results: `results/multiomics/cptac_miRNA_mRNA_vs_protein.csv`, covering miRNA -> target and
    TF -> target edges whose partners are measured (n = 101).
  - Null pairs from `cptac_null_pairs.csv`. Their `edge_row` indexes the miRNA -> target (A) or TF -> target (B)
    rows of that file, in order.
  - Strata (`21`'s): miRNA -> target (all and by tier), TF -> target Activation, Repression and all signed; each at
    mRNA and at protein level.
  - Test (`21`): concordance against the null concordance, `binom.test(k, n, p0, "greater")`.
  - New TRRUST-layer arcs whose partners are measured are computed from `results/multiomics/cptac_bundle.rds` with
    `21`'s `rho_pair` and `make_null` (after `set.seed(1234)`), and reported separately.
- **Paper values placed beside the results.** The whole-network rows of `results/sign_concordance_summary.csv` and
  `results/multiomics/cptac_concordance_summary.csv`. For example:
  - TCGA TF -> target Activation: 0.711 against a null of 0.603 (p 5.8e-5);
  - TCGA miRNA -> target: 0.479 against 0.483 (p 0.71).
- **Decision line.** It judges the paper's reading: TF activation edges are sign-concordant above matched nulls,
  while miRNA -> target edges are not, and the Bhat networks' edges behave like the whole network's.

## P. Prioritisation with Bhat FFL participation (written 2026-09-28 18:25)

- **Script.** `scripts/06_node_prioritisation/01_prioritisation.R`, run unchanged except for two points, following the F1 precedent
  (analyses/six_node_pattern):
  - its output folder;
  - one inserted line after the FFL count is computed. The count `fps`, which enters T_topology as
    rn(ffl_cores), becomes the node's number of instances in one Bhat network: the `n_ffl` column of the
    network's node file. A node absent from the network gets 0.
  - Everything else is unchanged: the six other domains; degree (from `data/canonical_edges.tsv`, as in the paper's
    run) and betweenness; the role-breakdown display columns; and the rule of at least 3 domains.
- **Runs.** Twelve versions, one per network.
  - A thirteenth run, with nothing inserted, must reproduce the stored `results/v5/node_prioritisation_full.csv`
    and `top10_per_class.csv` byte for byte. This is "the paper's run".
- **Ranks.** Among the 583 rankable nodes (priority not NA), overall and within class, with ties given the lowest
  rank ('min'), as `scripts/06_node_prioritisation/40_node_compendium_assemble.py` re-ranks.
  - Rank change = the paper's rank − the version's rank, so a positive value means the node moved up.
- **Outputs for each version (the authors' request of 28 Sep).**
  1. A table of all 583 ranked nodes with these columns:
     - node, class, in_this_network, the FFL count used, and the T_topology score;
     - score and rank, overall and within class, for the paper and for this version;
     - rank change overall and within class;
     - the version's seven domain scores.
  2. Summaries over all nodes:
     - Spearman correlation of the version's scores with the paper's, overall and per class;
     - median and 90th percentile of the absolute rank change, overall and within class;
     - the nodes that move 50 or more places overall, with the direction.
  3. The Table 4 comparison, against the paper's thirty prioritised nodes (the ten highest-ranked of each class,
     `top10_per_class.csv`): for each class, the nodes entering and leaving the ten.
     - As an extra: the overall top 10 and top 20, entering and leaving.
- **Also checked (reported, not part of any decision).** Table S4's `rank_within_type` against the ranks among
  rankable nodes.
- **Paper values placed beside the results.** Table 4's thirty nodes.
- **Decision line.** It judges the paper's Table 4 selection: whether the thirty prioritised nodes hold when typed-core
  participation is replaced by Bhat FFL participation.

## Added analyses (requested by the authors on 28 Sep; each section is written before its run)

## A4. KnockTF support for the TF–TF arcs of the six-node patterns (added; written 2026-09-28 18:33)

- **Purpose.** It checks the signs that A1 takes from TRRUST. The four- and five-node patterns have no TF–TF arc.
- **Arcs.** The TF -> TF arcs that the six-node networks use as their compulsory T1–T2 link: 519 distinct arcs, 75 from
  the analysed network and 444 from the TRRUST layer.
  - Each arc is grouped by how it is used: as T2 -> T1 (TF-FFL and composite; the arc the paper's model lacks), as
    T1 -> T2 (miRNA-FFL and composite), or both.
  - Arc sets tested: all used arcs; the arcs used as T2 -> T1; the arcs used as T1 -> T2. Each set is also split by
    layer: analysed network, or TRRUST layer.
  - An arc's sign is the `sign` column of the network edge files, which is TRRUST's Activation (+1) or Repression
    (-1); blank means no sign.
- **Data and rule: P3b's, unchanged** (analyses/perturbation_tests, "P3b datasets" and "P3b measures").
  - Datasets: KnockTF 2.0 human datasets for the source TFs, from the cached `P3/p3b_datasets.tsv` and KnockTF's
    processed differential expression (flagged "processed"). All 120 source TFs are network TFs with TF -> target
    edges, so the cache covers them and nothing is downloaded.
  - Inclusion: at least 2 control and 2 treated samples; the TF's own Log2FC ≤ −0.74; duplicates removed.
  - Expressed genes: Mean_Control ≥ the dataset's median.
  - |log2FC| test: the median |Log2FC| of the TF's targets in the set, minus that of expressed genes that are not TRRUST
    targets of the TF in any mode. Two-sided Wilcoxon test. SE from 1,000 bootstrap resamples, seeded as in P3b (20250908
    and the dataset key). A dataset with fewer than 3 targets in the set gives no test.
  - Sign agreement: an Activation arc is expected to fall on knockdown, a Repression arc to rise; Log2FC = 0 counts as
    disagreement. Arcs without a sign enter the |log2FC| test only.
  - Pooling: REML for the |log2FC| difference (metafor). Pooled agreement = agreements / signed arcs, with a two-sided
    binomial test against 0.5.
  - Decision per arc set: SUPPORTED if the pooled difference is > 0 with a 95 % CI excluding 0, AND the agreement is
    > 50 % with binomial p < 0.05. Otherwise NOT SUPPORTED.
- **Per arc (reported).**
  - TRRUST sign; the number of included datasets of its source TF; the target's Log2FC in each, and whether the target
    is expressed there; the agreement count.
  - "no perturbation data" when the source TF has no included dataset, or the target is absent from all of them.
- **Rule for A1, fixed now.**
  - A1 takes each TF–TF arc's sign from TRRUST. An arc without a sign is run both as activation and as repression.
  - If A4 finds the T2 -> T1 set NOT SUPPORTED on sign agreement (agreement not > 50 % with p < 0.05), A1 runs the
    T2 -> T1 arc both ways in every configuration, not only where TRRUST gives no sign.
- **Decision line.** It judges whether the TF–TF arcs that the six-node patterns need are supported by perturbation
  data, with the rule the paper applied to its TF -> target edges (P3b: SUPPORTED).

## A5. Survival replication in METABRIC (added; written 2026-09-28 18:33)

- **Condition.** Item 5 gives circuits with q < 0.05 in TCGA (section S), so A5 applies.
- **Data check, done before writing this section.**
  - `data/metabric.rds` (cBioPortal brca_metabric) holds Illumina HT-12 mRNA expression: 20,384 genes × 1,980
    tumours, with clinical data.
  - It has no miRNA assay. `scripts/10_survival_and_clinical/cox_models/15_external_metabric.R` also states "METABRIC has NO miRNA assay".
  - METABRIC miRNA profiles (Dvinge et al. 2013) are controlled access (EGA), which the rules exclude.
- **Outcome: NOT POSSIBLE.** Every Bhat circuit's core contains a miRNA (M1), and the five- and six-node circuits
  contain two. So the core, circuit and extension scores of section S cannot be built in METABRIC, and nothing is run.
- **Proposed, not run.** A protein-coding-only version of the model, fitted in TCGA and METABRIC alike, as
  `15_external_metabric.R` does for modules. It runs only if the authors confirm it.

## Changes and decisions after sections were written

- **2026-09-28 18:37, A5.** The authors declined the proposed protein-coding-only replication, so A5 stays NOT POSSIBLE and nothing is
  run.
- **2026-09-28 19:59, A1 nested check: two outputs added, rule unchanged.** The authors asked for the nested result of
  every six-node behaviour that Results 3.7 reports, certified sustained oscillation included.
  - The sets that the nested run flags as sustained oscillation are certified by a1_certify.py's procedure
    (`scripts/a1_nested_certify.py`). They are compared set by set with the paper's certification: COMP_C2_toggle, 83
    of 108 (03f, stored); COMP_I1_negfeedback, 31 of 37 (analyses/six_node_pattern S2).
  - `scripts/a1_nested_table.py` writes one row per module and behaviour (analyses/six_node_pattern S2's list): the paper's
    prevalence (S2's `s2_prevalence_gain_loss.csv`, column pct), the nested prevalence, the absolute difference, and the
    number of sets whose flag differs.
  - The tolerance is unchanged: every per-set flag must be equal, now including the certification. That is 0 sets
    differing, and therefore a difference of 0.
  - This entry was written before the COMP_I1_negfeedback nested run finished and before the flag comparison ran.
    - One exception: testing the table code read the finished COMP_C2_toggle nested file (written 19:49). It had 0
      sets differing on every flag.
    - Nothing was changed after that.
- **2026-09-28 20:00, E, A2 and A3: ORA term overlap added (descriptive).** The authors asked, for each ORA rerun, for
  its overlap with the paper's significant terms.
  - `scripts/ora_overlap.py` writes `E3_ora_overlap.csv` in each enrichment folder (section E and the A2a, A3a and A3b
    reruns). It gives one row per network, universe and ontology, with:
    - the significant terms here;
    - the paper's significant terms for the matching motif set (`results/enrichment_all_motifs_FULL.csv`);
    - the terms shared by ID, and the Jaccard index;
    - for the reruns, the same counts against the same network in section E.
  - Matching sets: 3node_miRNA_FFL to 3-miR, 3node_TF_FFL to 3-TF, 3node_composite_FFL to 3-Comp; the four- to six-node
    networks to the paper's 4-node, 5-node and 6-node sets.
  - KEGG overlaps include release drift, because KEGG is downloaded at run time.
  - The ORA stays part of no decision line, as section E states. A2b is not rerun (section A2), so its values are
    section E's.
  - Written after section E's ORA was seen and before any rerun's ORA finished.
- **2026-09-28 20:17, A1: the verdict rule of the decision line.** Section A1 named what the decision line judges but
  gave no verdict rule. This one is written before any configuration module has finished; no chunk of the 41 modules
  had completed.
  - **Statement judged.** Results 3.7's six-node statement as analyses/six_node_pattern S2 left it. The paper's six-node
    module's memory AND pulse exceeds every smaller module (McNemar p < 0.05); its memory AND oscillation does not.
  - **Reproduced on a Bhat class.** A Bhat class (miRNA-FFL, TF-FFL, composite) reproduces a joint behaviour if at
    least one of its six-node configurations passes the claim guard (`a1_claim_guard.csv`). Passing the claim guard
    means exceeding the class's core, GG and GG+MM patterns with the same s_TM and s_TY = +1, each at McNemar p < 0.05.
  - **TF-FFL and composite.** Judged twice, with TF2 -> TF1 as activation and as repression. The two are reported side
    by side, and neither is chosen, because A4 found the arc not supported. The runs with the arc absent are for
    attribution only and do not count.
  - **Verdict, on memory AND pulse.**
    - supports: reproduced in every class, under each TF2 -> TF1 sign the class needs;
    - qualifies: reproduced in some of these cells but not all;
    - qualifies (not reproduced on the Bhat patterns): reproduced in none.
    - The paper's own module is covered by the nested check, so this line cannot contradict the paper's statement
      about that module.
  - **Memory AND oscillation.** Reported under the same rule, as S2 reports it. It is not part of the verdict.
  - **Also reported, not part of the verdict.** The comparisons with the paper's module, and the TF2 -> TF1
    attribution.
- **2026-09-28 21:03, A1 verdict rule: the authors accepted it, with one label change.**
  - **The logic is unchanged.** The five cells, the claim guard, activation and repression side by side, and memory AND
    oscillation reported but not judged. "supports" (all five cells) and "qualifies" (some cells) stay as written.
  - **The third verdict is renamed** from "qualifies (not reproduced on the Bhat patterns)" to "not reproduced on the
    Bhat patterns".
    - Its decision line adds that the paper's own module is reproduced exactly by the nested check, so the Results 3.7
      statement is unchanged.
  - **No module output had been read** when this change was made: no per-set file, certification or analysis output
    of the 41 modules.
    - The only lines seen were chunk run times from the run's progress log, used for the finish estimate at 20:30.
- **2026-09-29 03:47, A1: a cross-check added after the runs, and a detail added to the decision line.** No setting,
  rule or verdict changes. No change to the authors' rule arrived before 03:00.
  - **Cross-check.** One configuration, Bhat_comp_TM-_TY+ GG+MM+TT_T12+_T210_T2Y-, is by construction the model of
    analyses/six_node_pattern S2's variant S3b (TF2 -| G1). This was noticed on reading the certification log, where both
    have 226 flagged sets and 182 certified.
    - `scripts/a1_s3b_crosscheck.py` compares the two set by set: every flag, every other numeric column, and the
      certification.
    - It checks the sign handling on a module that the nested check does not cover.
  - **Decision-line detail.** The decision line now also gives:
    - for each cell, the configurations that pass the claim guard, with their prevalence and number of sets (in one
      cell the pass rests on a handful of sets, and a reader should see that);
    - the direction of the significant TF2 -> TF1 effects;
    - the count of configurations above and below the paper's module;
    - a scope sentence. The smaller patterns are the Bhat chain (core, GG, GG+MM) fixed in section A1, and the
      TF2-containing sub-modules of S2's factorial were not run.
- **2026-09-29 11:10, VALUES.md: rows added for traceability.** The authors asked that every A1 number they receive
  trace to VALUES.md. Two sets of rows are appended to the end of the A1 table, so every earlier A1 line keeps its
  number:
  - the McNemar tests behind the claim guard, for each six-node configuration against its class's core, GG and GG+MM
    patterns, for both joint behaviours (`dynamics/a1_mcnemar.csv`);
  - the per-module certification lines (`dynamics/a1_certification.log`).
  - No value, rule or verdict changes.
- **2026-09-29 11:24, A1: the instance coverage of each configuration, counted in full.** (This entry first carried the
  time 11:45 by mistake; the typo was corrected after the run of 11:24.)
  - `a1_configurations.csv` gives the greedy cover's newly-covered count per configuration: instances compatible with
    it and with none chosen before it. That count is a lower bound on the instances compatible with the
    configuration.
  - The 11:25 reply's caveat on minority configurations used these counts.
  - `scripts/a1_config_coverage.py` now counts, for each chosen six-node configuration:
    - the instances compatible with it (an unsigned arc matches either sign);
    - the instances with exactly its signs.
  - It reads the instances and signs as `a1_configurations.py` does. No setting, rule or verdict changes.

## A2. Cluster sensitivity (added; written 2026-09-28 18:41)

- **Question.** What, if anything, stands beyond one miRNA cluster? The authors asked for items 2-5 and 7 to be repeated
  for every network in two ways.
- **(a) miR-17~92 instances left out.**
  - Every instance that contains a member of miR-17~92 is dropped: hsa-miR-17, -18a, -19a, -20a, -19b and -92a, the
    network's six names for the cluster. Note that miR-19b and miR-92a are also produced by the paralogous
    miR-106a~363 cluster.
  - This is applied the same way to the observed count and to every replicate's count.
- **(b) Co-transcribed clusters merged.**
  - Clusters are the connected components of the 10 kb miRNA-pair layer: 480 pairs, 125 miRNAs, 35 clusters of 2 to
    28 members. Because the names miR-19b and miR-92a are shared, miR-17~92 and miR-106a~363 form one 10-member
    component.
  - The cluster map is fixed from the observed graph and used in every replicate.
  - **Merged graph.** Each cluster becomes one miRNA node that carries the union of its members' arcs; self-loops are
    dropped.
    - The three- and four-node patterns are counted on the merged graph.
    - Every miRNA pair of the five- and six-node patterns lies inside one cluster, so after the merge these patterns
      have no instance on the observed graph. This zero is reported, and not tested against the nulls.
  - **Cluster-level node sets.** Every instance's node set, with each miRNA replaced by its cluster, counted once. This
    is the main count for the five- and six-node patterns under (b), and is also given for three and four nodes.
- **Item 2 (nulls).**
  - C, B and L as in section O: the same program, seed, replicates and settings.
  - The program writes each replicate's arc list, and the counts of both variants are made from those lists by
    `make_ffl_networks.py`'s enumeration.
  - Check: on every replicate the unmodified counts of all twelve networks must equal the C program's.
  - Statistics and the decision rule are section O's.
- **Items 3, 4, 5 and 7 under (a).** They are rerun as in sections E, X, S and P on the networks rebuilt from the kept
  instances.
- **Items 3, 4, 5 and 7 under (b).** A cluster node stands for its member miRNAs, so:
  - Item 3 (GSEA, background test, ORA) uses protein-coding nodes only. These do not change, so the result is identical
    to section E by construction, and it is not rerun.
  - Item 4 (sign concordance): a cluster node's arcs are its members' arcs, so the tested edges do not change. The
    result is identical to section X by construction, and it is not rerun.
  - Item 5 (survival): the circuits are the distinct cluster-level node sets, each represented by the first instance
    listed. Scores, models and correction are section S's, with BH over these circuits.
  - Item 7 (prioritisation): a node's FFL count is the number of distinct cluster-level instances containing it; a miRNA
    counts those containing its cluster. The full rankings are as in section P.
- **Decision line.** It judges whether the section O-S, X and P results stand without the miR-17~92 instances, and
  after co-transcribed clusters are merged.

## A3. Evidence sensitivity (added; written 2026-09-28 18:41)

- **Item 2 (nulls), all twelve networks, on two restricted graphs.** Both are analyses/six_node_pattern S7's inputs, copied
  unchanged (md5 in MANIFEST.tsv).
  - **(a) Physical STRING only.** Gene-gene links kept only where STRING v12.0 physical-subnetwork score >= 0.900
    (S7(b): 391 gene-gene arcs kept, 443 dropped). Files: `graph_physical_ge_0.900.txt` and
    `labels_physical_ge_0.900.txt`.
  - **(b) Validated miRNA targets only.** miRNA -> target edges of the strong and weak tiers only; the 31
    exemplar-circuit edges removed; the census layers added as `build_graph` does. This is S7(c), the definition
    behind Note S2's validated count. Files: `null_validated.txt` and `labels_validated.txt`.
  - Nulls B, C and L with section O's program and settings: 1,000 replicates, seed 20250908, 100 swaps per item. The
    decision rule is section O's.
  - Check: the observed counts must equal S7's where S7 reports them.
    - (a): six-node composite 906 / 604, miRNA-FFL 3,633 / 2,264, TF-FFL 4,153 / 2,688.
    - (b): three-node composite 783; six-node composite 3,012 / 1,926, miRNA-FFL 10,142 / 6,352, TF-FFL 13,252 / 8,292.
- **Other items with a positive result, fixed now.**
  - Item 4: TF -> target Activation, and miRNA -> target in 7 networks.
  - Item 5: GSVA circuits at q < 0.05.
  - Item 3: GSEA, if its result (pending when this was written) is significant.
  - These are rerun on the networks rebuilt from each restricted graph's observed instances, as in sections E, X and S.
- **Decision line.** It judges whether the section O results, and the other positive results, stand when gene-gene
  links are physical only, and when miRNA targets are validated only.

## A1. Extended ODE model for the Bhat patterns (added; written 2026-09-28 18:54)

- **Base: the paper's module model and sweep, unchanged.** No project file is edited; the extension is a copy of the
  model file in this folder (`scripts/dyn_models_bhat.py`).
  - `scripts/11_dynamics/dyn_models.py`: equations, AND gate and state slots.
  - `dyn_sample.py`: 35 parameters and their ranges, scrambled Sobol, seed 20260908, 16,384 sets.
  - `03_higher_order_sweep.py` `run_module()`: initial conditions, RK4 with dt 0.01, operating input S = 2.0, and
    every behaviour classifier and threshold.
  - analyses/six_node_pattern S2: the chunked runner, and the limit-cycle certification (1,000 tau).
- **Pattern sizes follow the paper's layer factorial.**
  - Bhat three-node = core.
  - Four-node = GG: G2 added, with mutual protein stabilisation.
  - Five-node = GG+MM: a co-transcribed miRNA2 added, sharing the RISC pool.
  - Six-node = GG+MM+TT: TF2 added.
- **Classes: the compulsory miRNA-TF arcs of each Bhat pattern.**
  - Composite: TF1 -> miR1 and miR1 -| TF1, as in the paper's COMP families.
  - TF-FFL: TF1 -> miR1 only, as in the paper's non-composite family.
  - miRNA-FFL: miR1 -| TF1 only (s_TM = 0). The input drives the miRNA (input_node = MIR1), as in the paper's
    miRNA-led three-node topologies.
- **Added arcs: the only structural additions.** Each has the kinetic form and ranges of the model's existing TF-TF
  arc, TF1 -> TF2: a Hill activation or repression with threshold in [0.05, 1.5] (log) and Hill coefficient in
  [1, 4] (linear).
  1. **TF2 -> TF1** (s_T21: +1 activation, -1 repression, 0 absent). It enters TF1's production through the model's
     gate, as a second input to a gene enters elsewhere in the model: u_x1 = G(u_x1, f(x2; Kx21, nx21)). The six-node
     TF-FFL and composite patterns need it.
  2. **TF1 -> TF2 made switchable** (s_T12 = 0). When absent, TF2 is transcribed constitutively, as the model already
     does for an absent TF1 -> miRNA arc. The six-node TF-FFL pattern needs this, because TF1 -> TF2 is not
     compulsory in it.
  - Kx21 and nx21 are drawn from a separate scrambled Sobol sequence (2 dimensions, seed 20260909), with the ranges
    of Kx12 and nx12, and paired to the main sets by index. The 35-parameter design is not changed, so every existing
    parameter set stays identical.
  - Arcs of the paper's modules that are not compulsory in a Bhat pattern are kept as the paper has them, because Bhat
    patterns allow extra edges. These are TF1 -> G1/G2 at four to six nodes, the co-regulation of G2, and
    TF2 -> miRNA.
- **Signs, from the real instances** (the network files' `sign` column; miRNA arcs are -1 by construction).
  - Three nodes: TF1 -> miR1 (s_TM) and TF1 -> G1 (s_TY).
  - Four and five nodes: TF1 -> miR1 (s_TM) and TF1 -> G2 (s_TY). In the model, s_TY governs TF1's action on G1 and
    G2 together.
  - Six nodes:
    - TF1 -> miR1 (s_TM), for TF-FFL and composite;
    - TF1 -> TF2 (s_T12), for miRNA-FFL and composite;
    - TF2 -> TF1 (s_T21), for TF-FFL and composite;
    - TF2 -> G2 (s_T2Y). In the model, s_T2Y governs TF2's action on G1 and G2 together.
    - s_TY stays at the paper's +1, because no TF1 -> gene arc is compulsory at six nodes.
  - An unsigned arc is compatible with either sign.
- **Configurations: fully signed assignments of these arcs, chosen by a greedy cover.**
  - At each step, add the configuration that is compatible with the most instances not yet covered. Ties go to +1
    before -1, taking the arcs in the order listed.
  - Stop when at least 90 % of the network's instances are compatible with a chosen configuration.
  - TF2 -> TF1 is not part of the cover. By the rule fixed in A4 (NOT SUPPORTED), every chosen six-node TF-FFL or
    composite configuration is run with TF2 -> TF1 as activation, as repression, and absent (s_T21 = 0, for the
    attribution).
  - The list is written by `scripts/a1_configurations.py` to `dynamics/a1_configurations.csv`, with the number of
    modules and a time estimate. It is sent to the authors before the configuration runs start.
- **Nested check: first, and reported before any other A1 result.**
  - The extended code, with TF2 -> TF1 absent and TF1 -> TF2 present, reruns the paper's stored six-node modules:
    COMP_C2_toggle and COMP_I1_negfeedback, GG+MM+TT, 16,384 sets, seed 20260908.
  - Every per-set behaviour flag must equal the stored one.
- **Output, per pattern (class × size) and configuration.** Prevalence as % of the 16,384 sets, with 95 %
  Clopper-Pearson CI, of:
  - memory;
  - pulse;
  - oscillation (damped or certified sustained), each kind also given separately;
  - memory AND pulse;
  - memory AND oscillation.
- **Comparisons.** Exact two-sided McNemar tests, paired by parameter set, since every module uses the same 16,384
  sets.
  - Each six-node configuration against the paper's six-node module (COMP_C2_toggle GG+MM+TT).
  - Each six-node configuration against the same class's three-, four- and five-node configurations with the same
    s_TM and s_TY = +1.
  - TF2 -> TF1 attribution: each six-node configuration with the arc as activation, as repression, and absent, as in
    the layer factorial.
- **Claim guard, fixed.**
  - No "only the six-node pattern ..." statement unless its prevalence exceeds that of every smaller pattern with
    McNemar p < 0.05.
  - Pulses are never attributed to the negative loop.
- **Decision line.** It judges Results 3.7's six-node dynamics statements on the Bhat patterns.

## A1b. The claim guard against the full layer factorial (added; plan written 2026-09-29 11:20; nothing runs until the authors confirm)

- **Why.** A1's claim guard compared each six-node configuration only with the Bhat chain: core, GG and GG+MM. None of
  these contains TF2. Results 3.7 and Discussion 4.2 compare the paper's six-node module with every smaller module of
  S2's layer factorial: core, GG, MM, TT, GG+MM, GG+TT and MM+TT. A1b adds the four missing smaller modules, so that
  A1's pass can be set beside the paper's statement.
- **Configurations.** The nine that passed A1's guard for memory AND pulse (`dynamics/a1_claim_guard.csv`):
  - miRNA-FFL: Bhat_miR_TM0_TY+ with T12-, T21 absent, T2Y+.
  - TF-FFL, TF2 -> TF1 activation: Bhat_TF_TM+_TY+ with T12 absent, T21+, T2Y-.
  - TF-FFL, TF2 -> TF1 repression: Bhat_TF_TM-_TY+ with T12 absent, T21-, T2Y+.
  - Composite, activation:
    - Bhat_comp_TM+_TY+ with T12-, T21+, T2Y+;
    - Bhat_comp_TM-_TY+ with T12-, T21+, T2Y+.
  - Composite, repression:
    - Bhat_comp_TM+_TY+ with T12+, T21-, T2Y-;
    - Bhat_comp_TM+_TY+ with T12-, T21-, T2Y+;
    - Bhat_comp_TM-_TY+ with T12+, T21-, T2Y-;
    - Bhat_comp_TM-_TY+ with T12-, T21-, T2Y+.
  - Only these nine are run. The extended guard contains the Bhat chain, so a configuration that failed A1's guard
    cannot pass it.
- **Modules to add.** The layers are those of S2's factorial. The model, parameter sample, seeds and classifiers are
  A1's.
  - **TT, GG+TT and MM+TT** (TT = core + TF2), with each configuration's own TF2 arcs (T12, T21, T2Y). TF2 -> miRNA is
    +1, as in A1. That is 9 × 3 = 27 modules.
  - **MM** (core + miRNA2). MM contains no TF2, so it depends only on the class and the TF1 signs: one MM per family.
    - Three families reuse S2's stored MM modules: Bhat_comp_TM-_TY+ (COMP_C2_toggle MM), Bhat_comp_TM+_TY+
      (COMP_I1_negfeedback MM) and Bhat_TF_TM+_TY+ (I1_miRNA_FFL MM). These are the same models whose core, GG and
      GG+MM A1 used.
    - New: Bhat_miR_TM0_TY+ MM and Bhat_TF_TM-_TY+ MM, 2 modules.
  - **Unchanged.** Core, GG, GG+MM and the six-node modules are A1's, run or stored.
  - **Not included: the arc-absent versions** (TF2 -> TF1 absent) of the new TF2-containing modules. The verdict does
    not need them, and A1 already reports the TF2 -> TF1 attribution at six nodes.
    - Optional (a): +24 modules, if the authors want the attribution at every TF2-containing layer.
  - **Optional (b): +3 modules.** The TF2-containing modules of Bhat_miR_TM0_TY+ with T12+, T2Y+. This configuration had
    A1's highest memory AND oscillation (3.27 %) but failed its guard for memory AND pulse. It is needed only if memory
    AND oscillation is to be reported for every configuration that passed A1's guard for it.
- **Nested check first, as in A1. It gates the runs.**
  - The new layers are built through `a1_run.py`'s topology path, extended with MM, TT, GG+TT and MM+TT in
    `scripts/a1b_run.py`.
  - The check reruns S2's stored COMP_I1_negfeedback TT, GG+TT, MM and MM+TT through that path, as Bhat_comp_TM+_TY+
    with T12+, T21 absent, T2Y+ (the same model): 4 modules.
  - Every per-set flag must equal the stored one.
  - The certification must match S2's set by set: TT 21 of 26, GG+TT 20 of 25, MM+TT 34 of 38, and no set flagged in
    MM.
  - The 29 modules start only if the check passes.
- **Settings are A1's.**
  - The 16,384 sets of seed 20260908, with Kx21 and nx21 from seed 20260909.
  - RK4, dt 0.01, S = 2.0.
  - The classifiers and thresholds of `run_module()`.
  - Certification of sustained oscillation by `a1_certify.py`'s procedure.
- **Comparisons.** Exact two-sided McNemar tests, paired by parameter set. A configuration "exceeds every smaller
  pattern" only if its prevalence is higher than each of core, GG, MM, TT, GG+MM, GG+TT and MM+TT of the same family
  and signs, each at p < 0.05.
- **Verdict, fixed now before any run.** It is on memory AND pulse, cell by cell as in A1.
  - The five cells: miRNA-FFL; TF-FFL with TF2 -> TF1 as activation and as repression; composite with TF2 -> TF1 as
    activation and as repression.
  - A cell passes if at least one of its configurations passes the extended guard.
  - The labels:
    - supports: all five cells pass;
    - qualifies: some cells pass;
    - not reproduced on the Bhat patterns: no cell passes. This decision line adds that the paper's own module is
      reproduced exactly by A1's nested check and passed S2's full-factorial guard for memory AND pulse, so the Results
      3.7 statement about the paper's module is unchanged.
  - For each cell, the decision line names the passing configurations, with their prevalence and number of sets. It
    also gives the instances their signs cover (`dynamics/a1_configurations.csv`) against the network's total, so a
    pass on a minority of instances is visible.
- **Reported, not judged.**
  - Memory AND oscillation, under the same guard, for the configurations that are run.
  - The comparisons with the paper's module.
- **Claim guard and wording.**
  - No "only ..." statement unless the extended guard passes.
  - Pulses are never attributed to the negative loop.
  - TF2 -> TF1 activation and repression are reported side by side, and neither is preferred (A4).
  - Bhat objects are called patterns, configurations or instances.
- **Estimate, with 12 workers.** Chunk times come from A1's run and S2's log, scaled to A1's load. Per chunk of 2,048
  sets: TT about 450 s, MM about 450 s, MM+TT about 700 s, GG+TT about 1,100 s.

  | step | modules | chunks | CPU-h | wall time |
  |---|---|---|---|---|
  | nested check | 4 | 32 | 6 | about 45 min |
  | main | 29 | 232 | 47 | about 4-4.5 h |
  | certification and analysis | | | | about 30 min |
  | **total** | | | | **about 5.5-6 h from the start** |
  | optional (a) | +24 | | +40 | +3.3 h |
  | optional (b) | +3 | | +5 | +0.4 h |

  - Storage: about 2.1 MB per module, so about 70 MB for the main runs and the nested check. The zip grows to about
    215 MB (8 parts).
- **Scripts.** They are written and tested on stand-in data before the start; nothing runs on the real modules before
  the authors confirm.
  - `a1b_run.py`: the extended layers.
  - `a1b_modules.py`: the module list.
  - `a1b_nested_check.py`: the nested check (named `a1b_nested_compare.py` in the 11:20 plan; corrected 16:34).
  - Certification: `a1_certify.py` with the extended layers.
  - `a1b_analyse.py`: the extended guard, the verdict and the coverage.
  - make_values.py: section A1b.
- **Decision line.** It judges Results 3.7's six-node statement on the Bhat patterns, against the full layer
  factorial.
- **2026-09-29 11:25: the authors confirmed this plan, with option (b) and without option (a).**
  - **Main runs: 32 modules.** The 29 above, plus option (b)'s 3: TT, GG+TT and MM+TT of Bhat_miR_TM0_TY+ with T12+,
    T21 absent, T2Y+.
    - The option (b) configuration is reported for memory AND oscillation only. It is not part of the verdict, because
      it failed A1's guard for memory AND pulse.
  - **One reporting change, made before any A1b run.** The per-cell instance coverage uses
    `dynamics/a1_config_coverage.csv`, not the greedy cover's newly-covered count (Changes and decisions, 11:24). Two
    counts are given: instances compatible with the configuration (an unsigned arc matches either sign), and
    instances with exactly its signs.
  - **Everything else is as written above.** The nested check gates the runs, and the verdict rule is unchanged.
  - **Handover files**, next to the zip and not in it: a status file and a reply note.
