# SETTINGS: tier-stratified miRNA–target sign concordance after stromal adjustment

The design of the analysis.

## Design

- **Base design.** The five-draw within-network matched null of
  `scripts/07_expression_validation/edge_correlations/04_caf_adjusted_vs_matched_null.R`, copied unchanged.
  - Samples: TCGA-BRCA primary tumours with both assays (n = 1,066).
  - Spearman correlation as rank-z dot products.
  - Null: for every edge, 10 pairs of network nodes of the same kind (gene or miRNA) in the same mean-expression
    decile, excluding self-pairs and real edges. Five draws, seeds 501–505.
  - Concordance = the fraction of edges whose correlation sign equals the predicted sign (−1 for miRNA -> target).
  - Test: exact two-sided binomial against each draw's null concordance; p is the mean over the five draws, and the
    largest p is also given.
- **Control.** The script first reproduces the published table
  `results/expression_validation/caf_adjusted_vs_network_null.csv` (raw, CAF-A and CAF-B columns), within a relative
  tolerance of 1e-9. Nothing else is reported unless it does.

## Adjustment

- **Method.** First-order partial Spearman correlation, the formula of 04's `mk()`:
  - (r − a·b) / √((1 − a²)(1 − b²)), where a and b are the rank correlations of source and target with the covariate;
  - applied identically to the real edges and to the null pairs of the same five draws.
- **Covariates, each alone.** The 45 stromal estimates of the mediation analysis
  (`results/systems_and_dynamics/deconv_fibroblast_estimates_extended.csv`). Among them:
  - ABSOLUTE purity, as `ABSOLUTE_nontumour` = 1 − purity; 1,023 of the 1,066 tumours;
  - the ESTIMATE stromal score (`ESTIMATE_StromalScore`);
  - the fibroblast and CAF estimates (xCell, MCP-counter, EPIC, ConsensusTME, CIBERSORT-type and other deconvolutions;
    signature scores; consensus scores; pathology and EpiDISH methylation estimates).
- **Missing values.** For a covariate with missing values (ABSOLUTE 43, a pathology score 1, the three EpiDISH columns
  299), the ranks are recomputed on the tumours that have it, for the edges and the null pairs alike. The unadjusted
  concordance on the same tumours is reported beside it.

## Edge sets

- miRNA -> target: all, strong, weak and predicted_only (the tiers of `data/network/edge_evidence_tier.tsv`).
- The TF sets of 04 are kept in the tables for completeness.

## Leave-stroma-enriched-targets-out sensitivity

- **The Farmer-correlated gene set.** The rule of `scripts/10_survival_and_clinical/reactive_stroma_and_neoadjuvant/43_farmer_module_overlap.R` (the rho > 0.5 fraction, lines
  145–150): genes whose
  Spearman rho with the Farmer stroma-related ssGSEA score is above 0.5, from
  `results/prioritisation_and_stroma/farmer_pergene_rho.csv`. That is 698 genes; 67 of them are network miRNA targets.
- **The sensitivity.** Every miRNA -> target edge whose target is in the set is dropped. The unadjusted and all
  adjusted concordances are then recomputed on the remaining edges, with the same null pairs of those edges.

## Reported

- **Per tier and edge set:**
  - unadjusted: concordance, null and p (mean and maximum over the 5 draws);
  - adjusted for ABSOLUTE alone and for ESTIMATE alone: the same, with the n of tumours;
  - for the 45 estimates: the median and range of the adjusted concordance and of the mean p, and the number of
    estimates with mean p < 0.05;
  - every per-estimate row in the full table.
- **Direction.** Concordance above the null = agreement with repression. The p is two-sided, as in 04.

No decision rule is attached; the numbers are reported whichever way they fall.

## Software

R, with matrixStats. Seeds 501–505 come from 04 and are unchanged. No other random numbers.

## Outputs

- `concordance_adjusted_by_estimate.csv`: every estimate × edge set × all or leave-out.
- `concordance_summary.csv`
- `control_reproduction.csv`

No sample-level values are written.

