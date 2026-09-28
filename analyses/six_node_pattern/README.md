# The six-node composite pattern (Supplementary Note S2)

Whether the six-node composite FFL pattern is over-represented, what it does in dynamical models and whether it carries
clinical signal. The folder names are the items of the analysis plan.

| Folder | What it holds |
|---|---|
| `S0/` | Reconciles the four six-node objects of the project (census modules, the composite pattern, the edge-addition sets and the modelled architecture) |
| `S1/` | Over-representation of six-node patterns against randomised networks (NULL-B, NULL-C and the core-conditional NULL-L); `jobs.txt` lists the runs |
| `S2/` | Dynamics of the layer-factorial module variants and the TF–TF sign controls; per-parameter-set results in `perset/` |
| `S4/` | How oscillatory the damped oscillations are |
| `S5/` | The signed real six-node circuits, two literature positive controls and MCODE without miRNA–miRNA links (`s5e_mcode/`) |
| `S6/` | Models of the fully signed six-node configurations found in the network |
| `S7/` | Robustness to the STRING and evidence filters: physical STRING links and validated-only graphs |
| `S8/` | Survival test of each six-node circuit against its three-node core. `six_node_pattern_survival_all_circuits.csv` has one row per circuit, scoring method and endpoint, with cohort-level counts only |
| `E/` | `E5`: re-typed graphs, coherence and the miRNA–TF pair filter; `E6`: the two miR-29a–*COL1A1* correlations; `E7`: mediation intervals across the 42 stromal estimates |
| `F/` | Re-runs on the analysed network: prioritisation (`F1`), survival hub list (`F2`), compendium topology (`F3`), enrichment gene sets from the exact census (`F4`), motif significance (`F5`) and the target rule of the named FFL classes (`f6_target_rule.py`) |
| `G/` | Benjamini–Hochberg correction over the eleven reported miRNAs (`g1`, Table S8) and the degree-matched driver test (`g2`) |

The C programs (`S1/s1_null6.c`, `S7/ffl_census_composition.c`, `S7/v2_null.c`, `F/F4/f4_node_union.c`) are given as
source; the repository README gives their compile commands.
