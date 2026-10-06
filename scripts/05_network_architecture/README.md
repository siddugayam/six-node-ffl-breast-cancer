# Network architecture, controllability and information flow (Supplementary Note S4)

`netlib.py` (shared library) → `01_architecture.py`, `01b_architecture_nulls.py`, `02_powerlaw_csn.R` (degree
distribution) → `03_controllability.py` → `04_information_flow.py` → `05_perturbation.py`, `06_knockout_landscape.py`,
`09_combination_stats.py` (with `perturblib.py`) → `07_dynamics.py`, `07b_dynamics_nulls.py` (linear response of the
network) → `10_figures.py`, `11_master_summary.py`. `12`–`16`: controls and cross-checks; `96_powerlaw_figure.R`: the
degree-distribution figure. The order above is the run order (`09` runs before `07`). Re-runs without the exemplar
edges: `analyses/analysed_network_reruns/v3/`.

## Scripts in this folder

| Script | What it does |
|---|---|
| `01_architecture.py` | (A) ARCHITECTURE of the canonical miRNA-TF-gene network. |
| `01b_architecture_nulls.py` | (A cont.) Are the architectural statistics anything other than a consequence of the degree sequence? |
| `02_powerlaw_csn.R` | (A) Is the degree distribution scale-free? |
| `03_controllability.py` | (B) STRUCTURAL CONTROLLABILITY of the canonical network. |
| `04_information_flow.py` | (C) INFORMATION FLOW through the regulatory network. |
| `05_perturbation.py` | (D) PERTURBATION AND CONTROL: in-silico single- and double-node knockouts. |
| `06_knockout_landscape.py` | (D) Complete single- and double-knockout landscape of the canonical network. |
| `07_dynamics.py` | (E) DYNAMIC RESPONSE OF THE REAL NETWORK. |
| `07b_dynamics_nulls.py` | (E cont.) Nulls and self-consistency checks for the dynamic-response analysis. |
| `09_combination_stats.py` | (D cont.) Does multi-target perturbation actually beat single-target inhibition here? |
| `10_figures.py` | Figures for the systems-pharmacology / network-control analysis. |
| `11_master_summary.py` | Consolidate every headline number produced by the v3 systems-pharmacology analysis into one table. |
| `12_controllability_degree_control.py` | (B cont.) Is the driver-node depletion among FFL hubs anything more than a degree effect? |
| `13_matching_crosscheck.py` | (B cont.) Independent cross-check of the maximum-matching / driver-node computation. |
| `14_regulator_subnetwork.py` | (A cont.) Is the "hierarchical regulatory structure" real, or an artefact of how the network was assembled? |
| `15_igraph_crosscheck.R` | (A cont.) Independent cross-check, in R/igraph, of the topology numbers that the Python analysis (01_architecture.py, 14_regulator_subnetwork.py) reports with networkx. |
| `16_flow_vs_degree_sinkcontrol.py` | (C cont.) The flow-vs-degree comparison, controlled for the structural sinks. |
| `96_powerlaw_figure.R` | Total-degree distribution with discrete power-law fits (xmin 33 and 24) and a log-normal fit (Supplementary Note S4); writes Fig_v3_A2_powerlaw_xmin_sensitivity.pdf and .png. |
| `netlib.py` | Shared network utilities for the v3 systems-pharmacology / network-control analysis. |
| `perturblib.py` | Fast knockout-damage kernel shared by the single- and double-knockout analyses. |
