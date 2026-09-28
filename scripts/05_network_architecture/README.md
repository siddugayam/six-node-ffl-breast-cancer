# Network architecture, controllability and information flow (Supplementary Note S4)

`netlib.py` (shared library) → `01_architecture.py`, `01b_architecture_nulls.py`, `02_powerlaw_csn.R` (degree
distribution) → `03_controllability.py` → `04_information_flow.py` → `05_perturbation.py`, `06_knockout_landscape.py`,
`09_combination_stats.py` (with `perturblib.py`) → `07_dynamics.py`, `07b_dynamics_nulls.py` (linear response of the
network) → `10_figures.py`, `11_master_summary.py`. `12`–`16`: controls and cross-checks; `96_powerlaw_figure.R`: the
degree-distribution figure. Re-runs without the exemplar edges: `analyses/analysed_network_reruns/v3/`.
