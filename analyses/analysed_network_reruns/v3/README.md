# V3

## Scripts in this folder

| Script | What it does |
|---|---|
| `02_powerlaw_csn_full.R` | (A) Is the degree distribution scale-free? |
| `02_powerlaw_csn_nolegacy.R` | (A) Is the degree distribution scale-free? |
| `exact_control_roles.py` | Exact (matching-invariant) structural-controllability roles, replacing the single-configuration driver flag behind Note S5's hub comparison. |
| `fig_s5_information_flow.py` | Fig. S5 (Fig_v3_C_information_flow) only: the Fig 3 block of scripts/05_network_architecture/10_figures.py, copied verbatim, with RES/FIG pointed at out_<mode>_hash0/ (the other v3 figures need analyses not re-run here). usage: … |
| `run_v3.py` | Re-run the v3 architecture / controllability / information-flow scripts (Supplementary Note S5) unchanged, with only two things swapped in scripts/05_network_architecture/netlib.py at run time: RES/FIG/LOG -> out_<mode>/ in this … |
