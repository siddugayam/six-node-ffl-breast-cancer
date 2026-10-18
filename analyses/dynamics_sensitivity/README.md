# Dynamics of the higher-order modules: intervals, PRCC, excursions and parameter table

This folder holds four analyses of the 16,384-point Sobol sweep behind Table 4 and Fig. 5 (families I1_miRNA_FFL,
COMP_C2_toggle and COMP_I1_negfeedback; modules n3 to n6 and n4d):
- **(a)** exact Clopper–Pearson intervals for every behaviour frequency;
- **(b)** PRCC of every acting parameter against memory, pulse and oscillation;
- **(c)** one-at-a-time ×0.5 and ×2 excursions of the top-5 parameters, re-running the sweep;
- **(d)** the table of parameter ranges and sampling distributions.

The design is in `SETTINGS.md`. The values, with file, line and column, are in `VALUES.md`.

| file | what it is |
|---|---|
| `behaviour_intervals_and_prcc.py` | (a), (b), (d), and the excursion plan |
| `run_excursions.py` | (c): runs `run_module()` of `scripts/11_dynamics/03_higher_order_sweep.py` unchanged, with the scaled parameter; first checks that an unscaled chunk reproduces the stored flags (`excursion_runs/gate.txt`) |
| `summarise_excursions.py` | (c): frequencies, intervals and matched changes against the published sweep |
| `write_values.py` | Writes `VALUES.md` |
| `behaviour_frequencies_clopper_pearson.csv` | (a) |
| `prcc_binary.csv`, `prcc_top5.csv`, `prcc_continuous.csv`, `parameters_acting_by_module.csv` | (b) |
| `excursion_plan.csv`, `excursion_frequencies.csv`, `excursion_range_by_module.csv`, `excursion_runs/` | (c), with per-point flags of every excursion run; the range table is descriptive |
| `parameter_table.csv`, `integration_settings.csv` | (d) |

**Run:** set `FFL_ROOT=<repository checkout>`, then:
1. `python3 behaviour_intervals_and_prcc.py`
2. `OMP_NUM_THREADS=1 python3 run_excursions.py 10`
3. `python3 summarise_excursions.py`
4. `python3 write_values.py`

## Scripts in this folder

| Script | What it does |
|---|---|
| `behaviour_intervals_and_prcc.py` | Behaviour frequencies of the higher-order module sweep (Table 4, Fig. 5): exact Clopper-Pearson intervals, PRCC of every acting parameter against memory, pulse and oscillation, the parameters chosen for the one-at-a-time … |
| `run_excursions.py` | One-at-a-time parameter excursions of the higher-order module sweep. |
| `summarise_excursions.py` | Behaviour frequencies of the one-at-a-time parameter excursions against the published sweep. |
| `write_values.py` | Writes VALUES.md for the dynamics sensitivity analysis from the tables in this folder: every number with its file, line (header = line 1) and column. usage (from this folder): python3 write_values.py |
