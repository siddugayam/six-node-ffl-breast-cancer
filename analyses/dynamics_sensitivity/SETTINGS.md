# SETTINGS: dynamics of the higher-order modules (Fig. 5 / Table 4)

The design of the analysis.

## Inputs (public repository; FFL_ROOT = a checkout of it)

- **Per-set behaviour flags for all 16,384 Sobol points:** `results/systems_and_dynamics/dynamics_higher_order_persetset.csv.gz`
  (families I1_miRNA_FFL and COMP_C2_toggle) and `dynamics_higher_order_persetset_compI1.csv.gz` (COMP_I1_negfeedback).
  - Modules: n3 (the three-node core), n4, n5, n6 and n4d.
- **Stored summaries:** `dynamics_higher_order_gain_loss.csv` and `_compI1.csv` (Table 4), and
  `dynamics_higher_order_fractions*.csv`.
- **Limit-cycle certification:**
  - `results/systems_and_dynamics/dynamics_n6_limit_cycle_certified.csv` (COMP_C2_toggle n6);
  - `analyses/six_node_pattern/S2/s2_limit_cycle_certification.csv` (COMP_I1 n6 = GG+MM+TT).
- **Code, imported unchanged:** `scripts/11_dynamics/dyn_sample.py`, `dyn_models.py`, `03_higher_order_sweep.py` (whose
  `run_module` is identical to that of `03d_higher_order_compI1.py`).
  - The parameter matrix is regenerated with `dyn_sample.sample(16384, seed=20260908)`.
  - It is checked against `dynamics_3node_sobol_params.npz` (columns 1–1,024, which must equal points 0–1,023) before
    use.

## (a) Exact 95 % Clopper–Pearson intervals

- **Rows:** every family (3) × module (n3, n4, n5, n6, n4d) × binary behaviour of the sweep:
  - is_memory, is_pulse, is_sustained_osc, is_damped_osc, is_bistable, bistability_unresolved, is_ultrasensitive,
    is_very_ultrasensitive, is_noise_rejecting, nonmonotone_dose.
- **Count:** k = the number of points whose flag is 1, out of n = 16,384.
  - It is checked against the stored percentages, and for n4–n4d against k = n_gain + n_both, and for n3 against
    k = n_loss + n_both.
- **Interval:** exact two-sided 95 % Clopper–Pearson: Beta(0.025; k, n − k + 1) and Beta(0.975; k + 1, n − k), with 0
  and 1 at the bounds.
- **Also given:** the certified limit-cycle counts (COMP_C2_toggle n6; COMP_I1 n6), with their intervals.

## (b) PRCC

- **Modules:** n3, n4, n5, n6 and n4d of the three families.
- **Behaviours (primary):** the binary flags of
  - memory (is_memory),
  - pulse (is_pulse),
  - oscillation (is_sustained_osc, the uncertified flag, as in Table 4).
- **Estimable only if** the module has at least 10 points with the behaviour and 10 without; otherwise the row is
  "not estimable".
- **Parameters:** those that act in the module, found empirically.
  - A parameter acts if multiplying it by 2 changes the right-hand side (`dyn_models.make_rhs`) at any of 64 random state
    vectors × 16 sampled points (numpy seed 1).
  - Parameters that do not act are excluded and listed.
- **PRCC.** Each parameter and the output are rank-transformed (mid-ranks for ties).
  - The parameter's rank and the output's rank are each regressed (ordinary least squares, with intercept) on the ranks
    of the other acting parameters.
  - PRCC = the Pearson correlation of the two residual vectors.
  - p from t = r √((N − 2 − m) / (1 − r²)) with N − 2 − m degrees of freedom, m = the number of other parameters. BH q
    is taken within module × behaviour.
- **Top 5:** the five parameters with the largest |PRCC| for each module × behaviour.
- **Secondary (continuous proxies; not used to choose parameters for (c)):**
  - memory: T_half_OFF, where "not reached in 50 τ" is ranked highest;
  - pulse: overshoot_ratio (points where it is undefined are dropped);
  - oscillation: rel_peak_to_peak.

## (c) One-at-a-time excursions

- **Modules:** n3, n4, n5 and n6 of the three families. n4d is not run.
- **Parameters per module:** the union over its estimable behaviours of the top-5 parameters from the primary PRCC.
- **Run.** For each such parameter and factor f ∈ {0.5, 2}, every one of the 16,384 sampled values of that parameter is
  multiplied by f; all other parameters keep their sampled values. `run_module` of `03_higher_order_sweep.py` is run
  unchanged.
  - The sample is injected in chunks of 2,048 points, as `analyses/six_node_pattern/S2/s2_run_modules.py` does.
- **Reported:** for every binary behaviour,
  - the excursion frequency with its Clopper–Pearson interval;
  - the baseline frequency;
  - the matched change (points gained and lost against the same baseline point, with an exact McNemar p).
- **Range note.** ×2 takes every parameter past the upper end of its sampled range, and ×0.5 below the lower end. theta
  and kappa start at 0, so ×0.5 keeps them within [0, 5]; ×2 extends them to [0, 20].
- **Gate, before any excursion run.** An unscaled run of points 0–2,047 of COMP_C2_toggle n6 must reproduce the stored
  flags of those points exactly (all 10 binary flags). If it does not, nothing is run and the failure is reported.
- **Order of the runs:**
  1. n6 of COMP_C2_toggle, COMP_I1_negfeedback, I1_miRNA_FFL;
  2. n5 of the same three;
  3. n3 of the same three;
  4. n4 of the same three.
- **Compute.** 10 worker processes, each single-threaded (OMP/OpenBLAS/MKL threads = 1).

## (d) Parameter table

One row per parameter of `dyn_sample.RANGES`, with:
- the symbol, quantity and meaning (`dynamics_parameter_definitions.csv`);
- low, high and scale; the distribution (log-uniform for "log", uniform for "lin");
- the justification;
- the modules in which it acts (from (b));
- the sampler: a scrambled Sobol sequence, d = 35, seed 20260908, the first 16,384 points via `random_base2(14)`;
- the module integration settings of `03_higher_order_sweep.py`: fixed-step RK4, dt 0.01, the input levels and run
  lengths.

## Outputs

- `behaviour_frequencies_clopper_pearson.csv`
- `parameters_acting_by_module.csv`
- `prcc_binary.csv`
- `prcc_continuous.csv`
- `prcc_top5.csv`
- `excursion_runs/`: per-run per-point flags, not zipped if large
- `excursion_frequencies.csv`
- `parameter_table.csv`

## Acting parameters

- The state vectors used to find the acting parameters are set to 0 in the module's inactive states, as in every
  integration of the model; a parameter acts if doubling it changes the right-hand side.
