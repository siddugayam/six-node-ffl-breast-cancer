# Dynamics of three-node and higher-order modules (Methods 2.5; Results 3.7)

Shared model code: `dyn_models.py`, `dyn_sample.py`, `dyn_metrics.py`, `dyn_style.py`. Three-node characterisation
(`01`, `02*`) → higher-order sweep and checks (`03*`) → miR-29 circuit fit (`04*`) → composite module (`06*`, `09b`,
`09c`) → consolidation (`07_consolidate.py`) and figures (`05_figures.py`, `10b_*`). The numbers group the steps;
`99*_driver.sh` run them in their exact order.
Controls and variants: `analyses/dynamics_controls/`, `analyses/six_node_pattern/S2`, `S4`, `S6`.

## Scripts in this folder

| Script | What it does |
|---|---|
| `01_integrator_validation.py` | Cross-check the vectorised fixed-step RK4 used everywhere else against scipy.integrate.solve_ivp (LSODA, rtol=1e-9, atol=1e-12) on every topology, and against itself at half the step size. |
| `02_characterise_3node.py` | Simulate every 3-node topology (Alon transcriptional references, miRNA-mediated cores, composite cores, and the regulated-cascade control) at (i) a stated reference parameter set and (ii) a matched scrambled-Sobol ensemble, and … |
| `02b_summary_table.py` | Build the paper's 3-node comparison table from the saved per-parameter-set metrics (dynamics_3node_metrics_ensemble.csv). |
| `02c_characterise_3node_N16384.py` | Re-run of 02_characterise_3node.py with the ensemble enlarged from 1,024 to 16,384 scrambled-Sobol parameter sets per topology, matching the size used in the higher-order sweep (03_higher_order_sweep.py) so that Parts B and C of … |
| `02d_adaptive_solver_crosscheck.py` | INDEPENDENT NUMERICAL VERIFICATION of the dynamical battery. |
| `02e_paired_tests_and_textbook_validation.py` | Paired tests of each three-node topology against the regulated cascade at the same Sobol parameter sets, and validation of the model against textbook FFL behaviour. |
| `02f_pulse_longwindow.py` | FIX FOR A REPORTING ARTEFACT IN THE PULSE STATISTICS. |
| `02g_fcd_recomputed.py` | CORRECTED FOLD-CHANGE-DETECTION TEST. |
| `02h_stochastic_noise_buffering.py` | DOES THE miRNA-MEDIATED INCOHERENT FFL BUFFER NOISE? |
| `02i_paper_table1.py` | A compact, print-ready version of the 3-node comparison table: the columns a journal table would actually carry, formatted as text, with the full 46-column version left in dynamics_3node_comparison_table.csv for the supplement. |
| `02j_noise_at_matched_mean.py` | NOISE COMPARED AT MATCHED MEAN OUTPUT. |
| `03_higher_order_sweep.py` | DOES ORDER ADD FUNCTION? |
| `03c_wellposedness_check.py` | The 4-, 5- and 6-node modules add an undirected STRING gene-gene edge, which we read generously as mutual protein stabilisation (complex formation protects both partners from degradation). |
| `03d_higher_order_compI1.py` | The three- to six-node escalation of 03_higher_order_sweep.py for the composite I1 circuit with negative feedback. |
| `03e_higher_order_convergence_audit.py` | HOW MUCH OF THE HIGHER-ORDER "GAIN" IS JUST AN UNCONVERGED INTEGRATION? |
| `03f_n6_limit_cycle_certification.py` | CERTIFY THE ONE QUALITATIVE GAIN OF FUNCTION FROM HIGHER ORDER. |
| `03g_higher_order_empirical_support.py` | HOW MUCH OF THE NETWORK ACTUALLY CARRIES THE EDGES THAT THE 4-, 5- AND 6-NODE ARCHITECTURES ARE BUILT FROM? |
| `04a_extract_mir29_data.R` | Measured quantities that constrain the miR-29–collagen circuit: TCGA-BRCA expression of miR-29a/b/c, NFKB1, SP1, RELA, COL1A1 and COL3A1, and matched CPTAC-BRCA mRNA and protein. |
| `04b_mir29_circuit_fit.py` | Parameterise one concrete circuit from the data. |
| `04c_mir29_modelfree_prediction.py` | The circuit fit in 04b gives "a 2.02-fold rise in miR-29 halves COL1A1 mRNA". |
| `05_figures.py` | Publication figures for the dynamical analysis. |
| `06_empirical_grounding.py` | Tie every modelled topology to its count in the BRCA network, and record the two structural facts that motivate the modelling choices: (i) no core in the network is Alon's C1 (all-activating): the miRNA->target edge is repressive … |
| `06b_composite_showcase.py` | Concrete, certified examples of the two behaviours that the COMPOSITE FFL has and no acyclic (transcriptional or non-composite miRNA) FFL can have: |
| `06c_composite_bistability_export.py` | CORRECTION TO THE BISTABILITY FIGURE. |
| `06d_composite_ringing_export.py` | CORRECTION TO THE RINGING FIGURE. |
| `07_consolidate.py` | Pull every headline number of the dynamical analysis into one table with its provenance (which script and which file produced it). |
| `07c_consolidate_additions.py` | Append the results of the later dynamics analyses to dynamics_key_results.csv, each with the file it came from. |
| `08_mir_vs_txn_paired.py` | (a) PAIRED comparison, at identical parameter vectors, of each miRNA-mediated FFL against the transcriptional FFL of the same sign structure. |
| `09_analytic_bounds.py` | Three analytic consequences of putting the repressing arm POST-transcriptionally, each derived in closed form and then checked numerically against the simulator. |
| `09b_composite_oscillation.py` | WHY A COMPOSITE FFL CAN DO SOMETHING NO TRANSCRIPTIONAL FFL CAN. |
| `09c_relaxation_rate_theorem.py` | THE EXACT STATEMENT THAT SEPARATES A miRNA FFL FROM A TRANSCRIPTIONAL FFL. |
| `10b_dynamics_extra_figures.py` | Figures for the dynamics analyses that 07c_consolidate_additions.py adds. |
| `99_dynamics_driver.sh` | Runs the remaining dynamical-analysis steps in order once the enlarged (N = 16,384) 3-node ensemble has been written by 02c_characterise_3node_N16384.py. |
| `99b_fcd_driver.sh` | Waits for the enlarged 3-node ensemble, then recomputes fold-change detection |
| `99c_dynamics_driver2.sh` | Runs the analyses that were written but never executed, plus the additions from this pass, on the existing 1,024-set 3-node ensemble. |
| `dyn_metrics.py` | Dynamical characterisation of the FFL topologies defined in dyn_models.py. |
| `dyn_models.py` | Non-dimensional ODE models of miRNA-TF-gene feed-forward loops, in the style of Mangan & Alon (2003 PNAS, doi:10.1073/pnas.2133841100) and Alon (2006, doi:10.1201/9781420011432), ADAPTED so that the miRNA arm acts POST- … |
| `dyn_sample.py` | scrambled-Sobol sampling of the dimensionless parameter space. |
| `dyn_style.py` | Shared figure style. Palette validated with the dataviz validator: node scripts/validate_palette.js "#2a78d6,#eb6834,#1baf7a,#4a3aa7" --mode light -> all checks PASS (aqua carries a contrast WARN, relieved by direct labels). |
