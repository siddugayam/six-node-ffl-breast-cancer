# S2

## Scripts in this folder

| Script | What it does |
|---|---|
| `s2_analyse.py` | S2 / S3 / S6 analysis (analyses/six_node_pattern). |
| `s2_certify_limit_cycles.py` | Certification of sustained oscillation for every S2/S3 module, exactly as scripts/11_dynamics/03f_n6_limit_cycle_certification.py does it for COMP_C2_toggle n6: every set flagged is_sustained_osc by the sweep is re-integrated for … |
| `s2_common.py` | Shared loader for S2/S3/S4 (analyses/six_node_pattern): every per-set result of the 16,384-set Sobol design (seed 20260908), stored and new, under the factorial module labels. |
| `s2_run_modules.py` | S2 (layer-factorial dynamics) and S3 (TF-TF sign controls) for analyses/six_node_pattern. |
