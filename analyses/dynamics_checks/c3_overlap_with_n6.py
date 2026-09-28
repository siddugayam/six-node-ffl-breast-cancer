#!/usr/bin/env python3
"""C3 follow-up: do the same parameter sets show the behaviour in the four-node TF2-only module and in
the stored six-node module?  Reads results/v3/dynamics_higher_order_persetset.csv.gz (read-only) and
c3_persetset.csv.gz; writes c3_overlap_with_n6.csv."""
import pandas as pd
S = pd.read_csv('/path/to/revision/results/v3/dynamics_higher_order_persetset.csv.gz')
n6 = S[(S.family == 'COMP_C2_toggle') & (S.module == 'n6')].set_index('param_set').sort_index()
C = pd.read_csv('c3_persetset.csv.gz')
tf = C[C.module == 'n4tf'].set_index('param_set').sort_index()
rows = []
for b in ['is_pulse', 'is_damped_osc', 'is_sustained_osc', 'is_ultrasensitive', 'is_bistable', 'is_memory']:
    a = tf[b] > 0.5; c = n6[b] > 0.5
    rows.append(dict(behaviour=b, n_sets_4node_TF2only=int(a.sum()), n_sets_6node=int(c.sum()), n_both=int((a & c).sum()),
                     pct_of_6node_sets_also_in_4node=round(100 * (a & c).sum() / max(c.sum(), 1), 1)))
R = pd.DataFrame(rows); R.to_csv('c3_overlap_with_n6.csv', index=False); print(R.to_string(index=False))
