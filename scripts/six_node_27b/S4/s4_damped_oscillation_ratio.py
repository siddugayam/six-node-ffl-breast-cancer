#!/usr/bin/env python3
"""S4: how oscillatory are the damped oscillations?  For every module with >= 50 damped-oscillation
sets (always including COMP_C2_toggle GG+MM+TT, 862 sets, and its TT control, 1,404 sets), the leading
eigenvalue at the operating point (lead_eig_re, lead_eig_im = |Im| of the eigenvalue with the largest
real part, stored per set by run_module() of 03_higher_order_sweep.py) gives |Im l| / |Re l|.
Reported: median and IQR of the ratio, and the number of sets with |Im l| >= |Re l| (ratio >= 1).
Reads the stored and new per-set files through S2/s2_common.py.  Read-only on the project."""
import sys, os
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'S2'))
import numpy as np, pandas as pd, s2_common as C
D = C.load_all()
rows = []
for (f, m), g in D.groupby(['family', 'module']):
    d = g[g.is_damped_osc > 0.5]
    must = (f, m) in (('COMP_C2_toggle', 'GG+MM+TT'), ('COMP_C2_toggle', 'TT'))
    if len(d) < 50 and not must: continue
    assert (d.lead_eig_re < 0).all() and (d.lead_eig_im > 1e-3).all()
    r = d.lead_eig_im.abs() / d.lead_eig_re.abs()
    q1, med, q3 = np.percentile(r, [25, 50, 75])
    rows.append(dict(family=f, module=m, n_damped_sets=len(d), median_ratio=med, iqr_q1=q1, iqr_q3=q3,
                     n_ratio_ge_1=int((r >= 1).sum()), pct_ratio_ge_1=100 * (r >= 1).mean(),
                     max_ratio=r.max(), source=g.source.iloc[0]))
R = pd.DataFrame(rows)
R.to_csv(os.path.join(HERE, 's4_damped_oscillation_ratio.csv'), index=False)
pd.set_option('display.width', 200); print(R.drop(columns='source').to_string(index=False))
