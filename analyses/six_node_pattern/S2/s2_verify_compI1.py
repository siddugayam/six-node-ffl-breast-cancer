#!/usr/bin/env python3
"""Verifies the stored COMP_I1 sweep (results/v3/dynamics_higher_order_persetset_compI1.csv.gz, written by
scripts/11_dynamics/03d_higher_order_compI1.py on 2026-09-09) before reusing it: sets 0-2047 of the core (n3) and of the
six-node module (n6) were re-run with the current code (perset/COMP_I1_negfeedback__{core,GG+MM+TT}@verify.csv.gz)
and are compared set by set with the stored values, as C3 of the first re-runs did for COMP_C2_toggle."""
import os, sys
sys.dont_write_bytecode = True
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REV = '/path/to/revision'
S = pd.read_csv(f'{REV}/results/v3/dynamics_higher_order_persetset_compI1.csv.gz')
rows = []
for mod, stored in (('core', 'n3'), ('GG+MM+TT', 'n6')):
    fn = f'{HERE}/perset/COMP_I1_negfeedback__{mod}@verify.csv.gz'
    if not os.path.exists(fn): print('not yet run:', fn); continue
    N = pd.read_csv(fn).set_index('param_set').sort_index()
    O = S[S.module == stored].set_index('param_set').loc[N.index]
    for c in N.columns:
        if c in ('family', 'module'): continue
        a = O[c].to_numpy(float); b = N[c].to_numpy(float)
        d = np.where(np.isnan(a) & np.isnan(b), 0.0, np.abs(a - b))
        rows.append(dict(module=mod, stored_module=stored, column=c, n_sets=len(a), n_differing_1e9=int((d > 1e-9).sum()), max_abs_diff=float(np.nanmax(d))))
R = pd.DataFrame(rows); R.to_csv(f'{HERE}/s2_verify_compI1.csv', index=False)
flags = [c for c in R.column.unique() if c.startswith('is_') or c in ('bistability_unresolved', 'nonmonotone_dose')]
for mod, g in R.groupby('module'):
    f = g[g.column.isin(flags)]
    print(f'{mod}: {g.n_sets.iloc[0]} sets; behaviour flags differing in {int(f.n_differing_1e9.sum())} set-flags; '
          f'max |difference| over all columns {g.max_abs_diff.max():.3g} (column {g.loc[g.max_abs_diff.idxmax(), "column"]})')
