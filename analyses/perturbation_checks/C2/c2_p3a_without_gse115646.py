#!/usr/bin/env python3
"""C2 (analyses/perturbation_checks): P3a at contrast level without series GSE115646.  Descriptive, added after the primary result;
the 27d decisions are unchanged (SETTINGS.md, C2).  Every GSE115646 row of the 27d per-dataset table is removed, the other
contrasts are kept as they are (no replacement is drawn), and each tier is pooled by REML as in 27d.
Input: analyses/perturbation_tests/P3/p3a_per_dataset.tsv and p3a_pooled.tsv (read only).
Writes c2_pool_in.tsv, c2_pooled.tsv and c2_summary.txt (numbered lines)."""
import os, subprocess
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); INB = os.path.dirname(HERE); REV = os.path.dirname(os.path.dirname(INB))
SRC = f'{REV}/analyses/perturbation_tests/P3'
R = pd.read_csv(f'{SRC}/p3a_per_dataset.tsv', sep='\t')
R = R[np.isfinite(R['shift']) & np.isfinite(R['se']) & (R['se'] > 0)]
drop = R.gse == 'GSE115646'
X = R[~drop]
X.rename(columns={'tier': 'group', 'shift': 'yi', 'se': 'sei'})[['group', 'yi', 'sei']].to_csv(f'{HERE}/c2_pool_in.tsv', sep='\t', index=False)
subprocess.run(['Rscript', f'{INB}/_rma.R', f'{HERE}/c2_pool_in.tsv', f'{HERE}/c2_pooled.tsv'], check=True)
PT = pd.read_csv(f'{HERE}/c2_pooled.tsv', sep='\t'); P0 = pd.read_csv(f'{SRC}/p3a_pooled.tsv', sep='\t')
L = []; P = L.append
for tier in ('strong', 'weak', 'predicted_only'):
    x = PT[PT.group == tier].iloc[0]; x0 = P0[P0.group == tier].iloc[0]
    nrm = int((drop & (R.tier == tier)).sum()); nser = X[X.tier == tier].gse.nunique()
    P(f'{tier} without GSE115646 (descriptive, added after the primary result): k = {int(x.k)} contrasts from {nser} series '
      f'({nrm} GSE115646 contrasts removed); pooled shift {x.est:+.4f} (95 % CI {x.ci_lb:+.4f} to {x.ci_ub:+.4f}), p {x.p:.3g}; tau2 {x.tau2:.3g}; '
      f'CI excludes 0 above 0: {"yes" if x.ci_lb > 0 else "no"} | all 40 contrasts (27d): k = {int(x0.k)}, {x0.est:+.4f} '
      f'(95 % CI {x0.ci_lb:+.4f} to {x0.ci_ub:+.4f}), p {x0.p:.3g}')
open(f'{HERE}/c2_summary.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
