#!/usr/bin/env python3
"""A1 nested check (SETTINGS.md): the extended model (dyn_models_bhat.py) with TF2 -> TF1 absent reruns the paper's
stored six-node modules (COMP_C2_toggle and COMP_I1_negfeedback, n6 = GG+MM+TT, 16,384 sets, seed 20260908). Every
per-set behaviour flag must equal the stored one; every other stored column is compared too.
usage: python3 a1_nested_compare.py <analysis root> <nested perset folder> <output file>
Exit status 0 only if every flag of every set is equal."""
import os, sys
import numpy as np, pandas as pd

FLAGS = ['is_memory', 'is_pulse', 'is_damped_osc', 'is_sustained_osc', 'is_bistable', 'bistability_unresolved', 'is_ultrasensitive',
         'is_very_ultrasensitive', 'is_noise_rejecting']
STORED = {'COMP_C2_toggle': 'results/v3/dynamics_higher_order_persetset.csv.gz',
          'COMP_I1_negfeedback': 'results/v3/dynamics_higher_order_persetset_compI1.csv.gz'}


def main(root, folder, out):
    L = []; ok = True
    for fam, rel in STORED.items():
        S = pd.read_csv(os.path.join(root, rel)); S = S[(S.family == fam) & (S.module == 'n6')].sort_values('param_set').reset_index(drop=True)
        X = pd.read_csv(os.path.join(folder, f'{fam}__GG+MM+TT@nested.csv.gz')).sort_values('param_set').reset_index(drop=True)
        assert len(S) == len(X) == 16384 and (S.param_set.values == X.param_set.values).all()
        nd = {f: int((S[f].fillna(-9).values != X[f].fillna(-9).values).sum()) for f in FLAGS}
        num = [c for c in S.columns if c not in ('family', 'module', 'param_set') and c in X.columns and np.issubdtype(S[c].dtype, np.number)]
        mx = max(((float(np.nanmax(np.abs(S[c].values - X[c].values))) if np.isfinite(S[c].values - X[c].values).any() else 0.0), c) for c in num)
        nan_mismatch = sum(int((S[c].isna().values != X[c].isna().values).sum()) for c in num)
        prev = ', '.join(f'{f} {100 * S[f].mean():.3f} % vs {100 * X[f].mean():.3f} %' for f in ('is_memory', 'is_pulse', 'is_damped_osc', 'is_sustained_osc'))
        good = all(v == 0 for v in nd.values())
        ok &= good
        L.append(f'NESTED CHECK {fam} GG+MM+TT: 16,384 sets; behaviour flags differing: ' + ', '.join(f'{f} {v}' for f, v in nd.items())
                 + f'; largest |difference| over the other numeric columns {mx[0]:.3g} ({mx[1]}); missing-value mismatches {nan_mismatch}; '
                 f'prevalence stored vs rerun: {prev} -> {"PASS" if good else "FAIL"}')
        print(L[-1], flush=True)
    open(out, 'w').write('\n'.join(L) + '\n')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main(*sys.argv[1:4])
