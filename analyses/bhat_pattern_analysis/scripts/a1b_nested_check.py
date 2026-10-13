#!/usr/bin/env python3
"""A1b nested check (SETTINGS.md A1b; it gates the A1b runs). The four nested modules (a1b_modules.py: S2's stored
COMP_I1_negfeedback TT, GG+TT, MM and MM+TT, rebuilt through a1b_run.py's layers) must equal S2's per-set files:
  - every behaviour flag of every set;
  - the certification of the flagged sets (a1_certify.py's procedure: 1,000 tau from two initial conditions) must
    equal S2's s2_limit_cycle_certification.csv set by set (S2: TT 21 of 26, GG+TT 20 of 25, MM+TT 34 of 38, MM none).
The largest difference over the other numeric columns and the missing-value mismatches are reported, not gated.
usage: python3 a1b_nested_check.py <analysis root> <nested module list csv> <nested perset folder> <S2 folder>
                                   <output txt> [workers]
Outputs: <output txt>; per-set certification in <output stem>_certification.csv. Exit status 0 only if all pass."""
import sys, os
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from multiprocessing import Pool
import a1b_run                       # the extended layers (reads the analysis root from argv[1])
import a1_certify as C

FLAGS = ['is_memory', 'is_pulse', 'is_damped_osc', 'is_sustained_osc', 'is_bistable', 'bistability_unresolved', 'is_ultrasensitive',
         'is_very_ultrasensitive', 'is_noise_rejecting']


def main(root, mlist, folder, s2, out, workers=4):
    M = pd.read_csv(mlist, dtype=str).fillna('')
    S2C = pd.read_csv(os.path.join(s2, 's2_limit_cycle_certification.csv'))
    res, tasks = {}, []
    for row in M.to_dict('records'):
        lay = row['size_key']
        X = pd.read_csv(os.path.join(folder, f"{row['family']}__{row['module']}.csv.gz")).sort_values('param_set').reset_index(drop=True)
        S = pd.read_csv(os.path.join(s2, 'perset', f'COMP_I1_negfeedback__{lay}.csv.gz')).sort_values('param_set').reset_index(drop=True)
        assert len(X) == len(S) == 16384 and (X.param_set.values == S.param_set.values).all()
        nd = {f: int((S[f].fillna(-9).values != X[f].fillna(-9).values).sum()) for f in FLAGS}
        num = [c for c in S.columns if c not in ('family', 'module', 'param_set') and c in X.columns and np.issubdtype(S[c].dtype, np.number)]
        d = [np.abs(S[c].values - X[c].values) for c in num]
        mx = max((float(np.nanmax(v)) if np.isfinite(v).any() else 0.0) for v in d)
        nanm = sum(int((S[c].isna().values != X[c].isna().values).sum()) for c in num)
        fl = np.sort(X.loc[X.is_sustained_osc > 0.5, 'param_set'].to_numpy(int))
        for k in range(0, len(fl), 64): tasks.append((row, fl[k:k + 64]))
        res[row['module']] = (lay, nd, mx, nanm)
    cols = ['family', 'module', 'param_set', 'certified_limit_cycle']
    if tasks:
        with Pool(int(workers)) as pool:
            cert = pd.concat(list(pool.imap_unordered(C.certify, tasks)), ignore_index=True).sort_values(['module', 'param_set'])
    else:
        cert = pd.DataFrame(columns=cols)
    cert.to_csv(os.path.splitext(out)[0] + '_certification.csv', index=False)
    L, ok = [], True
    for mod, (lay, nd, mx, nanm) in res.items():
        c = cert[cert.module == mod]
        s = S2C[(S2C.family == 'COMP_I1_negfeedback') & (S2C.module == lay)]
        m = c[['param_set', 'certified_limit_cycle']].merge(s[['param_set', 'certified_limit_cycle']], on='param_set', how='outer', indicator=True)
        same = bool((m._merge == 'both').all())
        agree = same and bool((m.certified_limit_cycle_x.astype(bool) == m.certified_limit_cycle_y.astype(bool)).all())
        good = all(v == 0 for v in nd.values()) and agree
        ok &= good
        L.append(f"A1B NESTED CHECK COMP_I1_negfeedback {lay}: 16,384 sets; behaviour flags differing: "
                 + ', '.join(f'{f} {v}' for f, v in nd.items())
                 + f"; largest |difference| over the other numeric columns {mx:.3g}; missing-value mismatches {nanm}; certification "
                   f"{int(c.certified_limit_cycle.astype(bool).sum())} of {len(c)} flagged here vs {int(s.certified_limit_cycle.astype(bool).sum())} of "
                   f"{len(s)} in S2, same flagged sets {same}, set-by-set agreement {agree} -> {'PASS' if good else 'FAIL'}")
        print(L[-1], flush=True)
    open(out, 'w').write('\n'.join(L) + '\n')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main(*sys.argv[1:7])
