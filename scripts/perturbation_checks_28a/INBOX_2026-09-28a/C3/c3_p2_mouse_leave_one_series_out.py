#!/usr/bin/env python3
"""C3 (INBOX_2026-09-28a): P2 mouse, leaving one series out.  Descriptive, added after the primary result; the 27d
decisions are unchanged (SETTINGS.md, C3).  For the six-node class test and the site-adjusted TargetScan coefficient, every
series (GSE) with a contrast in the test is removed in turn, with all of its contrasts, and the rest are pooled by REML.
Input: INBOX_2026-09-27d/P2/p2_datasets.tsv (read only).
Writes c3_pool_in.tsv, c3_pooled.tsv, c3_leave_one_series_out.tsv and c3_summary.txt (numbered lines)."""
import os, subprocess
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); INB = os.path.dirname(HERE); REV = os.path.dirname(INB)
D = pd.read_csv(f'{REV}/INBOX_2026-09-27d/P2/p2_datasets.tsv', sep='\t')
TESTS = [('six-node class', 'six-node class vs 1-member (network edges)'),
         ('site-adjusted TargetScan', 'site-adjusted: >=2 members coefficient (TargetScan)')]
M = D[(D.cls == 'mouse') & D.test.isin([t for _, t in TESTS])]
M = M[np.isfinite(M['diff']) & np.isfinite(M['se']) & (M['se'] > 0)]
pin = []
for lab, t in TESTS:
    x = M[M.test == t]
    pin += [dict(group=f'{lab}|all', yi=v, sei=s) for v, s in zip(x['diff'], x['se'])]
    for g in sorted(x.gse.unique()):
        y = x[x.gse != g]
        pin += [dict(group=f'{lab}|{g}', yi=v, sei=s) for v, s in zip(y['diff'], y['se'])]
pd.DataFrame(pin).to_csv(f'{HERE}/c3_pool_in.tsv', sep='\t', index=False)
subprocess.run(['Rscript', f'{INB}/_rma.R', f'{HERE}/c3_pool_in.tsv', f'{HERE}/c3_pooled.tsv'], check=True)
PT = pd.read_csv(f'{HERE}/c3_pooled.tsv', sep='\t')
PT['test'] = PT.group.str.split('|').str[0]; PT['left_out'] = PT.group.str.split('|').str[1]
rows, L = [], []; P = L.append
for lab, t in TESTS:
    x = M[M.test == t]; a = PT[(PT.test == lab) & (PT.left_out == 'all')].iloc[0]
    P(f'{lab}, mouse, all series (27d): k = {int(a.k)} contrasts from {x.gse.nunique()} series; {a.est:+.4f} (95 % CI {a.ci_lb:+.4f} to {a.ci_ub:+.4f}), p {a.p:.3g}')
    lo = PT[(PT.test == lab) & (PT.left_out != 'all')].copy()
    lo['n_removed'] = [int((x.gse == g).sum()) for g in lo.left_out]; lo['change'] = lo.est - a.est
    lo = lo.reindex(lo.change.abs().sort_values(ascending=False).index)
    for r in lo.itertuples():
        rows.append(dict(test=lab, left_out=r.left_out, n_removed=r.n_removed, k=int(r.k), est=r.est, ci_lb=r.ci_lb, ci_ub=r.ci_ub, p=r.p, change=r.change))
    top = lo.iloc[0]; inc0 = lo[(lo.ci_lb <= 0) & (lo.ci_ub >= 0)]
    P(f'{lab}, leave one series out (descriptive, added after the primary result): {len(lo)} runs; estimates {lo.est.min():+.4f} to {lo.est.max():+.4f}; '
      f'largest change without {top.left_out} ({top.n_removed} contrast{"" if top.n_removed == 1 else "s"}): {top.est:+.4f} (95 % CI {top.ci_lb:+.4f} to {top.ci_ub:+.4f}), p {top.p:.3g}, '
      f'change {top.change:+.4f}; runs whose 95 % CI includes 0: {len(inc0)}' + (f' ({", ".join(inc0.left_out)})' if len(inc0) else ''))
    for g in ('GSE118698', 'GSE63813'):
        z = lo[lo.left_out == g]
        P(f'{lab}, without {g}: ' + (f'{int(z.n_removed.iloc[0])} contrast{"" if int(z.n_removed.iloc[0]) == 1 else "s"} removed; k = {int(z.k.iloc[0])}; {z.est.iloc[0]:+.4f} (95 % CI {z.ci_lb.iloc[0]:+.4f} to '
          f'{z.ci_ub.iloc[0]:+.4f}), p {z.p.iloc[0]:.3g}; change {z.change.iloc[0]:+.4f}' if len(z) else 'no contrast of this series in the test'))
pd.DataFrame(rows).to_csv(f'{HERE}/c3_leave_one_series_out.tsv', sep='\t', index=False)
open(f'{HERE}/c3_summary.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
