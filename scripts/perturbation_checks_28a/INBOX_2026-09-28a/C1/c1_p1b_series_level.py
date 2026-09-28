#!/usr/bin/env python3
"""C1 (INBOX_2026-09-28a): P1b fibroblast lineage vs carcinoma at series level.  Descriptive, added after the primary
result; the 27d decision is unchanged (SETTINGS.md, C1).
The contrasts of one series are averaged: the sign-aligned value is the unweighted mean over the series' contrasts with a
value, and its SE is the mean of their SEs (the rule of the 27d P3a series check).  The series are then compared by REML
meta-regression on class (fibroblast lineage 1, carcinoma 0), for the collagen pseudo-gene, COL1A1 alone and COL3A1 alone.
Input: INBOX_2026-09-27d/P1b/p1b_contrasts.tsv and p1b_meta_regression.tsv (read only).
Writes c1_series_values.tsv, c1_mr_in.tsv, c1_meta_regression.tsv and c1_summary.txt (numbered lines)."""
import os, subprocess
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); INB = os.path.dirname(HERE); REV = os.path.dirname(INB)
SRC = f'{REV}/INBOX_2026-09-27d/P1b'
C = pd.read_csv(f'{SRC}/p1b_contrasts.tsv', sep='\t')
C = C[(C.mirna_set == 'miR-29') & (C.tag == 'primary') & C.cls.isin(['fibroblast lineage', 'carcinoma'])]
OUT = [('collagen pseudo-gene', 'pair_repression', 'pair_SE', 'collagen pseudo-gene (primary)'),
       ('COL1A1 alone', 'COL1A1_repression', 'COL1A1_SE', 'COL1A1 alone'),
       ('COL3A1 alone', 'COL3A1_repression', 'COL3A1_SE', 'COL3A1 alone')]
rows, mr = [], []
for lab, y, s, _ in OUT:
    x = C[np.isfinite(C[y]) & np.isfinite(C[s]) & (C[s] > 0)]
    for (gse, cls), g in x.groupby(['gse', 'cls']):
        r = dict(outcome=lab, gse=gse, cls=cls, n_contrasts=len(g), cells=';'.join(g.cell), yi=float(g[y].mean()), sei=float(g[s].mean()))
        rows.append(r); mr.append(dict(group=lab, yi=r['yi'], sei=r['sei'], mod=1 if cls == 'fibroblast lineage' else 0))
S = pd.DataFrame(rows); S.to_csv(f'{HERE}/c1_series_values.tsv', sep='\t', index=False)
pd.DataFrame(mr).to_csv(f'{HERE}/c1_mr_in.tsv', sep='\t', index=False)
subprocess.run(['Rscript', f'{INB}/_rma.R', f'{HERE}/c1_mr_in.tsv', f'{HERE}/c1_meta_regression.tsv'], check=True)
MR = pd.read_csv(f'{HERE}/c1_meta_regression.tsv', sep='\t')
M0 = pd.read_csv(f'{SRC}/p1b_meta_regression.tsv', sep='\t')
L = []; P = L.append
for lab, y, s, lab27 in OUT:
    x = MR[MR.group == lab].set_index('term'); s_ = S[S.outcome == lab]
    nf, nc = int((s_.cls == 'fibroblast lineage').sum()), int((s_.cls == 'carcinoma').sum())
    x0 = M0[(M0.group == lab27)].set_index('term')
    ref = (f' | contrast level (27d): {x0.loc["mod", "est"]:+.3f} (95 % CI {x0.loc["mod", "ci_lb"]:+.3f} to {x0.loc["mod", "ci_ub"]:+.3f}), '
           f'p {x0.loc["mod", "p"]:.3g}, k = {int(x0.loc["mod", "k"])} contrasts') if 'mod' in x0.index else ''
    if 'mod' not in x.index or not np.isfinite(x.loc['mod', 'est']):
        P(f'{lab}, series level (descriptive, added after the primary result): not estimable (fibroblast {nf}, carcinoma {nc} series){ref}'); continue
    b, p = x.loc['mod', 'est'], x.loc['mod', 'p']
    P(f'{lab}, series level (descriptive, added after the primary result): carcinoma mean {x.loc["intercept", "est"]:+.3f}; fibroblast minus carcinoma '
      f'{b:+.3f} (95 % CI {x.loc["mod", "ci_lb"]:+.3f} to {x.loc["mod", "ci_ub"]:+.3f}), z-test p {p:.3g}; k = {int(x.loc["mod", "k"])} series '
      f'(fibroblast {nf}, carcinoma {nc}); tau2 {x.loc["mod", "tau2"]:.3g}; 27d criterion (> 0 with p < 0.05) met: {"yes" if (b > 0 and p < 0.05) else "no"}{ref}')
for r in S.itertuples():
    P(f'series value, {r.outcome}: {r.gse} [{r.cls}; {r.n_contrasts} contrast(s): {r.cells}] {r.yi:+.3f} (SE {r.sei:.3f})')
open(f'{HERE}/c1_summary.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
