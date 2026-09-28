#!/usr/bin/env python3
"""P3a sensitivity analysis, descriptive, added on request after the primary result had been seen (2026-09-28 01:55):
25 of the 40 contrasts come from one series (GSE115646), so the contrasts are first averaged within each series and the
series are then pooled by REML, per tier.  Series value = unweighted mean of its contrasts' shifts; its SE = the mean of
their SEs (treats a series' contrasts as fully correlated, which is conservative for contrasts that share controls).
The decisions remain those of the fixed rule (p3a_summary.txt); this shows only whether a tier's result depends on one
series.  Input: p3a_per_dataset.tsv.  Writes p3a_series_sensitivity.tsv and p3a_series_sensitivity.txt (numbered lines)."""
import os, sys, subprocess
import pandas as pd
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = f'{INB}/P3'
R = pd.read_csv(f'{D}/p3a_per_dataset.tsv', sep='\t').dropna(subset=['shift', 'se'])
S = R.groupby(['tier', 'gse']).agg(yi=('shift', 'mean'), sei=('se', 'mean'), n_contrasts=('shift', 'size')).reset_index()
S.to_csv(f'{D}/p3a_series_sensitivity.tsv', sep='\t', index=False)
S.rename(columns={'tier': 'group'})[['group', 'yi', 'sei']].to_csv(f'{D}/p3a_series_pool_in.tsv', sep='\t', index=False)
subprocess.run(['Rscript', f'{INB}/_rma.R', f'{D}/p3a_series_pool_in.tsv', f'{D}/p3a_series_pooled.tsv'], check=True)
PT = pd.read_csv(f'{D}/p3a_series_pooled.tsv', sep='\t')
L = ['sensitivity analysis added after the primary result (descriptive; the decisions stay those of p3a_summary.txt): contrasts averaged within '
     'each series, then series pooled by REML']
for tier in ('strong', 'weak', 'predicted_only'):
    x = PT[PT.group == tier]
    if not len(x): L.append(f'{tier}: no series'); continue
    x = x.iloc[0]; s = S[S.tier == tier]
    g = s[s.gse == 'GSE115646']
    L.append(f'{tier}: k = {int(x.k)} series ({int(s.n_contrasts.sum())} contrasts); pooled shift {x.est:+.4f} (95 % CI {x.ci_lb:+.4f} to {x.ci_ub:+.4f}), '
             f'p {x.p:.3g}; tau2 {x.tau2:.3g}; CI excludes 0 above 0: {"yes" if x.ci_lb > 0 else "no"}; GSE115646 enters as one series'
             + (f' (mean of {int(g.n_contrasts.iloc[0])} contrasts: {g.yi.iloc[0]:+.4f}, SE {g.sei.iloc[0]:.4f})' if len(g) else ' (no GSE115646 contrast with >= 3 targets)'))
open(f'{D}/p3a_series_sensitivity.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
