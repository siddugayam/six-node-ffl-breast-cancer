#!/usr/bin/env python3
"""
S1 summary (INBOX_2026-09-27b).  Reads runs/null_{B,C,L,LT,LG,LM}.tsv (s1_null6) and obs/ (observed).

Checks first
  * NULL-B and NULL-C: the three-node columns (total, comp, tf, mir, other, recip) of every replicate
    must equal the stored runs INBOX_2026-09-26b/nulls/null_{B,C}_pub_nolegacy.tsv row for row.
  * observed three-node count 5,833 and the stored census-graph result (-13.9 % / -10.3 %).
Statistics per pattern and null: observed, null mean and SD, observed/mean, Z, empirical p in both
tails (k+1)/(R+1), replicates; and instances per three-node composite core (count3 'comp', the D1-D4
three-node FFLs whose TF->miRNA arc is a reciprocal pair), observed and null mean of the per-replicate
ratio; for BHAT6 also per Bhat three-node composite instance.
Decision rule (fixed in the analysis plan): 'over-represented' if p_upper < 0.05 under both NULL-C and NULL-L;
'explained by the three-node cores' if p_upper < 0.05 under NULL-C but not under NULL-L; 'depleted' if
p_lower < 0.05 under both; otherwise 'none of the three'.
"""
import os, csv
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
REV = '/path/to/revision'
lines = []
P = lambda *a: (lines.append(' '.join(str(x) for x in a)), print(*a))
O = pd.read_csv(f'{HERE}/obs_nolegacy.tsv', sep='\t').iloc[0]
PAT = [('three-node FFLs (reference)', 'total3'), ('three-node composite cores', 'comp3'),
       ('BHAT6 composite, instances', 'bhat6_comp_inst'), ('BHAT6 composite, node sets', 'bhat6_comp_sets'),
       ('BHAT6 miRNA-FFL, instances', 'bhat6_mirFFL_inst'), ('BHAT6 miRNA-FFL, node sets', 'bhat6_mirFFL_sets'),
       ('BHAT6 TF-FFL, instances', 'bhat6_TFFFL_inst'), ('BHAT6 TF-FFL, node sets', 'bhat6_TFFFL_sets'),
       ('MODEL6-architecture (FULL), instances', 'model6_full_inst'), ('MODEL6-architecture (FULL), node sets', 'model6_full_sets'),
       ('MODEL6 defining edges (DEF), instances', 'model6_def_inst'), ('MODEL6 defining edges (DEF), node sets', 'model6_def_sets'),
       ('CENSUS6 (c)-with-core modules', 'census6c')]
NULLS = [('NULL-B', 'B'), ('NULL-C', 'C'), ('NULL-L', 'L'), ('NULL-L TRRUST layer only', 'LT'),
         ('NULL-L STRING layer only', 'LG'), ('NULL-L 10 kb layer only', 'LM')]
runs = {}
for nm, m in NULLS:
    fn = f'{HERE}/runs/null_{m}.tsv'
    if os.path.exists(fn):
        X = pd.read_csv(fn, sep='\t'); X = X[X.rep.astype(str) != 'obs']; X['rep'] = X.rep.astype(int)
        runs[m] = X.sort_values('rep').reset_index(drop=True)
# ---- checks
for m in ('B', 'C'):
    if m not in runs: continue
    S = pd.read_csv(f'{REV}/INBOX_2026-09-26b/nulls/null_{m}_pub_nolegacy.tsv', sep='\t')
    so = S[S.rep == 'obs'].iloc[0]; S = S[S.rep != 'obs'].copy(); S['rep'] = S.rep.astype(int)
    X = runs[m].merge(S, on='rep')
    same = all((X[a] == X[b]).all() for a, b in (('total3', 'total'), ('comp3', 'comp'), ('tf3', 'tf'), ('mir3', 'mir'),
                                               ('other3', 'other'), ('recip_x', 'recip_y')))
    P(f'CHECK NULL-{m}: {len(X)} replicates compared with the stored run; three-node columns identical in every replicate: {same}; '
      f'stored observed total {so.total}, this observed total {O.total3}; stored-run mean {S.total.mean():.1f} '
      f'-> observed/mean - 1 = {100 * (int(so.total) / S.total.mean() - 1):+.1f} %')
rows = []
for nm, m in NULLS:
    if m not in runs: continue
    X = runs[m]
    for label, col in PAT:
        x = X[col].to_numpy(float); ok = x >= 0; x = x[ok]; R = len(x)
        if R == 0: continue
        o = float(O[col]); mu = x.mean(); sd = x.std(ddof=1)
        r = dict(pattern=label, column=col, null=nm, observed=o, null_mean=mu, null_sd=sd,
                 obs_over_mean=o / mu if mu > 0 else np.nan, Z=(o - mu) / sd if sd > 0 else np.nan,
                 p_upper=(np.sum(x >= o) + 1) / (R + 1), p_lower=(np.sum(x <= o) + 1) / (R + 1), replicates=R)
        if col not in ('total3', 'comp3'):
            c3 = X.loc[ok, 'comp3'].to_numpy(float)
            r['per_3node_composite_core_observed'] = o / float(O['comp3'])
            r['per_3node_composite_core_null_mean'] = np.mean(x / c3)
            if col.startswith('bhat6'):
                b3 = X.loc[ok, 'bhat3_comp'].to_numpy(float)
                r['per_Bhat3_composite_observed'] = o / float(O['bhat3_comp'])
                r['per_Bhat3_composite_null_mean'] = np.mean(x / b3)
        rows.append(r)
T = pd.DataFrame(rows)
T.to_csv(f'{HERE}/s1_null_summary.csv', index=False)
# ---- decision rule
dec = []
for label, col in PAT:
    c = T[(T.column == col) & (T.null == 'NULL-C')]; l = T[(T.column == col) & (T.null == 'NULL-L')]
    if c.empty or l.empty: continue
    c, l = c.iloc[0], l.iloc[0]
    if c.p_upper < 0.05 and l.p_upper < 0.05: d = 'over-represented'
    elif c.p_upper < 0.05: d = 'explained by the three-node cores'
    elif c.p_lower < 0.05 and l.p_lower < 0.05: d = 'depleted'
    else: d = 'none of the three (not over-represented, not depleted)'
    dec.append(dict(pattern=label, observed=c.observed, NULL_C_ratio=c.obs_over_mean, NULL_C_p_upper=c.p_upper,
                    NULL_C_p_lower=c.p_lower, NULL_L_ratio=l.obs_over_mean, NULL_L_p_upper=l.p_upper,
                    NULL_L_p_lower=l.p_lower, decision=d))
pd.DataFrame(dec).to_csv(f'{HERE}/s1_decision.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(T[['pattern', 'null', 'observed', 'null_mean', 'null_sd', 'obs_over_mean', 'Z', 'p_upper', 'p_lower', 'replicates']].to_string(index=False))
print(pd.DataFrame(dec).to_string(index=False))
open(f'{HERE}/s1_checks.txt', 'w').write('\n'.join(lines) + '\n')
