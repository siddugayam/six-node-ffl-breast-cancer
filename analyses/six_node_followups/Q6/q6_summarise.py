#!/usr/bin/env python3
"""Q6 summary: BHAT6 composite (instances, node sets) and MODEL6-architecture (instances, node sets) on the three variant
graphs, NULL-C and NULL-L, 1,000 replicates each (Q6/runs).  Statistics and the decision rule are those of
analyses/six_node_pattern/S1/s1_summarise.py, verbatim: null mean, SD (ddof 1), observed/mean, Z, empirical p in both tails
(count + 1)/(R + 1); 'over-represented' if p_upper < 0.05 under both NULL-C and NULL-L; 'explained by the three-node
cores' if p_upper < 0.05 under NULL-C only; 'depleted' if p_lower < 0.05 under both; otherwise 'none of the three'."""
import os
import numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
GR = [('retyped22', '(i) 22 co-regulators typed as genes'), ('validated', '(ii) validated-only miRNA-target edges, with STRING'),
      ('physical900', '(iii) gene-gene links = STRING physical >= 0.900')]
PAT = [('BHAT6 composite, instances', 'bhat6_comp_inst'), ('BHAT6 composite, node sets', 'bhat6_comp_sets'),
       ('MODEL6-architecture (FULL), instances', 'model6_full_inst'), ('MODEL6-architecture (FULL), node sets', 'model6_full_sets')]
rows, dec = [], []
for g, glab in GR:
    O = pd.read_csv(f'{H}/runs/obs_{g}.tsv', sep='\t').iloc[0]
    T = {}
    for nm, m in (('NULL-C', 'C'), ('NULL-L', 'L')):
        fn = f'{H}/runs/null_{m}_{g}.tsv'
        if not os.path.exists(fn): continue
        X = pd.read_csv(fn, sep='\t'); X = X[X.rep.astype(str) != 'obs']
        for label, col in PAT:
            x = X[col].to_numpy(float); x = x[x >= 0]; R = len(x)
            if R == 0: continue
            o = float(O[col]); mu = x.mean(); sd = x.std(ddof=1)
            r = dict(graph=g, graph_label=glab, pattern=label, column=col, null=nm, observed=o, null_mean=mu, null_sd=sd,
                     obs_over_mean=o / mu if mu > 0 else np.nan, Z=(o - mu) / sd if sd > 0 else np.nan,
                     p_upper=(np.sum(x >= o) + 1) / (R + 1), p_lower=(np.sum(x <= o) + 1) / (R + 1), replicates=R)
            rows.append(r); T[(col, nm)] = r
    for label, col in PAT:
        if (col, 'NULL-C') not in T or (col, 'NULL-L') not in T: continue
        c, l = T[(col, 'NULL-C')], T[(col, 'NULL-L')]
        if c['p_upper'] < 0.05 and l['p_upper'] < 0.05: d = 'over-represented'
        elif c['p_upper'] < 0.05: d = 'explained by the three-node cores'
        elif c['p_lower'] < 0.05 and l['p_lower'] < 0.05: d = 'depleted'
        else: d = 'none of the three (not over-represented, not depleted)'
        dec.append(dict(graph=g, graph_label=glab, pattern=label, observed=c['observed'], NULL_C_ratio=c['obs_over_mean'], NULL_C_p_upper=c['p_upper'],
                        NULL_C_p_lower=c['p_lower'], NULL_L_ratio=l['obs_over_mean'], NULL_L_p_upper=l['p_upper'], NULL_L_p_lower=l['p_lower'], decision=d))
pd.DataFrame(rows).to_csv(f'{H}/q6_null_summary.csv', index=False); pd.DataFrame(dec).to_csv(f'{H}/q6_decision.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(pd.DataFrame(rows)[['graph', 'pattern', 'null', 'observed', 'null_mean', 'null_sd', 'obs_over_mean', 'Z', 'p_upper', 'p_lower', 'replicates']].to_string(index=False))
print(pd.DataFrame(dec)[['graph', 'pattern', 'observed', 'NULL_C_ratio', 'NULL_C_p_upper', 'NULL_L_ratio', 'NULL_L_p_upper', 'decision']].to_string(index=False))
