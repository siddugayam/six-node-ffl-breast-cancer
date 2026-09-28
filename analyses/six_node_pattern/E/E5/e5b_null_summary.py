#!/usr/bin/env python3
"""E5(b) nulls: NULL-A, NULL-B, NULL-C (scripts/v2/v2_null.c, 1,000 randomisations, seed 20250908,
100 swaps per edge) on the analysed network (6,829 edges) with the 22 non-Lambert TF-typed nodes re-typed
as genes (null_dep_nolegacy_retyped22.txt: node-type line changed, fine classes recomputed), against the
stored original-typing runs analyses/census_and_motif_nulls/nulls/null_{A,B,C}_dep_nolegacy.tsv."""
import pandas as pd, numpy as np, os
H = os.path.dirname(os.path.abspath(__file__)); REV = '/path/to/revision'
rows = []
for typing, tmpl in (('original typing (stored)', f'{REV}/analyses/census_and_motif_nulls/nulls/null_{{m}}_dep_nolegacy.tsv'),
                     ('22 re-typed as genes', f'{H}/null_{{m}}_dep_nolegacy_retyped22.tsv')):
    for m in 'ABC':
        X = pd.read_csv(tmpl.format(m=m), sep='\t'); o = X[X.rep == 'obs'].iloc[0]; X = X[X.rep != 'obs']
        X = X.assign(named=X.comp + X.tf + X.mir); o = o.copy(); o['named'] = o.comp + o.tf + o.mir
        for col, lab in (('total', 'all three-node FFLs'), ('named', 'named-class (composite + TF-FFL + miRNA-FFL)'), ('comp', 'composite'),
                         ('tf', 'TF-FFL'), ('mir', 'miRNA-FFL')):
            x = X[col].astype(float).values; ob = float(o[col])
            rows.append(dict(typing=typing, null=f'NULL-{m}', motif=lab, observed=ob, null_mean=x.mean(), null_sd=x.std(ddof=1),
                             excess_pct=100 * (ob / x.mean() - 1) if x.mean() > 0 else np.nan,
                             p_upper=(np.sum(x >= ob) + 1) / (len(x) + 1), p_lower=(np.sum(x <= ob) + 1) / (len(x) + 1), replicates=len(x)))
T = pd.DataFrame(rows); T.to_csv(f'{H}/e5b_null_summary.csv', index=False)
pd.set_option('display.width', 220); print(T.to_string(index=False))
