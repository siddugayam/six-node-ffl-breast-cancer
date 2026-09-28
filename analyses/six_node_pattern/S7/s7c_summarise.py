#!/usr/bin/env python3
"""S7(c) summary: exhaustive census (ffl_census_composition, max-flow D4) at n = 3-6 on the validated-only
graphs with and without STRING; NULL-B / NULL-C (v2_null.c, 1,000 randomisations, seed 20250908, 100 swaps
per edge) for three-node FFLs; BHAT6 and MODEL6 counts (s1_null6 OBS).  Reference: the nolegacy graph."""
import glob, os, csv
import numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
REF = {'nolegacy (reference, all tiers)': dict(n3=5833, n4=116505, n5=1506263, n6=19401840),
       'nolegacy_nostring (reference)': dict(n3=3432, n4=81104, n5=791085, n6=8751860)}
rows = []
for v in ('validated', 'validated_nostring'):
    r = dict(graph=v)
    for k in (3, 4, 5):
        L = open(f'census/{v}/n{k}.txt').read().split('\n'); r[f'n{k}'] = int([l for l in L if l.startswith('found')][0].split()[1])
        r[f'max_classes_n{k}'] = max(bin(int(l.split()[1])).count('1') for l in L if l.startswith('mask'))
    tot = 0; masks = set()
    for f in glob.glob(f'census/{v}/n6_parts/part_*.txt'):
        for l in open(f):
            if l.startswith('found'): tot += int(l.split()[1])
            if l.startswith('mask'): masks.add(int(l.split()[1]))
    r['n6'] = tot; r['n6_parts'] = len(glob.glob(f'census/{v}/n6_parts/part_*.txt')); r['max_classes_n6'] = max(bin(m).count('1') for m in masks)
    o = pd.read_csv(f'obs_{v}.tsv', sep='\t').iloc[0]
    r['count3_check'] = int(o.total3) == r['n3']
    for c in ('bhat3_comp', 'bhat6_comp_inst', 'bhat6_comp_sets', 'bhat6_mirFFL_inst', 'bhat6_mirFFL_sets', 'bhat6_TFFFL_inst',
              'bhat6_TFFFL_sets', 'model6_full_inst', 'model6_full_sets', 'census6c'):
        r[c] = int(o[c])
    for m in 'BC':
        X = pd.read_csv(f'nulls/null_{m}_{v}.tsv', sep='\t'); ob = X[X.rep == 'obs'].iloc[0]; X = X[X.rep != 'obs']
        for col, lab in (('total', 'all'), ('comp', 'composite')):
            x = X[col].astype(float).values; o_ = float(ob[col])
            r[f'NULL{m}_{lab}_observed'] = o_; r[f'NULL{m}_{lab}_null_mean'] = x.mean(); r[f'NULL{m}_{lab}_null_sd'] = x.std(ddof=1)
            r[f'NULL{m}_{lab}_excess_pct'] = 100 * (o_ / x.mean() - 1); r[f'NULL{m}_{lab}_p_upper'] = (np.sum(x >= o_) + 1) / (len(x) + 1)
            r[f'NULL{m}_{lab}_p_lower'] = (np.sum(x <= o_) + 1) / (len(x) + 1); r[f'NULL{m}_replicates'] = len(x)
    rows.append(r)
T = pd.DataFrame(rows); T.to_csv(f'{H}/s7c_validated_only_summary.csv', index=False)
pd.set_option('display.width', 250); print(T.T.to_string())
