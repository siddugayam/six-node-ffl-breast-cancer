#!/usr/bin/env python3
"""Summarise the eight n=7 RAND-ESU seeds (same seeds and q-profiles as scripts/03_ffl_census/run_n7*.sh)."""
import json, statistics as st, collections
seeds = [11, 22, 33, 44, 55, 71, 72, 73]
rows = []
for s in seeds:
    d = json.load(open(f'n7_seeds/census_n7_seed_{s}.json'))
    r = dict(seed=s, q='x'.join(str(1/d['prod_q'])[:8] for _ in [0]), prod_q=d['prod_q'], found=d['found'],
             est_total=d['estimated_modules'])
    for k, v in d['flags_found'].items():
        r['found_' + k] = v; r['est_' + k] = v / d['prod_q']; r['pct_' + k] = 100 * v / d['found']
    for k, v in d['n_classes_hist'].items():
        r[f'found_{k}_classes'] = v
    rows.append(r)
keys = [k for k in rows[0] if k.startswith(('est_', 'pct_'))]
out = dict(per_seed=rows, mean={k: st.mean(r[k] for r in rows) for k in keys},
           sd={k: st.stdev(r[k] for r in rows) for k in keys},
           max_classes_any_seed=max(int(k.split('_')[1]) for r in rows for k in r if k.endswith('_classes') and k.startswith('found_')),
           seven_class_found_per_seed={r['seed']: r.get('found_7_classes', 0) for r in rows})
json.dump(out, open('census_n7_randesu_8seeds_summary.json', 'w'), indent=1)
for k in keys:
    print(f"{k:60s} mean {out['mean'][k]:>16,.2f}   sd {out['sd'][k]:>14,.2f}")
print('max classes seen in any seed:', out['max_classes_any_seed'], ' 7-class modules sampled per seed:', out['seven_class_found_per_seed'])
