#!/usr/bin/env python3
"""B2: modules with >=2 TF, >=2 miRNA and >=2 Gene, n = 3-7.
Source 1: the stored 03_ffl_census.py outputs (composition fields of results/ffl_census_k{4,6,7}.json;
          for k3 and k5 the JSONs now on disk were overwritten on 2026-09-12 by the v7 no-legacy-miRNA
          script, so the Table 3 graph composition is read from logs/classdiv_k3.log and logs/classdiv_k5.log).
Source 2: this re-run (exhaustive n=3-6; n=7 RAND-ESU, same 8 seeds as the v2 census)."""
import json, re, csv
REV = '/path/to/revision'

def qual(key):
    a, b, c = map(int, re.match(r'TF(\d+)_miR(\d+)_G(\d+)', key).groups()); return a >= 2 and b >= 2 and c >= 2

def from_log(path):
    comp, total = {}, None
    for line in open(path):
        m = re.search(r'TF=(\d+) miRNA=(\d+) Gene=(\d+)\s+->\s+([\d,]+)\s+\(sampled (\d+)\)', line)
        if m: comp[f"TF{m[1]}_miR{m[2]}_G{m[3]}"] = (float(m[4].replace(',', '')), int(m[5]))
        m = re.search(r'estimated total \d-node FFLs:\s+([\d,]+)', line)
        if m: total = float(m[1].replace(',', ''))
        m = re.search(r'FFLs among them=(\d+)', line)
        if m: sampled = int(m[1])
    return comp, total, sampled

rows = []
for k in (3, 4, 5, 6, 7):
    if k in (3, 5):
        comp, total, sampled = from_log(f'{REV}/logs/classdiv_k{k}.log')
        listed = sum(v[1] for v in comp.values())
        est = sum(v[0] for kk, v in comp.items() if qual(kk)); ns = sum(v[1] for kk, v in comp.items() if qual(kk))
        src = f'logs/classdiv_k{k}.log (Table 3 graph, 8,031 arcs; JSON overwritten 2026-09-12)'
        note = f'log lists the top 12 compositions ({listed} of {sampled} sampled modules)' if listed < sampled else 'log lists every composition'
    else:
        d = json.load(open(f'{REV}/results/ffl_census_k{k}.json'))
        total, sampled = d['est_ffl'], d['sampled_ffl']
        est = sum(v for kk, v in d['composition'].items() if qual(kk))
        ns = round(est / d['scale'])
        src = f'results/ffl_census_k{k}.json (03_ffl_census.py, Table 3 graph, 8,031 arcs)'
        note = f"RAND-ESU probs {d['probs']}, scale {d['scale']:.6g}"
    rows.append(dict(n=k, source=src, total_modules=total, sampled_modules=sampled,
                     modules_2TF_2miR_2Gene=est, sampled_2TF_2miR_2Gene=ns,
                     pct=100 * est / total, note=note))
for k in (3, 4, 5, 6):
    fn = f'val_agg_n{k}.json' if k < 6 else 'census_n6_exhaustive.json'
    d = json.load(open(fn)); v = d['flags_found']['a_2TF_2miR_2Gene']
    rows.append(dict(n=k, source=f'this re-run, exhaustive ({fn})', total_modules=d['found'], sampled_modules=d['found'],
                     modules_2TF_2miR_2Gene=v, sampled_2TF_2miR_2Gene=v, pct=100 * v / d['found'], note='exact count'))
s = json.load(open('census_n7_randesu_8seeds_summary.json'))
rows.append(dict(n=7, source='this re-run, RAND-ESU 8 seeds (census_n7_randesu_8seeds_summary.json)',
                 total_modules=s['mean']['est_total'], sampled_modules=sum(r['found'] for r in s['per_seed']),
                 modules_2TF_2miR_2Gene=s['mean']['est_a_2TF_2miR_2Gene'],
                 sampled_2TF_2miR_2Gene=sum(r['found_a_2TF_2miR_2Gene'] for r in s['per_seed']),
                 pct=s['mean']['pct_a_2TF_2miR_2Gene'],
                 note=f"mean of 8 seeds; SD {s['sd']['est_a_2TF_2miR_2Gene']:,.0f} modules, {s['sd']['pct_a_2TF_2miR_2Gene']:.2f} percentage points"))
with open('B2_2TF_2miR_2Gene_by_n.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    print(f"n={r['n']}  {r['modules_2TF_2miR_2Gene']:>14,.1f} of {r['total_modules']:>16,.1f}  ({r['pct']:.3f} %)  sampled {r['sampled_2TF_2miR_2Gene']}/{r['sampled_modules']}  | {r['source']} | {r['note']}")
