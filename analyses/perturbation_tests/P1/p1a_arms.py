#!/usr/bin/env python3
"""P1a: assign each candidate's samples to the perturbed and control arms from the regexes of p1a_triage.tsv
(matched against the sample title, then characteristics), and write p1a_arms.tsv (one row per dataset and sample)."""
import os, sys, re, csv
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, INB)
import _geo
SRC = sys.argv[1] if len(sys.argv) > 1 else f'{INB}/P1a/p1a_triage.tsv'
DST = sys.argv[2] if len(sys.argv) > 2 else f'{INB}/P1a/p1a_arms.tsv'
T = list(csv.DictReader(open(SRC), delimiter='\t'))
rows = []
for t in T:
    if t['decision'] not in ('candidate', 'secondary'): continue
    S = _geo.samples(t['gse'], 'P1a')
    pr, cr = re.compile(t['pert_title_regex']), re.compile(t['ctrl_title_regex'])
    for s in S:
        txt = s['title']
        arm = 'pert' if pr.search(txt) else ('ctrl' if cr.search(txt) else '')
        if arm: rows.append(dict(gse=t['gse'], tf=t['tf'], gsm=s['gsm'], arm=arm, title=s['title'], platform=s['platform'], strategy=s['strategy'],
                                 cell=t['cell'], cls=t['class']))
with open(DST, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
import collections
c = collections.OrderedDict()
for r in rows: c.setdefault((r['gse'], r['tf'], r['cell']), collections.Counter())[r['arm']] += 1
for (g, tf, cell), v in c.items(): print(g, tf, cell[:20], dict(v), 'OK' if v['pert'] >= 2 and v['ctrl'] >= 2 else 'CHECK')
cand = [(t['gse'], t['tf'], t['cell']) for t in T if t['decision'] == 'candidate']
print('candidates without any assigned sample:', [x for x in cand if x not in c])
