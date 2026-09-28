#!/usr/bin/env python3
"""P3a search: the P1b search pattern applied to every miRNA node of the analysed network (223), Homo sapiens,
expression profiling; then the first-pass screen (samples naming a network miRNA near a perturbation word)."""
import os, sys, re
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); REV = os.path.dirname(os.path.dirname(INB)); sys.path.insert(0, INB)
import pandas as pd, _screen
N = pd.read_csv(f'{REV}/data/canonical_nodes.tsv', sep='\t'); MIRS = sorted(N[N.type == 'miRNA'].name)
ACT = '(mimic OR overexpression OR transfection OR inhibitor OR antagomir OR antagomiR OR knockout OR sponge OR precursor OR "pre-miR" OR agomir)'
Q = []
for m in MIRS:
    s = m.replace('hsa-', '')                         # e.g. miR-29a, let-7a
    alts = {s, s.replace('miR-', 'miR'), s.replace('miR-', 'microRNA-'), s.replace('miR-', 'mir-')}
    Q.append((f'P3a_{s}', '(' + ' OR '.join(f'"{x}"' for x in sorted(alts)) + f') AND {ACT}'))
hits = _screen.search('P3a', Q)
MIRRX = r'(mi-?R-?\d+[a-z]?|let-?7[a-z]?|microRNA-?\d+)'
PERT = r'(mimic|inhibit|anti|antago|agomi|over-?express|\bOE\b|transfect|precursor|\bpre-?mi|sponge|knock|\bKO\b|lenti)'
rows = _screen.screen('P3a', hits, r'(' + MIRRX + r'.{0,30}' + PERT + r'|' + PERT + r'.{0,30}' + MIRRX + r')',
                      r'(control|ctrl|scrambl|non[- ]?target|\bNT\b|\bNC\b|negative|mock|miR-?NC|cel-miR-67|luciferase|\bGFP\b|empty|vector)')
print(len(Q), 'queries;', len(hits), 'GSE hits; >=2 auto-perturbed and >=2 auto-control:', sum(1 for r in rows if r['auto_perturbed'] >= 2 and r['auto_control'] >= 2))
