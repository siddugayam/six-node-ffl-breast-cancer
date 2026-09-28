#!/usr/bin/env python3
"""P1b search: GEO (the analysis plan's terms: miR-29 mimic / overexpression / transfection / inhibitor / antagomir / knockout;
and miR-101 gain for the positive control), Homo sapiens, expression profiling; then the first-pass screen."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _screen
M29 = '("miR-29" OR "miR-29a" OR "miR-29b" OR "miR-29c" OR "mir-29" OR "microRNA-29" OR miR29 OR miR29a OR miR29b OR miR29c)'
M101 = '("miR-101" OR "microRNA-101" OR miR101 OR "miR-101-3p")'
ACT = '(mimic OR overexpression OR transfection OR inhibitor OR antagomir OR antagomiR OR knockout OR sponge OR precursor OR "pre-miR" OR agomir)'
Q = [('P1b_miR29', f'{M29} AND {ACT}'), ('P1b_miR101', f'{M101} AND (mimic OR overexpression OR transfection OR precursor OR "pre-miR" OR agomir)')]
hits = _screen.search('P1b', Q)
MIR = r'(mi-?R-?(29|101)|mir-?(29|101)|microRNA-?(29|101)|pre-?(29|101))'
PERT = r'(mimic|inhibit|anti|antago|agomi|over-?express|\bOE\b|transfect|precursor|\bpre\b|sponge|knock|\bKO\b|lenti)'
rows = _screen.screen('P1b', hits, r'(' + MIR + r'.{0,30}' + PERT + r'|' + PERT + r'.{0,30}' + MIR + r'|' + MIR + r'[a-c]?(-3p|-5p)?\b)',
                      r'(control|ctrl|scrambl|non[- ]?target|\bNT\b|\bNC\b|negative|mock|miR-?NC|NC mimic|cel-miR-67|luciferase|\bGFP\b|empty|untreated|vector)')
print(len(hits), 'GSE hits; >=2 auto-perturbed and >=2 auto-control:', sum(1 for r in rows if r['auto_perturbed'] >= 2 and r['auto_control'] >= 2))
