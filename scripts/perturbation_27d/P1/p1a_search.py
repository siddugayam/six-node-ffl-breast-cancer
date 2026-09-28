#!/usr/bin/env python3
"""P1a search: GEO (the analysis plan's terms per TF) plus the KnockTF 2.0 series for the four TFs, then the first-pass screen."""
import os, sys, csv
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _screen
PERT = r'(knock|\bKD\b|\bKO\b|\bsi[-_ ]?|\bsh[-_ ]?|\bsg[-_ ]?|CRISPR|deplet|silenc|-/-|null)'
TFS = {'ETS1': 'ETS1', 'NFKB1': 'NFKB1', 'RELA': '(RELA OR p65)', 'SP1': 'SP1'}
Q = [(f'P1a_{t}', f'{q} AND (knockdown OR siRNA OR shRNA OR knockout OR CRISPR OR CRISPRi)') for t, q in TFS.items()]
hits = _screen.search('P1a', Q)
kt = [r for r in csv.DictReader(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'knocktf_datasets.tsv')), delimiter='\t')
      if r['tf'].upper() in TFS and r['profile_id'].startswith('GSE')]
for r in kt:
    h = hits.setdefault(r['profile_id'], dict(gse=r['profile_id'], title='', summary='', gdstype='', n_samples='', gpl=r['platform'], taxon='', pubmed=r['pubmed'], tags=set()))
    h['tags'].add(f"KnockTF_{r['tf'].upper()}")
TFRX = r'(ETS1|ETS-1|NFKB1|NF-?kB1|p50|p105|RELA|p65|SP1)'
rows = _screen.screen('P1a', hits, r'(' + PERT + r'.{0,25}' + TFRX + r'|' + TFRX + r'.{0,15}' + PERT + r')',
                      r'(control|ctrl|scrambl|non[- ]?target|\bNT\b|\bNTC\b|mock|siNC|shNC|sgNC|luciferase|\bGFP\b|empty|negative)')
print(len(hits), 'GSE hits; with >=2 auto-perturbed and >=2 auto-control samples:', sum(1 for r in rows if r['auto_perturbed'] >= 2 and r['auto_control'] >= 2))
