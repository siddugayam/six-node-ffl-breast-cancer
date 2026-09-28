#!/usr/bin/env python3
"""Build a GSM-labelled count matrix from a submitter's GEO supplementary count file.
usage: prep_counts.py <part> <gse> <suppl file name> <out.tsv> [--sep ,] [--title-key bracket|title|regex:<pattern>]
Columns are matched to GSMs by the bracketed library name in the GSM title (e.g. '[M45NC1]'), or by the whole title."""
import os, sys, re, gzip, argparse
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); import _geo
ap = argparse.ArgumentParser(); ap.add_argument('part'); ap.add_argument('gse'); ap.add_argument('file'); ap.add_argument('out')
ap.add_argument('--sep', default=None); ap.add_argument('--title-key', default='bracket'); ap.add_argument('--skip', type=int, default=0)
ap.add_argument('--gene-col', default=None)
a = ap.parse_args()
p = _geo.suppl_fetch(a.gse, a.file, a.part)
t = pd.read_csv(p, sep=a.sep, engine='python' if a.sep is None else 'c', skiprows=a.skip, comment=None)
gc = a.gene_col or t.columns[0]
S = _geo.samples(a.gse, a.part); col2gsm = {}
for s in S:
    if a.title_key == 'bracket':
        m = re.search(r'\[([^\]]+)\]', s['title']); key = m.group(1) if m else None
    elif a.title_key.startswith('regex:'):
        m = re.search(a.title_key[6:], s['title']); key = m.group(1) if m else None
    else: key = s['title']
    if key is None: continue
    for c in t.columns:
        if c == key or c.split('/')[-1].replace('.bam', '').replace('Aligned.sortedByCoord.out', '').strip('._') == key: col2gsm[c] = s['gsm']
out = t[[gc] + list(col2gsm)].rename(columns=col2gsm).rename(columns={gc: 'gene'})
out.to_csv(a.out, sep='\t', index=False)
print('matched', len(col2gsm), 'columns of', len(t.columns) - 1, '; genes', len(out))
