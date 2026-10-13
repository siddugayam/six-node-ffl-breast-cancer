#!/usr/bin/env python3
"""P3a: GSM-labelled count matrices from the submitters' files for the P3a series without NCBI-generated counts.
GSE246676: counts.xlsx (sheet exprFiltered), columns named '<n>_<sample title with ', rep' -> '_rep'>'.
GSE252885: count_matrix.csv, columns miR_122_0k_count / nc_0k_count = the 'duplication k' sample of each arm.
GSE196161: per-sample CLC 'GE' tables in the RAW archive (GSM in the file name); count = 'Total gene reads'.
GSE277200: the submitter's DESeq2 tables carry per-sample raw counts ('SampleN_count'); in each table the first block of
           samples is the control group and the second the MIR155 group (their means equal the tables' group-mean
           columns, checked below); columns are assigned to the group's GSMs in GEO order (identity within a group does
           not enter a two-group model).
Writes <cache>/P3/work/<GSE>[_<stage>]_counts.tsv (first column gene id)."""
import os, sys, io, re, gzip, tarfile
import numpy as np, pandas as pd
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, INB)
import _geo, _dl
W = os.path.join(_dl.CACHE, 'P3', 'work'); os.makedirs(W, exist_ok=True)
def out(name, M):
    p = os.path.join(W, f'{name}_counts.tsv'); M.to_csv(p, sep='\t', index=False); print(name, M.shape); return p
S = {s['title']: s['gsm'] for s in _geo.samples('GSE246676', 'P3')}
x = pd.read_excel(_geo.suppl_fetch('GSE246676', 'GSE246676_counts.xlsx', 'P3'), sheet_name='exprFiltered')
ren = {}
for c in x.columns[1:]:
    t = re.sub(r'^\d+_', '', c); t = re.sub(r'_rep(\d+)$', r', rep\1', t)
    if t in S: ren[c] = S[t]
assert len(ren) == 30, len(ren)
out('GSE246676', x.rename(columns=ren).rename(columns={'Name': 'gene'})[['gene'] + list(ren.values())])
S = _geo.samples('GSE252885', 'P3'); t = pd.read_csv(_geo.suppl_fetch('GSE252885', 'GSE252885_count_matrix.csv.gz', 'P3'))
ren = {}
for s in S:
    k = int(re.search(r'duplication (\d+)', s['title']).group(1)); arm = 'miR_122' if 'miR-122-5p' in s['title'] else 'nc'
    ren[f'{arm}_{k:02d}_count'] = s['gsm']
assert len(ren) == 12
out('GSE252885', t.rename(columns=ren).rename(columns={'gene_id': 'gene'})[['gene'] + list(ren.values())])
tf = tarfile.open(_geo.suppl_fetch('GSE196161', 'GSE196161_RAW.tar', 'P3')); M = None
for m in tf.getmembers():
    raw = tf.extractfile(m).read(); raw = gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw
    d = pd.read_csv(io.BytesIO(raw), sep='\t', usecols=['Name', 'Total gene reads'])
    s = d.groupby('Name')['Total gene reads'].sum().round().astype(int).rename(m.name.split('_')[0])
    M = s.to_frame() if M is None else M.join(s, how='outer')
M.index.name = 'gene'; out('GSE196161', M.fillna(0).astype(int).reset_index())
S = _geo.samples('GSE277200', 'P3')
for stage, fn in (('NSC', 'GSE277200_CONTR_NSCvsMIR155_NSC_deg.xls.gz'), ('NEURO', 'GSE277200_CONTR_NEUROvsMIR155_NEURO_deg.xls.gz')):
    t = pd.read_csv(io.BytesIO(gzip.decompress(open(_geo.suppl_fetch('GSE277200', fn, 'P3'), 'rb').read())), sep='\t')
    cols = [c for c in t.columns[1:] if re.fullmatch(r'Sample\d+', c)]
    ctl, mir = f'CONTR_{stage}', f'MIR155_{stage}'
    n_ctl = sum(1 for s in S if s['title'].startswith(ctl)); a, b = cols[:n_ctl], cols[n_ctl:]
    assert np.allclose(t[a].mean(axis=1), t[ctl]) and np.allclose(t[b].mean(axis=1), t[mir]), stage
    g_ctl = [s['gsm'] for s in S if s['title'].startswith(ctl)]; g_mir = [s['gsm'] for s in S if s['title'].startswith(mir)]
    assert len(g_ctl) == len(a) and len(g_mir) == len(b), (stage, len(g_ctl), len(a), len(g_mir), len(b))
    C = t[['gene_id'] + [f'{c}_count' for c in a + b]].copy(); C.columns = ['gene'] + g_ctl + g_mir
    out(f'GSE277200_{stage}', C)
