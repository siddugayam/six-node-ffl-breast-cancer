#!/usr/bin/env python3
"""P1a: COL1A1 expression in the controls of every candidate dataset (SETTINGS change 3: this check may use the NCBI
TPM table or the series matrix before raw data are fetched).
RNA-seq: mean TPM of COL1A1 over the control samples (NCBI-generated TPM, GRCh38.p13); pass if >= 10.
Arrays: series-matrix values; gene value = the COL1A1 probe with the highest mean over all samples of the contrast;
percentile of its control mean among all genes (same probe-to-gene rule); pass if >= 50.
Writes p1a_col1a1_screen.tsv (one row per dataset; optional arguments: arms file, output file) and records the data route (NCBI counts / series matrix)."""
import os, sys, csv, gzip, re, io, collections
import numpy as np, pandas as pd
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, INB)
import _geo, _dl, _arr
A = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else f'{INB}/P1a/p1a_arms.tsv', sep='\t')
annot = pd.read_csv(_geo.ncbi_annot('P1a'), sep='\t', usecols=['GeneID', 'Symbol'], dtype=str)
SYM = dict(zip(annot.GeneID, annot.Symbol))

def gpl_symbols(gpl):
    return _arr.gpl_symbols(gpl, 'P1a')

def series_matrix_values(gse, gsms):
    frames = []
    for p in _geo.series_matrix(gse, 'P1a'):
        lines = gzip.open(p, 'rt', errors='replace').read().split('\n')
        try: i = lines.index('!series_matrix_table_begin') + 1
        except ValueError: continue
        j = lines.index('!series_matrix_table_end')
        t = pd.read_csv(io.StringIO('\n'.join(lines[i:j])), sep='\t', index_col=0, low_memory=False)
        keep = [c for c in t.columns if c in gsms]
        if keep: frames.append(t[keep].apply(pd.to_numeric, errors='coerce'))
    return pd.concat(frames, axis=1) if frames else pd.DataFrame()

rows = []
for (gse, tf, cell), d in A.groupby(['gse', 'tf', 'cell'], sort=False):
    ctrl = d[d.arm == 'ctrl'].gsm.tolist(); pert = d[d.arm == 'pert'].gsm.tolist()
    rnaseq = (d.strategy == 'RNA-Seq').all()
    r = dict(gse=gse, tf=tf, cell=cell, cls=d.cls.iloc[0], n_pert=len(pert), n_ctrl=len(ctrl), platform=';'.join(sorted(set(d.platform))),
             data='', col1a1_ctrl='', col1a1_rule='', col1a1_pass='', note='')
    try:
        if rnaseq:
            ok, _ = _geo.ncbi_counts_available(gse, 'P1a')
            if not ok:
                r.update(data='RNA-seq: no NCBI-generated counts', col1a1_pass='NOT CHECKED', note='submitter files needed')
            else:
                tp = pd.read_csv(_geo.ncbi_counts(gse, 'P1a', 'norm_counts_TPM'), sep='\t', index_col=0, dtype={0: str})
                _geo.ncbi_counts(gse, 'P1a', 'raw_counts')
                tp.index = tp.index.astype(str)
                gid = next((g for g, s in SYM.items() if s == 'COL1A1'), None)
                v = float(tp.loc[gid, ctrl].mean()) if gid in tp.index else float('nan')
                r.update(data='RNA-seq: NCBI raw counts + TPM', col1a1_ctrl=f'{v:.3g} TPM', col1a1_rule='mean control TPM >= 10', col1a1_pass='yes' if v >= 10 else 'no')
        else:
            m = series_matrix_values(gse, set(ctrl + pert))
            if m.empty:
                r.update(data='array: no series-matrix values', col1a1_pass='NOT CHECKED'); rows.append(r); continue
            sym, how = gpl_symbols(d.platform.iloc[0])
            m.index = m.index.astype(str)
            if not sym: sym, how = _arr.symbol_map(m.index, d.platform.iloc[0], 'P1a')     # fallback: key column chosen by overlap
            s = m.index.map(lambda x: sym.get(x, ''))
            m = m[(s != '') & m.notna().all(axis=1)]; s = m.index.map(lambda x: sym.get(x, ''))
            if m.max().max() > 50: m = np.log2(m.clip(lower=0) + 1)
            mu = m.mean(axis=1); best = mu.groupby(s).idxmax()
            g = m.loc[best.values]; g.index = best.index
            cm = g[ctrl].mean(axis=1)
            if 'COL1A1' in cm.index:
                pct = 100 * (cm < cm['COL1A1']).mean()
                r.update(data=f'array: series matrix ({how})', col1a1_ctrl=f'percentile {pct:.1f}', col1a1_rule='control mean at or above the median gene', col1a1_pass='yes' if pct >= 50 else 'no')
            else:
                why = ('platform annotation could not be mapped to gene symbols' if not sym else
                       'COL1A1 probe on the platform but absent from the series-matrix values (filtered probe set)' if 'COL1A1' in set(sym.values()) else
                       'COL1A1 not on the platform annotation')
                r.update(data=f'array: series matrix ({how})', col1a1_pass='NOT FOUND', note=why)
    except Exception as e:
        r.update(col1a1_pass='ERROR', note=str(e)[:200])
    rows.append(r); print(gse, tf, cell[:18], r['data'][:40], r['col1a1_ctrl'], r['col1a1_pass'], flush=True)
pd.DataFrame(rows).to_csv(sys.argv[2] if len(sys.argv) > 2 else f'{INB}/P1a/p1a_col1a1_screen.tsv', sep='\t', index=False)
