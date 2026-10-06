#!/usr/bin/env python3
"""P1: differential expression for one contrast, from raw data where available (SETTINGS.md).
usage: p1_run_dataset.py <part> <arms.tsv> <gse> <tf> <cell> [--route ncbi|cel|matrix] [--tag primary|secondary] [--pseudo]
Routes:
  ncbi   RNA-seq: NCBI-generated raw counts (GeneID rows) -> _de.R counts (filterByExpr, TMM, voom); TPM from NCBI's table.
  cel    Affymetrix: per-GSM CEL files -> oligo RMA (_cel.R) over the arrays of the contrast -> _de.R log_matrix.
  matrix processed series-matrix values -> _de.R (log2 if max > 50); flagged 'processed'.
Output: <part>/de/<gse>__<tf>__<cell>__<tag>.tsv (one row per gene) and .meta.json (route, arm sizes); the
sample-level matrices and the R spec stay in the cache (data_cache/<part>/work/)."""
import os, sys, re, json, gzip, io, argparse, subprocess
import pandas as pd, numpy as np
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, INB)
import _geo, _arr, _dl
ap = argparse.ArgumentParser(); ap.add_argument('part'); ap.add_argument('arms'); ap.add_argument('gse'); ap.add_argument('tf'); ap.add_argument('cell')
ap.add_argument('--route', default='auto'); ap.add_argument('--tag', default='primary'); ap.add_argument('--pseudo', action='store_true')
ap.add_argument('--log10', action='store_true')
ap.add_argument('--force-keep', action='store_true', help='P1b sensitivity: keep COL1A1/COL3A1 even if filterByExpr drops them')
ap.add_argument('--mouse', action='store_true', help='mouse data: gene ids / symbols mapped to human symbols by one-to-one orthologues (_mouse.py)')
ap.add_argument('--matrix', default=None, help='counts_local: count matrix (first column gene symbol or Ensembl id; columns GSM ids)')
a = ap.parse_args()
A = pd.read_csv(a.arms, sep='\t'); d = A[(A.gse == a.gse) & (A.tf == a.tf) & (A.cell == a.cell)]
assert len(d), 'no arms'
slug = re.sub(r'[^A-Za-z0-9]+', '_', a.cell)[:40]
od = os.path.join(INB, a.part, 'de'); os.makedirs(od, exist_ok=True)
wd = os.path.join(_dl.CACHE, a.part, 'work'); os.makedirs(wd, exist_ok=True)     # sample-level matrices stay in the cache
name = f'{a.gse}__{re.sub(r"[^A-Za-z0-9.-]+", "_", a.tf)}__{slug}__{a.tag}'
base = os.path.join(wd, name); res = os.path.join(od, name)
samples = dict(zip(d.gsm, d.arm))
route = a.route
if route == 'auto':
    route = 'ncbi' if (d.strategy == 'RNA-Seq').all() else ('cel' if d.platform.iloc[0] in ('GPL570', 'GPL17586', 'GPL6244', 'GPL16686', 'GPL11532', 'GPL17692', 'GPL201', 'GPL96', 'GPL571', 'GPL23126') else 'matrix')
spec = dict(samples=samples)
if a.pseudo: spec['extra_pseudo'] = ['COL1A1', 'COL3A1']
if a.force_keep: spec['force_keep'] = ['COL1A1', 'COL3A1']
if route == 'ncbi':
    rc = _geo.ncbi_counts(a.gse, a.part, 'raw_counts'); tp = _geo.ncbi_counts(a.gse, a.part, 'norm_counts_TPM')
    gi = pd.read_csv(_geo.ncbi_annot(a.part), sep='\t', usecols=['GeneID', 'Symbol'], dtype=str)
    gm = base + '.genemap.tsv'; gi.rename(columns={'GeneID': 'id', 'Symbol': 'symbol'}).to_csv(gm, sep='\t', index=False)
    m = base + '.counts.tsv'
    c = pd.read_csv(rc, sep='\t', dtype={0: str})
    miss = [g for g in samples if g not in c.columns]              # a sample NCBI did not process is left out, and noted
    if miss:
        samples = {g: v for g, v in samples.items() if g not in miss}; d = d[~d.gsm.isin(miss)]; spec['samples'] = samples
        spec['note'] = f'{len(miss)} sample(s) absent from the NCBI-generated counts and left out: ' + ','.join(miss)
        assert (d.arm == 'pert').sum() >= 2 and (d.arm == 'ctrl').sum() >= 2, 'fewer than 2 samples per arm'
    c[[c.columns[0]] + list(samples)].to_csv(m, sep='\t', index=False)
    spec.update(kind='counts', matrix=m, gene_map=gm, tpm=tp)
elif route == 'counts_local':
    # submitter count matrix already re-labelled with GSM ids (P1/prep_counts.py); gene ids are symbols or Ensembl ids
    c = pd.read_csv(a.matrix, sep='\t', dtype={0: str}); ids = c.iloc[:, 0].astype(str)
    gm = base + '.genemap.tsv'
    if ids.str.fullmatch(r'\d+').mean() > 0.5:          # Entrez Gene ids
        gi = pd.read_csv(_geo.ncbi_annot(a.part), sep='\t', usecols=['GeneID', 'Symbol'], dtype=str)
        pd.DataFrame({'id': ids, 'symbol': ids.map(dict(zip(gi.GeneID, gi.Symbol)))}).to_csv(gm, sep='\t', index=False)
    elif ids.str.startswith('ENSG').mean() > 0.5:
        gi = pd.read_csv(_geo.ncbi_annot(a.part), sep='\t', usecols=['Symbol', 'dbXrefs'], dtype=str)
        ens = gi.dbXrefs.str.extract(r'Ensembl:(ENSG\d+)')[0]
        mp = dict(zip(ens.dropna(), gi.Symbol[ens.notna()]))
        pd.DataFrame({'id': ids, 'symbol': ids.str.replace(r'\.\d+$', '', regex=True).map(mp)}).to_csv(gm, sep='\t', index=False)
    else:
        pd.DataFrame({'id': ids, 'symbol': ids}).to_csv(gm, sep='\t', index=False)
    m = base + '.counts.tsv'; c[[c.columns[0]] + list(samples)].to_csv(m, sep='\t', index=False)
    spec.update(kind='counts', matrix=m, gene_map=gm); spec['note'] = 'submitter raw count matrix (no NCBI-generated counts); baseline as log-CPM percentile'
elif route == 'cel':
    m = base + '.rma.tsv'
    info = _arr.rma(list(samples), a.part, m)
    sym, how = _arr.symbol_map(pd.read_csv(m, sep='\t', usecols=[0]).iloc[:, 0].astype(str), d.platform.iloc[0], a.part)
    gm = base + '.genemap.tsv'; pd.DataFrame({'id': list(sym), 'symbol': list(sym.values())}).to_csv(gm, sep='\t', index=False)
    spec.update(kind='log_matrix', matrix=m, gene_map=gm); spec['note'] = f'RMA: {info}; annotation {how}'
elif route in ('agilent1', 'agilent2', 'matrix2'):
    part_dir = f'{a.part}/agilent/{a.gse}'
    def fe_file(gsm):
        stub = re.sub(r'\d{3}$', 'nnn', gsm); base_url = f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{stub}/{gsm}/suppl/'
        import urllib.request
        html = urllib.request.urlopen(urllib.request.Request(base_url, headers={'User-Agent': 'curl/8.5.0'}), timeout=60).read().decode()
        n = [x for x in re.findall(r'href="([^"?/][^"]*)"', html) if x.endswith('.txt.gz')]
        if not n: raise RuntimeError(f'no Feature Extraction file for {gsm}')
        return _dl.fetch(base_url + n[0], part_dir, n[0], a.part, 'NCBI GEO', gsm, _geo.LIC_GEO)
    gm = base + '.genemap.tsv'
    def write_map(ids):
        sym, how = _arr.symbol_map(ids, d.platform.iloc[0], a.part)
        pd.DataFrame({'id': list(sym), 'symbol': list(sym.values())}).to_csv(gm, sep='\t', index=False); return how
    if route == 'agilent1':
        m = base + '.agilent.tsv'
        r0 = subprocess.run(['Rscript', os.path.join(INB, '_agilent.R'), 'single', m] + [fe_file(g) for g in samples], capture_output=True, text=True)
        if r0.returncode: sys.exit('agilent failed: ' + r0.stderr[-800:])
        how = write_map(pd.read_csv(m, sep='\t', usecols=[0]).iloc[:, 0].astype(str))
        spec.update(kind='log_matrix', matrix=m, gene_map=gm); spec['note'] = f'Agilent raw (normexp + quantile): {r0.stdout.strip()}; annotation {how}'
    else:
        # two-colour: orientation from the channel annotations (the channel naming the perturbation is the numerator)
        soft = open(os.path.join(_dl.CACHE, a.part, 'geo_samples', f'{a.gse}_samples_brief.soft'), encoding='utf-8', errors='replace').read()
        orient = {}
        for blk in soft.split('^SAMPLE = ')[1:]:
            g = blk.split('\n', 1)[0].strip()
            if g not in samples: continue
            ch1 = ' '.join(re.findall(r'!Sample_(?:source_name|characteristics)_ch1 = (.*)', blk)); ch2 = ' '.join(re.findall(r'!Sample_(?:source_name|characteristics)_ch2 = (.*)', blk))
            kw = re.compile(r'miR-?29|miR-?101|mimic|pre-?miR|inhibitor', re.I); ctl = re.compile(r'control|NTC|negative|scram|miR-C\b', re.I)
            orient[g] = 1 if (kw.search(ch1) and not ctl.search(ch1)) else (-1 if (kw.search(ch2) and not ctl.search(ch2)) else 0)
        if 0 in orient.values(): sys.exit(f'orientation not determined: {orient}')
        spec2 = dict(gene_map=gm, orient=orient)
        if a.pseudo: spec2['extra_pseudo'] = ['COL1A1', 'COL3A1']
        if route == 'agilent2':
            r0 = subprocess.run(['Rscript', os.path.join(INB, '_agilent.R'), 'two', base] + [fe_file(g) for g in samples], capture_output=True, text=True)
            if r0.returncode: sys.exit('agilent failed: ' + r0.stderr[-800:])
            spec2.update(M=base + '.M.tsv', A=base + '.A.tsv'); note = f'Agilent two-colour raw (normexp, loess, Aquantile): {r0.stdout.strip()}'
            how = write_map(pd.read_csv(base + '.M.tsv', sep='\t', usecols=[0]).iloc[:, 0].astype(str))
        else:
            frames = []
            for p in _geo.series_matrix(a.gse, a.part):
                lines = gzip.open(p, 'rt', errors='replace').read().split('\n')
                i = lines.index('!series_matrix_table_begin') + 1; j = lines.index('!series_matrix_table_end')
                t = pd.read_csv(io.StringIO('\n'.join(lines[i:j])), sep='\t', index_col=0, low_memory=False)
                keep = [c for c in t.columns if c in samples]
                if keep: frames.append(t[keep].apply(pd.to_numeric, errors='coerce'))
            mat = pd.concat(frames, axis=1)[list(samples)]
            if a.log10: mat = mat * np.log2(10)
            mp = base + '.M.tsv'; mat.reset_index().to_csv(mp, sep='\t', index=False)
            how = write_map(mat.index.astype(str))
            spec2.update(M=mp, A=None); note = 'processed series-matrix log ratios (two-colour)' + (' converted from log10' if a.log10 else '')
        sp = base + '.spec.json'; json.dump(spec2, open(sp, 'w'), indent=1)
        r = subprocess.run(['Rscript', os.path.join(INB, '_de2.R'), sp, res + '.tsv'], capture_output=True, text=True)
        if r.returncode: sys.exit('DE failed: ' + r.stderr[-1500:])
        meta = dict(gse=a.gse, tf=a.tf, cell=a.cell, cls=d.cls.iloc[0], tag=a.tag, route=route, n_pert=len(samples), n_ctrl=len(samples),
                    platform=d.platform.iloc[0], note=f'{note}; orientation {orient}; annotation {how}', de=r.stdout.strip(),
                    direction=d.direction.iloc[0] if 'direction' in d else '', mirna=d.mirna.iloc[0] if 'mirna' in d else '')
        json.dump(meta, open(res + '.meta.json', 'w')); print(json.dumps(meta)); sys.exit(0)
else:
    frames = []
    for p in _geo.series_matrix(a.gse, a.part):
        lines = gzip.open(p, 'rt', errors='replace').read().split('\n')
        i = lines.index('!series_matrix_table_begin') + 1; j = lines.index('!series_matrix_table_end')
        t = pd.read_csv(io.StringIO('\n'.join(lines[i:j])), sep='\t', index_col=0, low_memory=False)
        keep = [c for c in t.columns if c in samples]
        if keep: frames.append(t[keep])
    mat = pd.concat(frames, axis=1)[list(samples)]
    m = base + '.matrix.tsv'; mat.reset_index().to_csv(m, sep='\t', index=False)
    sym, how = _arr.symbol_map(mat.index.astype(str), d.platform.iloc[0], a.part)
    gm = base + '.genemap.tsv'; pd.DataFrame({'id': list(sym), 'symbol': list(sym.values())}).to_csv(gm, sep='\t', index=False)
    spec.update(kind='log_matrix', matrix=m, gene_map=gm); spec['note'] = f'processed series-matrix values; annotation {how}'
if a.mouse:                                           # mouse -> human one-to-one orthologues; unmapped genes are dropped
    import _mouse
    g = pd.read_csv(spec['gene_map'], sep='\t', dtype=str).fillna('')
    g['symbol'] = _mouse.to_human(g['id'], g['symbol']); g = g[g.symbol != '']
    g.to_csv(spec['gene_map'], sep='\t', index=False)
    spec['note'] = (spec.get('note', '') + '; mouse -> human one-to-one orthologues (MGI)').lstrip('; ')
sp = base + '.spec.json'; json.dump(spec, open(sp, 'w'), indent=1)
r = subprocess.run(['Rscript', os.path.join(INB, '_de.R'), sp, res + '.tsv'], capture_output=True, text=True)
if r.returncode: sys.exit('DE failed: ' + r.stderr[-1500:])
meta = dict(gse=a.gse, tf=a.tf, cell=a.cell, cls=d.cls.iloc[0], tag=a.tag, route=route, n_pert=int((d.arm == 'pert').sum()),
            n_ctrl=int((d.arm == 'ctrl').sum()), platform=d.platform.iloc[0], note=spec.get('note', ''), de=r.stdout.strip(),
            direction=d.direction.iloc[0] if 'direction' in d else '', mirna=d.mirna.iloc[0] if 'mirna' in d else '')
json.dump(meta, open(res + '.meta.json', 'w'))
print(json.dumps(meta))
