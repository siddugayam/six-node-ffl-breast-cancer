#!/usr/bin/env python3
"""Array helpers for analyses/perturbation_tests: GPL probe -> gene symbol, and raw CEL download plus RMA (via _cel.R)."""
import os, re, io, gzip, subprocess, urllib.request
import pandas as pd
import _dl, _geo
INB = os.path.dirname(os.path.abspath(__file__))

def gpl_symbols(gpl, part):
    """probe id -> gene symbol from the GPL's annotation (GEO annot.gz if present, else the full SOFT table)."""
    stub = re.sub(r'\d{3}$', 'nnn', gpl)
    # the family SOFT is never used: it carries every sample of every series on the platform (4 GB for GPL17586);
    # the platform-only view (acc.cgi view=data) has the same annotation table
    for url, fn, kind in ((f'https://ftp.ncbi.nlm.nih.gov/geo/platforms/{stub}/{gpl}/annot/{gpl}.annot.gz', f'{gpl}.annot.gz', 'annot'),
                          (f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=data', f'{gpl}_platform_table.soft', 'platform table')):
        try: p = _dl.fetch(url, 'common/gpl', fn, part, 'NCBI GEO', gpl, _geo.LIC_GEO)
        except RuntimeError: continue
        op = gzip.open if fn.endswith('.gz') else open
        lines = op(p, 'rt', errors='replace').read().split('\n')
        try: i = next(k for k, l in enumerate(lines) if l.startswith('!platform_table_begin') or l.startswith('!Annotation_table_begin')) + 1
        except StopIteration: continue
        j = next(k for k in range(i, len(lines)) if lines[k].startswith('!platform_table_end') or lines[k].startswith('!Annotation_table_end'))
        t = pd.read_csv(io.StringIO('\n'.join(lines[i:j])), sep='\t', dtype=str, low_memory=False)
        cols = {c.lower(): c for c in t.columns}
        for c in ('gene symbol', 'gene_symbol', 'symbol', 'ilmn_gene', 'genesymbol'):
            if c in cols:
                s = t[cols[c]].fillna('').str.split(r'\s*///\s*|\s*//\s*|;').str[0].str.strip()
                return dict(zip(t[t.columns[0]].astype(str), s)), f'{kind}:{cols[c]}'
        if 'gene_assignment' in cols:
            s = t[cols['gene_assignment']].fillna('').str.split('//').str[1].fillna('').str.strip()
            return dict(zip(t[t.columns[0]].astype(str), s)), f'{kind}:gene_assignment'
    return {}, 'no annotation'

def gsm_cel(gsm, part):
    stub = re.sub(r'\d{3}$', 'nnn', gsm)
    base = f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{stub}/{gsm}/suppl/'
    html = urllib.request.urlopen(urllib.request.Request(base, headers={'User-Agent': 'curl/8.5.0'}), timeout=60).read().decode()
    names = [n for n in re.findall(r'href="([^"]+)"', html) if re.search(r'\.cel(\.gz)?$', n, re.I)]
    if not names: raise RuntimeError(f'no CEL file for {gsm}')
    return _dl.fetch(base + names[0], f'{part}/cel/{gsm[:7]}', names[0], part, 'NCBI GEO', gsm, _geo.LIC_GEO)

def rma(gsms, part, out):
    cels = [gsm_cel(g, part) for g in gsms]
    r = subprocess.run(['Rscript', os.path.join(INB, '_cel.R'), out] + cels, capture_output=True, text=True)
    if r.returncode: raise RuntimeError(r.stderr[-800:])
    return r.stdout.strip().splitlines()[-1]          # the summary line only (the rest lists local file paths)

def gpl_table(gpl, part):
    """The GPL annotation table (GEO annot.gz if present, else the platform-only table)."""
    stub = re.sub(r'\d{3}$', 'nnn', gpl)
    for url, fn in ((f'https://ftp.ncbi.nlm.nih.gov/geo/platforms/{stub}/{gpl}/annot/{gpl}.annot.gz', f'{gpl}.annot.gz'),
                    (f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=data', f'{gpl}_platform_table.soft')):
        try: p = _dl.fetch(url, 'common/gpl', fn, part, 'NCBI GEO', gpl, _geo.LIC_GEO)
        except RuntimeError: continue
        op = gzip.open if fn.endswith('.gz') else open
        lines = op(p, 'rt', errors='replace').read().split('\n')
        try: i = next(k for k, l in enumerate(lines) if l.startswith('!platform_table_begin') or l.startswith('!Annotation_table_begin')) + 1
        except StopIteration: continue
        j = next(k for k in range(i, len(lines)) if lines[k].startswith('!platform_table_end') or lines[k].startswith('!Annotation_table_end'))
        return pd.read_csv(io.StringIO('\n'.join(lines[i:j])), sep='\t', dtype=str, low_memory=False), fn
    return None, 'no annotation'

ANNODB = {'GPL23159': 'clariomshumantranscriptcluster.db', 'GPL16686': 'hugene20sttranscriptcluster.db'}

def symbol_map(ids, gpl, part):
    """id -> gene symbol, keyed on the GPL column that shares the most values with the matrix row ids."""
    t, how = gpl_table(gpl, part)
    if t is None: return {}, how
    ids = set(map(str, ids))
    key = max(t.columns, key=lambda c: len(ids & set(t[c].dropna().astype(str))))
    cols = {c.lower(): c for c in t.columns}
    sym = None
    for c in ('gene symbol', 'gene_symbol', 'symbol', 'ilmn_gene', 'genesymbol'):
        if c in cols: sym = t[cols[c]].fillna('').str.split(r'\s*///\s*|\s*//\s*|;').str[0].str.strip(); how += f' key {key}, symbol {cols[c]}'; break
    if sym is None and 'gene_assignment' in cols:
        sym = t[cols['gene_assignment']].fillna('').str.split('//').str[1].fillna('').str.strip(); how += f' key {key}, symbol gene_assignment'
    if sym is None:
        for c in ('name', 'spot_id', 'description'):
            if c in cols and t[cols[c]].fillna('').str.contains(';').mean() > 0.3:      # "description; SYMBOL; accessions"
                sym = t[cols[c]].fillna('').str.split(';').str[1].fillna('').str.strip(); how += f' key {key}, symbol 2nd field of {cols[c]}'; break
    if sym is None and gpl in ANNODB:                     # platform table without symbols: Bioconductor annotation package
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            open(f'{td}/ids.txt', 'w').write('\n'.join(sorted(ids)))
            r = subprocess.run(['Rscript', os.path.join(INB, '_annodb.R'), ANNODB[gpl], f'{td}/ids.txt', f'{td}/map.tsv'], capture_output=True, text=True)
            if r.returncode: return {}, how + ' (annotation package failed)'
            mp = pd.read_csv(f'{td}/map.tsv', sep='\t', dtype=str).fillna('')
        return dict(zip(mp.id, mp.symbol)), f'Bioconductor {ANNODB[gpl]} ({r.stdout.strip()})'
    if sym is None: return {}, how + ' (no symbol column)'
    return dict(zip(t[key].astype(str), sym)), how
