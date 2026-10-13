#!/usr/bin/env python3
"""GEO helpers for analyses/perturbation_tests.  Searches go through NCBI E-utilities (db=gds).  Every response saved to disk goes
through _dl.fetch, so it is recorded in DOWNLOADS.tsv.  Sample metadata come from GEO's brief SOFT view of the series'
samples (acc.cgi?targ=gsm&view=brief&form=text).  NCBI-generated RNA-seq counts are looked up on the GEO download page."""
import os, re, json, time, urllib.parse, urllib.request
import _dl

EU = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
LIC_GEO = 'NCBI GEO: NCBI places no restrictions on use or distribution; submitters may claim rights (GEO disclaimer)'
LIC_EU = 'NCBI E-utilities response (public metadata; NCBI disclaimer)'
_last = [0.0]

def _pace():
    dt = time.time() - _last[0]
    if dt < 0.4: time.sleep(0.4 - dt)
    _last[0] = time.time()

def esearch(term, part, tag, retmax=5000):
    _pace()
    url = EU + 'esearch.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'term': term, 'retmax': retmax, 'retmode': 'json'})
    p = _dl.fetch(url, f'{part}/geo_search', f'esearch_{tag}.json', part, 'NCBI E-utilities', f'esearch {tag}', LIC_EU)
    r = json.load(open(p))['esearchresult']
    return int(r['count']), r['idlist']

def esummary(ids, part, tag):
    out = {}
    for k in range(0, len(ids), 200):
        _pace()
        url = EU + 'esummary.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'id': ','.join(ids[k:k + 200]), 'retmode': 'json'})
        p = _dl.fetch(url, f'{part}/geo_search', f'esummary_{tag}_{k}.json', part, 'NCBI E-utilities', f'esummary {tag} {k}', LIC_EU)
        res = json.load(open(p))['result']
        for i in res.get('uids', []): out[i] = res[i]
    return out

def samples(gse, part):
    """All samples of a series: list of dicts (gsm, title, source, characteristics, description, platform, strategy)."""
    _pace()
    url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=gsm&view=brief&form=text'
    p = _dl.fetch(url, f'{part}/geo_samples', f'{gse}_samples_brief.soft', part, 'NCBI GEO', gse, LIC_GEO)
    S, cur = [], None
    for line in open(p, encoding='utf-8', errors='replace'):
        line = line.rstrip('\n')
        if line.startswith('^SAMPLE'):
            cur = dict(gsm=line.split('=', 1)[1].strip(), title='', source='', characteristics=[], description='', platform='',
                       strategy='', molecule='', organism=''); S.append(cur); continue
        if cur is None or ' = ' not in line: continue
        k, v = line.split(' = ', 1)
        if k == '!Sample_title': cur['title'] = v
        elif k == '!Sample_source_name_ch1': cur['source'] = v
        elif k.startswith('!Sample_characteristics_ch1'): cur['characteristics'].append(v)
        elif k == '!Sample_description': cur['description'] += (' | ' if cur['description'] else '') + v
        elif k == '!Sample_platform_id': cur['platform'] = v
        elif k == '!Sample_library_strategy': cur['strategy'] = v
        elif k == '!Sample_molecule_ch1': cur['molecule'] = v
        elif k == '!Sample_organism_ch1': cur['organism'] = v
    return S

def ncbi_counts_available(gse, part):
    """True if GEO lists NCBI-generated RNA-seq raw counts for the series."""
    _pace()
    url = f'https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse}'
    p = _dl.fetch(url, f'{part}/geo_download_pages', f'{gse}_download_page.html', part, 'NCBI GEO', gse, LIC_GEO)
    t = open(p, encoding='utf-8', errors='replace').read()
    return f'{gse}_raw_counts_GRCh38.p13_NCBI.tsv.gz' in t, t

def ncbi_counts(gse, part, kind='raw_counts'):
    """kind: raw_counts or norm_counts_TPM."""
    fn = f'{gse}_{kind}_GRCh38.p13_NCBI.tsv.gz'
    url = f'https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc={gse}&format=file&file={fn}'
    return _dl.fetch(url, f'{part}/geo_counts', fn, part, 'NCBI GEO (NCBI-generated RNA-seq counts)', gse, LIC_GEO)

def ncbi_annot(part):
    """GeneID -> Symbol for the NCBI-generated counts.  GEO's own annotation file (Human.GRCh38.p13.annot.tsv.gz) is
    served only behind a reCAPTCHA page, so NCBI Gene's gene_info table (FTP) is used instead (columns GeneID, Symbol)."""
    fn = 'Homo_sapiens.gene_info.gz'
    url = f'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/{fn}'
    return _dl.fetch(url, 'common', fn, part, 'NCBI Gene (gene_info)', 'Homo_sapiens.gene_info',
                     'NCBI: no restrictions on use or distribution (NCBI website and data usage policies)')

def series_matrix(gse, part):
    """Series matrix file(s) of a series (processed values).  Returns local paths."""
    stub = re.sub(r'\d{3}$', 'nnn', gse)
    base = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/'
    _pace()
    html = urllib.request.urlopen(urllib.request.Request(base, headers={'User-Agent': 'Mozilla/5.0'}), timeout=60).read().decode()
    names = sorted(set(re.findall(r'href="([^"]+_series_matrix\.txt\.gz)"', html)))
    return [_dl.fetch(base + n, f'{part}/geo_matrix', n, part, 'NCBI GEO', gse, LIC_GEO) for n in names]

def suppl_list(gse):
    stub = re.sub(r'\d{3}$', 'nnn', gse)
    base = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/'
    _pace()
    try:
        html = urllib.request.urlopen(urllib.request.Request(base, headers={'User-Agent': 'Mozilla/5.0'}), timeout=60).read().decode()
    except Exception:
        return base, []
    items = re.findall(r'href="([^"?/][^"]*)"[^>]*>[^<]*</a>\s+(\S+ \S+)\s+(\S+)', html)
    return base, [(n, sz) for n, _, sz in items]

def suppl_fetch(gse, name, part):
    stub = re.sub(r'\d{3}$', 'nnn', gse)
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/{name}'
    return _dl.fetch(url, f'{part}/geo_suppl/{gse}', name, part, 'NCBI GEO', gse, LIC_GEO)
