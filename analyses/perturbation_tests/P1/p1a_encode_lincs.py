#!/usr/bin/env python3
"""P1a: (1) ENCODE portal search for knockdown/knockout RNA-seq of ETS1, NFKB1, RELA and SP1 (human);
(2) LINCS L1000 landmark status of COL1A1 and COL3A1 (GSE92742 gene table, pr_is_lm).  Writes p1a_encode_lincs.txt."""
import os, sys, json, gzip, csv, urllib.parse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _dl
OUT = os.path.dirname(os.path.abspath(__file__)); L = []
LIC_ENC = 'ENCODE: data released under the ENCODE data use policy (no restrictions on use; cite)'
for tf in ('ETS1', 'NFKB1', 'RELA', 'SP1'):
    q = {'type': 'Experiment', 'target.label': tf,
         'assay_title': ['shRNA RNA-seq', 'CRISPRi RNA-seq', 'siRNA RNA-seq', 'CRISPR RNA-seq'], 'format': 'json', 'limit': 'all'}
    url = 'https://www.encodeproject.org/search/?' + urllib.parse.urlencode(q, doseq=True)
    # ENCODE answers an empty search with HTTP 404 and a JSON body ("No results found"); keep that body as the record
    import urllib.request, urllib.error
    p = os.path.join(_dl.CACHE, 'P1/encode', f'encode_search_{tf}.json'); os.makedirs(os.path.dirname(p), exist_ok=True)
    if not os.path.exists(p):
        try: body = urllib.request.urlopen(urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'curl/8.5.0'}), timeout=60).read()
        except urllib.error.HTTPError as e: body = e.read()
        open(p, 'wb').write(body)
        _dl.register_existing(p, 'P1a', 'ENCODE portal', f'search {tf} knockdown RNA-seq', url, LIC_ENC)
    J = json.load(open(p)); g = J.get('@graph', [])
    if not g: L.append(f"ENCODE {tf}: portal answer: {J.get('notification', '')}")
    for e in g:
        L.append(f"ENCODE {tf}: {e['accession']} | {e.get('assay_title')} | {e.get('biosample_summary', '')} | status {e.get('status')}")
    if not g: L.append(f'ENCODE {tf}: no knockdown/knockout RNA-seq experiment (search: target.label={tf}, assay_title in shRNA/CRISPRi/siRNA/CRISPR RNA-seq)')
url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz'
p = _dl.fetch(url, 'P1/lincs', 'GSE92742_Broad_LINCS_gene_info.txt.gz', 'P1a', 'NCBI GEO (LINCS L1000, GSE92742)', 'GSE92742 gene info',
              'NCBI GEO: no NCBI restrictions (submitter terms may apply)')
rows = list(csv.DictReader(gzip.open(p, 'rt'), delimiter='\t'))
lm = {r['pr_gene_symbol']: r for r in rows}
L.append(f'L1000 gene table: {len(rows)} genes; landmarks (pr_is_lm = 1): {sum(1 for r in rows if r["pr_is_lm"] == "1")}')
for g in ('COL1A1', 'COL3A1', 'ETS1', 'NFKB1', 'RELA', 'SP1'):
    r = lm.get(g)
    L.append(f'L1000 {g}: ' + (f'pr_is_lm = {r["pr_is_lm"]}; pr_is_bing = {r.get("pr_is_bing")}' if r else 'not in the gene table'))
open(os.path.join(OUT, 'p1a_encode_lincs.txt'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
