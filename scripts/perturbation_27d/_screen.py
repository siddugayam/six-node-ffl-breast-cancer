#!/usr/bin/env python3
"""Search and first-pass screen shared by P1a, P1b, P2 and P3a.
search(part, queries): runs each GEO query (db=gds, GSE, Homo sapiens unless stated, expression profiling), logs the
query, date and count to <part>/<part>_search_log.tsv, and returns {gse: record}.
screen(part, hits, pert_rx, ctrl_rx): fetches every hit's sample metadata and counts samples whose text (title, source,
characteristics, description) matches the perturbation pattern or the control pattern.  The counts only rank the hits
for reading; inclusion is decided by reading each series (DATASETS.tsv gives the reason for every hit)."""
import os, re, csv, json, time, datetime
import _dl, _geo

INB = os.path.dirname(os.path.abspath(__file__))
EXPR = '("expression profiling by array"[DataSet Type] OR "expression profiling by high throughput sequencing"[DataSet Type])'

def search(part, queries, organism='Homo sapiens'):
    d = os.path.join(INB, part.split('_')[0]); os.makedirs(d, exist_ok=True)
    logp = os.path.join(d, f'{part}_search_log.tsv'); new = not os.path.exists(logp)
    hits = {}
    with open(logp, 'a', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        if new: w.writerow(['date_utc', 'tag', 'query', 'count'])
        for tag, q in queries:
            term = f'({q}) AND "{organism}"[Organism] AND "gse"[Entry Type] AND {EXPR}'
            n, ids = _geo.esearch(term, part, tag)
            w.writerow([datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), tag, term, n])
            summ = _geo.esummary(ids, part, tag) if ids else {}
            for i, r in summ.items():
                acc = r.get('accession', '')
                if not acc.startswith('GSE'): continue
                h = hits.setdefault(acc, dict(gse=acc, title=r.get('title', ''), summary=r.get('summary', ''), gdstype=r.get('gdstype', ''),
                                              n_samples=r.get('n_samples', ''), gpl='GPL' + str(r.get('gpl', '')).replace(';', ';GPL'),
                                              taxon=r.get('taxon', ''), pubmed=';'.join(str(x) for x in r.get('pubmedids', [])), tags=set()))
                h['tags'].add(tag)
    return hits

def sample_text(s):
    return ' '.join([s['title'], s['source'], ' '.join(s['characteristics']), s['description']])

def cell_fields(S):
    out = set()
    for s in S:
        for c in s['characteristics']:
            if re.match(r'\s*(cell ?type|cell ?line|cell|tissue|cells)\s*:', c, flags=re.I): out.add(c.split(':', 1)[1].strip())
        if s['source']: out.add(s['source'])
    return '; '.join(sorted(out))[:300]

def screen(part, hits, pert_rx, ctrl_rx):
    rows = []
    for g in sorted(hits, key=lambda x: int(x[3:])):
        h = hits[g]
        try: S = _geo.samples(g, part)
        except Exception as e: S = []; h['error'] = str(e)
        hs = [s for s in S if ('Homo sapiens' in s['organism'] or 'Mus musculus' in s['organism'] or not s['organism'])] if part.startswith('P2') else \
             [s for s in S if 'Homo sapiens' in s['organism'] or not s['organism']]
        npert = sum(1 for s in hs if re.search(pert_rx, sample_text(s), flags=re.I))
        nctrl = sum(1 for s in hs if re.search(ctrl_rx, sample_text(s), flags=re.I) and not re.search(pert_rx, sample_text(s), flags=re.I))
        rows.append(dict(gse=g, tags=';'.join(sorted(h['tags'])), gdstype=h['gdstype'], n_samples=h['n_samples'], n_human_samples=len(hs),
                         gpl=h['gpl'], auto_perturbed=npert, auto_control=nctrl, strategies=';'.join(sorted({s['strategy'] for s in hs if s['strategy']})),
                         cells=cell_fields(hs), title=h['title'], pubmed=h['pubmed'], summary=h['summary'][:600].replace('\t', ' ').replace('\n', ' ')))
    d = os.path.join(INB, part.split('_')[0])
    with open(os.path.join(d, f'{part}_hits.tsv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ['gse'], delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
    return rows
