#!/usr/bin/env python3
"""
S7(b): BHAT6 (all classes) and MODEL6 with the G1-G2 link restricted to
  (i)  STRING v12.0 physical-subnetwork links with score >= 0.900
  (ii) STRING v12.0 co-expression-channel scores >= 0.900
Downloaded public files (not in the project; kept outside it, URL and md5 recorded):
  https://stringdb-downloads.org/download/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz
      md5 112d22188a72885c00e7baec9dc7d74d
  https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz
      md5 1fbb6885583541bf0b3d82919ed0b649
Symbols map to STRING proteins through the project's own exact map (cache/string_symbol_map_exact.tsv,
scripts/01_network_assembly/05a_string_map_and_fetch.R).  The downloaded files are cross-checked against the cached API
output behind the gene-gene layer (cache/string_physical_900_ensp.tsv, string_functional_900_ensp.tsv).
The G1-G2 link of BHAT6 and MODEL6 is a gene_gene arc of the census graph (labels 4 and 5 of
S1/graph_nolegacy_labels.txt; the deposit's COL1A1-COL3A1 edge and the 833 STRING arcs).  Every such
arc whose pair fails the criterion is removed; nothing else changes (gene_gene arcs enter BHAT6 and
MODEL6 only as the G1-G2 link).  Counts: S1/s1_null6 in OBS mode on each restricted graph.
usage: s7b_restricted_gene_gene.py <dir with the two downloaded files>
"""
import sys, os, gzip, csv, subprocess, collections
HERE = os.path.dirname(os.path.abspath(__file__))
REV = '/path/to/revision'
S1 = os.path.join(HERE, '..', 'S1')
DL = sys.argv[1]
smap = {r['ensp']: r['symbol'] for r in csv.DictReader(open(f'{REV}/cache/string_symbol_map_exact.tsv'), delimiter='\t') if r['ensp']}
U = set(smap)
phys, comb, coex = {}, {}, {}
with gzip.open(f'{DL}/9606.protein.physical.links.v12.0.txt.gz', 'rt') as f:
    next(f)
    for l in f:
        a, b, s = l.split()
        if a in U and b in U and a != b: phys[frozenset((smap[a], smap[b]))] = int(s)
with gzip.open(f'{DL}/9606.protein.links.detailed.v12.0.txt.gz', 'rt') as f:
    next(f)
    for l in f:
        x = l.split()
        if x[0] in U and x[1] in U and x[0] != x[1]:
            p = frozenset((smap[x[0]], smap[x[1]])); coex[p] = int(x[5]); comb[p] = int(x[9])
P_phys = {p for p, s in phys.items() if s >= 900}
P_coex = {p for p, s in coex.items() if s >= 900}
P_comb = {p for p, s in comb.items() if s >= 900}
out = []
Pr = lambda *a: (out.append(' '.join(str(x) for x in a)), print(*a, flush=True))
# cross-check with the cached API output behind the gene-gene layer
def cached(fn):
    d = {}
    for r in csv.DictReader(open(f'{REV}/cache/{fn}'), delimiter='\t'):
        if r['stringId_A'] in U and r['stringId_B'] in U and r['stringId_A'] != r['stringId_B']:
            d[frozenset((smap[r['stringId_A']], smap[r['stringId_B']]))] = r
    return d
cp, cf = cached('string_physical_900_ensp.tsv'), cached('string_functional_900_ensp.tsv')
Pr(f'network proteins mapped to STRING: {len(U)}')
Pr(f'physical >= 900: downloaded {len(P_phys)} pairs; cached API {len(cp)}; identical sets: {set(cp) == P_phys}; '
   f'max |score diff| {max(abs(round(float(cp[p]["score"]) * 1000) - phys[p]) for p in cp if p in phys)}')
Pr(f'combined >= 900: downloaded {len(P_comb)} pairs; cached API {len(cf)}; identical sets: {set(cf) == P_comb}')
cx = {p: round(float(r['ascore']) * 1000) for p, r in cf.items()}
Pr(f'co-expression >= 900: downloaded {len(P_coex)} pairs; cached API ascore >= 0.9: {sum(1 for v in cx.values() if v >= 900)}; '
   f'max |ascore - coexpression| over cached pairs {max(abs(cx[p] - coex.get(p, 0)) for p in cx)}; '
   f'co-expression pairs outside the combined >= 900 set: {len(P_coex - P_comb)}')
# restricted graphs
names = open(f'{S1}/node_names.txt').read().split('\n')[:587]
L = open(f'{S1}/graph_nolegacy_null.txt').read().rstrip('\n').split('\n')
nv, ne, ncls = map(int, L[0].split()); labs = open(f'{S1}/graph_nolegacy_labels.txt').read().split()
rows = []
for tag, keepset in (('all_gene_gene_links', None), ('physical_ge_0.900', P_phys), ('coexpression_ge_0.900', P_coex)):
    body, lb, kept, dropped = [], [], collections.Counter(), collections.Counter()
    for l, x in zip(L[2:2 + ne], labs):
        u, v = map(int, l.split()[:2])
        if x in ('4', '5'):
            p = frozenset((names[u], names[v])); t = 'Gene-Gene' if (names[u][:4] != 'hsa-') else '?'
            if keepset is not None and p not in keepset:
                dropped[x] += 1; continue
            kept[x] += 1
        body.append(l); lb.append(x)
    g = f'{HERE}/graph_{tag}.txt'; lf = f'{HERE}/labels_{tag}.txt'
    open(g, 'w').write(f'{nv} {len(body)} {ncls}\n{L[1]}\n' + '\n'.join(body) + '\n'); open(lf, 'w').write('\n'.join(lb) + '\n')
    r = subprocess.run([f'{S1}/s1_null6', g, lf, 'OBS', '0', '20250908', '100', '0', '1', '0'], capture_output=True, text=True)
    hdr, val = r.stdout.strip().split('\n')[:2]
    d = dict(zip(hdr.split('\t'), val.split('\t')))
    gg_gene = sum(1 for l, x in zip(body, lb) if x in ('4', '5') and all(not names[int(z)].startswith('hsa-') for z in l.split()[:2]))
    row = dict(restriction=tag, gene_gene_arcs_kept=sum(kept.values()), deposit_link_kept=kept['4'], STRING_arcs_kept=kept['5'],
               gene_gene_arcs_dropped=sum(dropped.values()),
               **{k: d[k] for k in ('bhat6_comp_inst', 'bhat6_comp_sets', 'bhat6_mirFFL_inst', 'bhat6_mirFFL_sets',
                                    'bhat6_TFFFL_inst', 'bhat6_TFFFL_sets', 'model6_full_inst', 'model6_full_sets',
                                    'model6_def_inst', 'model6_def_sets')})
    rows.append(row); Pr(row)
# Gene-Gene (both Gene-typed) links kept, for the BHAT6 reading
with open(f'{HERE}/s7b_restricted_gene_gene_counts.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
ggpairs = lambda S: sorted(tuple(sorted(p)) for p in S)
ntype = {r['name']: r['type'] for r in csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'), delimiter='\t')}
Pr('co-expression >= 0.900 pairs among network proteins: ' + '; '.join(f"{'-'.join(sorted(p))} ({'/'.join(ntype[x] for x in sorted(p))}) coexpression {coex[p]} combined {comb[p]}" for p in sorted(P_coex, key=sorted)))
Pr('COL1A1-COL3A1: physical score', phys.get(frozenset(('COL1A1', 'COL3A1'))), '; co-expression score', coex.get(frozenset(('COL1A1', 'COL3A1'))),
   '; combined', comb.get(frozenset(('COL1A1', 'COL3A1'))))
open(f'{HERE}/s7b_restricted_gene_gene.txt', 'w').write('\n'.join(out) + '\n')
