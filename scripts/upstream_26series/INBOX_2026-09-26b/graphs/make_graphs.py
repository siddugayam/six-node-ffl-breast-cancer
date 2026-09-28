#!/usr/bin/env python3
"""Graph variants for the INBOX_2026-09-26b re-runs.  Reads the project read-only.

Census format (for ffl_census_composition.c), built with the census scripts' own build_graph():
  table3             scripts/03_ffl_census.py                       (legacy miRNA-miRNA arcs kept, STRING kept)
  nolegacy           scripts/v7/03_ffl_census_nolegacymirna.py      (legacy dropped, STRING kept)
  table3_nostring    table3 minus STRING arcs (definition of scripts/v2/R9_no_string_census.py)
  nolegacy_nostring  nolegacy minus STRING arcs
Null format (for scripts/v2/v2_null.c), made by deleting lines from the v2 input files so that the
order of the remaining arcs, and hence the random stream, is unchanged:
  dep, pub                     = results/v2/orig_dep.txt, orig_pub.txt (copied unchanged)
  dep_nolegacy                 = orig_dep minus the 30 canonical miRNA_miRNA edges
  pub_nolegacy                 = orig_pub minus the 28 legacy arcs not re-added by the 10 kb layer
  pub_nostring, pub_nolegacy_nostring = the same minus the 833 STRING arcs
Every variant is checked for consistency between the two formats."""
import sys, os, csv, json, importlib.util, collections, shutil
sys.dont_write_bytecode = True
REV = '/path/to/revision'
CLS = ['TF->Gene', 'TF->TF', 'TF->miRNA', 'miRNA->Gene', 'miRNA->TF', 'Gene->Gene', 'miRNA-miRNA']
TMAP = {'miRNA': 0, 'TF': 1, 'Gene': 2}


def load_mod(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); sys.argv = ['x', '3']
    spec.loader.exec_module(m); return m


def paper_cls(nodes, edges, a, b):
    t = edges[(a, b)]
    if t == 'TF_target': return 'TF->TF' if nodes[b] == 'TF' else 'TF->Gene'
    if t == 'miRNA_target': return 'miRNA->TF' if nodes[b] == 'TF' else 'miRNA->Gene'
    return {'TF_miRNA': 'TF->miRNA', 'gene_gene': 'Gene->Gene', 'miRNA_miRNA': 'miRNA-miRNA'}[t]


dep_pairs = {(r['source'], r['target']) for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'), delimiter='\t')}
legacy = {(r['source'], r['target']) for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'), delimiter='\t')
          if r['edge_type'] == 'miRNA_miRNA'}


def string_arcs(nodes, edges):
    sp = set()
    for r in csv.DictReader(open(f'{REV}/data/layer_gene_gene.tsv'), delimiter='\t'):
        if r['evidence'] == 'STRING' and r['source'] in nodes and r['target'] in nodes and r['source'] != r['target']:
            sp.add((r['source'], r['target']))
    return {e for e in edges if edges[e] == 'gene_gene' and e in sp and e not in dep_pairs}


def export_census(tag, nodes, edges):
    V = sorted(nodes); idx = {v: i for i, v in enumerate(V)}
    orig = set(edges)
    for (a, b) in dep_pairs:     # contraction removed miRNA->TF arcs of reciprocal pairs; restore them
        if nodes.get(a) == 'miRNA' and nodes.get(b) == 'TF' and (b, a) in edges and (a, b) not in edges:
            orig.add((a, b))
    with open(f'census_{tag}_fine.txt', 'w') as f:
        f.write(f"{len(V)} {len(edges)} 1\n" + ' '.join(str(TMAP[nodes[v]]) for v in V) + '\n')
        for (a, b) in sorted(edges): f.write(f"{idx[a]} {idx[b]} 0\n")
    with open(f'census_{tag}_R4.txt', 'w') as f:
        f.write(f"{len(V)} {len(edges)}\n")
        for v in V: f.write(f"{v}\t{nodes[v]}\n")
        for (a, b) in sorted(edges): f.write(f"{idx[a]} {idx[b]} {CLS.index(paper_cls(nodes, edges, a, b))}\n")
    with open(f'census_{tag}_orig.txt', 'w') as f:
        f.write(f"{len(V)} {len(orig)} 1\n" + ' '.join(str(TMAP[nodes[v]]) for v in V) + '\n')
        for (a, b) in sorted(orig): f.write(f"{idx[a]} {idx[b]} 0\n")
    cc = collections.Counter(paper_cls(nodes, edges, a, b) for (a, b) in edges)
    print(f"census {tag:18s} contracted arcs {len(edges):5d}  uncontracted {len(orig):5d}  classes {dict(cc)}")
    return {(a, b) for (a, b) in orig}


census_orig = {}
for tag, path, env, nostr in (('table3', f'{REV}/scripts/03_ffl_census.py', None, False),
                              ('nolegacy', f'{REV}/scripts/v7/03_ffl_census_nolegacymirna.py', '1', False),
                              ('table3_nostring', f'{REV}/scripts/03_ffl_census.py', None, True),
                              ('nolegacy_nostring', f'{REV}/scripts/v7/03_ffl_census_nolegacymirna.py', '1', True)):
    if env: os.environ['DROP_LEGACY_MIRNA'] = env
    m = load_mod(path, 'm_' + tag)
    nodes, edges, nrecip = m.build_graph(True)
    assert nrecip == 1223
    if nostr:
        drop = string_arcs(nodes, edges); assert len(drop) == 833, len(drop)
        for e in drop: del edges[e]
    census_orig[tag] = export_census(tag, nodes, edges)
    NODES = nodes


# ---- null-format files by line deletion -------------------------------------------------------
def read_null(p):
    L = open(p).read().rstrip('\n').split('\n')
    meta = json.load(open(p + '.meta.json'))
    return L, meta


def write_null(tag, L, meta, keep):
    hdr = L[0].split(); body = [l for l, k in zip(L[2:], keep) if k]
    with open(f'null_{tag}.txt', 'w') as f:
        f.write(f"{hdr[0]} {len(body)} {hdr[2]}\n{L[1]}\n" + '\n'.join(body) + '\n')
    json.dump(meta, open(f'null_{tag}.txt.meta.json', 'w'))
    names = meta['nodes']
    arcs = {(names[int(l.split()[0])], names[int(l.split()[1])]) for l in body}
    print(f"null   {tag:22s} arcs {len(body)}")
    return arcs


for g in ('dep', 'pub'):
    for ext in ('', '.meta.json'):
        shutil.copy2(f'{REV}/results/v2/orig_{g}.txt{ext}', f'null_{g}.txt{ext}')
Ld, md = read_null(f'{REV}/results/v2/orig_dep.txt'); Lp, mp = read_null(f'{REV}/results/v2/orig_pub.txt')
assert md['nodes'] == mp['nodes'] == sorted(NODES)
pair = lambda L, m: [(m['nodes'][int(l.split()[0])], m['nodes'][int(l.split()[1])]) for l in L[2:]]
Pd, Pp = pair(Ld, md), pair(Lp, mp)
# legacy arcs that the 10 kb layer re-adds (co-clustered pairs) stay in the augmented graph
layer = set()
for r in csv.DictReader(open(f'{REV}/data/layer_miRNA_miRNA.tsv'), delimiter='\t'):
    if str(r['threshold_10kb']).upper() in ('TRUE', '1'):
        layer.add((r['miRNA_1'], r['miRNA_2'])); layer.add((r['miRNA_2'], r['miRNA_1']))
legacy_pub = {e for e in legacy if e not in layer}
print('legacy canonical arcs', len(legacy), ' of which not re-added by the 10 kb layer', len(legacy_pub))
nodes_all = {v: t for v, t in NODES.items()}
string_pub = {e for e in string_arcs(nodes_all, {e: 'gene_gene' for e in Pp}) }  # candidates by pair
# restrict to arcs that are gene_gene-layer arcs in the census graph
string_pub &= {e for e in census_orig['table3']} - census_orig['table3_nostring']
assert len(string_pub) == 833, len(string_pub)
out = {}
out['dep'] = set(Pd); out['pub'] = set(Pp)
out['dep_nolegacy'] = write_null('dep_nolegacy', Ld, md, [e not in legacy for e in Pd])
out['pub_nolegacy'] = write_null('pub_nolegacy', Lp, mp, [e not in legacy_pub for e in Pp])
out['pub_nostring'] = write_null('pub_nostring', Lp, mp, [e not in string_pub for e in Pp])
out['pub_nolegacy_nostring'] = write_null('pub_nolegacy_nostring', Lp, mp, [(e not in legacy_pub) and (e not in string_pub) for e in Pp])
# consistency: null-format arc sets == census uncontracted arc sets
for nt, ct in (('pub', 'table3'), ('pub_nolegacy', 'nolegacy'), ('pub_nostring', 'table3_nostring'),
               ('pub_nolegacy_nostring', 'nolegacy_nostring')):
    print(f"check null_{nt} == census_{ct} uncontracted:", out[nt] == census_orig[ct])
v2 = lambda p: {tuple(map(int, l.split()[:2])) for l in open(p).read().split('\n')[2:] if l.strip()}
print('check census_table3_fine == results/v2/graph_pub_fine.txt:',
      v2('census_table3_fine.txt') == v2(f'{REV}/results/v2/graph_pub_fine.txt'))
print('deposited minus legacy arcs:', len(out['dep_nolegacy']), ' (6,859 - 30 =', 6859 - 30, ')')
