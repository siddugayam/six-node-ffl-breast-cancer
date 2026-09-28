#!/usr/bin/env python3
"""
S7(c) validated-only census graphs (analyses/six_node_pattern).  Read-only on the project.

Definition (read from the analysis plan; stated here so it can be checked):
  * miRNA -> target edges: strong and weak tier only (data/edge_evidence_tier.tsv, = Table S1 'tier');
    the 1,803 predicted_only edges are removed
  * the 31 exemplar-circuit edges are removed: the deposit's 30 miRNA_miRNA edges and its one gene_gene
    edge (COL1A1 -> COL3A1).  These are the deposit's only edges that come solely from the authors'
    exemplar SIF files (4-TF, 5-TF, 6-TF.sif, 'interacts with') and not from a regulatory database
    (the analysis plan calls the 30 'the miRNA-miRNA edges of the deposited exemplar circuits').
  * everything else is built exactly as scripts/v7/03_ffl_census_nolegacymirna.py build_graph(True) does:
    the census layers (TRRUST TF_target, the gene-gene layer file, 10 kb miRNA pairs in both
    directions) are added where the arc is not already present, then every reciprocal TF<->miRNA pair is
    contracted to its TF->miRNA arc.  Because the gene-gene layer lists COL1A1-COL3A1 as a STRING link
    (combined 0.999), the builder re-adds it as a STRING arc when STRING is kept.
  * no-STRING variant: every gene_gene arc from a STRING row that is not a deposit edge is removed
    (the rule of analyses/census_and_motif_nulls/graphs/make_graphs.py).
Validation: the same builder on the unfiltered deposit (legacy edges dropped only) must reproduce the
nolegacy census graph arc for arc (analyses/census_and_motif_nulls/graphs/census_nolegacy_{fine,orig}.txt).
Outputs per variant: census format (fine, R4, orig) for ffl_census_composition, null format (+ meta)
for v2_null / s1_null6, and a layer label file for s1_null6.
"""
import os, csv, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
REV = '/path/to/revision'
CLS = ['TF->Gene', 'TF->TF', 'TF->miRNA', 'miRNA->Gene', 'miRNA->TF', 'Gene->Gene', 'miRNA-miRNA']
TMAP = {'miRNA': 0, 'TF': 1, 'Gene': 2}
nodes = {r['name']: r['type'] for r in csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'), delimiter='\t')}
DEP = list(csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'), delimiter='\t'))
tier = {(r['source'], r['target']): r['tier'] for r in csv.DictReader(open(f'{REV}/data/edge_evidence_tier.tsv'), delimiter='\t')}
string_pairs = {(r['source'], r['target']) for r in csv.DictReader(open(f'{REV}/data/layer_gene_gene.tsv'), delimiter='\t')
                if r['evidence'] == 'STRING'}


def build(validated_only):
    edges, dep = {}, set()
    removed = collections.Counter()
    for r in DEP:
        e = (r['source'], r['target'])
        if r['edge_type'] == 'miRNA_miRNA': removed['exemplar miRNA_miRNA (legacy)'] += 1; continue
        if validated_only and r['edge_type'] == 'gene_gene': removed['exemplar gene_gene'] += 1; continue
        if validated_only and r['edge_type'] == 'miRNA_target' and tier[e] == 'predicted_only':
            removed['miRNA_target predicted_only'] += 1; continue
        edges[e] = r['edge_type']; dep.add(e)

    def add(path, sc, tc, et, undirected=False, filt=None):
        for r in csv.DictReader(open(path), delimiter='\t'):
            if filt and not filt(r): continue
            a, b = r.get(sc), r.get(tc)
            if a not in nodes or b not in nodes or a == b: continue
            for e in ([(a, b), (b, a)] if undirected else [(a, b)]):
                if e not in edges: edges[e] = et
    add(f'{REV}/data/layer_TF_target.tsv', 'source', 'target', 'TF_target')
    add(f'{REV}/data/layer_gene_gene.tsv', 'source', 'target', 'gene_gene')
    add(f'{REV}/data/layer_miRNA_miRNA.tsv', 'miRNA_1', 'miRNA_2', 'miRNA_miRNA', undirected=True,
        filt=lambda r: str(r.get('threshold_10kb', 'TRUE')).upper() in ('TRUE', '1'))
    orig = dict(edges)
    nrecip = 0
    for (a, b) in list(edges):
        if nodes.get(a) == 'TF' and nodes.get(b) == 'miRNA' and (b, a) in edges:
            if edges.pop((b, a), None) is not None: nrecip += 1
    return edges, orig, dep, nrecip, removed


def nostring(edges, orig, dep):
    drop = {e for e, t in orig.items() if t == 'gene_gene' and e in string_pairs and e not in dep}
    return {e: t for e, t in edges.items() if e not in drop}, {e: t for e, t in orig.items() if e not in drop}, len(drop)


def paper_cls(t, b):
    if t == 'TF_target': return 'TF->TF' if nodes[b] == 'TF' else 'TF->Gene'
    if t == 'miRNA_target': return 'miRNA->TF' if nodes[b] == 'TF' else 'miRNA->Gene'
    return {'TF_miRNA': 'TF->miRNA', 'gene_gene': 'Gene->Gene', 'miRNA_miRNA': 'miRNA-miRNA'}[t]


V = sorted(nodes); ix = {v: i for i, v in enumerate(V)}


def export(tag, edges, orig, dep):
    with open(f'{HERE}/census_{tag}_fine.txt', 'w') as f:
        f.write(f"{len(V)} {len(edges)} 1\n" + ' '.join(str(TMAP[nodes[v]]) for v in V) + '\n')
        for (a, b) in sorted(edges, key=lambda e: (ix[e[0]], ix[e[1]])): f.write(f"{ix[a]} {ix[b]} 0\n")
    with open(f'{HERE}/census_{tag}_R4.txt', 'w') as f:
        f.write(f"{len(V)} {len(edges)}\n")
        for v in V: f.write(f"{v}\t{nodes[v]}\n")
        for (a, b) in sorted(edges, key=lambda e: (ix[e[0]], ix[e[1]])): f.write(f"{ix[a]} {ix[b]} {CLS.index(paper_cls(edges[(a, b)], b))}\n")
    with open(f'{HERE}/census_{tag}_orig.txt', 'w') as f:
        f.write(f"{len(V)} {len(orig)} 1\n" + ' '.join(str(TMAP[nodes[v]]) for v in V) + '\n')
        for (a, b) in sorted(orig, key=lambda e: (ix[e[0]], ix[e[1]])): f.write(f"{ix[a]} {ix[b]} 0\n")
    # null format (v2_null.c): uncontracted arcs, class = source type -> target type (sorted names)
    fine = sorted({f"{nodes[a]}->{nodes[b]}" for (a, b) in orig})
    with open(f'{HERE}/null_{tag}.txt', 'w') as f:
        f.write(f"{len(V)} {len(orig)} {len(fine)}\n" + ' '.join(str(TMAP[nodes[v]]) for v in V) + '\n')
        for (a, b) in sorted(orig, key=lambda e: (ix[e[0]], ix[e[1]])):
            f.write(f"{ix[a]} {ix[b]} {fine.index(nodes[a] + '->' + nodes[b])}\n")
    json.dump({'nodes': V, 'classes': fine}, open(f'{HERE}/null_{tag}.txt.meta.json', 'w'))
    # s1_null6 layer labels (same codes as S1/make_s1_inputs.py)
    lab = []
    for (a, b) in sorted(orig, key=lambda e: (ix[e[0]], ix[e[1]])):
        t = orig[(a, b)]; ind = (a, b) in dep
        lab.append({'TF_miRNA': 0, 'miRNA_target': 1}.get(t) if t in ('TF_miRNA', 'miRNA_target') else
                   (2 if ind else 3) if t == 'TF_target' else (4 if ind else 5) if t == 'gene_gene' else 6)
    open(f'{HERE}/labels_{tag}.txt', 'w').write('\n'.join(map(str, lab)) + '\n')
    cc = collections.Counter(paper_cls(t, b) for (a, b), t in edges.items())
    print(f"{tag:28s} contracted {len(edges):5d} uncontracted {len(orig):5d} classes {dict(cc)}")


# ---- validation: unfiltered deposit reproduces the nolegacy census graph
e0, o0, d0, n0, r0 = build(False)
G = f'{REV}/analyses/census_and_motif_nulls/graphs'
arcset = lambda p: {tuple(map(int, l.split()[:2])) for l in open(p).read().split('\n')[2:] if l.strip()}
ok = ({(ix[a], ix[b]) for (a, b) in e0} == arcset(f'{G}/census_nolegacy_fine.txt') and
      {(ix[a], ix[b]) for (a, b) in o0} == arcset(f'{G}/census_nolegacy_orig.txt') and n0 == 1223)
print('VALIDATION builder reproduces the nolegacy census graph (fine and orig arc sets):', ok, len(e0), len(o0))
assert ok
e1, o1, d1, n1, r1 = build(True)
print('validated-only: removed from the deposit', dict(r1), '; reciprocal TF<->miRNA pairs contracted', n1)
export('validated', e1, o1, d1)
e2, o2, k2 = nostring(e1, o1, d1)
print('   no-STRING variant: STRING arcs removed', k2)
export('validated_nostring', e2, o2, d1)
