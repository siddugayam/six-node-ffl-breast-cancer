#!/usr/bin/env python3
"""Export the census graphs built by scripts/03_ffl_census.py (Table 3 graph) and by
scripts/v7/03_ffl_census_nolegacymirna.py (DROP_LEGACY_MIRNA=1) in the formats read by
../ffl_census_composition.c.  Both scripts are imported read-only; nothing in the project is written.
The Table 3 export is compared arc-for-arc with results/v2/graph_pub_fine.txt as a check."""
import sys, os, importlib.util, collections
sys.dont_write_bytecode = True
REV = '/path/to/revision'
CLS = ['TF->Gene', 'TF->TF', 'TF->miRNA', 'miRNA->Gene', 'miRNA->TF', 'Gene->Gene', 'miRNA-miRNA']
TMAP = {'miRNA': 0, 'TF': 1, 'Gene': 2}

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); sys.argv = ['x', '3']
    spec.loader.exec_module(m); return m

def cls(nodes, edges, a, b):
    t = edges[(a, b)]
    if t == 'TF_target': return 'TF->TF' if nodes[b] == 'TF' else 'TF->Gene'
    if t == 'miRNA_target': return 'miRNA->TF' if nodes[b] == 'TF' else 'miRNA->Gene'
    return {'TF_miRNA': 'TF->miRNA', 'gene_gene': 'Gene->Gene', 'miRNA_miRNA': 'miRNA-miRNA'}[t]

def export(tag, nodes, edges, orig):
    V = sorted(nodes); idx = {v: i for i, v in enumerate(V)}
    with open(f'graph_{tag}_fine.txt', 'w') as f:            # structure + node types (class column unused)
        f.write(f"{len(V)} {len(edges)} 1\n" + ' '.join(str(TMAP[nodes[v]]) for v in V) + '\n')
        for (a, b) in sorted(edges): f.write(f"{idx[a]} {idx[b]} 0\n")
    with open(f'graph_{tag}_R4.txt', 'w') as f:              # paper 7-class codes
        f.write(f"{len(V)} {len(edges)}\n")
        for v in V: f.write(f"{v}\t{nodes[v]}\n")
        for (a, b) in sorted(edges): f.write(f"{idx[a]} {idx[b]} {CLS.index(cls(nodes, edges, a, b))}\n")
    with open(f'graph_{tag}_orig.txt', 'w') as f:            # uncontracted
        f.write(f"{len(V)} {len(orig)} 1\n" + ' '.join(str(TMAP[nodes[v]]) for v in V) + '\n')
        for (a, b) in sorted(orig): f.write(f"{idx[a]} {idx[b]} 0\n")
    print(tag, 'nodes', len(V), 'arcs (contracted)', len(edges), 'arcs (uncontracted)', len(orig),
          'classes', dict(collections.Counter(cls(nodes, edges, a, b) for (a, b) in edges)))

for tag, path, env in (('table3', f'{REV}/scripts/03_ffl_census.py', None),
                       ('nolegacy', f'{REV}/scripts/v7/03_ffl_census_nolegacymirna.py', '1')):
    if env is not None: os.environ['DROP_LEGACY_MIRNA'] = env
    m = load(path, 'cen_' + tag)
    nodes, edges, nrecip = m.build_graph(True)
    # the uncontracted arc set: contracted arcs plus every miRNA->TF arc removed by contraction
    # (a removed arc (m,t) is one whose reverse (t,m) is a TF->miRNA arc and which the source files contain)
    import csv
    orig_graph = set(edges)
    for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'), delimiter='\t'):
        a, b = r['source'], r['target']
        if nodes.get(a) == 'miRNA' and nodes.get(b) == 'TF' and (b, a) in edges and (a, b) not in edges:
            orig_graph.add((a, b))
    print(tag, 'reciprocal pairs contracted by build_graph:', nrecip,
          ' restored miRNA->TF arcs:', len(orig_graph) - len(edges))
    export(tag, nodes, edges, orig_graph)
