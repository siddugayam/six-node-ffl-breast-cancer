#!/usr/bin/env python3
"""
S1 inputs.  Reads the project read-only; writes only into this folder.

graph_nolegacy_null.txt    = analyses/census_and_motif_nulls/graphs/null_pub_nolegacy.txt, copied unchanged: the input
                             of the stored NULL-B / NULL-C runs (uncontracted census graph without the
                             legacy miRNA-miRNA arcs; 9,226 arcs, fine classes = source type -> target type)
graph_nolegacy_labels.txt  one layer label per arc line of that file, in the same order:
    0 TF_miRNA                 (analysed network)
    1 miRNA_target             (analysed network; miRNA->Gene and miRNA->TF, incl. the reciprocal arcs)
    2 TF_target                (analysed network)
    3 TF_target, TRRUST layer  (census layer, NOT in the analysed network)
    4 gene_gene, deposit       (analysed network; one edge)
    5 gene_gene, STRING        (census layer, NOT in the analysed network; oriented alphabetically)
    6 miRNA_miRNA, 10 kb       (census layer, NOT in the analysed network; both directions)
The layer of each arc comes from the census graph builder itself
(scripts/03_ffl_census/03_ffl_census_nolegacymirna.py, build_graph(True), before contraction), and 'analysed
network' = data/canonical_edges.tsv minus its 30 miRNA_miRNA rows (6,829 edges).
node_names.txt             node index -> name (the order of the null-format meta file)
"""
import sys, os, csv, json, shutil, importlib.util, collections
sys.dont_write_bytecode = True
REV = '/path/to/revision'
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = f'{REV}/analyses/census_and_motif_nulls/graphs/null_pub_nolegacy.txt'

os.environ['DROP_LEGACY_MIRNA'] = '1'
spec = importlib.util.spec_from_file_location('m', f'{REV}/scripts/03_ffl_census/03_ffl_census_nolegacymirna.py')
m = importlib.util.module_from_spec(spec); sys.argv = ['x', '3']; spec.loader.exec_module(m)
nodes, edges, nrecip = m.build_graph(True)            # contracted
assert nrecip == 1223 and len(edges) == 8003
dep = {(r['source'], r['target']): r['edge_type'] for r in
       csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'), delimiter='\t') if r['edge_type'] != 'miRNA_miRNA'}
assert len(dep) == 6829
layer = dict(edges)
for (a, b), t in dep.items():                            # restore the contracted miRNA->TF arcs
    if (a, b) not in layer:
        assert nodes[a] == 'miRNA' and nodes[b] == 'TF' and (b, a) in edges and t == 'miRNA_target'
        layer[(a, b)] = t
assert len(layer) == 9226

shutil.copy2(SRC, f'{HERE}/graph_nolegacy_null.txt')
meta = json.load(open(SRC + '.meta.json'))
names = meta['nodes']
L = open(SRC).read().rstrip('\n').split('\n')
nv, ne, ncls = map(int, L[0].split())
types = list(map(int, L[1].split()))
TMAP = {'miRNA': 0, 'TF': 1, 'Gene': 2}
assert all(TMAP[nodes[names[i]]] == types[i] for i in range(nv))
labs = []
cnt = collections.Counter()
for l in L[2:2 + ne]:
    u, v, c = map(int, l.split()); a, b = names[u], names[v]
    t = layer[(a, b)]; ind = (a, b) in dep
    if t == 'TF_miRNA': lab = 0; assert ind
    elif t == 'miRNA_target': lab = 1; assert ind
    elif t == 'TF_target': lab = 2 if ind else 3
    elif t == 'gene_gene': lab = 4 if ind else 5
    elif t == 'miRNA_miRNA': lab = 6; assert not ind
    else: raise ValueError(t)
    if lab == 5: assert a < b
    labs.append(lab); cnt[(lab, nodes[a], nodes[b])] += 1
assert len(labs) == ne == 9226
open(f'{HERE}/graph_nolegacy_labels.txt', 'w').write('\n'.join(map(str, labs)) + '\n')
open(f'{HERE}/node_names.txt', 'w').write('\n'.join(names) + '\n')
for k in sorted(cnt): print(k, cnt[k])
print('analysed-network arcs', sum(1 for x in labs if x in (0, 1, 2, 4)), ' layer arcs', sum(1 for x in labs if x in (3, 5, 6)))
