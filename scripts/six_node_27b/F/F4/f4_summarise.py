#!/usr/bin/env python3
"""F4 summary.  Exhaustive D1-D4 modules of the Table 3 census graph (8,031 arcs, 28 legacy one-direction
miRNA-miRNA arcs) at n = 4, 5, 6 (f4_node_union.c = ffl_census_composition.c + node marks).  For each node:
'any' = in at least one module; 'free' = in at least one module without a legacy arc.  Compared with the
stored enrichment gene sets results/motif_node_sets.tsv (08_motif_node_sets.py; RAND-ESU samples on the
same graph; protein-coding members 352 / 340 / 314 symbols, 349 / 337 / 311 Entrez-mapped).
The stored sample cannot be regenerated exactly: 08_motif_node_sets.py pops candidates from Python sets
of node names, whose order depends on the per-process string-hash seed, which was not fixed."""
import glob, os, collections, csv
H = os.path.dirname(os.path.abspath(__file__)); REV = '/path/to/revision'
names = [l.split('\t')[0] for l in open(f'{REV}/INBOX_2026-09-26b/graphs/census_table3_R4.txt').read().split('\n')[1:588]]
ntype = {r['name']: r['type'] for r in csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'), delimiter='\t')}
out = []
P = lambda *a: (out.append(' '.join(str(x) for x in a)), print(*a))
stored = collections.defaultdict(set)
for r in csv.DictReader(open(f'{REV}/results/motif_node_sets.tsv'), delimiter='\t'): stored[r['motif_set']].add(r['node'])
EXPECT = {4: 117967, 5: 1534131, 6: 19838598}
rows = []
for k, files in ((4, ['parts/n4.txt']), (5, ['parts/n5.txt']), (6, sorted(glob.glob('parts/n6_part_*.txt')))):
    anym, freem, found = set(), set(), 0
    for f in files:
        for l in open(f):
            x = l.split()
            if x[0] == 'found': found += int(x[1])
            elif x[0] == 'node':
                if x[2] == '1': anym.add(names[int(x[1])])
                if x[3] == '1': freem.add(names[int(x[1])])
    S = {v for v in stored[f'{k}-node'] if ntype[v] != 'miRNA'}
    only_leg = sorted(v for v in S if v not in freem)
    not_any = sorted(v for v in S if v not in anym)
    P(f'n = {k}: modules {found:,} (Table 3 exhaustive {EXPECT[k]:,}: {found == EXPECT[k]}); nodes in any module {len(anym)}, '
      f'in a legacy-free module {len(freem)}; stored protein-coding set {len(S)}: members found in no module at all {not_any}; '
      f'members found ONLY in modules with a legacy arc {only_leg}')
    for v in sorted(S): rows.append(dict(n=k, gene=v, type=ntype[v], in_any_module=v in anym, in_legacy_free_module=v in freem))
    pc_free = {v for v in freem if ntype[v] != 'miRNA'}; pc_any = {v for v in anym if ntype[v] != 'miRNA'}
    P(f'   protein-coding nodes in a legacy-free module: {len(pc_free)}; in any module: {len(pc_any)}; '
      f'stored set is a subset of the legacy-free union: {S <= pc_free}')
with open(f'{H}/f4_gene_membership.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
open(f'{H}/f4_summary.txt', 'w').write('\n'.join(out) + '\n')
