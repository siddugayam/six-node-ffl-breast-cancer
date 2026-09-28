#!/usr/bin/env python3
"""
S0(a): trace the EA6 objects (analyses/six_node_pattern).  Read-only on the project; writes into this folder.

EA6 objects and the scripts that built them
  scripts/11_ffl_module_membership.R  (2026-09-08 16:56) -> data/ffl_module_sets.rds,
      results/ffl_module_membership.csv: N3..N6 node sets and HIGHER_ONLY = N6 minus N3
  scripts/13_circuit_level_3node_vs_6node.R (2026-09-08 17:07) -> results/circuit_level_3node_vs_6node.csv:
      2,000 six-node circuits sampled (seed 1234), 1,989 distinct, each with its own three-node core
Both scripts use the same link sets, re-derived here exactly as they build them:
  three-node core  TF->miRNA  data/layer_TF_miRNA.tsv; TF->Gene data/layer_TF_target.tsv (TF-typed source,
                   Gene-typed target); miRNA->Gene data/canonical_edges.tsv edge_type miRNA_target
  G2 (gene-gene)   data/layer_gene_gene.tsv, ALL rows (STRING and TRRUST), undirected, both Gene-typed
  M2 (miRNA-miRNA) data/layer_miRNA_miRNA.tsv, ALL rows, undirected
  T2 (TF-TF)       data/layer_TF_target.tsv, TF-typed both ends, undirected
"""
import csv, collections, os
REV = '/path/to/revision'
OUT = os.path.dirname(os.path.abspath(__file__))
rd = lambda p: list(csv.DictReader(open(f'{REV}/{p}'), delimiter='\t'))
ntype = {r['name']: r['type'] for r in rd('data/canonical_nodes.tsv')}
ced = rd('data/canonical_edges.tsv')
legacy = {frozenset((r['source'], r['target'])) for r in ced if r['edge_type'] == 'miRNA_miRNA'}
assert len(legacy) == 30
dep_gg = {frozenset((r['source'], r['target'])) for r in ced if r['edge_type'] == 'gene_gene'}

# ---- link sets exactly as 11_/13_ build them
gg = rd('data/layer_gene_gene.tsv')
ggu = collections.defaultdict(set)          # pair -> evidence
for r in gg:
    a, b = r['source'], r['target']
    if a != b and ntype.get(a) == 'Gene' and ntype.get(b) == 'Gene':
        ggu[frozenset((a, b))].add(r['evidence'])
mm = rd('data/layer_miRNA_miRNA.tsv')
mmu = {}
for r in mm:
    if r['miRNA_1'] != r['miRNA_2']:
        mmu[frozenset((r['miRNA_1'], r['miRNA_2']))] = r
tfl = rd('data/layer_TF_target.tsv')
tftf = {frozenset((r['source'], r['target'])) for r in tfl
        if ntype.get(r['source']) == 'TF' and ntype.get(r['target']) == 'TF' and r['source'] != r['target']}

lines = []
P = lambda *a: lines.append(' '.join(str(x) for x in a))
P('## link sets used by 11_ffl_module_membership.R and 13_circuit_level_3node_vs_6node.R')
ev = collections.Counter('+'.join(sorted(v)) for v in ggu.values())
P('G2 gene-gene links (undirected, Gene-Gene):', len(ggu), 'pairs; evidence', dict(ev),
  '; of which the one deposit gene_gene edge:', sum(1 for p in ggu if p in dep_gg))
thr = collections.Counter((r['threshold_3kb'], r['threshold_10kb'], r['threshold_50kb']) for r in mmu.values())
P('M2 miRNA-miRNA pairs (undirected):', len(mmu), 'pairs; (3kb,10kb,50kb) flags', dict(thr),
  '; max distance bp', max(int(float(r['distance_bp'])) for r in mmu.values()))
P('T2 TF-TF pairs (undirected, TRRUST):', len(tftf))
leg_in = sorted(tuple(sorted(p)) for p in legacy if p in mmu)
P('legacy GeneMANIA pairs that are ALSO genomic pairs in layer_miRNA_miRNA.tsv:', len(leg_in), leg_in,
  [(mmu[frozenset(p)]['distance_bp'], mmu[frozenset(p)]['threshold_10kb']) for p in leg_in])

# ---- node sets
mem = rd('results/ffl_module_membership.csv') if False else list(csv.DictReader(open(f'{REV}/results/ffl_module_membership.csv')))
T = lambda col: {r['node'] for r in mem if r[col] == 'TRUE'}
N3, N4, N5, N6, HO = T('in_3node'), T('in_4node_set'), T('in_5node_set'), T('in_6node_set'), T('higher_order_only')
pc = lambda S: {v for v in S if ntype[v] != 'miRNA'}
P('## node sets (results/ffl_module_membership.csv)')
P('N6 =', len(N6), 'nodes; protein-coding', len(pc(N6)), '; HIGHER_ONLY =', len(HO), 'nodes; protein-coding', len(pc(HO)))
# would N5/N6 change if the legacy pairs were dropped from the M2 link set?
ffl3 = list(csv.DictReader(open(f'{REV}/results/ffl3_composite_enumerated.csv')))
mir3 = {r['miR'] for r in ffl3}
M2_all = {m for p in mmu for m in p for o in p if o != m and o in mir3}
M2_noleg = {m for p in mmu if p not in legacy for m in p for o in p if o != m and o in mir3}
P('miRNA2 partners of three-node miRNAs: with the legacy-coinciding pairs', len(M2_all),
  '; without them', len(M2_noleg), '; lost:', sorted(M2_all - M2_noleg))
lost_nodes = (M2_all - M2_noleg) - N4 if (M2_all - M2_noleg) else set()
P('nodes that would leave N5/N6 without those pairs (not in N4 and not a TF2):', sorted(lost_nodes))

# ---- the 1,989 circuits
C = list(csv.DictReader(open(f'{REV}/results/circuit_level_3node_vs_6node.csv')))
six = [r for r in C if r['circuit_class'] == '6-node']
circ = sorted({r['members'] for r in six})
P('## circuits (results/circuit_level_3node_vs_6node.csv)')
P('distinct six-node circuits:', len(circ), '; rows per endpoint', dict(collections.Counter(r['endpoint'] for r in six)))
n_leg = 0; used_leg = collections.Counter(); ggev = collections.Counter(); mmthr = collections.Counter()
rows = []
for m in circ:
    TF, miR, G1, T2, M2, G2 = m.split(';')
    pm = frozenset((miR, M2)); pg = frozenset((G1, G2)); pt = frozenset((TF, T2))
    assert pm in mmu and pg in ggu and pt in tftf, m
    leg = pm in legacy
    n_leg += leg
    if leg: used_leg[tuple(sorted(pm))] += 1
    ggev['+'.join(sorted(ggu[pg]))] += 1
    mmthr['10kb' if mmu[pm]['threshold_10kb'] == 'TRUE' else ('50kb only' if mmu[pm]['threshold_50kb'] == 'TRUE' else 'other')] += 1
    rows.append(dict(TF1=TF, miR1=miR, G1=G1, TF2=T2, miR2=M2, G2=G2, miR_pair_is_legacy=leg,
                     miR_pair_distance_bp=mmu[pm]['distance_bp'], miR_pair_10kb=mmu[pm]['threshold_10kb'],
                     G1_G2_evidence='+'.join(sorted(ggu[pg])), G1_G2_is_deposit_edge=pg in dep_gg))
P('circuits whose miRNA pair is one of the 30 legacy pairs:', n_leg, dict(used_leg))
P('circuits by G1-G2 evidence:', dict(ggev))
P('circuits by miRNA-pair distance class:', dict(mmthr))
with open(f'{OUT}/ea6_1989_circuits_links.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
open(f'{OUT}/s0a_ea6_trace.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
