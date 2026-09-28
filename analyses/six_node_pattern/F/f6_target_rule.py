#!/usr/bin/env python3
"""
F6: target rule of the named-class D1-D4 three-node FFLs (analyses/six_node_pattern).  Read-only.

Graphs = the null-format inputs of the stored motif runs (analyses/census_and_motif_nulls/graphs):
  null_dep.txt           deposited network, 6,859 edges (30 legacy miRNA-miRNA edges included)
  null_dep_nolegacy.txt  analysed network, 6,829 edges
Enumeration = count3() of scripts/04_motif_significance/v2_null.c, written out instance by instance: contract every
reciprocal TF<->miRNA pair to its TF->miRNA arc; for every arc u->v of the contracted graph without a
reverse arc, every T with u->T and v->T and no arc from T to u or v is one FFL.  Class by the u->v arc:
composite (reciprocal TF->miRNA), TF-FFL (non-reciprocal TF->miRNA), miRNA-FFL (miRNA->TF), other.
"""
import json, collections, os
G = '/path/to/revision/analyses/census_and_motif_nulls/graphs'
OUT = os.path.dirname(os.path.abspath(__file__))
TN = {0: 'miRNA', 1: 'TF', 2: 'Gene'}
lines = []
for tag in ('dep', 'dep_nolegacy'):
    L = open(f'{G}/null_{tag}.txt').read().rstrip('\n').split('\n')
    nv, ne, _ = map(int, L[0].split()); ty = list(map(int, L[1].split()))
    names = json.load(open(f'{G}/null_{tag}.txt.meta.json'))['nodes']
    A = {tuple(map(int, l.split()[:2])) for l in L[2:2 + ne]}
    red, comp = set(), set()
    for (u, v) in A:
        rec = (v, u) in A and {ty[u], ty[v]} == {0, 1}
        if rec and ty[u] == 0: continue
        red.add((u, v))
        if rec: comp.add((u, v))
    out = collections.defaultdict(set); inn = collections.defaultdict(set)
    for u, v in red: out[u].add(v); inn[v].add(u)
    tab = collections.Counter(); via_mm = collections.Counter()
    total = 0
    for (u, v) in red:
        if u in out[v]: continue
        Ts = {t for t in out[u] & out[v] if t not in inn[u] and t not in inn[v]}
        total += len(Ts)
        if (u, v) in comp: k = 'composite'
        elif ty[u] == 1 and ty[v] == 0: k = 'TF-FFL'
        elif ty[u] == 0 and ty[v] == 1: k = 'miRNA-FFL'
        else: k = 'other'
        for t in Ts:
            tab[(k, TN[ty[t]])] += 1
            if ty[t] == 0 and k != 'other': via_mm[(names[u], names[v], names[t])] += 1
    named = sum(c for (k, _), c in tab.items() if k != 'other')
    lines.append(f'== {tag}: arcs {ne}; all D1-D4 three-node FFLs {total}; named-class {named}')
    for k in ('composite', 'TF-FFL', 'miRNA-FFL', 'other'):
        lines.append(f'   {k:10s} ' + ', '.join(f'target {tt}: {tab[(k, tt)]}' for tt in ('Gene', 'TF', 'miRNA')))
    mir_named = sum(c for (k, tt), c in tab.items() if k != 'other' and tt == 'miRNA')
    lines.append(f'   named-class FFLs with a miRNA target: {mir_named}')
    if mir_named:
        mm_arcs = collections.Counter()
        for (u, v, t), c in via_mm.items(): mm_arcs[(v, t)] += c
        lines.append(f'   their miRNA->miRNA arcs (v -> target): {dict(mm_arcs)}')
open(f'{OUT}/f6_target_rule.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
