#!/usr/bin/env python3
"""E5(b)/(c) inputs: the 22 TF-typed nodes absent from the Lambert census re-typed as genes (type code 1 -> 2) in
  null_dep_nolegacy_retyped22.txt   = INBOX_2026-09-26b/graphs/null_dep_nolegacy.txt (analysed network, 6,829 arcs)
                                      with the node-type line changed and the fine classes (source type -> target
                                      type, sorted names) recomputed; input of E/E5/v2_null (NULL-A/B/C)
  census_nolegacy_retyped22_null.txt = S1/graph_nolegacy_null.txt (census graph, 9,226 arcs) with the node-type line
                                      changed (classes not used); input of S1/s1_null6 OBS mode for BHAT6 (E5c)
Arc order is unchanged in both.  Read-only on the project."""
import json, os
H = os.path.dirname(os.path.abspath(__file__))
G = '/path/to/revision/INBOX_2026-09-26b/graphs'
T22 = 'APEX1 BMI1 BRCA1 CREBBP CTNNB1 EP300 EZH2 HDAC1 HDAC2 HDAC3 HDAC4 HDAC9 ILF3 MEN1 MKL1 MTA1 NCOR1 NF1 RB1 SIRT1 SUZ12 VHL'.split()
def retype(src, dst, reclass):
    L = open(src).read().rstrip('\n').split('\n'); nv, ne, nc = map(int, L[0].split()); ty = list(map(int, L[1].split()))
    names = json.load(open(src + '.meta.json'))['nodes'] if os.path.exists(src + '.meta.json') else \
        open(os.path.join(H, '..', '..', 'S1', 'node_names.txt')).read().split('\n')[:nv]
    ix = {n: i for i, n in enumerate(names)}
    for x in T22: assert ty[ix[x]] == 1; ty[ix[x]] = 2
    arcs = [tuple(map(int, l.split())) for l in L[2:2 + ne]]
    if reclass:
        TN = {0: 'miRNA', 1: 'TF', 2: 'Gene'}; cls = sorted({f'{TN[ty[u]]}->{TN[ty[v]]}' for u, v, _ in arcs})
        arcs = [(u, v, cls.index(f'{TN[ty[u]]}->{TN[ty[v]]}')) for u, v, _ in arcs]; nc = len(cls)
        json.dump({'nodes': names, 'classes': cls}, open(dst + '.meta.json', 'w'))
    open(dst, 'w').write(f'{nv} {ne} {nc}\n' + ' '.join(map(str, ty)) + '\n' + '\n'.join(f'{u} {v} {c}' for u, v, c in arcs) + '\n')
    print(dst, 'nodes', nv, 'arcs', ne, 'classes', nc)
retype(f'{G}/null_dep_nolegacy.txt', os.path.join(H, 'null_dep_nolegacy_retyped22.txt'), True)
retype(os.path.join(H, '..', '..', 'S1', 'graph_nolegacy_null.txt'), os.path.join(H, 'census_nolegacy_retyped22_null.txt'), False)
