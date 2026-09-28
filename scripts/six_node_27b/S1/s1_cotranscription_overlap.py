#!/usr/bin/env python3
"""S1 context: do the co-transcribed miRNA pairs (10 kb layer) share targets and TF inputs more than other miRNA
pairs?  Jaccard overlap of miRNA_target target sets and of TF_miRNA regulator sets, for the 480 co-transcribed
pairs versus all other pairs of network miRNAs (nolegacy census graph, S1 inputs).  Read-only."""
import os, itertools, statistics, collections
H = os.path.dirname(os.path.abspath(__file__))
names = open(f'{H}/node_names.txt').read().split('\n')[:587]
L = open(f'{H}/graph_nolegacy_null.txt').read().rstrip('\n').split('\n'); ty = list(map(int, L[1].split()))
labs = list(map(int, open(f'{H}/graph_nolegacy_labels.txt').read().split()))
arcs = [tuple(map(int, l.split()[:2])) for l in L[2:]]
tg = collections.defaultdict(set); rg = collections.defaultdict(set); pairs = set()
for (u, v), lb in zip(arcs, labs):
    if lb == 1: tg[u].add(v)
    elif lb == 0: rg[v].add(u)
    elif lb == 6: pairs.add(frozenset((u, v)))
mirs = [i for i in range(len(ty)) if ty[i] == 0]
jac = lambda a, b: len(a & b) / len(a | b) if (a | b) else 0.0
co_t = [jac(tg[a], tg[b]) for a, b in map(tuple, pairs)]; co_r = [jac(rg[a], rg[b]) for a, b in map(tuple, pairs)]
ot_t, ot_r = [], []
for a, b in itertools.combinations(mirs, 2):
    if frozenset((a, b)) in pairs: continue
    ot_t.append(jac(tg[a], tg[b])); ot_r.append(jac(rg[a], rg[b]))
out = [f'co-transcribed pairs {len(pairs)}; other miRNA pairs {len(ot_t)}',
       f'target-set Jaccard: co-transcribed mean {statistics.mean(co_t):.3f} (median {statistics.median(co_t):.3f}); other pairs mean {statistics.mean(ot_t):.3f} (median {statistics.median(ot_t):.3f})',
       f'TF-regulator-set Jaccard: co-transcribed mean {statistics.mean(co_r):.3f} (median {statistics.median(co_r):.3f}); other pairs mean {statistics.mean(ot_r):.3f} (median {statistics.median(ot_r):.3f})']
open(f'{H}/s1_cotranscription_overlap.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
