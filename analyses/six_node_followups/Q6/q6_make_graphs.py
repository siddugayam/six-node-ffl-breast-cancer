#!/usr/bin/env python3
"""Q6 inputs.  The three variant census graphs for S1's s1_null6 (analyses/six_node_pattern/S1/s1_null6, unchanged):
  (i)   retyped22: analyses/six_node_pattern/E/E5/census_nolegacy_retyped22_null.txt (the E5 re-typing of the 22 non-Lambert
        TF-typed nodes as genes on S1's census graph) with its fine arc classes RECOMPUTED from the new node types
        (sorted 'source->target' type names, the convention of S1/graph_nolegacy_null.txt, whose 8 classes are
        Gene->Gene ... miRNA->miRNA).  The E5 file kept the old classes because E5(c) used it only in OBS mode, where
        classes are not read; NULL-B/C read them.  Labels: S1/graph_nolegacy_labels.txt (same arc order).
  (ii)  validated: analyses/six_node_pattern/S7/null_validated.txt + labels_validated.txt (the S7c graph, with STRING).
  (iii) physical900: analyses/six_node_pattern/S7/graph_physical_ge_0.900.txt + labels_physical_ge_0.900.txt (S7b).
Writes Q6/graph_retyped22_reclassed.txt and Q6/q6_inputs.txt (md5 of every input).  Read-only on the project."""
import os, hashlib
H = os.path.dirname(os.path.abspath(__file__))
B = '/path/to/revision/analyses/six_node_pattern'
src = f'{B}/E/E5/census_nolegacy_retyped22_null.txt'
L = open(src).read().rstrip('\n').split('\n'); nv, ne, nc = map(int, L[0].split()); ty = list(map(int, L[1].split()))
TN = {0: 'miRNA', 1: 'TF', 2: 'Gene'}
arcs = [tuple(map(int, l.split()[:2])) for l in L[2:2 + ne]]
cls = sorted({f'{TN[ty[u]]}->{TN[ty[v]]}' for u, v in arcs})
out = f'{H}/graph_retyped22_reclassed.txt'
open(out, 'w').write(f'{nv} {ne} {len(cls)}\n{L[1]}\n' + '\n'.join(f'{u} {v} {cls.index(TN[ty[u]] + "->" + TN[ty[v]])}' for u, v in arcs) + '\n')
md5 = lambda p: hashlib.md5(open(p, 'rb').read()).hexdigest()
lines = [f'retyped22 classes ({len(cls)}): {cls}; arcs {ne}; nodes {nv}; TF-typed nodes {ty.count(1)}, genes {ty.count(2)}, miRNAs {ty.count(0)}']
for p in [src, out, f'{B}/S1/graph_nolegacy_labels.txt', f'{B}/S7/null_validated.txt', f'{B}/S7/labels_validated.txt',
          f'{B}/S7/graph_physical_ge_0.900.txt', f'{B}/S7/labels_physical_ge_0.900.txt', f'{B}/S1/s1_null6', f'{B}/S1/s1_null6.c']:
    lines.append(f'md5 {md5(p)}  {p}')
open(f'{H}/q6_inputs.txt', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
