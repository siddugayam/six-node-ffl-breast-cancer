#!/usr/bin/env python3
"""F1(b): compare the prioritisation re-run on the analysed network (6,829 edges) with the stored run.
Stored run = out_reproduce/ (byte-identical to results/v5/node_prioritisation_full.csv).
The thirty prioritised nodes = top 10 RANKABLE nodes per class (priority not NA), as written to
top10_per_class.csv by 01_prioritisation.R (Table 4).  Within-class ranks are recomputed among rankable
nodes (ties 'min'), as scripts/06_node_prioritisation/40_node_compendium_assemble.py does: the stored rank_within_type column
sorts the 4 unrankable miRNAs first and so offsets every miRNA rank by 4."""
import csv, os
H = os.path.dirname(os.path.abspath(__file__))
rd = lambda p: {r['name']: r for r in csv.DictReader(open(p))}
A = rd(f'{H}/out_reproduce/node_prioritisation_full.csv'); B = rd(f'{H}/out_analysed_network/node_prioritisation_full.csv')
S4 = {r['name']: r for r in csv.DictReader(open('/path/to/revision/results/v5/tables/TableS4_node_prioritisation_full.csv'))}
def ranks(X):
    out = {}
    for ty in ('TF', 'Gene', 'miRNA'):
        pr = sorted(((float(r['priority']), n) for n, r in X.items() if r['type'] == ty and r['priority'] not in ('', 'NA')), reverse=True)
        for n in X:
            if X[n]['type'] == ty and X[n]['priority'] not in ('', 'NA'):
                out[n] = 1 + sum(1 for p, _ in pr if p > float(X[n]['priority']))
    return out
RA, RB = ranks(A), ranks(B)
out = []
P = lambda *a: (out.append(' '.join(str(x) for x in a)), print(*a))
deg_changed = sorted((n, A[n]['degree'], B[n]['degree']) for n in A if A[n]['degree'] != B[n]['degree'])
P('nodes whose degree changes (stored -> analysed network):', deg_changed)
P('Table S4 degree column equals the stored (deposited-network) degree for these nodes:', all(S4[n]['degree'] == A[n]['degree'] for n, _, _ in deg_changed))
ta = [r['name'] for r in csv.DictReader(open(f'{H}/out_reproduce/top10_per_class.csv'))]
tb = [r['name'] for r in csv.DictReader(open(f'{H}/out_analysed_network/top10_per_class.csv'))]
st = [r['name'] for r in csv.DictReader(open('/path/to/revision/results/v5/top10_per_class.csv'))]
P('stored top10_per_class.csv reproduced:', ta == st)
P('thirty prioritised nodes identical after re-run:', set(ta) == set(tb), '; identical order:', ta == tb)
chg = [(n, A[n]['type'], RA[n], RB[n]) for n in ta if RA[n] != RB[n]]
P('within-class rank changes among the thirty (name, class, stored, analysed):', chg)
for n in ('SP1', 'RELA', 'NFKB1'):
    P(f"{n}: rank among the 157 TFs stored {RA[n]}, analysed network {RB[n]}")
for n in sorted({x for x, _, _ in deg_changed}):
    P(f"{n}: rank among rankable miRNAs stored {RA.get(n)}, analysed network {RB.get(n)}")
open(f'{H}/f1_compare.txt', 'w').write('\n'.join(out) + '\n')
