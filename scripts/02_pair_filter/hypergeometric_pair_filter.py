#!/usr/bin/env python3
"""The miRNA-TF pair filter (Methods 2.1), recomputed from its input block.

Each row of data/pair_filter/value.txt gives a tested miRNA-TF pair with Nmir (targets of the miRNA), Ntf (targets of
the TF) and i (targets shared). The filter computed the lower tail of the hypergeometric distribution, P(X <= i), in a
universe of N = 256 genes, and retained the rows whose Benjamini-Hochberg-adjusted P, taken over all 2,676 rows, was
below 0.05. Rows with Ntf > 256 have no P value ("#NUM!" in the workbook) and are not retained. The upper tail,
P(X >= i), is added for comparison.

Writes results/pair_filter/hypergeometric_pair_filter.tsv and prints the counts given in data/pair_filter/README.md.
Needs only the Python standard library.  Usage: python3 scripts/02_pair_filter/hypergeometric_pair_filter.py
"""
import csv
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
N = 256        # the filter's universe (data/pair_filter/universe_256_genes.txt)
FDR = 0.05


def tail(m, t, i, upper):
    xs = range(i, min(m, t) + 1) if upper else range(i + 1)
    return sum(comb(m, x) * comb(N - m, t - x) for x in xs) / comb(N, t)


rows = []
with open(ROOT / 'data/pair_filter/value.txt') as fh:
    for mirna, tf, m, t, i in csv.reader(fh, delimiter='\t'):
        m, t, i = int(m), int(t), int(i)
        valid = m <= N and t <= N
        rows.append(dict(mirna=mirna, tf=tf, Nmir=m, Ntf=t, i=i,
                         p_lower=tail(m, t, i, False) if valid else None,
                         p_upper=tail(m, t, i, True) if valid else None))

# Benjamini-Hochberg over all rows: rank the defined P values (ties in file order), keep ranks up to the largest k
# with P_(k) <= k / n * FDR, where n counts every row
ranked = sorted((r for r in rows if r['p_lower'] is not None), key=lambda r: r['p_lower'])
k = max((j for j, r in enumerate(ranked, 1) if r['p_lower'] <= j / len(rows) * FDR), default=0)
for j, r in enumerate(ranked, 1):
    r['bh_rank'], r['retained'] = j, j <= k

out = ROOT / 'results/pair_filter/hypergeometric_pair_filter.tsv'
out.parent.mkdir(parents=True, exist_ok=True)
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mirna', 'tf', 'Nmir', 'Ntf', 'i', 'p_lower', 'p_upper', 'bh_rank', 'retained'])
    for r in rows:
        w.writerow([r['mirna'], r['tf'], r['Nmir'], r['Ntf'], r['i'],
                    '' if r['p_lower'] is None else f"{r['p_lower']:.6g}",
                    '' if r['p_upper'] is None else f"{r['p_upper']:.6g}",
                    r.get('bh_rank', ''), 'yes' if r.get('retained') else 'no'])

pairs = {(r['mirna'], r['tf']) for r in rows}
kept = {(r['mirna'], r['tf']): r['i'] for r in rows if r.get('retained')}
defined = {(r['mirna'], r['tf']): r['p_upper'] for r in rows if r['p_upper'] is not None}
print(f"{len(rows)} rows, {len(pairs)} distinct pairs; {sum(1 for r in rows if r.get('retained'))} rows retained, "
      f"{len(kept)} distinct pairs, of which {sum(1 for i in kept.values() if i == 0)} share no target")
print(f"{len(defined)} distinct pairs with a defined P; upper-tail P < 0.05: {sum(1 for p in defined.values() if p < 0.05)}")
print(f"wrote {out.relative_to(ROOT)}")
