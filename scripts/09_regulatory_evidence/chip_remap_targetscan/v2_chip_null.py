#!/usr/bin/env python3
"""
Is the ChIP-seq support of the 770 published TF->target edges better than chance?

For every network TF we read only columns 1-2 of its ChIP-Atlas target table
(gene symbol, mean MACS2 score over all experiments for that antigen).  A gene absent from
the table has no peak in any experiment at that TSS distance.

Null models (edge count per TF preserved, TF identity preserved, target resampled):
  NULL-G  target drawn uniformly from the ChIP-Atlas gene universe
  NULL-N  target drawn uniformly from the 223 target genes used by the network
  NULL-H  target drawn from the same decile of "promoter hotness" (= number of the 148
          profiled network TFs that bind that gene) as the real target
Self-edges and duplicate (TF,target) pairs inside a replicate are allowed (they are rare and
excluded would bias the degree sequence); real edges are NOT excluded from the draw.
"""
import os, sys, csv, random
from collections import defaultdict
import numpy as np

REV = "/path/to/revision"
CHIP = "/path/to/scratch/chip"
OUT = os.path.join(REV, "results", "v2")
NPERM = 1000
SEED = 20260908
THRS = ["1", "5", "10"]

edges = []
with open(os.path.join(REV, "data", "canonical_edges.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["edge_type"] == "TF_target":
            edges.append((r["source"], r["target"]))
tfs = sorted({s for s, t in edges})
net_targets = sorted({t for s, t in edges})

res = []
for thr in THRS:
    avg = {}          # tf -> {gene: mean MACS2}
    for tf in tfs:
        p = os.path.join(CHIP, thr, f"{tf}.{thr}.tsv")
        if not os.path.exists(p):
            continue
        d = {}
        with open(p) as fh:
            fh.readline()
            for line in fh:
                i = line.find("\t"); j = line.find("\t", i + 1)
                g = line[:i]
                try: d[g] = float(line[i + 1:j])
                except ValueError: d[g] = 0.0
        avg[tf] = d
    universe = sorted(set().union(*[set(d) for d in avg.values()]))
    hot = {g: sum(1 for tf in avg if g in avg[tf]) for g in universe}
    hv = np.array([hot[g] for g in universe])
    dec = np.digitize(hv, np.quantile(hv, [.1, .2, .3, .4, .5, .6, .7, .8, .9]))
    bydec = defaultdict(list)
    for g, dd in zip(universe, dec):
        bydec[int(dd)].append(g)
    gene_dec = {g: int(dd) for g, dd in zip(universe, dec)}

    usable = [(s, t) for s, t in edges if s in avg]
    def score(pairs, crit):
        if crit == "any":
            return sum(1 for s, t in pairs if avg[s].get(t, 0.0) > 0)
        return sum(1 for s, t in pairs if avg[s].get(t, 0.0) >= 50)

    rng = random.Random(SEED)
    nulls = {("G", "any"): [], ("G", "strong"): [], ("N", "any"): [], ("N", "strong"): [],
             ("H", "any"): [], ("H", "strong"): []}
    npermdone = 0
    for it in range(NPERM):
        pg = [(s, universe[rng.randrange(len(universe))]) for s, t in usable]
        pn = [(s, net_targets[rng.randrange(len(net_targets))]) for s, t in usable]
        ph = []
        for s, t in usable:
            pool = bydec[gene_dec.get(t, 0)]
            ph.append((s, pool[rng.randrange(len(pool))]))
        for tag, pp in (("G", pg), ("N", pn), ("H", ph)):
            for crit in ("any", "strong"):
                nulls[(tag, crit)].append(score(pp, crit))
        npermdone += 1
    assert npermdone == NPERM, npermdone

    for crit in ("any", "strong"):
        obs = score(usable, crit)
        for tag in ("G", "N", "H"):
            arr = np.array(nulls[(tag, crit)], dtype=float)
            z = (obs - arr.mean()) / arr.std(ddof=1)
            p = (np.sum(arr >= obs) + 1) / (len(arr) + 1)
            res.append(dict(threshold_kb=thr, criterion=f"bound_{crit}", null=f"NULL-{tag}",
                            n_edges_testable=len(usable), observed=obs,
                            observed_frac=round(obs / len(usable), 4),
                            null_mean=round(float(arr.mean()), 2),
                            null_sd=round(float(arr.std(ddof=1)), 2),
                            null_frac=round(float(arr.mean()) / len(usable), 4),
                            Z=round(float(z), 2), fold=round(obs / arr.mean(), 3),
                            p_one_sided=round(float(p), 5), n_permutations=npermdone))
            print(res[-1])
    print(f"[thr={thr}] universe genes={len(universe)}  TFs with data={len(avg)}  "
          f"edges usable={len(usable)}  hotness median={np.median(hv):.0f} max={hv.max()}")

with open(os.path.join(OUT, "chip_edge_null.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(res[0].keys())); w.writeheader(); w.writerows(res)
print("wrote chip_edge_null.csv;  permutations per cell:", NPERM)
