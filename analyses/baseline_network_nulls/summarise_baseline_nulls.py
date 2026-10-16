#!/usr/bin/env python3
"""Summary of the three-node FFL motif test on the baseline networks, beside the analysed network (see SETTINGS.md).

Reads runs/null_<null model>_<network>_<typing>.tsv (this folder) and, for the analysed network, the published
replicate files of the same program (root: environment variable FFL_ROOT), which carry the program's model argument in
their names: analyses/census_and_motif_nulls/nulls/null_<argument>_dep_noexemplar.tsv (published typing) and
analyses/six_node_pattern/E/E5/null_<argument>_dep_noexemplar_retyped22.tsv (22 co-regulators as genes).
Writes null_model_summary.csv and decision.csv.
usage: FFL_ROOT=<repository checkout> python3 summarise_baseline_nulls.py"""
import csv, os
import numpy as np

ROOT = os.environ["FFL_ROOT"]
HERE = os.path.dirname(os.path.abspath(__file__))
# program argument -> null model
NULLS = {"A": "degree-preserving", "B": "reciprocity-preserving", "C": "bipartite-layer"}
METRICS = {"total": "all three-node FFLs", "comp": "composite FFLs", "tf": "TF-FFLs", "mir": "miRNA-FFLs"}
TYPING = {"published": "published typing (157 TFs)", "retyped": "22 co-regulators typed as genes"}
SRC = {}
for t in TYPING:
    for n in ("validated", "validated_plus_predicted"):
        SRC[(n, t)] = {m: os.path.join(HERE, "runs", f"null_{NULLS[m].replace('-', '_')}_{n}_{t}.tsv") for m in NULLS}
SRC[("analysed network", "published")] = {m: os.path.join(ROOT, "analyses/census_and_motif_nulls/nulls", f"null_{m}_dep_noexemplar.tsv") for m in NULLS}
SRC[("analysed network", "retyped")] = {m: os.path.join(ROOT, "analyses/six_node_pattern/E/E5", f"null_{m}_dep_noexemplar_retyped22.tsv") for m in NULLS}

rows = []
for (net, t), files in SRC.items():
    for m, fn in files.items():
        L = [l.rstrip("\n").split("\t") for l in open(fn)]
        hdr = L[0]; obs = dict(zip(hdr[1:], map(int, L[1][1:])))
        arr = np.array([[int(x) for x in r[1:]] for r in L[2:]], float); R = arr.shape[0]; assert R == 1000
        col = {c: arr[:, i] for i, c in enumerate(hdr[1:])}
        for k, lab in METRICS.items():
            v = col[k]; o = obs[k]; mu, sd = v.mean(), v.std(ddof=1)
            rows.append(dict(network=net, typing=TYPING[t], null_model=NULLS[m], metric=lab, observed=o, expected=mu, sd=sd,
                             fold=o / mu if mu > 0 else float("nan"), excess_pct=100 * (o / mu - 1) if mu > 0 else float("nan"),
                             p_upper=(1 + int((v >= o).sum())) / (R + 1), p_lower=(1 + int((v <= o).sum())) / (R + 1),
                             reciprocal_pairs=obs["recip"], replicates=R))
with open(os.path.join(HERE, "null_model_summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
dec = []
for (net, t) in SRC:
    for lab in METRICS.values():
        b = [r for r in rows if r["network"] == net and r["typing"] == TYPING[t] and r["metric"] == lab and r["null_model"] == "reciprocity-preserving"][0]
        c = [r for r in rows if r["network"] == net and r["typing"] == TYPING[t] and r["metric"] == lab and r["null_model"] == "bipartite-layer"][0]
        d = "persists" if b["p_upper"] < 0.05 and c["p_upper"] < 0.05 else "reverses" if b["p_lower"] < 0.05 and c["p_lower"] < 0.05 else "vanishes"
        dec.append(dict(network=net, typing=TYPING[t], metric=lab, observed=b["observed"],
                        excess_reciprocity_preserving_pct=b["excess_pct"], p_upper_reciprocity_preserving=b["p_upper"],
                        p_lower_reciprocity_preserving=b["p_lower"], excess_bipartite_layer_pct=c["excess_pct"],
                        p_upper_bipartite_layer=c["p_upper"], p_lower_bipartite_layer=c["p_lower"],
                        decision=d if net != "analysed network" else d + " (reference)"))
with open(os.path.join(HERE, "decision.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(dec[0])); w.writeheader(); w.writerows(dec)
for r in dec:
    print(f"{r['network']:26s} {r['typing'][:22]:22s} {r['metric']:20s} obs {r['observed']:6d} "
          f"reciprocity-preserving {r['excess_reciprocity_preserving_pct']:+6.1f}% (p_up {r['p_upper_reciprocity_preserving']:.3f}) "
          f"bipartite-layer {r['excess_bipartite_layer_pct']:+6.1f}% (p_up {r['p_upper_bipartite_layer']:.3f}) -> {r['decision']}")
