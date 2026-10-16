#!/usr/bin/env python3
"""Summary of the tier- and sign-aware null models (see SETTINGS.md).

Controls: the NULL-A/B/C runs reproduce the published replicate files line for line
(analyses/census_and_motif_nulls/nulls/null_{A,B,C}_{dep,pub}_noexemplar.tsv, published typing;
analyses/six_node_pattern/E/E5/null_{A,B,C}_dep_noexemplar_retyped22.tsv, retyped analysed network).
Writes control_reproduction.csv, null_model_summary.csv (observed, expected, SD, fold, excess %, empirical p for every
run and metric) and waterfall.csv (excess % from NULL-A to NULL-C + tier + sign, with the step changes).
usage: FFL_ROOT=<repository checkout> python3 summarise_tier_sign_nulls.py"""
import csv, os
import numpy as np

ROOT = os.environ["FFL_ROOT"]
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
MODELS = {"A": "NULL-A", "B": "NULL-B", "C": "NULL-C", "BT": "(i) NULL-B + tier", "BS": "(ii) NULL-B + sign",
          "BTS": "(iii) NULL-B + tier + sign", "CT": "NULL-C + tier", "CTS": "NULL-C + tier + sign"}
GRAPHS = {"analysed": "analysed network (6,829 arcs)", "census": "census graph (9,226 arcs)"}
TYPINGS = {"published": "published typing (157 TFs)", "retyped": "22 co-regulators retyped as non-TF (135 TFs)"}
METRICS = {"total": "all three-node FFLs", "comp": "composite FFLs", "tf": "TF-FFLs", "mir": "miRNA-FFLs"}


def read(fn):
    L = [l.rstrip("\n").split("\t") for l in open(fn)]
    hdr = L[0]; obs = dict(zip(hdr[1:], map(int, L[1][1:])))
    arr = np.array([[int(x) for x in r[1:]] for r in L[2:]], dtype=float)
    return hdr, obs, arr


ctrl = []
pub = {("analysed", "published"): ("analyses/census_and_motif_nulls/nulls", "null_{m}_dep_noexemplar.tsv"),
       ("census", "published"): ("analyses/census_and_motif_nulls/nulls", "null_{m}_pub_noexemplar.tsv"),
       ("analysed", "retyped"): ("analyses/six_node_pattern/E/E5", "null_{m}_dep_noexemplar_retyped22.tsv")}
for (g, t), (d, pat) in pub.items():
    for m in "ABC":
        a = open(os.path.join(RUNS, f"null_{m}_{g}_{t}.tsv")).read()
        b = open(os.path.join(ROOT, d, pat.format(m=m))).read()
        ctrl.append(dict(graph=g, typing=t, null=MODELS[m], published_file=f"{d}/{pat.format(m=m)}",
                         replicate_lines=a.count("\n"), identical=a == b))
assert all(c["identical"] for c in ctrl), "a control run differs from the published file"
with open(os.path.join(HERE, "control_reproduction.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(ctrl[0])); w.writeheader(); w.writerows(ctrl)

rows, ex = [], {}
for g in GRAPHS:
    for t in TYPINGS:
        for m in MODELS:
            hdr, obs, arr = read(os.path.join(RUNS, f"null_{m}_{g}_{t}.tsv"))
            col = {c: arr[:, i] for i, c in enumerate(hdr[1:])}
            R = arr.shape[0]; assert R == 1000
            for k, lab in METRICS.items():
                v = col[k]; o = obs[k]; mu, sd = v.mean(), v.std(ddof=1)
                fold = o / mu if mu > 0 else float("nan")
                r = dict(graph=g, typing=t, null=MODELS[m], null_code=m, metric=lab, observed=o, expected=mu, sd=sd,
                         fold=fold, excess_pct=100 * (fold - 1) if mu > 0 else float("nan"),
                         z=(o - mu) / sd if sd > 0 else float("nan"),
                         p_upper=(1 + int((v >= o).sum())) / (R + 1), p_lower=(1 + int((v <= o).sum())) / (R + 1),
                         replicates=R, reciprocal_pairs_mean=col["recip"].mean())
                rows.append(r); ex[(g, t, m, k)] = r["excess_pct"]
with open(os.path.join(HERE, "null_model_summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

STEPS = [("A", "NULL-A"), ("B", "+ reciprocity (NULL-B)"), ("C", "+ bipartite curveball (NULL-C)"),
         ("CT", "+ tier (NULL-C + tier)"), ("CTS", "+ sign (NULL-C + tier + sign)")]
wf = []
for g in GRAPHS:
    for t in TYPINGS:
        for k, lab in METRICS.items():
            prev = None
            for code, step in STEPS:
                e = ex[(g, t, code, k)]
                wf.append(dict(graph=g, typing=t, metric=lab, step=step, null_code=code, excess_pct=e,
                               change_pp=(e - prev) if prev is not None else float("nan")))
                prev = e
with open(os.path.join(HERE, "waterfall.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(wf[0])); w.writeheader(); w.writerows(wf)
for r in rows:
    if r["metric"] in ("all three-node FFLs", "composite FFLs"):
        print(f"{r['graph']:9s} {r['typing']:10s} {r['null']:28s} {r['metric']:20s} obs {r['observed']:5d} exp {r['expected']:8.1f} "
              f"fold {r['fold']:.3f} ({r['excess_pct']:+.1f} %) p_up {r['p_upper']:.4f} p_low {r['p_lower']:.4f}")
print("controls identical:", all(c["identical"] for c in ctrl))
