#!/usr/bin/env python3
"""Summary of the STRING orientation analysis (see SETTINGS.md), from counts/ written by run_string_orientation.sh.

Checks: the listed-orientation graph reproduces the published census counts at n = 3-6 and the published n = 7
RAND-ESU samples seed by seed (analyses/census_and_motif_nulls/census/noexemplar/), and the no-STRING graph reproduces the
published no-STRING counts (census/noexemplar_nostring/). Writes:
  gene_gene_layer_variants.csv           counts per n for the census graph, TRRUST-only gene-gene links, no gene-gene layer
  orientation_robust_core.csv            per n: modules with no STRING arc / no STRING pair (n = 7: RAND-ESU estimate)
  string_bootstrap_counts.csv            per replicate and n: modules, modules also in the census, Jaccard
  string_bootstrap_summary.csv           per n: median, 2.5 and 97.5 % percentiles, mean, SD, CV, percentile of the
                                         published count, and the recovery of the published modules
usage: FFL_ROOT=<repository checkout> python3 summarise_string_orientation.py"""
import csv, glob, os, statistics as st
import numpy as np

ROOT = os.environ["FFL_ROOT"]
HERE = os.path.dirname(os.path.abspath(__file__))
C = os.path.join(HERE, "counts")
PUBC = os.path.join(ROOT, "analyses", "census_and_motif_nulls", "census")
SEEDS = (11, 22, 33, 44, 55, 71, 72, 73)


def parse(fn):
    d = {}
    for line in open(fn):
        w = line.split()
        if len(w) >= 2:
            try: d[w[0]] = float(w[1]) if w[0] == "prod_q" else int(w[1])
            except ValueError: d[w[0]] = " ".join(w[1:])
    return d


def total(g, n):
    fs = [f"{C}/{g}/n{n}.txt"] if n < 6 else sorted(glob.glob(f"{C}/{g}/n6_part_*.txt"))
    if n == 6: assert len(fs) == 64
    out = {}
    for fn in fs:
        for k, v in parse(fn).items():
            if isinstance(v, int) and k in ("found", "visited", "no_string_arc", "no_string_pair"): out[k] = out.get(k, 0) + v
    return out


def published(variant, n):
    fs = [f"{PUBC}/{variant}/n{n}.txt"] if n < 6 else sorted(glob.glob(f"{PUBC}/{variant}/n6_parts/part_*.txt"))
    return sum(parse(f)["found"] for f in fs)


rows, core = [], []
for g, desc in (("listed", "census graph (STRING in the listed orientation)"),
                ("nostring", "TRRUST-only gene-gene links (census graph minus STRING arcs)"),
                ("nogenegene", "no gene-gene layer (minus STRING and the 604 TRRUST-layer arcs)")):
    for n in (3, 4, 5, 6):
        t = total(g, n)
        if g == "listed": assert t["found"] == published("noexemplar", n), (g, n)
        if g == "nostring": assert t["found"] == published("noexemplar_nostring", n), (g, n)
        rows.append(dict(graph=g, description=desc, n=n, method="exhaustive", modules=t["found"], sd=""))
        if g == "listed":
            core.append(dict(n=n, method="exhaustive", modules=t["found"], no_string_arc=t["no_string_arc"],
                             fraction_no_string_arc=t["no_string_arc"] / t["found"], no_string_pair=t["no_string_pair"],
                             fraction_no_string_pair=t["no_string_pair"] / t["found"], sd_fraction_no_string_arc="",
                             sd_fraction_no_string_pair=""))
    ds = [parse(f"{C}/{g}/n7_seed_{s}.txt") for s in SEEDS]
    if g in ("listed", "nostring"):
        pv = "noexemplar" if g == "listed" else "noexemplar_nostring"
        for s, d in zip(SEEDS, ds):
            p = parse(f"{PUBC}/{pv}/n7_seeds/seed_{s}.txt"); assert (d["found"], d["visited"]) == (p["found"], p["visited"]), (g, s)
    est = [d["found"] / d["prod_q"] for d in ds]
    rows.append(dict(graph=g, description=desc, n=7, method="RAND-ESU, 8 seeds", modules=f"{st.mean(est):.1f}",
                     sd=f"{st.stdev(est):.1f}"))
    if g == "listed":
        fa = [d["no_string_arc"] / d["found"] for d in ds]; fp = [d["no_string_pair"] / d["found"] for d in ds]
        E = st.mean(est)
        core.append(dict(n=7, method="RAND-ESU, 8 seeds", modules=f"{E:.1f}", no_string_arc=f"{E * st.mean(fa):.1f}",
                         fraction_no_string_arc=st.mean(fa), no_string_pair=f"{E * st.mean(fp):.1f}",
                         fraction_no_string_pair=st.mean(fp), sd_fraction_no_string_arc=st.stdev(fa),
                         sd_fraction_no_string_pair=st.stdev(fp)))
for fn, R in (("gene_gene_layer_variants.csv", rows), ("orientation_robust_core.csv", core)):
    with open(os.path.join(HERE, fn), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(R[0])); w.writeheader(); w.writerows(R)

# ---- bootstrap -------------------------------------------------------------------------------------------------------
base = {n: total("listed", n)["found"] for n in (3, 4, 5, 6)}
graphs = {int(r["replicate"]): r for r in csv.DictReader(open(os.path.join(HERE, "string_bootstrap_graphs.csv")))}
B, S = [], []
for n in (3, 4, 5, 6):
    found, rec = [], []
    for r in sorted(graphs):
        fn = f"{C}/boot/boot_{r:03d}_n{n}.txt"
        if not os.path.exists(fn): continue
        d = parse(fn)
        if "found" not in d: continue
        assert d["baseline_keys"] == base[n]
        found.append(d["found"]); rec.append(d["in_baseline"])
        B.append(dict(replicate=r, n=n, string_arcs_added=graphs[r]["string_arcs_added"], modules=d["found"],
                      also_in_census=d["in_baseline"], recovery_of_census_modules=d["in_baseline"] / base[n],
                      jaccard=d["in_baseline"] / (base[n] + d["found"] - d["in_baseline"])))
    f = np.array(found, float); rc = np.array(rec, float) / base[n]
    jac = np.array(rec, float) / (base[n] + f - np.array(rec, float))
    S.append(dict(n=n, replicates=len(f), published_count=base[n], median=float(np.median(f)),
                  p2_5=float(np.percentile(f, 2.5)), p97_5=float(np.percentile(f, 97.5)), mean=float(f.mean()),
                  sd=float(f.std(ddof=1)), cv=float(f.std(ddof=1) / f.mean()),
                  published_percentile=float(100 * (f < base[n]).mean() + 50 * (f == base[n]).mean()),
                  recovery_median=float(np.median(rc)), recovery_p2_5=float(np.percentile(rc, 2.5)),
                  recovery_p97_5=float(np.percentile(rc, 97.5)), jaccard_median=float(np.median(jac)),
                  jaccard_p2_5=float(np.percentile(jac, 2.5)), jaccard_p97_5=float(np.percentile(jac, 97.5))))
with open(os.path.join(HERE, "string_bootstrap_counts.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(B[0])); w.writeheader(); w.writerows(B)
with open(os.path.join(HERE, "string_bootstrap_summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(S[0])); w.writeheader(); w.writerows(S)
for r in rows + core + S: print(r)
