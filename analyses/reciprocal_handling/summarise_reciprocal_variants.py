#!/usr/bin/env python3
"""Summary of the census under four treatments of reciprocal TF<->miRNA pairs (see SETTINGS.md).

Reads counts/{base,unc,mod,scc}/*.txt and the module keys in work/keys/ written by run_reciprocal_variants.sh, the SCC
file work/scc_uncontracted.txt and the census graphs (root: environment variable FFL_ROOT). Writes:
  census_reciprocal_variants.csv                        one row per n x variant
  lost_under_module_condensation_first_condition.csv    baseline modules lost under M, by first failing condition
  composite_three_node_under_module_condensation.csv    the composite three-node FFLs under M
usage: FFL_ROOT=<repository checkout> python3 summarise_reciprocal_variants.py"""
import csv, glob, os, statistics as st
import numpy as np

ROOT = os.environ["FFL_ROOT"]
HERE = os.path.dirname(os.path.abspath(__file__))
C = os.path.join(HERE, "counts"); K = os.path.join(HERE, "work", "keys")
G_DIR = os.path.join(ROOT, "analyses", "census_and_motif_nulls", "graphs")
PUBLISHED = {3: 5833, 4: 116505, 5: 1506263, 6: 19401840}
SEEDS = (11, 22, 33, 44, 55, 71, 72, 73)


def parse(fn):
    d = {}
    for line in open(fn):
        w = line.split()
        if len(w) >= 2:
            try: d[w[0]] = float(w[1]) if w[0] == "prod_q" else int(w[1])
            except ValueError: d[w[0]] = " ".join(w[1:])
    return d


def total(variant, n):
    files = [f"{C}/{variant}/n{n}.txt"] if n < 6 or variant == "scc" else sorted(glob.glob(f"{C}/{variant}/n6_part_*.txt"))
    if variant != "scc" and n == 6: assert len(files) == 64, (variant, len(files))
    out = {}
    for fn in files:
        for k, v in parse(fn).items():
            if isinstance(v, int) and k not in ("K", "nparts", "part", "seed", "rng_final"): out[k] = out.get(k, 0) + v
    return out


def keys(variant, n, ext="keys"):
    fs = [f"{K}/{variant}/n{n}.{ext}"] if n < 6 or variant == "scc" else \
        [f"{K}/{variant}/n6_part_{p}.{ext}" for p in range(64)]
    dt = np.uint8 if ext == "flags" else np.uint64
    return np.concatenate([np.fromfile(f, dtype=dt) for f in fs])


sccL = open(os.path.join(HERE, "work", "scc_uncontracted.txt")).read().split("\n")
scc_of = np.array(list(map(int, sccL[1].split())), dtype=np.uint64)

rows, lost_rows = [], []
for n in (3, 4, 5, 6):
    b = total("base", n); u = total("unc", n); m = total("mod", n); s = total("scc", n)
    assert b["found"] == PUBLISHED[n], f"baseline n={n} does not reproduce the published count"
    kb = keys("base", n); fl = keys("base", n, "flags"); ku = keys("unc", n); km = keys("mod", n); ks = keys("scc", n)
    assert len(kb) == b["found"] and len(np.unique(kb)) == len(kb)
    assert len(km) == m["found"] and len(np.unique(km)) == len(km)
    # U: exact comparison of instance sets
    shared_u = len(np.intersect1d(kb, ku, assume_unique=True))
    assert shared_u == len(ku) == int(((fl & 1) > 0).sum()) == u["found"] == b["valid_U"]
    assert set(kb[(fl & 1) > 0].tolist()) == set(ku.tolist())
    # M
    shared_m = len(np.intersect1d(kb, km, assume_unique=True))
    assert shared_m == b["valid_M"] == m["shared_with_baseline"] == int(((fl & 4) > 0).sum())
    # S: images of baseline modules
    ids = np.stack([(kb >> np.uint64(10 * (n - 1 - i))) & np.uint64(1023) for i in range(n)], axis=1)
    img = np.sort(scc_of[ids.astype(np.int64)], axis=1)
    distinct = (np.diff(img, axis=1) > 0).all(axis=1)
    ik = np.zeros(len(kb), dtype=np.uint64)
    for i in range(n): ik = (ik << np.uint64(10)) | img[:, i]
    imgs_all = np.unique(ik[distinct]); ks_u = np.unique(ks)
    shared_s = len(np.intersect1d(imgs_all, ks_u, assume_unique=True))
    retained = distinct & np.isin(ik, ks_u)
    assert int(retained.sum()) == b["retained_S"] == int(((fl & 2) > 0).sum())
    assert len(ks_u) == s["found"]
    B = b["found"]
    for var, X, sh, jac, lostf, gainf, meth in (
            ("U no contraction", u["found"], shared_u, shared_u / (B + u["found"] - shared_u), 1 - shared_u / B,
             (u["found"] - shared_u) / u["found"] if u["found"] else float("nan"), "exhaustive"),
            ("S graph-level SCC condensation", s["found"], shared_s,
             shared_s / (len(imgs_all) + s["found"] - shared_s), 1 - b["retained_S"] / B,
             (s["found"] - shared_s) / s["found"] if s["found"] else float("nan"), "exhaustive (images)"),
            ("M module-level condensation", m["found"], shared_m, shared_m / (B + m["found"] - shared_m),
             1 - shared_m / B, (m["found"] - shared_m) / m["found"], "exhaustive")):
        rows.append(dict(n=n, variant=var, method=meth, baseline_count=B, variant_count=X, shared=sh,
                         jaccard=jac, lost_fraction=lostf, gained_fraction=gainf, variant_sd="", baseline_sd="",
                         note=(f"baseline modules retained: {b['retained_S']}; distinct images of baseline modules "
                               f"with n distinct SCCs: {len(imgs_all)}" if var.startswith("S") else "")))
    lost_rows += [dict(n=n, baseline_modules=B, lost_under_M=B - shared_m, first_failing_D3=b["M_fail_D3"],
                       first_failing_D4=b["M_fail_D4"], lost_under_U_all_by_D2=B - shared_u,
                       lost_under_S_share_an_SCC=int((~distinct).sum()),
                       lost_under_S_image_not_a_module=int((distinct & ~retained).sum()))]
    print(f"n={n}: B {B}  U {u['found']}  M {m['found']} (shared {shared_m})  S {s['found']} (shared {shared_s})", flush=True)

# n = 7: RAND-ESU, 8 seeds
def seeds(variant):
    return [parse(f"{C}/{variant}/n7_seed_{sd}.txt") for sd in SEEDS]


bs, us, ms = seeds("base"), seeds("unc"), seeds("mod"); s7 = parse(f"{C}/scc/n7.txt")
PUB7 = os.path.join(ROOT, "analyses", "census_and_motif_nulls", "census", "noexemplar", "n7_seeds")
for sd, d in zip(SEEDS, bs):        # the baseline samples are the published ones
    p = parse(os.path.join(PUB7, f"seed_{sd}.txt"))
    assert (d["found"], d["visited"], d["rng_final"]) == (p["found"], p["visited"], p["rng_final"]), f"seed {sd} differs"
eb = [d["found"] / d["prod_q"] for d in bs]; eu = [d["found"] / d["prod_q"] for d in us]; em = [d["found"] / d["prod_q"] for d in ms]
fU = [d["valid_U"] / d["found"] for d in bs]; fS = [d["retained_S"] / d["found"] for d in bs]
fM = [d["valid_M"] / d["found"] for d in bs]; gM = [1 - d["shared_with_baseline"] / d["found"] for d in ms]
EB, EU, EM = st.mean(eb), st.mean(eu), st.mean(em)
for var, EX, sdX, f_ret, gain in (("U no contraction", EU, st.stdev(eu), st.mean(fU), 0.0),
                                  ("S graph-level SCC condensation", s7["found"], 0.0, st.mean(fS), float("nan")),
                                  ("M module-level condensation", EM, st.stdev(em), st.mean(fM), st.mean(gM))):
    sh = EB * f_ret
    rows.append(dict(n=7, variant=var, method="RAND-ESU, 8 seeds (S: exhaustive on the condensation DAG)",
                     baseline_count=f"{EB:.1f}", variant_count=f"{EX:.1f}", shared=f"{sh:.1f}",
                     jaccard=sh / (EB + EX - sh) if var[0] != "S" else float("nan"), lost_fraction=1 - f_ret,
                     gained_fraction=gain, variant_sd=f"{sdX:.1f}", baseline_sd=f"{st.stdev(eb):.1f}",
                     note=("Jaccard on images not estimable from the samples; shared = baseline estimate x retained "
                           "fraction" if var[0] == "S" else "shared = baseline estimate x fraction of sampled baseline "
                           "modules valid under the variant")))
with open(os.path.join(HERE, "census_reciprocal_variants.csv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
with open(os.path.join(HERE, "lost_under_module_condensation_first_condition.csv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(lost_rows[0])); w.writeheader(); w.writerows(lost_rows)

# composite three-node FFLs under M: the reciprocal TF<->miRNA pair with both regulating the third node, read in the
# contracted graph as count3() of motif_null_models.c does (1,375 on the census graph)
L = open(os.path.join(G_DIR, "census_noexemplar_uncontracted.txt")).read().split("\n")
nv, ne, _ = map(int, L[0].split()); ty = list(map(int, L[1].split()))
U = set(tuple(map(int, l.split()[:2])) for l in L[2:2 + ne])
LC = open(os.path.join(G_DIR, "census_noexemplar_fine.txt")).read().split("\n")
CA = set(tuple(map(int, l.split()[:2])) for l in LC[2:2 + int(LC[0].split()[1])])
kb3 = keys("base", 3); fl3 = keys("base", 3, "flags")
comp = fail4 = fail3 = valid = 0
for k, f in zip(kb3.tolist(), fl3.tolist()):
    a, b, c = (k >> 20) & 1023, (k >> 10) & 1023, k & 1023
    S = (a, b, c); iscomp = False
    for x in S:
        for y in S:
            if x != y and ty[x] == 1 and ty[y] == 0 and (x, y) in U and (y, x) in U:
                z = [w for w in S if w not in (x, y)][0]
                if (x, z) in CA and (y, z) in CA: iscomp = True      # TF and miRNA both regulate z (contracted graph)
    if iscomp:
        comp += 1
        if f & 4: valid += 1
        elif (f >> 4) == 1: fail4 += 1
        elif (f >> 4) == 3: fail3 += 1
with open(os.path.join(HERE, "composite_three_node_under_module_condensation.csv"), "w") as fh:
    w = csv.writer(fh); w.writerow(["composite_three_node_FFLs", "valid_under_M", "fail_D4_under_M", "fail_D3_under_M"])
    w.writerow([comp, valid, fail4, fail3])
print("composite three-node FFLs", comp, "valid under M", valid, "fail D4", fail4, "fail D3", fail3)
for r in rows: print(r)
