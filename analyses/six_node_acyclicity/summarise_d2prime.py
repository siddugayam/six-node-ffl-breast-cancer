#!/usr/bin/env python3
"""Part (b) summary: D2' six-node counts on the census graph against NULL-C and NULL-L (see SETTINGS.md).

Reads runs/obs.tsv and runs/null_{C,L}_part*.tsv written by run_d2prime_nulls.sh. Checks that the replicates are the
graphs of the S1 runs: for every column of s1_null6.c the mean and SD over the replicates equal those of
analyses/six_node_pattern/S1/s1_null_summary.csv (relative tolerance 1e-9), and, if the folder of the S1 replicate files
is given as the first argument, every replicate row is identical to the stored one. Writes d2prime_counts_observed.csv,
d2prime_null_summary.csv, d2prime_decision.csv and check_replicates_against_s1.txt.
usage: FFL_ROOT=<repository checkout> python3 summarise_d2prime.py [folder with the S1 files null_C.tsv, null_L.tsv]"""
import csv, glob, os, sys
import numpy as np

ROOT = os.environ["FFL_ROOT"]
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
S1 = os.path.join(ROOT, "analyses", "six_node_pattern", "S1")
S1RUNS = sys.argv[1] if len(sys.argv) > 1 else None

o = [l.rstrip("\n").split("\t") for l in open(os.path.join(RUNS, "obs.tsv"))]
hdr = o[0]; obs = dict(zip(hdr, o[[r[0] for r in o].index("obs")]))
od = o[[r[0] for r in o].index("obs_d2prime")]
OBS = dict(seeds=int(od[1]), primary=int(od[2]), secondary=int(od[3]), census6c=int(obs["census6c"]),
           census6c_seeds=int(obs["census6c_seeds"]), bhat6_comp_inst=int(obs["bhat6_comp_inst"]))
with open(os.path.join(HERE, "d2prime_counts_observed.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["quantity", "value"])
    for k, lab in (("seeds", "D2' seeds (composite core + second TF, uncontracted census graph)"),
                   ("primary", "D2' six-node modules, primary (D3-D4 on the module with its cycle)"),
                   ("secondary", "D2' six-node modules, secondary (D3-D4 with the cycle condensed)"),
                   ("census6c_seeds", "CENSUS6c seeds (contracted census graph, s1_null6.c)"),
                   ("census6c", "CENSUS6c six-node topological-definition modules with the composite core plus TF2"),
                   ("bhat6_comp_inst", "six-node composite pattern instances (check)")):
        w.writerow([lab, OBS[k]])

S1SUM = list(csv.DictReader(open(os.path.join(S1, "s1_null_summary.csv"))))
log, summ, dec = [], [], {}
for m in "CL":
    rows, d2 = {}, {}
    for fn in sorted(glob.glob(os.path.join(RUNS, f"null_{m}_part*.tsv"))):
        for line in open(fn):
            f = line.rstrip("\n").split("\t")
            if f[0] == "d2p":
                if f[1].isdigit(): d2[int(f[1])] = (int(f[2]), int(f[3]), int(f[4]))
            elif f[0].isdigit():
                rows[int(f[0])] = f
    R = len(rows); assert set(rows) == set(d2) == set(range(R)), "replicates missing"
    arr = {c: np.array([float(rows[r][i]) for r in range(R)]) for i, c in enumerate(hdr) if c not in ("rep", "sec")}
    nm = {"C": "NULL-C", "L": "NULL-L"}[m]
    ok = True
    for s in S1SUM:
        if s["null"] != nm or s["column"] not in arr: continue
        mu, sd = arr[s["column"]].mean(), arr[s["column"]].std(ddof=1)
        good = abs(mu - float(s["null_mean"])) <= 1e-9 * max(1, abs(float(s["null_mean"]))) and \
            abs(sd - float(s["null_sd"])) <= 1e-9 * max(1, abs(float(s["null_sd"])))
        ok &= good
        log.append(f"{nm} {s['column']}: mean {mu:.6f} vs S1 {float(s['null_mean']):.6f}; sd {sd:.6f} vs S1 "
                   f"{float(s['null_sd']):.6f}: {'equal' if good else 'DIFFERENT'}")
    if S1RUNS:
        st = [l.rstrip("\n").split("\t") for l in open(os.path.join(S1RUNS, f"null_{m}.tsv"))][1:]
        same = sum(1 for f in st if f[:-1] == rows[int(f[0])][:-1])
        log.append(f"{nm}: {same} of {len(st)} stored S1 replicate rows identical (all columns except the timing)")
        ok &= same == len(st) == R
    log.append(f"{nm}: replicates {R}; check {'passed' if ok else 'FAILED'}")
    assert ok
    for k, col in ((0, "seeds"), (1, "primary"), (2, "secondary")):
        v = np.array([d2[r][k] for r in range(R)], float); ob = OBS[col]
        summ.append(dict(null=nm, count=f"D2' {col}", observed=ob, null_mean=v.mean(), null_sd=v.std(ddof=1),
                         fold=ob / v.mean() if v.mean() > 0 else float("nan"),
                         p_upper=(1 + int((v >= ob).sum())) / (R + 1), p_lower=(1 + int((v <= ob).sum())) / (R + 1),
                         replicates=R))
    v = arr["census6c"]; ob = OBS["census6c"]
    summ.append(dict(null=nm, count="CENSUS6c (the topological definition, contracted; comparison)", observed=ob, null_mean=v.mean(),
                     null_sd=v.std(ddof=1), fold=ob / v.mean(), p_upper=(1 + int((v >= ob).sum())) / (R + 1),
                     p_lower=(1 + int((v <= ob).sum())) / (R + 1), replicates=R))
with open(os.path.join(HERE, "d2prime_null_summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summ[0])); w.writeheader(); w.writerows(summ)
out = []
for cnt in ("D2' primary", "D2' secondary", "CENSUS6c (the topological definition, contracted; comparison)"):
    c = [s for s in summ if s["count"] == cnt and s["null"] == "NULL-C"][0]; l = [s for s in summ if s["count"] == cnt and s["null"] == "NULL-L"][0]
    if c["p_upper"] < 0.05 and l["p_upper"] < 0.05: d = "over-represented"
    elif c["p_lower"] < 0.05 and l["p_lower"] < 0.05: d = "depleted"
    else: d = "none of the three (not over-represented, not depleted)"
    out.append(dict(count=cnt, observed=c["observed"], NULL_C_fold=c["fold"], NULL_C_p_upper=c["p_upper"],
                    NULL_C_p_lower=c["p_lower"], NULL_L_fold=l["fold"], NULL_L_p_upper=l["p_upper"],
                    NULL_L_p_lower=l["p_lower"], decision=d, rule_applies=("yes" if cnt == "D2' primary" else "reported")))
with open(os.path.join(HERE, "d2prime_decision.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
open(os.path.join(HERE, "check_replicates_against_s1.txt"), "w").write("\n".join(log) + "\n")
print("\n".join(log[-4:])); [print(r) for r in out]
