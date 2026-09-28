#!/usr/bin/env python3
"""
TargetScan 8.0 sequence-level evidence for the miRNA->target edges of the network.

Downloaded from https://www.targetscan.org/vert_80/vert_80_data_download/ and filtered to
human (Gene Tax ID / Species ID == 9606):
  Predicted_Targets_Context_Scores.default_predictions.txt  - TargetScan's DEFAULT prediction
        set (sites of conserved miRNA families)
  Conserved_Site_Context_Scores.txt                          - CONSERVED sites
  Nonconserved_Site_Context_Scores.txt                       - NON-conserved sites (optional;
        542 MB zip, fetched with a 12-way parallel range download)
  Predicted_Targets_Info.default_predictions.txt             - seed-match string + PCT
  miR_Family_Info.txt                                        - mature miRNA -> family

Site Type codes in the context-score files were verified against the "Seed match" column of
Predicted_Targets_Info by joining on (transcript, UTR end):
  1 = 7mer-A1, 2 = 7mer-m8, 3 = 8mer, negative codes = 3'-compensatory sites.

Three nested site definitions are scored for every miRNA-target pair:
  conserved : sites in Conserved_Site_Context_Scores
  default   : sites in TargetScan's default prediction set
  all       : conserved + default + non-conserved
"""
import os, re, sys, csv, json
from collections import defaultdict
import numpy as np
from scipy import stats

REV = "/path/to/revision"
TS = "/path/to/scratch/ts"
OUT = os.path.join(REV, "results", "v2")
SITE = {"1": "7mer-A1", "2": "7mer-m8", "3": "8mer"}

def load_ctx(path, tag, store, seen):
    n = 0
    if not os.path.exists(path) or os.path.getsize(path) < 1000:
        sys.stderr.write(f"WARNING: {path} missing or empty - skipped\n")
        return 0
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or p[3] != "9606":
                continue
            key = (p[4], p[1], p[2], p[6], p[7])
            r = store.get(key)
            if r is None:
                r = store[key] = dict(mirna=p[4], gene=p[1], tx=p[2], start=p[6], end=p[7],
                                      stype=SITE.get(p[5], "3comp"), ctx=p[8], ctxp=p[9],
                                      wctx=p[10], wctxp=p[11],
                                      conserved=False, default=False, nonconserved=False)
            r[tag] = True
            seen.add(p[4]); n += 1
    return n

sites, mirnas_seen = {}, set()
n_def = load_ctx(os.path.join(TS, "Predicted_Targets_Context_Scores.default_predictions.txt"), "default", sites, mirnas_seen)
n_con = load_ctx(os.path.join(TS, "Conserved_Site_Context_Scores.txt"), "conserved", sites, mirnas_seen)
n_non = load_ctx(os.path.join(TS, "Nonconserved_Site_Context_Scores.txt"), "nonconserved", sites, mirnas_seen)
print(f"human site rows read: default={n_def} conserved={n_con} nonconserved={n_non}; "
      f"unique sites={len(sites)}; distinct miRNAs={len(mirnas_seen)}")

pct, fam, famcons = {}, {}, {}
with open(os.path.join(TS, "Predicted_Targets_Info.default_predictions.txt")) as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if p[4] == "9606":
            pct[(p[0], p[3], p[6])] = p[10]
with open(os.path.join(TS, "miR_Family_Info.txt")) as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if p[2] == "9606":
            fam[p[3]] = p[0]; famcons[p[3]] = p[5]

by_pair = defaultdict(list)
for v in sites.values():
    by_pair[(v["mirna"], v["gene"])].append(v)

all_ts_mirnas = sorted(mirnas_seen | set(fam))
_armcache = {}
def arms_for(label):
    if label not in _armcache:
        pat = re.compile(r"^" + re.escape(label) + r"(-\d+)?(-(3p|5p))?(\.\d+)?$")
        _armcache[label] = [m for m in all_ts_mirnas if pat.match(m)]
    return _armcache[label]

def fnum(x):
    try: return float(x)
    except Exception: return float("nan")

def agg(recs, sel):
    rs = [r for r in recs if sel(r)]
    if not rs:
        return dict(n=0, n8=0, n7m8=0, n7a1=0, noth=0, best="", bestp="", tot="", bestw="", pct="")
    c = [fnum(r["ctx"]) for r in rs]
    i = int(np.nanargmin(c))
    ps = []
    for r in rs:
        v = pct.get((fam.get(r["mirna"], ""), r["tx"], r["end"]))
        if v not in (None, "NULL", ""):
            try: ps.append(float(v))
            except Exception: pass
    return dict(n=len(rs), n8=sum(1 for r in rs if r["stype"] == "8mer"),
                n7m8=sum(1 for r in rs if r["stype"] == "7mer-m8"),
                n7a1=sum(1 for r in rs if r["stype"] == "7mer-A1"),
                noth=sum(1 for r in rs if r["stype"] == "3comp"),
                best=rs[i]["ctx"], bestp=rs[i]["ctxp"], tot=round(float(np.nansum(c)), 4),
                bestw=min((fnum(r["wctx"]) for r in rs), default=""),
                pct=(max(ps) if ps else ""))

def summarise(label, gene):
    arms = arms_for(label)
    recs = []
    for a in arms:
        recs.extend(by_pair.get((a, gene), []))
    A = agg(recs, lambda r: True)
    D = agg(recs, lambda r: r["default"])
    C = agg(recs, lambda r: r["conserved"])
    detail = ";".join(f"{r['mirna']}|{r['stype']}|{r['start']}-{r['end']}|ctx={r['ctx']}|"
                      f"pctile={r['ctxp']}|cons={'Y' if r['conserved'] else 'N'}|"
                      f"default={'Y' if r['default'] else 'N'}"
                      for r in sorted(recs, key=lambda r: fnum(r["ctx"]))[:12])
    return dict(
        miRNA_in_targetscan=str(len(arms) > 0).upper(),
        miRNA_has_any_prediction=str(any((a, g) in by_pair for a in arms for g in [gene]) or
                                     any(a in mirnas_seen for a in arms)).upper(),
        arms_matched=";".join(arms),
        n_sites=A["n"], n_8mer=A["n8"], n_7mer_m8=A["n7m8"], n_7mer_A1=A["n7a1"], n_3comp=A["noth"],
        best_ctx=A["best"], best_ctx_pct=A["bestp"], sum_ctx=A["tot"], best_wctx=A["bestw"], max_PCT=A["pct"],
        n_sites_default=D["n"], best_ctx_default=D["best"], sum_ctx_default=D["tot"],
        n_sites_conserved=C["n"], n_8mer_conserved=C["n8"], best_ctx_conserved=C["best"],
        best_ctx_pct_conserved=C["bestp"], sum_ctx_conserved=C["tot"], max_PCT_conserved=C["pct"],
        conserved_family=";".join(sorted({famcons.get(r["mirna"], "") for r in recs})),
        arms_with_sites=";".join(sorted({r["mirna"] for r in recs})),
        site_detail=detail)

# ---------------- inputs ----------------
edges = [r for r in csv.DictReader(open(os.path.join(REV, "data", "canonical_edges.tsv")), delimiter="\t")]
mir_edges = [(r["source"], r["target"]) for r in edges if r["edge_type"] == "miRNA_target"]
mir_edge_set = set(mir_edges)
rho, tier = {}, {}
for r in csv.DictReader(open(os.path.join(REV, "results", "v2", "v2_edge_rho.csv"))):
    if r["edge_type"] == "miRNA_target":
        try: rho[(r["source"], r["target"])] = float(r["rho"])
        except Exception: pass
        tier[(r["source"], r["target"])] = r["tier"]

NAMED = [("hsa-miR-29a", "COL1A1"), ("hsa-miR-29b", "COL1A1"), ("hsa-miR-29c", "COL1A1"),
         ("hsa-miR-29a", "COL3A1"), ("hsa-miR-29b", "COL3A1"), ("hsa-miR-29c", "COL3A1"),
         ("hsa-let-7b", "COL3A1"), ("hsa-let-7e", "COL3A1"),
         ("hsa-let-7b", "COL1A1"), ("hsa-let-7e", "COL1A1"),
         ("hsa-miR-101", "EZH2"), ("hsa-miR-130a", "VEGFA"),
         ("hsa-miR-124", "STAT3"), ("hsa-let-7b", "HK2")]
named_set = set(NAMED)

rows = []
def rowfor(label, gene, analysis):
    s = summarise(label, gene)
    s.update(analysis=analysis, miRNA=label, target_gene=gene,
             edge_in_network=str((label, gene) in mir_edge_set).upper(),
             tier=tier.get((label, gene), ""), tcga_rho=rho.get((label, gene), ""),
             has_targetscan_site=str(s["n_sites"] > 0).upper())
    rows.append(s); return s

print("\n=== NAMED miRNA-TARGET PAIRS (TargetScan 8.0, human) ===")
for lab, g in NAMED:
    s = rowfor(lab, g, "named_pair")
    print(f"{lab} -> {g}: total sites={s['n_sites']} (8mer={s['n_8mer']}, 7mer-m8={s['n_7mer_m8']}, "
          f"7mer-A1={s['n_7mer_A1']}, 3'comp={s['n_3comp']}); conserved sites={s['n_sites_conserved']}; "
          f"best ctx++={s['best_ctx']} (pctile {s['best_ctx_pct']}); total ctx++={s['sum_ctx']}; "
          f"max PCT={s['max_PCT']}; in_network={s['edge_in_network']}; tier={s['tier']}; rho={s['tcga_rho']}")
    for d in (s["site_detail"].split(";") if s["site_detail"] else []):
        print("      ", d)

for lab, g in mir_edges:
    if (lab, g) not in named_set:
        rowfor(lab, g, "network_miRNA_target_edge")

flds = ["analysis", "miRNA", "target_gene", "edge_in_network", "tier", "tcga_rho",
        "has_targetscan_site", "miRNA_in_targetscan", "arms_matched", "arms_with_sites",
        "n_sites", "n_8mer", "n_7mer_m8", "n_7mer_A1", "n_3comp",
        "best_ctx", "best_ctx_pct", "sum_ctx", "best_wctx", "max_PCT",
        "n_sites_default", "best_ctx_default", "sum_ctx_default",
        "n_sites_conserved", "n_8mer_conserved", "best_ctx_conserved", "best_ctx_pct_conserved",
        "sum_ctx_conserved", "max_PCT_conserved", "conserved_family", "site_detail"]
with open(os.path.join(OUT, "targetscan_sites.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=flds, extrasaction="ignore")
    w.writeheader(); w.writerows(rows)
print("\ntargetscan_sites.csv rows:", len(rows))

net = [r for r in rows if r["edge_in_network"] == "TRUE"]
print("network miRNA_target edges:", len(mir_edges), "| rows for them:", len(net))
mirs = sorted({m for m, g in mir_edges})
print(f"network miRNAs matched to a TargetScan mature miRNA: {sum(1 for m in mirs if arms_for(m))}/{len(mirs)}")
print(f"network miRNAs with >=1 predicted site anywhere: "
      f"{len({r['miRNA'] for r in net if r['n_sites']>0})}/{len(mirs)}")
for lab, key in (("any site", "n_sites"), ("default-set site", "n_sites_default"),
                 ("conserved site", "n_sites_conserved")):
    k = sum(1 for r in net if r[key] > 0)
    print(f"edges with >=1 {lab}: {k} ({100*k/len(net):.1f}%)")

# ---------------- does site quality predict TCGA sign concordance? ----------------
res = []
def analyse(sub, name, bestkey="best_ctx", nkey="n_sites"):
    x = [r for r in sub if r[nkey] > 0 and r["tcga_rho"] != "" and r[bestkey] != ""]
    y = [r for r in sub if r["tcga_rho"] != ""]
    if len(x) < 20:
        return
    ctx = np.array([float(r[bestkey]) for r in x])
    rr = np.array([float(r["tcga_rho"]) for r in x])
    conc = rr < 0
    u = stats.mannwhitneyu(ctx[conc], ctx[~conc], alternative="two-sided")
    sp = stats.spearmanr(ctx, rr)
    q = np.quantile(ctx, [0, .25, .5, .75, 1.0])
    qs = []
    for i in range(4):
        m = (ctx >= q[i]) & (ctx <= q[i+1]) if i == 3 else (ctx >= q[i]) & (ctx < q[i+1])
        qs.append((int(m.sum()), float(conc[m].mean()) if m.sum() else float("nan"),
                   float(rr[m].mean()) if m.sum() else float("nan")))
    ws = np.array([float(r["tcga_rho"]) for r in y if r[nkey] > 0])
    ns = np.array([float(r["tcga_rho"]) for r in y if r[nkey] == 0])
    sv = stats.mannwhitneyu(ws, ns, alternative="two-sided") if len(ns) >= 20 else None
    pv = (stats.fisher_exact([[int((ws < 0).sum()), int((ws >= 0).sum())],
                              [int((ns < 0).sum()), int((ns >= 0).sum())]])
          if len(ns) >= 20 else None)
    d = dict(subset=name, site_definition=nkey, n_with_site=len(x),
             frac_concordant=round(float(conc.mean()), 4),
             mean_ctx_concordant=round(float(ctx[conc].mean()), 4),
             mean_ctx_discordant=round(float(ctx[~conc].mean()), 4),
             MWU_p=float(u.pvalue),
             spearman_ctx_vs_rho=round(float(sp.statistic), 4), spearman_p=float(sp.pvalue),
             n_q1=qs[0][0], conc_q1=round(qs[0][1], 4), rho_q1=round(qs[0][2], 4),
             n_q2=qs[1][0], conc_q2=round(qs[1][1], 4), rho_q2=round(qs[1][2], 4),
             n_q3=qs[2][0], conc_q3=round(qs[2][1], 4), rho_q3=round(qs[2][2], 4),
             n_q4=qs[3][0], conc_q4=round(qs[3][1], 4), rho_q4=round(qs[3][2], 4),
             n_no_site=len(ns),
             frac_conc_with_site=round(float((ws < 0).mean()), 4) if len(ws) else "",
             frac_conc_no_site=round(float((ns < 0).mean()), 4) if len(ns) else "",
             mean_rho_with_site=round(float(ws.mean()), 4) if len(ws) else "",
             mean_rho_no_site=round(float(ns.mean()), 4) if len(ns) else "",
             site_vs_nosite_MWU_p=(float(sv.pvalue) if sv else ""),
             site_vs_nosite_fisher_p=(float(pv[1]) if pv else ""),
             site_vs_nosite_OR=(round(float(pv[0]), 4) if pv else ""))
    res.append(d)
    print(f"\n[{name} | {nkey}] n_with_site={len(x)} concordant={100*conc.mean():.1f}%")
    print(f"   mean best ctx++: concordant {ctx[conc].mean():.4f} vs discordant {ctx[~conc].mean():.4f} (MWU p={u.pvalue:.3g})")
    print(f"   Spearman(best ctx++, rho) = {sp.statistic:.4f}, p={sp.pvalue:.3g}")
    print("   by ctx++ quartile (Q1 = strongest sites): " +
          ", ".join(f"Q{i+1} {qs[i][1]*100:.1f}% conc / mean rho {qs[i][2]:+.4f} (n={qs[i][0]})" for i in range(4)))
    if sv:
        print(f"   with site n={len(ws)}: {100*(ws<0).mean():.1f}% conc, mean rho {ws.mean():+.4f} | "
              f"no site n={len(ns)}: {100*(ns<0).mean():.1f}% conc, mean rho {ns.mean():+.4f} "
              f"(Fisher p={pv[1]:.3g}, OR={pv[0]:.3f}; MWU p={sv.pvalue:.3g})")

for nkey, bkey in (("n_sites", "best_ctx"), ("n_sites_default", "best_ctx_default"),
                   ("n_sites_conserved", "best_ctx_conserved")):
    analyse(net, "all_miRNA_target_edges", bkey, nkey)
    for t in ("strong", "weak", "predicted_only"):
        analyse([r for r in net if r["tier"] == t], f"tier_{t}", bkey, nkey)

with open(os.path.join(OUT, "targetscan_concordance_stats.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(res[0].keys())); w.writeheader(); w.writerows(res)
print("\nwrote targetscan_concordance_stats.csv rows:", len(res))
