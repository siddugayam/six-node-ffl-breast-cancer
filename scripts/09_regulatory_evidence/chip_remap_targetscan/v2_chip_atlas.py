#!/usr/bin/env python3
"""
ChIP-seq evidence for the TF->target edges of the breast-cancer miRNA-TF FFL network.

Source: ChIP-Atlas (https://chip-atlas.dbcls.jp) hg38 "Target Genes" tables,
        https://chip-atlas.dbcls.jp/data/hg38/target/{ANTIGEN}.{THR}.tsv
        THR in {1,5,10} = maximum distance (kb) between a MACS2 peak and the RefSeq TSS.
        Column 1 = gene symbol, column 2 = mean MACS2 score (-10*log10 q) averaged over ALL
        experiments for that antigen (experiments without a peak contribute 0), columns 3..N
        = one column per ChIP-seq experiment, header "SRXid|CellType", value = MACS2 score
        of the peak assigned to that gene in that experiment (0 = no peak).
        A gene absent from the table has no peak in any experiment at that threshold.

Outputs (results/v2/):
  chip_evidence.csv                  - tidy, one row per (TF, target, threshold, source)
  chip_edge_support_summary.csv      - fraction of the 770 TF->target edges with ChIP support
  chip_atlas_tf_inventory.csv        - per-TF experiment counts / genes bound
"""
import os, re, sys, csv, glob, json
from collections import defaultdict

REV = "/path/to/revision"
SCR = "/path/to/scratch"
CHIP = os.path.join(SCR, "chip")
OUT = os.path.join(REV, "results", "v2")
THRS = [t for t in ("1", "5", "10") if os.path.isdir(os.path.join(CHIP, t))]

# ---------- network ----------
edges = []
with open(os.path.join(REV, "data", "canonical_edges.tsv")) as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        edges.append(r)
tf_edges = [(r["source"], r["target"]) for r in edges if r["edge_type"] == "TF_target"]
tf_edge_set = set(tf_edges)
tfs = sorted({s for s, t in tf_edges})
targets = sorted({t for s, t in tf_edges})
sys.stderr.write(f"TF_target edges={len(tf_edges)} TFs={len(tfs)} targets={len(targets)}\n")

FOCUS_TF = ["NFKB1", "RELA", "SP1", "ETS1"]
FOCUS_GENES = ["COL1A1", "COL3A1"]
# genes we always extract, even if not a network target of that TF
extra_genes = set(FOCUS_GENES)
want_genes = set(targets) | extra_genes

BREAST = ("MCF-7","MCF7","MCF_7","MDA-MB","SUM_","SUM1","T-47D","T47D","ZR-75","ZR75","BT-474",
          "BT474","HMEC","breast","MCF_10A","MCF-10A","SK-BR","SKBR","CAMA","HCC70","HCC1937",
          "HCC1806","HCC1954","MCF-12A","mammary","Luminal","Basal")
def is_breast(ct):
    c = ct.lower()
    return any(k.lower() in c for k in BREAST)

rows = []            # chip_evidence.csv rows
inventory = []       # per TF per threshold
bound = defaultdict(dict)   # (thr) -> (tf,gene) -> dict

for thr in THRS:
    for tf in tfs:
        path = os.path.join(CHIP, thr, f"{tf}.{thr}.tsv")
        if not os.path.exists(path):
            inventory.append(dict(tf=tf, threshold_kb=thr, data_available="FALSE",
                                  n_experiments="NA", n_genes_in_table="NA", n_celltypes="NA"))
            continue
        with open(path) as fh:
            hdr = [h.strip() for h in fh.readline().rstrip("\n").split("\t")]
            # columns 3..N are one per ChIP-seq experiment EXCEPT the trailing "STRING"
            # column, which ChIP-Atlas appends (a STRING protein-interaction score, not a
            # ChIP experiment).  Experiment columns are "SRX/ERX/DRX<digits>|CellType".
            expcols = [i for i, h in enumerate(hdr)
                       if i >= 2 and re.match(r"^[SED]RX\d+\|", h)]
            assert all(h == "STRING" or re.match(r"^[SED]RX\d+\|", h)
                       for i, h in enumerate(hdr) if i >= 2), f"unexpected column in {path}"
            exps = [hdr[i] for i in expcols]
            cts = [e.split("|", 1)[1] if "|" in e else "NA" for e in exps]
            ngen = 0
            allavg = []
            for line in fh:
                p = line.rstrip("\n").split("\t")
                ngen += 1
                g = p[0]
                try:
                    allavg.append(float(p[1]))
                except ValueError:
                    allavg.append(0.0)
                if g not in want_genes:
                    continue
                try:
                    avg = float(p[1])
                except ValueError:
                    avg = float("nan")
                vals = []
                for i in expcols:
                    try:
                        vals.append(float(p[i]))
                    except (ValueError, IndexError):
                        vals.append(0.0)
                # consistency check: ChIP-Atlas "Average" must be the mean over experiments
                if vals and avg == avg:
                    assert abs(sum(vals) / len(vals) - avg) < 0.01, \
                        f"average mismatch {tf} {g}: {sum(vals)/len(vals)} vs {avg}"
                nz = [(exps[i], cts[i], vals[i]) for i in range(len(vals)) if vals[i] > 0]
                nz.sort(key=lambda z: -z[2])
                ct_nz = sorted({z[1] for z in nz})
                br = [z for z in nz if is_breast(z[1])]
                bound[thr][(tf, g)] = dict(
                    avg=avg, n_exp=len(exps), n_peak=len(nz),
                    n_ct=len(ct_nz), maxs=(nz[0][2] if nz else 0.0),
                    top=";".join(f"{z[1]}:{int(z[2])}" for z in nz[:6]),
                    cts=";".join(ct_nz[:12]),
                    n_breast_peak=len(br),
                    breast_cts=";".join(sorted({f"{z[1]}:{int(z[2])}" for z in br})[:8]))
            allavg.sort(reverse=True)
            import bisect
            desc = allavg
            for (tf2, g2), d in list(bound[thr].items()):
                if tf2 != tf:
                    continue
                # rank of this gene's average score among ALL genes in this TF's table
                lo = 0; hi = len(desc)
                while lo < hi:
                    mid = (lo + hi) // 2
                    if desc[mid] > d["avg"]: lo = mid + 1
                    else: hi = mid
                d["rank"] = lo + 1
                d["pctile"] = round(100.0 * (1.0 - lo / max(len(desc), 1)), 2)
                d["n_genes_table"] = len(desc)
            inventory.append(dict(tf=tf, threshold_kb=thr, data_available="TRUE",
                                  n_experiments=len(exps), n_genes_in_table=ngen,
                                  n_celltypes=len(set(cts)),
                                  n_breast_experiments=sum(1 for c in cts if is_breast(c))))
        sys.stderr.write(".")
    sys.stderr.write(f" thr={thr} done\n")

inv_idx = {(d["tf"], d["threshold_kb"]): d for d in inventory}

def emit(tf, gene, thr, analysis):
    d = bound[thr].get((tf, gene))
    inv = inv_idx.get((tf, thr), {})
    avail = inv.get("data_available", "FALSE")
    if avail != "TRUE":
        rows.append(dict(analysis=analysis, source_TF=tf, target_gene=gene,
                         edge_in_network=str((tf, gene) in tf_edge_set).upper(),
                         resource="ChIP-Atlas_hg38", threshold_kb=thr,
                         chip_data_available="FALSE", n_experiments="NA",
                         n_experiments_with_peak="NA", frac_experiments_with_peak="NA",
                         n_celltypes_with_peak="NA", avg_MACS2_score="NA",
                         max_MACS2_score="NA", bound_any="NA", bound_reproducible="NA",
                         bound_strong="NA", top_celltypes="", celltypes_with_peak="",
                         n_breast_experiments_with_peak="NA", breast_celltypes_with_peak="",
                         rank_among_TF_targets="NA", pctile_among_TF_targets="NA",
                         n_genes_bound_by_TF="NA"))
        return
    nexp = inv["n_experiments"]
    if d is None:
        rows.append(dict(analysis=analysis, source_TF=tf, target_gene=gene,
                         edge_in_network=str((tf, gene) in tf_edge_set).upper(),
                         resource="ChIP-Atlas_hg38", threshold_kb=thr,
                         chip_data_available="TRUE", n_experiments=nexp,
                         n_experiments_with_peak=0, frac_experiments_with_peak=0.0,
                         n_celltypes_with_peak=0, avg_MACS2_score=0.0, max_MACS2_score=0.0,
                         bound_any="FALSE", bound_reproducible="FALSE", bound_strong="FALSE",
                         top_celltypes="", celltypes_with_peak="",
                         n_breast_experiments_with_peak=0, breast_celltypes_with_peak="",
                         rank_among_TF_targets="", pctile_among_TF_targets="",
                         n_genes_bound_by_TF=inv["n_genes_in_table"]))
        return
    rows.append(dict(analysis=analysis, source_TF=tf, target_gene=gene,
                     edge_in_network=str((tf, gene) in tf_edge_set).upper(),
                     resource="ChIP-Atlas_hg38", threshold_kb=thr,
                     chip_data_available="TRUE", n_experiments=nexp,
                     n_experiments_with_peak=d["n_peak"],
                     frac_experiments_with_peak=round(d["n_peak"] / nexp, 4),
                     n_celltypes_with_peak=d["n_ct"],
                     avg_MACS2_score=round(d["avg"], 4), max_MACS2_score=d["maxs"],
                     bound_any=str(d["n_peak"] >= 1).upper(),
                     bound_reproducible=str(d["n_peak"] >= 2 and d["n_ct"] >= 2).upper(),
                     bound_strong=str(d["avg"] >= 50).upper(),
                     top_celltypes=d["top"], celltypes_with_peak=d["cts"],
                     n_breast_experiments_with_peak=d["n_breast_peak"],
                     breast_celltypes_with_peak=d["breast_cts"],
                     rank_among_TF_targets=d.get("rank", ""),
                     pctile_among_TF_targets=d.get("pctile", ""),
                     n_genes_bound_by_TF=inv["n_genes_in_table"]))

for thr in THRS:
    for tf, g in tf_edges:
        emit(tf, g, thr, "network_TF_target_edge")
    for tf in tfs:
        for g in FOCUS_GENES:
            if (tf, g) not in tf_edge_set:
                emit(tf, g, thr, "collagen_reverse_query")

# ---------- summary over the 770 edges ----------
summary = []
for thr in THRS:
    ev = [r for r in rows if r["analysis"] == "network_TF_target_edge" and r["threshold_kb"] == thr]
    testable = [r for r in ev if r["chip_data_available"] == "TRUE"]
    for crit in ("bound_any", "bound_reproducible", "bound_strong"):
        n = sum(1 for r in testable if r[crit] == "TRUE")
        summary.append(dict(threshold_kb=thr, criterion=crit, n_edges_total=len(ev),
                            n_edges_testable=len(testable), n_supported=n,
                            frac_of_testable=round(n / len(testable), 4) if testable else "NA",
                            frac_of_all_770=round(n / len(ev), 4) if ev else "NA"))

os.makedirs(OUT, exist_ok=True)
flds = ["analysis", "source_TF", "target_gene", "edge_in_network", "resource", "threshold_kb",
        "chip_data_available", "n_experiments", "n_experiments_with_peak",
        "frac_experiments_with_peak", "n_celltypes_with_peak", "avg_MACS2_score",
        "max_MACS2_score", "bound_any", "bound_reproducible", "bound_strong",
        "n_breast_experiments_with_peak", "rank_among_TF_targets", "pctile_among_TF_targets",
        "n_genes_bound_by_TF", "top_celltypes", "celltypes_with_peak",
        "breast_celltypes_with_peak"]
with open(os.path.join(OUT, "chip_evidence.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=flds); w.writeheader(); w.writerows(rows)
with open(os.path.join(OUT, "chip_edge_support_summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
with open(os.path.join(OUT, "chip_atlas_tf_inventory.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(inventory[0].keys())); w.writeheader(); w.writerows(inventory)

print("thresholds used:", THRS)
print("rows in chip_evidence.csv:", len(rows))
for s in summary:
    print(s)
print("\nFOCUS:")
for thr in THRS:
    for tf in FOCUS_TF:
        for g in FOCUS_GENES:
            d = bound[thr].get((tf, g))
            inv = inv_idx.get((tf, thr), {})
            if inv.get("data_available") != "TRUE":
                print(f"  {thr}kb {tf}->{g}: NO ChIP-Atlas DATA"); continue
            if d is None:
                print(f"  {thr}kb {tf}->{g}: NOT BOUND (0/{inv['n_experiments']} experiments)")
            else:
                print(f"  {thr}kb {tf}->{g}: avg={d['avg']:.2f} peak in {d['n_peak']}/{inv['n_experiments']} exp "
                      f"({d['n_breast_peak']} breast), {d['n_ct']} cell types, rank {d.get('rank')}/{d.get('n_genes_table')} "
                      f"of this TF's bound genes (top {100-d.get('pctile',0):.1f}%) | top: {d['top']}")
