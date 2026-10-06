#!/usr/bin/env python3
"""Consolidation / reporting for Parts A-C (fresh run 2026-09-09)."""
import csv, os, math, collections
import numpy as np
from scipy import stats
OUT="/path/to/revision/results/v2"
def rd(f):
    with open(os.path.join(OUT,f)) as fh: return list(csv.DictReader(fh))
def fl(v):
    try:
        x=float(v); return x if math.isfinite(x) else float('nan')
    except: return float('nan')

print("="*100); print("PART B  PAN-CANCER"); print("="*100)
pa = rd("pancancer_axis.csv"); cs = rd("pancancer_cohort_summary_v3.csv")
cohorts = [c["cohort"] for c in cs]
print(f"cohorts analysed: {len(cohorts)}; total axis rows: {len(pa)}")
idx = {(r["cohort"], r["axis"]): r for r in pa}

def table(axes, cols=("rho","p","rho_partial_CAFA")):
    hdr = f"{'cohort':38s}{'n':>6s}" + "".join(f"{a.replace('hsa-','')[:22]:>24s}" for a in axes)
    print(hdr); print("-"*len(hdr))
    for c in cohorts:
        line=f"{c[:37]:38s}"
        n=""
        cells=[]
        for a in axes:
            r=idx.get((c,a))
            if r is None: cells.append(f"{'-':>24s}"); continue
            n=r["n"]
            cells.append(f"{fl(r['rho']):+.3f}({fl(r['p']):.0e})".rjust(24))
        print(line+f"{n:>6s}"+"".join(cells))

print("\n--- miR-29 (arm-averaged) -> COL1A1 ---")
table(["hsa-miR-29a -> COL1A1","hsa-miR-29b -> COL1A1","hsa-miR-29c -> COL1A1"])
print("\n--- miR-29 (arm-averaged) -> COL3A1 ---")
table(["hsa-miR-29a -> COL3A1","hsa-miR-29b -> COL3A1","hsa-miR-29c -> COL3A1"])
print("\n--- miR-29-3p (dominant arm) -> COL1A1 / COL3A1 ---")
table(["hsa-miR-29a-3p -> COL1A1","hsa-miR-29b-3p -> COL3A1","hsa-miR-29c-3p -> COL3A1"])
print("\n--- ETS1/NFKB1 -> COL1A1, and COL1A1~COL3A1 ---")
table(["ETS1 -> COL1A1","NFKB1 -> COL1A1","COL1A1 -> COL3A1"])
print("\n--- CAF-adjusted partial rho (same axes) ---")
for a in ["hsa-miR-29a -> COL3A1","hsa-miR-29c -> COL3A1","ETS1 -> COL1A1","NFKB1 -> COL1A1"]:
    vals=[(c, fl(idx[(c,a)]["rho"]), fl(idx[(c,a)]["rho_partial_CAFA"])) for c in cohorts if (c,a) in idx]
    neg_raw=sum(1 for _,r,_ in vals if r<0); neg_p=sum(1 for _,_,p in vals if p<0)
    print(f"{a:26s} raw mean {np.nanmean([v[1] for v in vals]):+.3f}  partial mean {np.nanmean([v[2] for v in vals]):+.3f}"
          f"   raw negative in {neg_raw}/{len(vals)}   partial negative in {neg_p}/{len(vals)}")

print("\n--- consistency of sign across cohorts (raw rho) ---")
for a in sorted({r["axis"] for r in pa}):
    vals=[fl(idx[(c,a)]["rho"]) for c in cohorts if (c,a) in idx and math.isfinite(fl(idx[(c,a)]["rho"]))]
    ps  =[fl(idx[(c,a)]["p"]) for c in cohorts if (c,a) in idx]
    if len(vals)<8: continue
    neg=sum(1 for v in vals if v<0); sig_neg=sum(1 for v,p in zip(vals,ps) if v<0 and p<0.05)
    sig_pos=sum(1 for v,p in zip(vals,ps) if v>0 and p<0.05)
    bt=stats.binomtest(neg,len(vals),0.5)
    print(f"{a:28s} k={len(vals):2d} mean_rho={np.mean(vals):+.3f} median={np.median(vals):+.3f} "
          f"neg={neg:2d} sig_neg={sig_neg:2d} sig_pos={sig_pos:2d} binom_p={bt.pvalue:.2e}")

print("\n--- cohort stromal / collagen summary ---")
print(f"{'cohort':38s}{'n':>5s}{'stromIdx_med':>13s}{'stromIdx_sd':>12s}{'epiIdx_med':>11s}"
      f"{'COL1A1_med':>11s}{'COL1A1_sd':>10s}{'rCAF_COL1A1':>12s}{'rCOL1_COL3':>11s}")
for c in cs:
    print(f"{c['cohort'][:37]:38s}{c['n']:>5s}{fl(c['strom_abs_median']):>13.3f}{fl(c['strom_abs_sd']):>12.3f}"
          f"{fl(c['epi_abs_median']):>11.3f}{fl(c['COL1A1_median']):>11.2f}{fl(c['COL1A1_sd']):>10.2f}"
          f"{fl(c['rho_CAFA_COL1A1']):>12.3f}{fl(c['rho_COL1A1_COL3A1']):>11.3f}")

print("\n--- ECOLOGICAL TEST: does axis strength track cohort stromal content? ---")
eco = rd("pancancer_ecological_v3.csv")
keep_ax=["ETS1 -> COL1A1","NFKB1 -> COL1A1","SP1 -> COL1A1","RELA -> COL1A1","MYC -> COL1A1",
         "ETS1 -> COL3A1","NFKB1 -> COL3A1",
         "hsa-miR-29a -> COL1A1","hsa-miR-29b -> COL1A1","hsa-miR-29c -> COL1A1",
         "hsa-miR-29a -> COL3A1","hsa-miR-29b -> COL3A1","hsa-miR-29c -> COL3A1",
         "COL1A1 -> COL3A1","CAF_A -> COL1A1"]
print(f"{'axis':28s}{'cohort stat':20s}{'k':>4s}{'spearman':>10s}{'p':>10s}{'pearson_z':>11s}{'p':>10s}")
for e in eco:
    if e["axis"] in keep_ax and e["cohort_statistic"] in ("strom_abs_median","strom_abs_sd","rho_CAFA_COL1A1","HETEROGENEITY_Q"):
        print(f"{e['axis']:28s}{e['cohort_statistic']:20s}{e['k_cohorts']:>4s}"
              f"{fl(e['spearman_rho']):>10.3f}{fl(e['spearman_p']):>10.2e}"
              f"{fl(e['pearson_on_fisherz']):>11.3f}{fl(e['pearson_p']):>10.2e}")

print("\n"+"="*100); print("PART C  CPTAC PAN-CANCER PROTEOME"); print("="*100)
cp = rd("cptac_pancancer_mir29.csv")
inv = rd("cptac_pancancer_inventory_v3.csv")
print(f"{'cohort':8s}{'nMir':>6s}{'nProt':>7s}{'nRNA':>7s}{'ov(mir&prot)':>14s}{'ov(all3)':>10s}")
for i in inv:
    print(f"{i['cohort']:8s}{i['n_mirna_samples']:>6s}{i['n_protein_samples']:>7s}{i['n_rna_samples']:>7s}"
          f"{i['n_overlap_mirna_protein']:>14s}{i['n_overlap_all3']:>10s}")
cohC=[i["cohort"] for i in inv]
for tgt in ["COL1A1","COL3A1"]:
    print(f"\n--- miR-29a/b/c (arm-averaged) -> {tgt}: rho_protein (rho_mRNA) ---")
    print(f"{'cohort':8s}{'n':>5s}" + "".join(f"{m:>26s}" for m in ["miR-29a","miR-29b","miR-29c"]))
    for c in cohC:
        cells=[]; n=""
        for m in ["hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"]:
            r=[x for x in cp if x["cohort"]==c and x["miRNA"]==m and x["target"]==tgt]
            if not r: cells.append(f"{'-':>26s}"); continue
            r=r[0]; n=r["n_protein"]
            cells.append(f"{fl(r['rho_protein']):+.3f}[p{fl(r['p_protein']):.0e}]({fl(r['rho_mrna']):+.3f})".rjust(26))
        print(f"{c:8s}{n:>5s}"+"".join(cells))
print("\n--- protein-level sign consistency across cohorts (miR-29 x COL1A1/COL3A1, arm-avg) ---")
for m in ["hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"]:
    for tgt in ["COL1A1","COL3A1"]:
        rs=[fl(x["rho_protein"]) for x in cp if x["miRNA"]==m and x["target"]==tgt and math.isfinite(fl(x["rho_protein"]))]
        ps=[fl(x["p_protein"]) for x in cp if x["miRNA"]==m and x["target"]==tgt]
        neg=sum(1 for v in rs if v<0); sn=sum(1 for v,p in zip(rs,ps) if v<0 and p<0.05)
        print(f"{m} -> {tgt:7s} k={len(rs)} mean={np.mean(rs):+.3f} negative={neg}/{len(rs)} sig_neg={sn} "
              f"binom_p={stats.binomtest(neg,len(rs),0.5).pvalue:.3e}")
print("\n--- protein vs mRNA, per cohort ---")
pv = rd("cptac_pancancer_protein_vs_mrna_v3.csv")
print(f"{'cohort':8s}{'pairs':>6s}{'mean_rho_prot':>14s}{'mean_rho_mrna':>14s}{'prot_stronger':>14s}{'p':>9s}{'neg_prot':>10s}{'p':>10s}")
for r in pv:
    print(f"{r['cohort']:8s}{r['n_pairs']:>6s}{fl(r['mean_rho_protein']):>14.3f}{fl(r['mean_rho_mrna']):>14.3f}"
          f"{r['n_protein_stronger']:>14s}{fl(r['binom_p_stronger']):>9.3f}{r['n_negative_protein']:>10s}{fl(r['binom_p_negative']):>10.2e}")
print("\nDONE")
