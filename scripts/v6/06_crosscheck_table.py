#!/usr/bin/env python3
"""Part D: HPA Pathology-Atlas breast prognostic calls vs our own TCGA-BRCA and
METABRIC Cox results.

HPA calls a gene "prognostic" only at p < 0.001 after a best-cut-off scan (their
own stated threshold); we therefore compare on two axes:
  (i)  DIRECTION  - favourable vs unfavourable, taken for EVERY gene from the HPA
       pathology table (v23 columns carry the direction even when the call is
       'unprognostic'); compared with the sign of our Cox coefficient.
  (ii) CALL       - HPA prognostic (p<0.001) vs ours significant at p<0.05 and at
       the matched p<0.001.
"""
import csv, os, collections
R="/path/to/revision/results/v6"
D="/path/to/revision/data/hpa"

hp={r["gene"]:r for r in csv.DictReader(open(os.path.join(R,"hpa_prognostic_breast.csv")))}
order=[r["gene"] for r in csv.DictReader(open(os.path.join(R,"hpa_prognostic_breast.csv")))]
roles={r["gene"]:r["role"] for r in csv.DictReader(open(os.path.join(R,"hpa_prognostic_breast.csv")))}

# direction from v23 pathology.tsv (present for prognostic AND unprognostic genes)
v23={}
with open(os.path.join(D,"pathology.tsv"),newline="",encoding="utf-8") as f:
    for r in csv.DictReader(f,delimiter="\t"):
        if r["Cancer"]=="breast cancer" and r["Gene name"] in hp:
            for col,dirn,call in (("prognostic - favorable","favourable","prognostic"),
                                  ("unprognostic - favorable","favourable","unprognostic"),
                                  ("prognostic - unfavorable","unfavourable","prognostic"),
                                  ("unprognostic - unfavorable","unfavourable","unprognostic")):
                if r[col].strip():
                    v23[r["Gene name"]]=dict(direction=dirn,call=call,p=r[col].strip())
cox=list(csv.DictReader(open(os.path.join(R,"hpa_crosscheck_our_cox.csv"))))
by=collections.defaultdict(dict)
for r in cox: by[r["gene"]][(r["cohort"],r["endpoint"])]=r

rows=[]
for g in order:
    h=hp[g]; v=v23.get(g,{}); d=by[g]
    hpa_dir=v.get("direction","")
    hpa_p25=float(h["hpa_TCGA_p"]) if h["hpa_TCGA_p"] else None
    hpa_sig=(h["hpa_TCGA_call"].startswith("potential prognostic"))
    row=dict(gene=g, role=roles[g],
        hpa_v25_call=h["hpa_TCGA_call"], hpa_v25_p=h["hpa_TCGA_p"],
        hpa_v25_validation_call=h["hpa_validation_call"], hpa_v25_validation_p=h["hpa_validation_p"],
        hpa_v23_call=v.get("call",""), hpa_v23_direction=hpa_dir, hpa_v23_p=v.get("p",""))
    for label,key in (("TCGA_OS",("TCGA-BRCA","OS")),("TCGA_DSS",("TCGA-BRCA","DSS")),
                      ("METABRIC_OS",("METABRIC","OS")),("METABRIC_RFS",("METABRIC","RFS"))):
        r=d.get(key)
        row[f"our_{label}_HR_per_SD"]=round(float(r["HR_per_SD"]),3) if r else ""
        row[f"our_{label}_p"]=f'{float(r["p_cox"]):.3g}' if r else ""
        row[f"our_{label}_direction"]=r["cox_direction"] if r else ""
        row[f"dir_agree_{label}"]=("agree" if r and hpa_dir and r["cox_direction"]==hpa_dir
                                   else ("DISAGREE" if r and hpa_dir else ""))
    r=d.get(("TCGA-BRCA","OS"))
    ours_sig05 = bool(r) and float(r["p_cox"])<0.05
    ours_sig001= bool(r) and float(r["p_cox"])<0.001
    row["call_agree_TCGA_OS_alpha0.05"]=("agree: both prognostic" if hpa_sig and ours_sig05 else
        "agree: both non-prognostic" if not hpa_sig and not ours_sig05 else
        "HPA prognostic, ours not" if hpa_sig else "ours prognostic, HPA not")
    row["call_agree_TCGA_OS_matched_alpha0.001"]=("agree: both prognostic" if hpa_sig and ours_sig001 else
        "agree: both non-prognostic" if not hpa_sig and not ours_sig001 else
        "HPA prognostic, ours not" if hpa_sig else "ours prognostic, HPA not")
    rows.append(row)

with open(os.path.join(R,"hpa_vs_our_survival_crosscheck.csv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("hpa_vs_our_survival_crosscheck.csv",len(rows))
t=collections.Counter()
for r in rows:
    for k in ("dir_agree_TCGA_OS","dir_agree_TCGA_DSS","dir_agree_METABRIC_OS","dir_agree_METABRIC_RFS",
              "call_agree_TCGA_OS_alpha0.05","call_agree_TCGA_OS_matched_alpha0.001"):
        t[(k,r[k])]+=1
for k in sorted(t): print(k,t[k])
print()
print(f"{'gene':8s} {'HPAdir':12s} {'HPAv25 p':10s} {'TCGA-OS':14s} {'p':9s} {'MB-OS':14s} {'p':9s} {'dirTCGA':9s} {'dirMB':9s}")
for r in rows:
    print(f"{r['gene']:8s} {r['hpa_v23_direction']:12s} {r['hpa_v25_p']:10s} {r['our_TCGA_OS_direction']:14s} {r['our_TCGA_OS_p']:9s} {r['our_METABRIC_OS_direction']:14s} {r['our_METABRIC_OS_p']:9s} {r['dir_agree_TCGA_OS']:9s} {r['dir_agree_METABRIC_OS']:9s}")
