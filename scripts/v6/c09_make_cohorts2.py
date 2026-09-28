#!/usr/bin/env python3
import sys, os, json, collections, numpy as np
sys.path.insert(0,"/path/to/revision/scripts/v6")
from c07_build_mirna_cohorts import build, write_cohort

def num(v):
    try: return float(v)
    except Exception: return ""

report=json.load(open("cache/v6/cohorts/built/probe_report.json"))

def show(b):
    c=b["chars"]
    print("== %s n=%d gpl=%s"%(b["cohort"],b["n"],b["gpl"]))
    for k in c:
        vals=[x for x in c[k] if x]
        print("   %-45s uniq=%d e.g. %s"%(k[:45],len(set(vals)),sorted(set(vals))[:5]))
    print("   targets:",sorted(b["expr"].keys()))

# GSE19783 (existing miRNA cohort, Oslo MicMa) GPL8227
b=build("GSE19783","GSE19783-GPL8227_series_matrix.txt.gz","GPL8227"); show(b)
c=b["chars"]
clin=dict(dfs_time=[num(x) for x in c["disease free survival time (months)"]],
          death_status=[x if x else "" for x in c["death status"]],
          er=[x if x else "" for x in c["estrogen receptor status"]],
          subtype=[x if x else "" for x in c.get("breast cancer subtype",[""]*b["n"])])
write_cohort(b,clin); report["GSE19783"]=dict(n=b["n"],probe_used=b["probe_used"],probe_all=b["probe_all"])

for cohort,f,gpl in [("GSE26666","GSE26666-GPL8227_series_matrix.txt.gz","GPL8227"),
                     ("GSE57897","GSE57897_series_matrix.txt.gz","GPL18722"),
                     ("GSE103161","GSE103161_series_matrix.txt.gz","GPL23960"),
                     ("GSE126125","GSE126125_series_matrix.txt.gz","GPL23960"),
                     ("GSE97811","GSE97811_series_matrix.txt.gz","GPL21263"),
                     ("GSE28321","GSE28321_series_matrix.txt.gz","GPL5106"),
                     ("GSE40267","GSE40267_series_matrix.txt.gz","GPL10850")]:
    try:
        b=build(cohort,f,gpl); show(b)
        clin={k:[x if x is not None else "" for x in v] for k,v in b["chars"].items() if k!="_free"}
        clin={("c_"+k.replace(" ","_").replace(",","")[:40]):v for k,v in clin.items()}
        write_cohort(b,clin)
        report[cohort]=dict(n=b["n"],probe_used=b["probe_used"],probe_all=b["probe_all"])
    except Exception as e:
        print("FAILED",cohort,e)
json.dump(report, open("cache/v6/cohorts/built/probe_report.json","w"), indent=1)
