#!/usr/bin/env python3
import sys, os, json, numpy as np, csv
sys.path.insert(0,"/path/to/revision/scripts/10_survival_and_clinical/mirna_survival_meta_analysis")
from c07_build_mirna_cohorts import build, write_cohort, TARGETS

def num(v):
    try:
        f=float(v); return f
    except Exception: return ""

report={}
# ---------------- GSE22216: Illumina GPL8178, DRFS years + event -----------
b=build("GSE22216","GSE22216_series_matrix.txt.gz","GPL8178")
c=b["chars"]
print("GSE22216 keys:", list(c.keys()))
clin=dict(time=[num(x) for x in c["distant-relapse free survival"]],
          event=[num(x) for x in c["distant-relapse event"]],
          age=[num(x) for x in c["patient age"]],
          er=[num(x) for x in c["er status"]],
          grade=[num(x) for x in c["tumour grade"]],
          size=[num(x) for x in c["tumour size"]],
          nodes=[num(x) for x in c["nodes involved"]])
write_cohort(b,clin); report["GSE22216"]=dict(n=b["n"],probe_used=b["probe_used"],probe_all=b["probe_all"])

# ---------------- GSE59829: Illumina GPL8179, DFS months + distant met -----
b=build("GSE59829","GSE59829_series_matrix.txt.gz","GPL8179")
c=b["chars"]; print("GSE59829 keys:", list(c.keys()))
clin=dict(time=[num(x) for x in c["dfs (months)"]],
          event=[num(x) for x in c["distant_metastasis"]],
          age=[num(x) for x in c["age"]],
          er=[num(x) for x in c["esr1_status"]],
          her2=[num(x) for x in c["erbb2_status"]],
          size=[num(x) for x in c["tumor_size"]])
write_cohort(b,clin); report["GSE59829"]=dict(n=b["n"],probe_used=b["probe_used"],probe_all=b["probe_all"])

# ---------------- GSE78870: TLDA GPL20662, TTP months + progression -------
b=build("GSE78870","GSE78870_series_matrix.txt.gz","GPL20662")
c=b["chars"]; print("GSE78870 keys:", list(c.keys()))
clin=dict(time=[num(x) for x in c["time-to-progression (in months)"]],
          event=[num(x) for x in c["disease progression (event)"]],
          dfi_time=[num(x) for x in c["disease-free interval (in months)"]],
          dfi_event=[num(x) for x in c["disease relapse (event)"]],
          age=[num(x) for x in c["age at start therapy"]],
          grade=[x for x in c["grade"]], er=[x for x in c["estrogen receptor status"]])
write_cohort(b,clin); report["GSE78870"]=dict(n=b["n"],probe_used=b["probe_used"],probe_all=b["probe_all"])

# ---------------- GSE37405 (3 platforms): time.to.recurrence / death ------
for gpl in ["GPL13703","GPL14149","GPL15462"]:
    b=build("GSE37405_%s"%gpl,"GSE37405-%s_series_matrix.txt.gz"%gpl,gpl)
    c=b["chars"]; print("GSE37405 %s keys:"%gpl, list(c.keys()))
    n=b["n"]
    def g(k): return c.get(k,[""]*n)
    clin=dict(t_rec=[num(x) for x in g("time.to.recurrence (years)")],
              t_death=[num(x) for x in g("time.to.death (years)")],
              t_ok=[num(x) for x in g("time.to.ok")],
              age=[num(x) for x in g("age at op")],
              er=[num(x) for x in g("er-status")],
              size=[num(x) for x in g("size (mm)")],
              nodepos=[num(x) for x in g("nodal-status, positive")])
    write_cohort(b,clin); report["GSE37405_%s"%gpl]=dict(n=b["n"],probe_used=b["probe_used"],probe_all=b["probe_all"])

json.dump(report, open("/path/to/revision/cache/v6/cohorts/built/probe_report.json","w"), indent=1)
for k,v in report.items():
    print("==",k,"n=",v["n"])
    for t in sorted(v["probe_used"]): print("   %-9s -> %s   [all: %s]"%(t,v["probe_used"][t],v["probe_all"][t][:150]))
