#!/usr/bin/env python3
"""Auto-classify every probed breast miRNA series by the outcome fields actually
deposited in GEO: time+event (Cox-able), event-only (logistic), or none."""
import json,os,re,csv,collections
BASE="/path/to/revision"; C=BASE+"/cache/v6/cohorts"; S=C+"/soft"
TIME=re.compile(r"(time|months?|years?|days?|duration|interval|survival|follow.?up|dfs|rfs|os|dmfs|drfs|ttp)",re.I)
EVT =re.compile(r"(event|status|relapse|recurrenc|metasta|death|dead|alive|censor|progression|survival|pcr|response)",re.I)
NUM =re.compile(r"^-?\d+(\.\d+)?$")
out={}
for f in sorted(os.listdir(S)):
    if not f.endswith(".soft.txt"): continue
    acc=f[:-9]; t=open(os.path.join(S,f),encoding="utf-8",errors="replace").read()
    vals=collections.defaultdict(list)
    for line in t.splitlines():
        if line.startswith("!Sample_characteristics_ch1"):
            v=line.split("=",1)[1].strip()
            if ":" in v:
                k,vv=v.split(":",1); vals[k.strip().lower()].append(vv.strip())
    timek=[]; evtk=[]
    for k,vv in vals.items():
        numeric=[x for x in vv if NUM.match(x)]
        frac_num=len(numeric)/max(1,len(vv))
        uniq=set(vv)
        if TIME.search(k) and frac_num>0.5 and len(set(numeric))>10: timek.append(k)
        if EVT.search(k) and len(uniq)<=4 and len(uniq)>=2: evtk.append(k)
    plats=sorted(set(re.findall(r"!Sample_platform_id = (GPL\d+)",t)))
    n=len(re.findall(r"^\^SAMPLE",t,re.M))
    title=(re.findall(r"!Series_title = (.+)",t) or [""])[0]
    cls = "time+event (Cox possible)" if (timek and evtk) else \
          ("event only (logistic possible)" if evtk else "no outcome field deposited")
    out[acc]=dict(accession=acc,n_samples=n,platforms=";".join(plats),title=title[:200],
                  time_fields=";".join(timek), event_fields=";".join(evtk),
                  auto_class=cls, all_char_keys=";".join(sorted(vals.keys()))[:600])
json.dump(out, open(C+"/autoclass.json","w"), indent=1)
c=collections.Counter(v["auto_class"] for v in out.values())
print("probed series:",len(out))
for k,v in c.most_common(): print("%4d  %s"%(v,k))
print("\n--- time+event series ---")
for a,v in sorted(out.items(), key=lambda kv:-kv[1]["n_samples"]):
    if v["auto_class"].startswith("time+event"):
        print("%-11s n=%-5d %s | time=[%s] event=[%s]"%(a,v["n_samples"],v["title"][:60],v["time_fields"][:60],v["event_fields"][:60]))
