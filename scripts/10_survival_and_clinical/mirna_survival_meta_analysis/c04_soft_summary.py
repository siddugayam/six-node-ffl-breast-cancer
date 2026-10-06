#!/usr/bin/env python3
# Screens the GEO SOFT descriptions in cache/v6/cohorts/soft/ for outcome annotation (survival, relapse,
# response) and prints a summary for each series.
import os, re, sys, json, collections
CACHE="/path/to/revision/cache/v6/cohorts/soft"
OUTCOME_RE=re.compile(r"surviv|dfs|rfs|\bos\b|relapse|recur|event|follow.?up|death|dead|alive|"
  r"dmfs|drfs|\bmfs\b|\bpcr\b|response|outcome|metasta|censor|prognos|time|status", re.I)
accs=sys.argv[1:] or sorted(f[:-9] for f in os.listdir(CACHE) if f.endswith(".soft.txt"))
for acc in accs:
    p=os.path.join(CACHE,acc+".soft.txt")
    if not os.path.exists(p): print(acc,"MISSING"); continue
    t=open(p,encoding="utf-8",errors="replace").read()
    chars=collections.Counter(); ex={}
    for line in t.splitlines():
        if line.startswith("!Sample_characteristics_ch1"):
            v=line.split("=",1)[1].strip()
            k=v.split(":")[0].strip() if ":" in v else v[:40]
            chars[k]+=1; ex.setdefault(k, v[:80])
    plat=sorted(set(re.findall(r"!Sample_platform_id = (GPL\d+)",t)))
    per_plat=collections.Counter(re.findall(r"!Sample_platform_id = (GPL\d+)",t))
    n=len(re.findall(r"^\^SAMPLE",t,re.M))
    title=(re.findall(r"!Series_title = (.+)",t) or [""])[0]
    hits=[k for k in chars if OUTCOME_RE.search(k)]
    print("="*100)
    print("%s  n=%d  platforms=%s %s" % (acc,n,plat,dict(per_plat)))
    print("  title: %s" % title[:150])
    print("  OUTCOME-LIKE KEYS: %s" % hits)
    for k in list(chars.keys()):
        print("     [%3d] %-38s | e.g. %s" % (chars[k], k[:38], ex[k]))
