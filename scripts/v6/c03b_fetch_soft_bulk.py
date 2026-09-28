#!/usr/bin/env python3
import os, sys, time, urllib.request, json
CACHE="/path/to/revision/cache/v6/cohorts/soft"
os.makedirs(CACHE, exist_ok=True)
URL="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=all&form=text&view=brief"
accs=[l.strip() for l in open(sys.argv[1]) if l.strip()]
for acc in accs:
    p=os.path.join(CACHE, acc+".soft.txt")
    if os.path.exists(p) and os.path.getsize(p)>200: continue
    ok=False
    for i in range(4):
        try:
            req=urllib.request.Request(URL%acc, headers={"User-Agent":"mirna-ffl/1.0"})
            with urllib.request.urlopen(req, timeout=300) as r: t=r.read().decode("utf-8","replace")
            open(p,"w").write(t); ok=True; break
        except Exception as e:
            sys.stderr.write("retry %s %s\n"%(acc,e)); time.sleep(3+3*i)
    print(acc, "OK" if ok else "FAIL", flush=True)
    time.sleep(0.4)
