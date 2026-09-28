#!/usr/bin/env python3
"""Probe GEO series (SOFT brief) for outcome/survival fields in sample characteristics."""
import json, os, re, sys, time, urllib.request, gzip, io, collections

CACHE = "/path/to/revision/cache/v6/cohorts/soft"
os.makedirs(CACHE, exist_ok=True)
URL = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=all&form=text&view=brief")

def fetch(acc):
    p = os.path.join(CACHE, acc + ".soft.txt")
    if os.path.exists(p) and os.path.getsize(p) > 200:
        return open(p, encoding="utf-8", errors="replace").read()
    for i in range(4):
        try:
            req = urllib.request.Request(URL % acc, headers={"User-Agent":"mirna-ffl/1.0"})
            with urllib.request.urlopen(req, timeout=180) as r:
                t = r.read().decode("utf-8","replace")
            open(p,"w").write(t); time.sleep(0.5); return t
        except Exception as e:
            sys.stderr.write("retry %s %s\n"%(acc,e)); time.sleep(3+3*i)
    return ""

OUTCOME_RE = re.compile(r"surviv|dfs|rfs|os_|_os\b|relapse|recur|event|follow.?up|"
    r"death|dead|alive|status|time|month|year|dmfs|drfs|mfs|pcr|response|outcome|"
    r"metasta|censor|prognos", re.I)

accs = sys.argv[1:]
res = {}
for acc in accs:
    t = fetch(acc)
    if not t:
        res[acc] = dict(ok=False); print(acc, "FETCH FAILED"); continue
    chars = collections.Counter()
    for line in t.splitlines():
        if line.startswith("!Sample_characteristics_ch1"):
            v = line.split("=",1)[1].strip()
            key = v.split(":")[0].strip() if ":" in v else v[:40]
            chars[key] += 1
    plat = sorted(set(re.findall(r"!Sample_platform_id = (GPL\d+)", t)))
    nsamp = len(re.findall(r"^\^SAMPLE", t, re.M))
    src = collections.Counter(re.findall(r"!Sample_source_name_ch1 = (.+)", t))
    title = (re.findall(r"!Series_title = (.+)", t) or [""])[0]
    hits = [k for k in chars if OUTCOME_RE.search(k)]
    res[acc] = dict(ok=True, n=nsamp, platforms=plat, title=title,
                    char_keys=dict(chars), outcome_keys=hits,
                    top_sources=src.most_common(6))
    print("%-11s n=%-5d %s" % (acc, nsamp, ",".join(plat)))
    print("    title: %s" % title[:120])
    print("    OUTCOME KEYS: %s" % (hits if hits else "NONE"))
    print("    all keys: %s" % list(chars.keys())[:25])
json.dump(res, open("/path/to/revision/cache/v6/cohorts/probe_%s.json"%(accs[0]),"w"), indent=1)
