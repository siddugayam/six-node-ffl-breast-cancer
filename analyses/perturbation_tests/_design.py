#!/usr/bin/env python3
"""Group a series' samples by their title with replicate numbers removed, to read the design quickly.
usage: _design.py <part> <regex|-> <GSE> [<GSE> ...]  (groups not matching the regex are only counted)"""
import sys, os, re, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _geo
def norm(t):
    t = re.sub(r'(?i)(rep(licate)?|_r|-r|#|\bbio|\bsample)[ _-]?\d+\b', 'REP', t)
    t = re.sub(r'(?<![A-Za-z])\d{1,2}$', 'N', t.strip())
    return re.sub(r'\s+', ' ', t)
RX = None if sys.argv[2] == '-' else re.compile(sys.argv[2], re.I)
for g in sys.argv[3:]:
    S = _geo.samples(g, sys.argv[1])
    c = collections.OrderedDict()
    for s in S:
        k = norm(s['title']); c.setdefault(k, []).append(s)
    plat = collections.Counter((s['platform'], s['strategy'], s['organism']) for s in S)
    print(f'== {g}: {len(S)} samples; {dict(plat)}')
    hidden = 0
    for k, v in c.items():
        ch = '; '.join(x[:40] for x in v[0]['characteristics'])[:110]
        if RX and not RX.search(k + ' ' + ch): hidden += 1; continue
        print(f'  {len(v):3d} x {k[:70]} | {ch}')
    if hidden: print(f'  (+{hidden} other groups)')
