#!/usr/bin/env python3
"""Print a compact sample table of a GEO series (from the cached brief SOFT): gsm | title | characteristics | platform."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _geo
part = sys.argv[1]
for g in sys.argv[2:]:
    S = _geo.samples(g, part)
    print(f'== {g}: {len(S)} samples')
    for s in S:
        print(f"{s['gsm']} | {s['title'][:60]} | {'; '.join(c[:50] for c in s['characteristics'])[:150]} | {s['platform']} {s['strategy']}")
