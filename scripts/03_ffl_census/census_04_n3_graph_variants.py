#!/usr/bin/env python3
"""Exhaustive n=3 FFL census under the formal definition (D1-D4), plus the two
alternative 'core' counting conventions, for several graph-augmentation variants."""
import sys, os, csv
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_03_build_graph import build, summarise, NODE_TYPE

def reduce_graph(E):
    """D1 contraction: replace every reciprocal TF<->miRNA pair by its TF->miRNA arc,
    i.e. drop the miRNA->TF arc of the pair.  Returns (E2, set_of_composite_arcs)."""
    comp = set()
    drop = set()
    for (s, t) in E:
        if (t, s) in E and {NODE_TYPE[s], NODE_TYPE[t]} == {"TF", "miRNA"}:
            if NODE_TYPE[s] == "TF":
                comp.add((s, t))          # keep TF->miRNA
            else:
                drop.add((s, t))          # drop miRNA->TF
    E2 = {k: v for k, v in E.items() if k not in drop}
    return E2, comp

def bitsets(E, nodes):
    idx = {n: i for i, n in enumerate(nodes)}
    out = [0] * len(nodes); inn = [0] * len(nodes)
    for (s, t) in E:
        out[idx[s]] |= 1 << idx[t]
        inn[idx[t]] |= 1 << idx[s]
    return idx, out, inn

def bits(x):
    while x:
        b = x & -x
        yield b.bit_length() - 1
        x ^= b

def census3(E):
    """Return (induced_D1D4_count, class_counter, ordered_core_count_G, ordered_core_count_reduced)"""
    E2, comp = reduce_graph(E)
    nodes = sorted(NODE_TYPE)
    idx, out, inn = bitsets(E2, nodes)
    full = (1 << len(nodes)) - 1
    cls = Counter()
    total = 0
    for (s, m) in E2:
        i, j = idx[s], idx[m]
        cand = out[i] & out[j] & ~inn[i] & ~inn[j] & full
        # v->u absent is guaranteed in E2 unless a non-TF/miRNA reciprocal survived
        if (m, s) in E2:
            continue
        k = bin(cand).count("1")
        if k == 0:
            continue
        total += k
        ts, tm = NODE_TYPE[s], NODE_TYPE[m]
        if (s, m) in comp:
            c = "Composite-FFL"
        elif ts == "TF" and tm == "miRNA":
            c = "TF-FFL"
        elif ts == "miRNA" and tm == "TF":
            c = "miRNA-FFL"
        else:
            c = "other(%s->%s)" % (ts, tm)
        cls[c] += k
    # ordered cores, no induced restriction
    idxF, outF, innF = bitsets(E, nodes)
    ord_G = 0
    for (s, m) in E:
        ord_G += bin(outF[idx[s]] & outF[idx[m]]).count("1")
    ord_R = 0
    for (s, m) in E2:
        ord_R += bin(out[idx[s]] & out[idx[m]]).count("1")
    return total, cls, ord_G, ord_R

def brute3(E):
    """Independent brute force over all connected triples (validation)."""
    nodes = sorted(NODE_TYPE)
    E2, comp = reduce_graph(E)
    adj = {}
    for (s, t) in E2:
        adj.setdefault(s, set()).add(t)
    und = {}
    for (s, t) in E2:
        und.setdefault(s, set()).add(t); und.setdefault(t, set()).add(s)
    seen = set(); cnt = 0
    for a in nodes:
        for b in und.get(a, ()):
            for c in und.get(a, ()) | und.get(b, ()):
                if c in (a, b):
                    continue
                key = tuple(sorted((a, b, c)))
                if key in seen:
                    continue
                seen.add(key)
                V = list(key)
                arcs = [(x, y) for x in V for y in V if x != y and y in adj.get(x, ())]
                if len(arcs) != 3:
                    continue
                outd = {v: sum(1 for x, y in arcs if x == v) for v in V}
                ind = {v: sum(1 for x, y in arcs if y == v) for v in V}
                src = [v for v in V if ind[v] == 0]; snk = [v for v in V if outd[v] == 0]
                if len(src) != 1 or len(snk) != 1:
                    continue
                s0, t0 = src[0], snk[0]
                mid = [v for v in V if v not in (s0, t0)]
                if len(mid) != 1:
                    continue
                m0 = mid[0]
                if (s0, m0) in arcs and (m0, t0) in arcs and (s0, t0) in arcs:
                    cnt += 1
    return cnt

VARIANTS = {
 "V0_deposited_only":      dict(),
 "V1_TT+GGdir":            dict(tf_target=1, gg_directed=1),
 "V2_TT+GGdir+mm_asis":    dict(tf_target=1, gg_directed=1, mirna_mirna="asis"),
 "V3_TT+GGdir+mm_both":    dict(tf_target=1, gg_directed=1, mirna_mirna="both"),
 "V4_all+STRINGasis":      dict(tf_target=1, gg_directed=1, mirna_mirna="asis", gg_string="asis"),
 "V5_all+STRINGboth":      dict(tf_target=1, gg_directed=1, mirna_mirna="asis", gg_string="both"),
 "V6_all_both":            dict(tf_target=1, gg_directed=1, mirna_mirna="both", gg_string="both"),
 "V7_TTonly":              dict(tf_target=1),
 "V8_mm_asis_only":        dict(mirna_mirna="asis"),
}

if __name__ == "__main__":
    rows = []
    for name, v in VARIANTS.items():
        E, prov, added = build(v)
        E2, comp = reduce_graph(E)
        c, recip = summarise(E)
        tot, cls, ordG, ordR = census3(E)
        rows.append((name, len(E), len(E2), len(recip), tot, ordG, ordR, dict(cls)))
        print("%-24s |E|=%5d |E_reduced|=%5d recip=%4d  D1-D4_n3=%6d  ordered_cores_G=%6d ordered_cores_reduced=%6d"
              % (name, len(E), len(E2), len(recip), tot, ordG, ordR))
        print("      classes:", dict(cls))
    import json
    json.dump([{"variant": r[0], "n_edges": r[1], "n_edges_reduced": r[2], "recip": r[3],
                "ffl3": r[4], "ordered_cores_full": r[5], "ordered_cores_reduced": r[6],
                "classes": r[7]} for r in rows],
              open("/path/to/revision/results/v2/n3_variants.json", "w"), indent=1)
