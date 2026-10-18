#!/usr/bin/env python3
"""Benchmark of the topological-definition census enumerator on planted modules (see SETTINGS.md).

For each replicate:
- background: the analysed network (data/network/canonical_edges.tsv without the 30 exemplar miRNA-miRNA edges, each
  reciprocal TF<->miRNA pair contracted to its TF->miRNA arc) as an undirected simple graph, rewired by degree-preserving
  double-edge swaps (100 attempts per edge) and oriented by a random vertex order (a DAG);
- 10 planted modules per n = 3-6 on disjoint random vertex sets, drawn from all the topological definition topologies on n vertices;
- noise arcs (10 % of the edge count) between unlinked vertex pairs, oriented by the same order;
- census_benchmark (strict, d3relax, d4relax) at n = 3-6; an independent Python test of the planted modules and of a
  sample of 1,000 reported modules; an independent exhaustive enumeration at n = 3 (all replicates) and n = 4 (first 5).
Writes planted_topologies.csv, benchmark_per_replicate.csv and benchmark_summary.csv.
usage: FFL_ROOT=<repository checkout> python3 benchmark.py [replicates, default 200] [parallel processes, default 4]"""
import csv, itertools, json, os, subprocess, sys
from multiprocessing import Pool
import numpy as np
import networkx as nx

ROOT = os.environ["FFL_ROOT"]
HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.join(HERE, "work"); os.makedirs(W, exist_ok=True)
REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 200
PROCS = int(sys.argv[2]) if len(sys.argv) > 2 else 4
SEED0, PLANT, NOISE, SAMPLE = 20261015, 10, 0.10, 1000
EXE = os.path.join(W, "census_benchmark")


# ---- independent topological-definition test (networkx) --------------------------------------------------------------------------
def d1234(H):
    if not nx.is_weakly_connected(H): return False
    if not nx.is_directed_acyclic_graph(H): return False
    src = [v for v in H if H.in_degree(v) == 0]; snk = [v for v in H if H.out_degree(v) == 0]
    if len(src) != 1 or len(snk) != 1 or src[0] == snk[0]: return False
    s, t = src[0], snk[0]
    if (nx.descendants(H, s) | {s}) != set(H) or (nx.ancestors(H, t) | {t}) != set(H): return False
    if H.has_edge(s, t):
        H2 = H.copy(); H2.remove_edge(s, t); return nx.has_path(H2, s, t)
    return nx.node_connectivity(H, s, t) >= 2


def topologies(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    seen, out = set(), []
    for mask in range(1, 1 << len(pairs)):
        arcs = [pairs[k] for k in range(len(pairs)) if mask >> k & 1]
        H = nx.DiGraph(); H.add_nodes_from(range(n)); H.add_edges_from(arcs)
        if not d1234(H): continue
        canon = min(tuple(sorted((p[a], p[b]) for a, b in arcs)) for p in itertools.permutations(range(n)))
        if canon in seen: continue
        seen.add(canon); out.append(arcs)
    return out


def skeleton():
    nodes = {r["name"]: r["type"] for r in csv.DictReader(open(os.path.join(ROOT, "data", "network", "canonical_nodes.tsv")), delimiter="\t")}
    arcs = {(r["source"], r["target"]) for r in csv.DictReader(open(os.path.join(ROOT, "data", "network", "canonical_edges.tsv")), delimiter="\t")
            if r["edge_type"] != "miRNA_miRNA"}
    for a, b in list(arcs):
        if nodes[a] == "TF" and nodes[b] == "miRNA" and (b, a) in arcs: arcs.discard((b, a))
    V = sorted(nodes); ix = {v: i for i, v in enumerate(V)}
    E = sorted({tuple(sorted((ix[a], ix[b]))) for a, b in arcs if a != b})
    return len(V), E, len(arcs)


NV, EDGES, NARCS = skeleton()


def esu_all(n_nodes, und, k, cb):
    for v in range(n_nodes):
        _ext([v], {u for u in und[v] if u > v}, v, und, k, cb)


def _ext(sub, ext, v, und, k, cb):
    if len(sub) == k: cb(sub); return
    ext = set(ext)
    while ext:
        w = ext.pop()
        nb = set(sub) | {u for s in sub for u in und[s]}
        _ext(sub + [w], ext | {u for u in und[w] if u > v and u not in nb}, v, und, k, cb)


def key(S):
    k = 0
    for x in sorted(S): k = (k << 10) | int(x)
    return k


def one_replicate(r):
    done = os.path.join(W, f"result_{r:03d}.json")                 # finished replicates are kept, so a stopped run resumes
    if os.path.exists(done):
        return json.load(open(done))
    rng = np.random.default_rng(SEED0 + r)
    E = [list(e) for e in EDGES]; m = len(E); eset = {tuple(e) for e in EDGES}
    for _ in range(100 * m):                                     # degree-preserving double-edge swaps
        i, j = rng.integers(m), rng.integers(m)
        if i == j: continue
        a, b = E[i]; c, d = E[j]
        if rng.random() < 0.5: c, d = d, c
        if len({a, b, c, d}) < 4: continue
        e1, e2 = tuple(sorted((a, d))), tuple(sorted((c, b)))
        if e1 in eset or e2 in eset: continue
        eset.discard(tuple(sorted((a, b)))); eset.discard(tuple(sorted((c, d)))); eset.add(e1); eset.add(e2)
        E[i], E[j] = list(e1), list(e2)
    rank = np.empty(NV, int); rank[rng.permutation(NV)] = np.arange(NV)
    arcs = {(a, b) if rank[a] < rank[b] else (b, a) for a, b in eset}
    background_edges = len(arcs)
    free = list(rng.permutation(NV)); planted = []
    for n in (3, 4, 5, 6):
        for _ in range(PLANT):
            S = [int(free.pop()) for _ in range(n)]
            S.sort(key=lambda v: rank[v])
            arcs = {(a, b) for a, b in arcs if not (a in S and b in S)}
            T = TOPS[n][int(rng.integers(len(TOPS[n])))]
            arcs |= {(S[i], S[j]) for i, j in T}
            planted.append((n, tuple(S)))
    added, target = 0, int(np.rint(NOISE * background_edges))
    while added < target:
        a, b = int(rng.integers(NV)), int(rng.integers(NV))
        if a == b or (a, b) in arcs or (b, a) in arcs: continue
        arcs.add((a, b) if rank[a] < rank[b] else (b, a)); added += 1
    G = nx.DiGraph(); G.add_nodes_from(range(NV)); G.add_edges_from(arcs)
    assert nx.is_directed_acyclic_graph(G)
    fn = os.path.join(W, f"graph_{r:03d}.txt")
    with open(fn, "w") as fh:
        fh.write(f"{NV} {len(arcs)} 1\n" + " ".join("2" for _ in range(NV)) + "\n")
        for a, b in sorted(arcs): fh.write(f"{a} {b} 0\n")
    res = dict(replicate=r, background_edges=background_edges, noise_arcs=added, arcs=len(arcs))
    keys = {}
    for n in (3, 4, 5, 6):
        for mode in ("strict", "d3relax", "d4relax"):
            kf = os.path.join(W, f"keys_{r:03d}_{n}_{mode}.bin")
            subprocess.run([EXE, fn, str(n), mode, kf], check=True, capture_output=True)
            keys[(n, mode)] = np.fromfile(kf, dtype=np.uint64); os.remove(kf)
        ks = set(keys[(n, "strict")].tolist())
        P = [S for nn, S in planted if nn == n]
        intact = [S for S in P if d1234(G.subgraph(S).copy())]
        res[f"n{n}_planted"] = len(P); res[f"n{n}_intact"] = len(intact)
        res[f"n{n}_recall_intact"] = sum(key(S) in ks for S in intact) / len(intact) if intact else float("nan")
        res[f"n{n}_destroyed_reported"] = sum(key(S) in ks for S in P if S not in intact)
        allk = keys[(n, "strict")]
        smp = allk if len(allk) <= SAMPLE else rng.choice(allk, SAMPLE, replace=False)
        ok = 0
        for k in smp.tolist():
            S = [(k >> (10 * (n - 1 - i))) & 1023 for i in range(n)]
            ok += d1234(G.subgraph(S).copy())
        res[f"n{n}_reported"] = len(allk); res[f"n{n}_sample"] = len(smp)
        res[f"n{n}_precision_sample"] = ok / len(smp) if len(smp) else float("nan")
        for mode in ("d3relax", "d4relax"):
            km = set(keys[(n, mode)].tolist())
            res[f"n{n}_{mode}_reported"] = len(km)
            res[f"n{n}_{mode}_ratio"] = len(km) / len(allk) if len(allk) else float("nan")
            res[f"n{n}_{mode}_contains_strict"] = ks <= km
            res[f"n{n}_{mode}_recall_intact"] = sum(key(S) in km for S in intact) / len(intact) if intact else float("nan")
            res[f"n{n}_{mode}_share_planted"] = sum(key(S) in km for S in P) / len(km) if km else float("nan")
        res[f"n{n}_share_planted"] = sum(key(S) in ks for S in P) / len(ks) if ks else float("nan")
    und = [set() for _ in range(NV)]
    for a, b in arcs: und[a].add(b); und[b].add(a)
    for n in ((3, 4) if r < 5 else (3,)):
        truth = set()
        esu_all(NV, und, n, lambda sub: truth.add(key(sub)) if d1234(G.subgraph(sub).copy()) else None)
        ks = set(keys[(n, "strict")].tolist())
        res[f"n{n}_exhaustive_truth"] = len(truth)
        res[f"n{n}_exhaustive_recall"] = len(truth & ks) / len(truth) if truth else float("nan")
        res[f"n{n}_exhaustive_precision"] = len(truth & ks) / len(ks) if ks else float("nan")
    os.remove(fn)
    json.dump(res, open(done + ".part", "w")); os.replace(done + ".part", done)
    print(f"replicate {r} done", flush=True)
    return res


TOPS = {n: topologies(n) for n in (3, 4, 5, 6)}

if __name__ == "__main__":
    subprocess.run(["gcc", "-O3", "-march=native", "-o", EXE, os.path.join(HERE, "census_benchmark.c")], check=True)
    with open(os.path.join(HERE, "planted_topologies.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["n", "D1_D4_topologies_up_to_isomorphism"])
        for n in (3, 4, 5, 6): w.writerow([n, len(TOPS[n])])
    print({n: len(t) for n, t in TOPS.items()}, "skeleton edges", len(EDGES), "contracted arcs", NARCS, flush=True)
    with Pool(PROCS) as pool:
        R = pool.map(one_replicate, range(REPS), chunksize=1)
    cols = sorted({k for r in R for k in r}, key=lambda c: (c != "replicate", c))
    with open(os.path.join(HERE, "benchmark_per_replicate.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(R)
    summ = []
    for c in cols[1:]:
        v = np.array([float(r[c]) for r in R if c in r and r[c] == r[c]], float)
        if not len(v): continue
        summ.append(dict(measure=c, replicates=len(v), median=float(np.median(v)), min=float(v.min()), max=float(v.max()),
                         mean=float(v.mean())))
    with open(os.path.join(HERE, "benchmark_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0])); w.writeheader(); w.writerows(summ)
    print("done")
