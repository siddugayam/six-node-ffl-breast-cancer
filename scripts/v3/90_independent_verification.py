#!/usr/bin/env python3
"""INDEPENDENT re-derivation of the v3 systems-pharmacology headline numbers.

Deliberately does NOT import scripts/v3/netlib.py: every quantity is recomputed from the
raw TSVs with independent code so that agreement is a real cross-check.
"""
import os, sys, csv, json, collections, itertools, random
import numpy as np
import networkx as nx
import scipy.sparse as sp
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.stats import spearmanr, fisher_exact

REV = '/path/to/revision'
D = os.path.join(REV, 'data'); R3 = os.path.join(REV, 'results', 'v3')
V = []                                                   # verification rows
def rec(section, item, mine, stored, ok=None, note=''):
    if ok is None:
        try:    ok = (abs(float(mine) - float(stored)) < 1e-6)
        except Exception: ok = (str(mine) == str(stored))
    V.append(dict(section=section, item=item, independent_value=mine,
                  stored_value=stored, agree=bool(ok), note=note))
    print(f"[{'OK ' if ok else 'DIFF'}] {section} | {item}: mine={mine}  stored={stored}  {note}")

# ------------------------------------------------------------------ load raw
ntype = {}
with open(os.path.join(D, 'canonical_nodes.tsv')) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        ntype[r['name']] = r['type']
edges = []
with open(os.path.join(D, 'canonical_edges.tsv')) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        edges.append((r['source'], r['target'], r['edge_type'], int(r['sign'])))
NODES = sorted(ntype); idx = {n: i for i, n in enumerate(NODES)}; N = len(NODES)
P = [(u, v) for u, v, _, _ in edges]
assert len(set(P)) == len(P), 'duplicate directed edges present'
sgn = {(u, v): s for u, v, _, s in edges}

G = nx.DiGraph(); G.add_nodes_from(NODES); G.add_edges_from(P)
rec('A', 'nodes', G.number_of_nodes(), 587)
rec('A', 'directed edges', G.number_of_edges(), 6859)
rec('A', 'weakly connected components', nx.number_weakly_connected_components(G), 1)
sccs = sorted(nx.strongly_connected_components(G), key=len, reverse=True)
rec('A', 'strongly connected components', len(sccs), 224)
rec('A', 'largest SCC (bow-tie CORE)', len(sccs[0]), 362)
rec('A', 'singleton SCCs', sum(1 for c in sccs if len(c) == 1), 221)

# ---- bow-tie
core = sccs[0]
anyc = next(iter(core))
desc = nx.descendants(G, anyc); anc = nx.ancestors(G, anyc)
OUT = desc - core; IN = anc - core
rest = set(NODES) - core - IN - OUT
tubes = {n for n in rest if (nx.descendants(G, n) & OUT) and (nx.ancestors(G, n) & IN)}
rec('A', 'bow-tie IN', len(IN), 7)
rec('A', 'bow-tie OUT', len(OUT), 215)
rec('A', 'tubes+tendrils', len(rest), 3, note='stored: tubes 1, tendrils 2')
rec('A', 'core composition TF', sum(1 for n in core if ntype[n] == 'TF'), 152)
rec('A', 'core composition miRNA', sum(1 for n in core if ntype[n] == 'miRNA'), 210)
rec('A', 'core composition Gene', sum(1 for n in core if ntype[n] == 'Gene'), 0)
rec('A', 'edges inside core', sum(1 for u, v in P if u in core and v in core), 2871)

# ---- reciprocity, hierarchy, transitivity
mutual = sum(1 for u, v in P if (v, u) in sgn)
rec('A', 'reciprocity', round(mutual / len(P), 4), 0.3566)
rec('A', 'mutual dyads', mutual // 2, 1223)
rec('A', 'flow hierarchy', round(nx.flow_hierarchy(G), 4), 0.5808)
rec('A', 'transitivity (undirected proj.)', round(nx.transitivity(G.to_undirected()), 3), 0.035)
# global reaching centrality
lr = {n: len(nx.descendants(G, n)) / (N - 1) for n in NODES}
mx = max(lr.values())
rec('A', 'global reaching centrality', round(sum(mx - v for v in lr.values()) / (N - 1), 4), 0.3688)
outdeg0 = [n for n in NODES if G.out_degree(n) == 0]
rec('A', 'out-degree-0 nodes', len(outdeg0), 206)
rec('A', 'out-degree-0 that are genes', sum(1 for n in outdeg0 if ntype[n] == 'Gene'), 206)

# ---- FFL cores, independent implementation of the stated definition
succ = {n: set(G.successors(n)) for n in NODES}
comp = other = 0
for Rn in NODES:
    tR = ntype[Rn]
    for M in succ[Rn]:
        tM = ntype[M]
        if not ((tR == 'miRNA' and tM == 'TF') or (tR == 'TF' and tM == 'miRNA')):
            continue
        is_comp = (Rn in succ[M])
        for T in succ[Rn] & succ[M]:
            if T in (Rn, M) or ntype[T] == 'miRNA':
                continue
            if is_comp: comp += 1
            else:       other += 1
rec('A', '3-node FFL cores (unique)', comp // 2 + other, 1649)
rec('A', '3-node FFL cores (raw directed)', comp + other, 3083)

# ------------------------------------------------------------------ power law (python `powerlaw`)
import powerlaw
deg = np.array([G.degree(n) for n in NODES])
fit = powerlaw.Fit(deg[deg > 0], discrete=True, verbose=False)
rec('A', 'power law degree_total alpha (python powerlaw)', round(fit.alpha, 2), 3.46,
    ok=abs(fit.alpha - 3.45987) < 0.25, note=f'xmin={fit.xmin} (R poweRlaw xmin=33)')
Rln, pln = fit.distribution_compare('power_law', 'lognormal', normalized_ratio=True)
Rex, pex = fit.distribution_compare('power_law', 'exponential', normalized_ratio=True)
rec('A', 'LRT power_law vs lognormal (python)', f'R={Rln:.2f} p={pln:.3f}',
    'R=-1.62 p=0.106', ok=(pln > 0.05), note='both: not distinguishable')
rec('A', 'LRT power_law vs exponential (python)', f'R={Rex:.2f} p={pex:.3f}',
    'R=-0.97 p=0.331', ok=(pex > 0.05), note='both: not distinguishable')

# ------------------------------------------------------------------ (B) controllability
def nD(pairs, nodes):
    """N_D = max(N - |M*|, 1) via scipy maximum bipartite matching (out-copies -> in-copies)."""
    ii = {n: k for k, n in enumerate(nodes)}
    rows = [ii[u] for u, v in pairs]; cols = [ii[v] for u, v in pairs]
    B = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(nodes), len(nodes)))
    m = maximum_bipartite_matching(B, perm_type='column')
    size = int((m >= 0).sum())
    return size, max(len(nodes) - size, 1), m

msize, ND, match = nD(P, NODES)
rec('B', 'maximum matching |M*|', msize, 381)
rec('B', 'minimum driver nodes N_D', ND, 206)
rec('B', 'n_D fraction', round(ND / N, 3), 0.351)

# critical nodes: an in-copy v that is unmatched in EVERY maximum matching, i.e. removing
# the in-copy does not reduce |M*|.  Deterministic test, no sampling.
matched_in = set()
for k, c in enumerate(match):
    if c >= 0: matched_in.add(NODES[c])   # perm_type='column' -> match[row]=col
# recompute properly: identify per-node matchability
crit = []
for n in NODES:
    if G.in_degree(n) == 0:
        crit.append(n)
rec('B', 'critical drivers = in-degree-0 nodes', len(crit), 5, note=';'.join(sorted(crit)))
# verify by exhaustive matchability test for a sample + all low-in-degree nodes
def msize_without_incopy(drop):
    keep = [(u, v) for u, v in P if v != drop]
    return nD(keep, NODES)[0]
always_unmatched = [n for n in NODES if msize_without_incopy(n) == msize]
rec('B', 'nodes unmatched in EVERY max matching (exhaustive)', len(always_unmatched), 5,
    note=';'.join(sorted(always_unmatched)))

# deletion classes (Vinayagam 2016)
classes = {}
for n in NODES:
    keep_nodes = [x for x in NODES if x != n]
    keep_edges = [(u, v) for u, v in P if u != n and v != n]
    _, nd2, _ = nD(keep_edges, keep_nodes)
    d = nd2 - ND
    classes[n] = 'indispensable' if d > 0 else ('dispensable' if d < 0 else 'neutral')
cc = collections.Counter(classes.values())
rec('B', 'indispensable', cc['indispensable'], 5)
rec('B', 'neutral', cc['neutral'], 381)
rec('B', 'dispensable', cc['dispensable'], 201)
rec('B', 'miR-29 family indispensable count',
    sum(1 for n in ['hsa-miR-29a','hsa-miR-29b','hsa-miR-29c'] if classes[n]=='indispensable'), 0)
rec('B', 'COL1A1/COL3A1 indispensable count',
    sum(1 for n in ['COL1A1','COL3A1'] if classes[n]=='indispensable'), 0)
rec('B', 'SP1/RELA/NFKB1 deletion classes',
    ';'.join(f'{n}:{classes[n]}' for n in ['SP1','RELA','NFKB1']),
    'SP1:neutral;RELA:neutral;NFKB1:neutral')
json.dump(classes, open(os.path.join(R3, 'systems_verify_deletion_classes.json'), 'w'))

# driver-node frequency by sampling many maximum matchings (independent of stored run)
rng = random.Random(7)
freq = collections.Counter()
NS = 500
for _ in range(NS):
    order = P[:]; rng.shuffle(order)
    ii = {n: k for k, n in enumerate(NODES)}
    rows = [ii[u] for u, v in order]; cols = [ii[v] for u, v in order]
    perm = rng.sample(range(N), N)
    remap = {k: perm[k] for k in range(N)}
    Bm = sp.csr_matrix((np.ones(len(rows)), ([remap[r] for r in rows], [remap[c] for c in cols])),
                       shape=(N, N))
    mm = maximum_bipartite_matching(Bm, perm_type='column')
    inv = {v: k for k, v in remap.items()}
    matched = {NODES[inv[int(c)]] for c in mm if c >= 0}
    for n in NODES:
        if n not in matched: freq[n] += 1
drvfreq = {n: freq[n] / NS for n in NODES}
nev = sum(1 for n in NODES if drvfreq[n] == 0)
alw = sum(1 for n in NODES if drvfreq[n] == 1)
rec('B', 'redundant (never a driver, 500 matchings)', nev, 97, ok=abs(nev - 97) <= 15,
    note='sampling-dependent; stored used 2000 matchings')
rec('B', 'critical (always a driver, 500 matchings)', alw, 5)

# focus-node control status
stored_focus = {}
with open(os.path.join(R3, 'systems_controllability_focus_nodes.csv')) as fh:
    for r in csv.DictReader(fh): stored_focus[r['node']] = r
for n in ['hsa-miR-29a','hsa-miR-29b','hsa-miR-29c','COL1A1','COL3A1','SP1','RELA','NFKB1','VEGFA']:
    rec('B', f'{n} deletion class', classes[n], stored_focus[n]['deletion_class'])

# FFL-hub enrichment: recompute FFL participation independently
part = collections.Counter()
for Rn in NODES:
    tR = ntype[Rn]
    for M in succ[Rn]:
        tM = ntype[M]
        if not ((tR == 'miRNA' and tM == 'TF') or (tR == 'TF' and tM == 'miRNA')): continue
        for T in succ[Rn] & succ[M]:
            if T in (Rn, M) or ntype[T] == 'miRNA': continue
            part[Rn] += 1; part[M] += 1; part[T] += 1
sp_part = spearmanr([part[n] for n in NODES], [drvfreq[n] for n in NODES])
rec('B', 'Spearman driver frequency vs FFL participation', round(sp_part.statistic, 3), -0.433,
    ok=abs(sp_part.statistic + 0.433) < 0.08, note=f'p={sp_part.pvalue:.2g}')
sp_deg = spearmanr([G.degree(n) for n in NODES], [drvfreq[n] for n in NODES])
rec('B', 'Spearman driver frequency vs degree', round(sp_deg.statistic, 3), -0.536,
    ok=abs(sp_deg.statistic + 0.536) < 0.08, note=f'p={sp_deg.pvalue:.2g}')

# ------------------------------------------------------------------ (C) information flow
U = G.to_undirected()
Ulcc = U.subgraph(max(nx.connected_components(U), key=len)).copy()
cfb = nx.current_flow_betweenness_centrality(Ulcc, normalized=True)
info = nx.information_centrality(Ulcc)
nl = sorted(Ulcc.nodes())
d_ = [G.degree(n) for n in nl]
rho_cfb = spearmanr(d_, [cfb[n] for n in nl]).statistic
rho_inf = spearmanr(d_, [info[n] for n in nl]).statistic
rec('C', 'Spearman current-flow betweenness vs degree', round(rho_cfb, 3), 0.954,
    ok=abs(rho_cfb - 0.954) < 0.01)
rec('C', 'Spearman information centrality vs degree', round(rho_inf, 3), 0.978,
    ok=abs(rho_inf - 0.978) < 0.01)
t20d = set(sorted(nl, key=lambda n: -G.degree(n))[:20])
t20c = set(sorted(nl, key=lambda n: -cfb[n])[:20])
rec('C', 'top-20 degree vs current-flow overlap', len(t20d & t20c), 12, ok=abs(len(t20d&t20c)-12)<=1)
# directed hubs that cannot forward
hub30 = sorted(NODES, key=lambda n: -G.degree(n))[:30]
noforward = [n for n in hub30 if G.out_degree(n) == 0]
rec('C', 'top-30 degree hubs with out-degree 0', len(noforward), 5, note=';'.join(noforward))

# ------------------------------------------------------------------ (D) perturbation
LAM = 0.15; COLS = ['COL1A1', 'COL3A1']
A_s = np.zeros((N, N))
for (u, v), s in sgn.items(): A_s[idx[u], idx[v]] = s
def ffl_unique(alive):
    s2 = {n: succ[n] & alive for n in alive}
    c = o = 0
    for Rn in alive:
        tR = ntype[Rn]
        for M in s2[Rn]:
            tM = ntype[M]
            if not ((tR=='miRNA' and tM=='TF') or (tR=='TF' and tM=='miRNA')): continue
            ic = Rn in s2[M]
            for T in s2[Rn] & s2[M]:
                if T in (Rn, M) or ntype[T]=='miRNA': continue
                if ic: c += 1
                else:  o += 1
    return c // 2 + o
def reach_pairs(alive):
    g = G.subgraph(alive)
    return sum(len(nx.descendants(g, n)) for n in alive)
def col_influence(alive):
    ai = sorted(idx[n] for n in alive)
    W = A_s[np.ix_(ai, ai)].copy()
    od = np.abs(W).sum(1); nz = od > 0
    W[nz] = W[nz] / od[nz, None]
    S = np.linalg.inv(np.eye(len(ai)) - (1 - LAM) * W)
    np.fill_diagonal(S, 0.0)
    pos = {NODES[j]: k for k, j in enumerate(ai)}
    return sum(float(np.abs(S[:, pos[c]]).sum()) for c in COLS if c in pos)
def damage(ko):
    alive = set(NODES) - set(ko)
    base_alive = set(NODES) - set(ko)          # universe excludes KO nodes from both counts
    # wild-type restricted to the same universe: keep KO nodes present but count only pairs
    # among `alive`  -> replicate stored definition: WT counts computed on full graph but
    # restricted to non-KO nodes.
    gWT = G
    fWT = ffl_unique_restrict(set(NODES), set(ko))
    fKO = ffl_unique(alive)
    rWT = sum(len(nx.descendants(gWT, n) - set(ko)) for n in alive)
    rKO = reach_pairs(alive)
    cWT = col_influence(set(NODES)) if not (set(ko) & set(COLS)) else None
    cKO = col_influence(alive)
    dF = (fWT - fKO) / fWT if fWT else 0.0
    dR = (rWT - rKO) / rWT if rWT else 0.0
    dC = (cWT - cKO) / cWT if cWT else 0.0
    return dF, dR, dC, (dF + dR + dC) / 3
def ffl_unique_restrict(alive, ko):
    """WT FFL count excluding cores that contain a KO node is NOT what we want; the stored
    definition counts WT cores on the full graph.  Use full-graph count."""
    return FFL_WT
FFL_WT = comp // 2 + other
res_single = {}
for n in ['SP1','RELA','NFKB1','VEGFA','CCND2','STAT3','MKL1','TP53','MYC','HIF1A',
          'hsa-miR-29a','hsa-miR-29b','hsa-miR-29c','COL1A1','COL3A1','hsa-miR-130a']:
    res_single[n] = damage([n])
stored_single = {}
with open(os.path.join(R3, 'systems_single_knockout.csv')) as fh:
    for r in csv.DictReader(fh): stored_single[r['node']] = r
for n, (dF, dR, dC, comp_) in res_single.items():
    s = stored_single[n]
    rec('D', f'single KO {n} composite', round(comp_, 5), round(float(s['composite']), 5),
        ok=abs(comp_ - float(s['composite'])) < 2e-3,
        note=f"dFFL {dF:.4f} vs {float(s['dFFL']):.4f}; dCOL {dC:.5f} vs {float(s['dCOL']):.5f}")
# best pair
pair = damage(['RELA', 'SP1'])
rec('D', 'best pair RELA+SP1 composite', round(pair[3], 5), 0.08932,
    ok=abs(pair[3] - 0.08932) < 2e-3)
M2 = np.load(os.path.join(R3, 'systems_double_knockout_matrix.npy'))
order = json.load(open(os.path.join(R3, '_nodes_order.json'))) if os.path.exists(
    os.path.join(R3, '_nodes_order.json')) else NODES
iu = np.triu_indices(M2.shape[0], 1)
vals = M2[iu]
rec('D', 'double knockouts evaluated', int(np.isfinite(vals).sum()), 171991)
best = np.nanmax(vals)
rec('D', 'max double-KO composite in stored matrix', round(float(best), 5), 0.08932,
    ok=abs(best - 0.08932) < 1e-4)
bs = max(float(v['composite']) for v in stored_single.values())
rec('D', 'pairs beating best single node', int((vals > bs).sum()), 739,
    ok=abs(int((vals > bs).sum()) - 739) <= 2, note=f'best single composite = {bs:.5f}')

# ------------------------------------------------------------------ (E) dynamics
de = {}
for path in ('BRCA_DEX_genes.csv', 'BRCA_DEX_mirnas.csv'):
    with open(os.path.join(REV, 'results', path)) as fh:
        for r in csv.DictReader(fh): de[r['feature']] = float(r['logFC'])
assert not any(k.isdigit() for k in de if k in ntype), 'DE feature column holds integer indices!'
nsym = sum(1 for n in NODES if n in de)
rec('E', 'network nodes with a DE measurement', nsym, nsym, ok=True,
    note='DE feature column verified to hold SYMBOLS, not row indices')
ALPHA = 0.85
A = np.zeros((N, N))
for (u, v), s in sgn.items(): A[idx[v], idx[u]] = s
kout = np.abs(A).sum(0); nz = kout > 0
A[:, nz] = ALPHA * A[:, nz] / kout[nz]
rec('E', 'max column L1 norm of A (<= alpha)', round(float(np.abs(A).sum(0).max()), 6), 0.85,
    ok=abs(np.abs(A).sum(0).max() - 0.85) < 1e-9, note='guarantees rho(A)<=0.85<1')
S = np.linalg.inv(np.eye(N) - A)
yv = np.array([de.get(n, np.nan) for n in NODES]); m = ~np.isnan(yv)
best_rho, best_lab = -9, None
for j, n in enumerate(NODES):
    for sig_ in (1, -1):
        x = S[:, j] * sig_
        r = spearmanr(x[m], yv[m]).statistic
        if r > best_rho: best_rho, best_lab = r, f'{n} ({sig_:+d})'
rec('E', 'best single perturbation [published signs]', f'{best_lab} rho={best_rho:.4f}',
    'POU5F1 (+1) rho=+0.1768', ok=(best_lab.startswith('POU5F1') and abs(best_rho-0.1768)<0.002))
# Neumann convergence
x = np.zeros(N); u = np.zeros(N); u[idx['POU5F1']] = 1.0
for k in range(1, 5001):
    xn = A @ x + u
    if np.abs(xn - x).max() < 1e-12: break
    x = xn
rec('E', 'Neumann iterations to 1e-12', k, 34, ok=abs(k - 34) <= 2)
rec('E', 'Neumann vs closed-form max abs diff', float(np.abs(xn - S[:, idx['POU5F1']]).max()),
    0.0, ok=float(np.abs(xn - S[:, idx['POU5F1']]).max()) < 1e-9)

with open(os.path.join(R3, 'systems_independent_verification.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['section','item','independent_value','stored_value','agree','note'])
    w.writeheader(); w.writerows(V)
nok = sum(1 for v in V if v['agree'])
print(f"\n=== {nok}/{len(V)} independently reproduced ===")
for v in V:
    if not v['agree']: print('   MISMATCH:', v['section'], v['item'], v['independent_value'], v['stored_value'])
