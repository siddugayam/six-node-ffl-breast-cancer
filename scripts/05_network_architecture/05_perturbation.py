#!/usr/bin/env python3
"""(D) PERTURBATION AND CONTROL: in-silico single- and double-node knockouts.

Damage metrics, all computed on a node universe that EXCLUDES the knocked-out nodes from
both the wild-type and the perturbed count, so that the score measures loss of routing and
not merely the loss of the deleted node itself:

  1. dFFL     fraction of 3-node FFL cores destroyed (cores are 3-node, so a single KO
              destroys exactly the cores it belongs to; the informative quantity is how
              far a PAIR falls short of / exceeds additivity)
  2. dReach   fraction of ordered reachable pairs (u,v), u,v not in KO set, that are lost
  3. dCOL     fraction of total absolute signed influence arriving at COL1A1 + COL3A1 that
              is lost.  Influence S = (I - (1-lambda) W)^-1 with W[u,v] = sign(u,v)/outdeg(u),
              lambda = 0.15; spectral radius of (1-lambda)W <= 0.85 so the Neumann series
              converges and the closed-form inverse is exact.

Composite damage = mean of the three fractions.

Druggability from DGIdb (interactions.tsv / drugs.tsv).  Combinations are tested against a
null of random node pairs matched on the marginal single-KO damage, which is the correct
control for Iyengar's claim that multi-target perturbation beats single-target inhibition
(Azeloglu & Iyengar 2015 CSH Perspect Biol 7:a005934).
"""
import os, sys, csv, json, collections, itertools, random, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import load_nodes, load_edges, RES, FIG, ffl_cores, ffl_unique_count

DGI = '/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data'
rng = random.Random(20260908)
np.random.seed(20260908)
LAM = 0.15
COLS = ['COL1A1', 'COL3A1']

nt = load_nodes()
E = load_edges()
P = [(e['source'], e['target']) for e in E]
sgn = {(e['source'], e['target']): e['sign_published'] for e in E}
NODES = sorted(nt); idx = {n: i for i, n in enumerate(NODES)}; N = len(NODES)
out_adj = collections.defaultdict(set)
for u, v in P:
    out_adj[u].add(v)

# ------------------------------------------------------------------ reachability (bitsets)
import networkx as nx

def reach_bits(alive_mask):
    """reach[i] = python-int bitset of nodes reachable from i by a path of length >= 1
    through alive nodes.  Computed on the SCC condensation (a DAG) in reverse topological
    order, so it is exact, not a truncation."""
    alive = [n for n in NODES if alive_mask[idx[n]]]
    aset = set(alive)
    g = nx.DiGraph(); g.add_nodes_from(alive)
    g.add_edges_from((u, v) for u, v in P if u in aset and v in aset)
    C = nx.condensation(g)                       # C.nodes have 'members'
    memb = {c: C.nodes[c]['members'] for c in C.nodes}
    cself = {}
    for c, ms in memb.items():
        b = 0
        for m in ms:
            b |= (1 << idx[m])
        cself[c] = b
    cbits = {}
    for c in reversed(list(nx.topological_sort(C))):
        acc = 0
        for s2 in C.successors(c):
            acc |= cbits[s2] | cself[s2]
        cbits[c] = acc
    res = {}
    for c, ms in memb.items():
        internal = cself[c] if len(ms) > 1 else 0
        for m in ms:
            r = cbits[c] | internal
            if len(ms) == 1 and g.has_edge(m, m):
                r |= (1 << idx[m])
            res[m] = r
    return res


def n_reach_pairs(alive_mask, restrict_mask):
    r = reach_bits(alive_mask)
    rm = 0
    for i in range(N):
        if restrict_mask[i]:
            rm |= (1 << i)
    tot = 0
    for n, b in r.items():
        if restrict_mask[idx[n]]:
            tot += bin(b & rm & ~(1 << idx[n])).count('1')
    return tot


# ------------------------------------------------------------------ signed influence
A_sgn = np.zeros((N, N))
for (u, v), s in sgn.items():
    A_sgn[idx[u], idx[v]] = s
outdeg = np.abs(A_sgn).sum(1)


def influence_to_cols(alive_mask):
    ai = np.where(alive_mask)[0]
    W = A_sgn[np.ix_(ai, ai)]
    od = np.abs(W).sum(1)
    nz = od > 0
    W[nz] = W[nz] / od[nz, None]
    S = np.linalg.inv(np.eye(len(ai)) - (1 - LAM) * W)
    np.fill_diagonal(S, 0.0)
    pos = {NODES[j]: k for k, j in enumerate(ai)}
    tot = 0.0
    for c in COLS:
        if c in pos:
            tot += float(np.abs(S[:, pos[c]]).sum())
    return tot


allmask = np.ones(N, bool)
FFL_WT = ffl_unique_count(ffl_cores(set(P), nt))
t0 = time.time()
R_WT_full = n_reach_pairs(allmask, allmask)
print('reachable ordered pairs (WT):', R_WT_full, f'{time.time()-t0:.2f}s')
INF_WT = influence_to_cols(allmask)
print('WT influence into collagens:', round(INF_WT, 3))
spec = np.max(np.abs(np.linalg.eigvals((1 - LAM) * np.where(outdeg[:, None] > 0,
                                                            A_sgn / np.maximum(outdeg, 1)[:, None],
                                                            0))))
print('spectral radius (1-lambda)W =', round(float(spec), 4))


def damage(ko):
    ko = list(ko)
    m = allmask.copy()
    for n in ko:
        m[idx[n]] = False
    # FFL
    keep = {(u, v) for u, v in P if u not in ko and v not in ko}
    f = ffl_unique_count(ffl_cores(keep, nt))
    # reachability, restricted to surviving nodes in BOTH graphs
    r_wt = n_reach_pairs(allmask, m)
    r_ko = n_reach_pairs(m, m)
    inf_ko = influence_to_cols(m)
    inf_wt = INF_WT if not (set(ko) & set(COLS)) else None
    dF = (FFL_WT - f) / FFL_WT
    dR = (r_wt - r_ko) / r_wt if r_wt else 0.0
    if inf_wt is None or inf_wt == 0:
        dC = float('nan')
    else:
        dC = (inf_wt - inf_ko) / inf_wt
    return dict(ffl_after=f, dFFL=dF, reach_wt_restricted=r_wt, reach_after=r_ko, dReach=dR,
                influence_after=inf_ko, dCOL=dC)


# ------------------------------------------------------------------ single knockouts
t0 = time.time()
single = {}
for k, n in enumerate(NODES):
    single[n] = damage([n])
    if (k + 1) % 100 == 0:
        print(f'  single KO {k+1}/{N}  {time.time()-t0:.1f}s', flush=True)
print(f'single KOs done in {time.time()-t0:.1f}s')

# druggability
drug_targets = collections.defaultdict(set); appr_targets = collections.defaultdict(set)
anti_targets = collections.defaultdict(set)
with open(f'{DGI}/interactions.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        g = r['gene_name']
        if g in nt:
            drug_targets[g].add(r['drug_name'])
            if r['approved'] == 'TRUE':
                appr_targets[g].add(r['drug_name'])
            if r['anti_neoplastic'] == 'TRUE':
                anti_targets[g].add(r['drug_name'])
print('druggable network nodes:', len(drug_targets), 'with approved drug:', len(appr_targets))

cores = ffl_cores(set(P), nt)
part = collections.Counter()
for R, M, T, c in cores:
    part[R] += 1; part[M] += 1; part[T] += 1
deg = collections.Counter()
for u, v in P:
    deg[u] += 1; deg[v] += 1

def composite(d):
    vals = [d['dFFL'], d['dReach'], 0.0 if np.isnan(d['dCOL']) else d['dCOL']]
    return float(np.mean(vals))

rows = []
for n in NODES:
    d = single[n]
    rows.append(dict(node=n, type=nt[n], degree=deg[n], ffl_participation=part[n],
                     dFFL=round(d['dFFL'], 6), dReach=round(d['dReach'], 6),
                     dCOL=('' if np.isnan(d['dCOL']) else round(d['dCOL'], 6)),
                     composite_damage=round(composite(d), 6),
                     n_drugs=len(drug_targets.get(n, ())),
                     n_approved_drugs=len(appr_targets.get(n, ())),
                     n_antineoplastic_drugs=len(anti_targets.get(n, ())),
                     druggable=n in drug_targets,
                     example_approved_drugs=';'.join(sorted(appr_targets.get(n, ()))[:5])))
rows.sort(key=lambda r: -r['composite_damage'])
with open(f'{RES}/systems_single_knockout.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('\nTop 15 single knockouts by composite damage:')
for r in rows[:15]:
    print(f"  {r['node']:<14} {r['type']:<6} dFFL={r['dFFL']:.3f} dReach={r['dReach']:.4f} "
          f"dCOL={r['dCOL']} comp={r['composite_damage']:.4f} drugs={r['n_drugs']}")
print('\nTop 10 DRUGGABLE single knockouts:')
for r in [x for x in rows if x['druggable']][:10]:
    print(f"  {r['node']:<14} comp={r['composite_damage']:.4f} approved={r['n_approved_drugs']} "
          f"{r['example_approved_drugs'][:60]}")

json.dump(dict(FFL_WT=FFL_WT, R_WT=R_WT_full, INF_WT=INF_WT, lambda_damping=LAM,
               spectral_radius=float(spec), n_druggable=len(drug_targets)),
          open(f'{RES}/systems_perturbation_meta.json', 'w'), indent=2)
np.save(f'{RES}/_single_damage.npy',
        np.array([[single[n]['dFFL'], single[n]['dReach'],
                   0.0 if np.isnan(single[n]['dCOL']) else single[n]['dCOL']] for n in NODES]))
with open(f'{RES}/_nodes_order.json', 'w') as fh:
    json.dump(NODES, fh)
print('single-KO stage complete')
