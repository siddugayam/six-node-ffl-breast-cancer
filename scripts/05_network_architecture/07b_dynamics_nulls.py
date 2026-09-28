#!/usr/bin/env python3
"""(E cont.) Nulls and self-consistency checks for the dynamic-response analysis.

  1. Boolean threshold model: is the 71.1% sign agreement between the attractor and the
     observed phenotype better than (a) a permuted seed, (b) a degree-preserving rewired
     network with the true seed?
  2. Best single-node clamp restricted to clamps that actually reach >= 50 measured nodes,
     with a binomial test and a permutation test.
  3. Self-consistency: does the best-fitting linear-response perturbation act in the same
     direction as the driver node's own observed tumour-vs-normal change?
  4. Focused readout: which perturbation best reproduces the observed COL1A1/COL3A1 change?
"""
import os, sys, csv, json, collections
import numpy as np
from scipy.stats import binomtest, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import load_nodes, load_edges, RES

ALPHA = 0.85
rng = np.random.default_rng(4242)
nt = load_nodes(); E = load_edges()
NODES = sorted(nt); idx = {n: i for i, n in enumerate(NODES)}; N = len(NODES)

de = {}
for p in ('BRCA_DEX_genes.csv', 'BRCA_DEX_mirnas.csv'):
    for r in csv.DictReader(open(f'/path/to/revision/results/{p}')):
        de[r['feature']] = (float(r['logFC']), float(r['adj.P.Val']))
# --- Guard against the project's known limma failure mode (topTable returning integer row
# indices when the expression matrix has duplicated rownames).  Checked properly:
#   (a) no duplicated feature names,
#   (b) the network's own gene/miRNA symbols are present,
#   (c) the only all-digit features are the 28 unmapped Entrez IDs that TCGA HiSeqV2 carries
#       as "?|<entrez>", none of which is a network node.
_feat = list(de)
assert len(_feat) == len(set(_feat)), 'duplicated feature names in the DE tables'
_dig = [k for k in _feat if k.isdigit()]
assert len(_dig) <= 40, f'{len(_dig)} all-digit features: topTable index artefact suspected'
assert not any(n.isdigit() for n in NODES), 'network node names are digits'
for _s in ('COL1A1', 'COL3A1', 'hsa-miR-29a', 'SP1', 'RELA'):
    assert _s in de, f'DE table missing expected symbol {_s}'
print(f'DE feature-name check passed: {len(_feat)} unique symbols, {len(_dig)} unmapped '
      f'Entrez IDs (none a network node)')
y = np.full(N, np.nan); fdr = np.full(N, np.nan)
for n in NODES:
    if n in de:
        y[idx[n]], fdr[idx[n]] = de[n]
meas = ~np.isnan(y)
sig = meas & (fdr < 0.05) & (np.abs(y) > 0.5)
print(f'measured {meas.sum()}  DE-significant {sig.sum()}')

rows = np.array([idx[e['source']] for e in E])
cols = np.array([idx[e['target']] for e in E])
vals = np.array([float(e['sign_published']) for e in E])
kout = np.bincount(rows, weights=np.abs(vals), minlength=N)
w = vals / np.where(kout[rows] > 0, kout[rows], 1.0)


def boolean_attractor(seed, r_, c_, w_, clamp=None, maxsteps=300):
    x = seed.copy()
    if clamp:
        for k, v in clamp.items():
            x[idx[k]] = v
    indeg0 = np.bincount(c_, minlength=N) == 0
    hist = {}
    for t in range(maxsteps):
        key = x.tobytes()
        if key in hist:
            return x, t, t - hist[key]
        hist[key] = t
        ax = np.bincount(c_, weights=w_ * x[r_], minlength=N)
        xn = np.sign(ax)
        xn[indeg0] = x[indeg0]
        if clamp:
            for k, v in clamp.items():
                xn[idx[k]] = v
        x = xn
    return x, maxsteps, -1


def agreement(att):
    m = sig & (np.abs(att) > 0)
    if m.sum() == 0:
        return np.nan, 0
    return float((np.sign(att[m]) == np.sign(y[m])).mean()), int(m.sum())


seed = np.zeros(N); seed[sig] = np.sign(y[sig])
att, steps, period = boolean_attractor(seed, rows, cols, w)
obs_acc, obs_n = agreement(att)
print(f'observed: attractor after {steps} steps, period {period}, agreement {obs_acc:.4f} on n={obs_n}')

# --- null 1: permuted seed signs (same number of +1/-1, reassigned among DE nodes)
NP = 2000
null1 = np.empty(NP)
sidx = np.where(sig)[0]
sv = np.sign(y[sig])
for i in range(NP):
    s2 = np.zeros(N); s2[sidx] = rng.permutation(sv)
    a2, _, _ = boolean_attractor(s2, rows, cols, w, maxsteps=120)
    m = sig & (np.abs(a2) > 0)
    null1[i] = (np.sign(a2[m]) == np.sign(y[m])).mean() if m.sum() else np.nan
null1 = null1[~np.isnan(null1)]
p1 = (np.sum(null1 >= obs_acc) + 1) / (len(null1) + 1)
print(f'NULL 1 permuted seed : mean {null1.mean():.4f} sd {null1.std(ddof=1):.4f}  p = {p1:.4g}')

# --- null 2: rewired network (curveball within edge class), true seed
byc = collections.defaultdict(list)
for e in E:
    byc[e['edge_type']].append((idx[e['source']], idx[e['target']], float(e['sign_published'])))


def curveball(edges):
    adj = collections.defaultdict(set); sgn_ = {}
    for u, v, s in edges:
        adj[u].add(v); sgn_[(u, v)] = s
    keys = list(adj)
    if len(keys) < 2:
        return edges
    for _ in range(5 * len(keys)):
        a, b = rng.choice(len(keys), 2, replace=False)
        A_, B_ = keys[a], keys[b]
        sa, sb = adj[A_], adj[B_]
        com = sa & sb
        ea = list(sa - com); eb = list(sb - com)
        pool = ea + eb
        if not pool:
            continue
        rng.shuffle(pool)
        adj[A_] = com | set(pool[:len(ea)]); adj[B_] = com | set(pool[len(ea):])
    s0 = edges[0][2]
    return [(u, v, s0) for u in adj for v in adj[u]]


NR = 300
null2 = np.empty(NR)
for i in range(NR):
    ed = []
    for k, es in byc.items():
        ed += curveball(es)
    r2 = np.array([e[0] for e in ed]); c2 = np.array([e[1] for e in ed])
    v2 = np.array([e[2] for e in ed])
    k2 = np.bincount(r2, weights=np.abs(v2), minlength=N)
    w2 = v2 / np.where(k2[r2] > 0, k2[r2], 1.0)
    a2, _, _ = boolean_attractor(seed, r2, c2, w2, maxsteps=120)
    m = sig & (np.abs(a2) > 0)
    null2[i] = (np.sign(a2[m]) == np.sign(y[m])).mean() if m.sum() else np.nan
null2 = null2[~np.isnan(null2)]
p2 = (np.sum(null2 >= obs_acc) + 1) / (len(null2) + 1)
print(f'NULL 2 rewired net   : mean {null2.mean():.4f} sd {null2.std(ddof=1):.4f}  p = {p2:.4g}')

bt_obs = binomtest(int(round(obs_acc * obs_n)), obs_n, 0.5)
print(f'binomial vs 0.5: p = {bt_obs.pvalue:.3g}')

# --- best clamp with a usable readout size
cl = list(csv.DictReader(open(f'{RES}/systems_dynamics_boolean_clamp.csv')))
cl = [r for r in cl if int(r['n_nodes_reached']) >= 50]
cl.sort(key=lambda r: -float(r['sign_accuracy']))
best = cl[0]
nb = int(best['n_nodes_reached']); ab = float(best['sign_accuracy'])
btb = binomtest(int(round(ab * nb)), nb, 0.5)
print(f"best clamp with n>=50: {best['node']} clamp={best['clamp']} acc={ab:.4f} n={nb} "
      f"binom p={btb.pvalue:.3g}")

# --- self-consistency of the linear-response drivers
out = {}
for tag in ('published', 'curated'):
    c = list(csv.DictReader(open(f'{RES}/systems_dynamics_single_perturbations_{tag}.csv')))
    top = c[:20]
    cons = []
    for r in top:
        try:
            fc = float(r['own_logFC'])
        except ValueError:
            continue
        cons.append(int(r['direction']) * np.sign(fc) > 0)
    out[tag] = dict(n_top_with_logFC=len(cons), n_direction_consistent=int(sum(cons)),
                    frac_consistent=float(np.mean(cons)) if cons else float('nan'),
                    best=top[0]['node'], best_direction=int(top[0]['direction']),
                    best_rho=float(top[0]['spearman']), best_own_logFC=top[0]['own_logFC'])
    print(f"[{tag}] top-20 drivers whose perturbation direction matches their own observed "
          f"logFC: {sum(cons)}/{len(cons)}")

# --- focused readout: reproducing the observed COL1A1/COL3A1 change
R = np.load(f'{RES}/systems_response_matrix_published.npy')
tg = [idx['COL1A1'], idx['COL3A1']]
obs_col = y[tg]
frows = []
for s in range(N):
    if s in tg:
        continue
    r = R[tg, s]
    if np.all(np.abs(r) < 1e-14):
        continue
    for sigma in (1, -1):
        rr = sigma * r
        frows.append(dict(node=NODES[s], type=nt[NODES[s]], direction=sigma,
                          response_COL1A1=float(rr[0]), response_COL3A1=float(rr[1]),
                          observed_COL1A1=float(obs_col[0]), observed_COL3A1=float(obs_col[1]),
                          both_signs_match=bool(np.all(np.sign(rr) == np.sign(obs_col))),
                          magnitude=float(np.abs(rr).sum())))
frows = [f for f in frows if f['both_signs_match']]
frows.sort(key=lambda d: -d['magnitude'])
with open(f'{RES}/systems_dynamics_collagen_drivers.csv', 'w', newline='') as fh:
    wtr = csv.DictWriter(fh, fieldnames=list(frows[0].keys())); wtr.writeheader(); wtr.writerows(frows)
print(f'\nperturbations reproducing the OBSERVED sign of both collagens '
      f'(COL1A1 {obs_col[0]:+.3f}, COL3A1 {obs_col[1]:+.3f}), ranked by magnitude:')
for f in frows[:12]:
    print(f"  {f['node']:<14} {f['type']:<6} dir={f['direction']:+d} "
          f"COL1A1={f['response_COL1A1']:+.5f} COL3A1={f['response_COL3A1']:+.5f}")

json.dump(dict(boolean=dict(observed_agreement=obs_acc, n=obs_n, steps=int(steps),
                            period=int(period),
                            null_permuted_seed_mean=float(null1.mean()),
                            null_permuted_seed_sd=float(null1.std(ddof=1)),
                            p_permuted_seed=float(p1), n_perm=int(len(null1)),
                            null_rewired_mean=float(null2.mean()),
                            null_rewired_sd=float(null2.std(ddof=1)),
                            p_rewired=float(p2), n_rewire=int(len(null2)),
                            binomial_p=float(bt_obs.pvalue),
                            best_clamp_n50=dict(node=best['node'], clamp=int(best['clamp']),
                                                accuracy=ab, n=nb, binom_p=float(btb.pvalue))),
               direction_consistency=out,
               collagen_top=[f['node'] for f in frows[:10]]),
          open(f'{RES}/systems_dynamics_nulls.json', 'w'), indent=2)
print('\nwrote systems_dynamics_nulls.json')
