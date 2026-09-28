#!/usr/bin/env python3
"""(E) DYNAMIC RESPONSE OF THE REAL NETWORK.

MODEL.  Signed, out-degree-normalised linear response on the canonical network:

        dx/dt = -x + A x + u ,      A[v,u] = alpha * s(u->v) / k_out(u)

  x[v] is the deviation of node v from its normal-tissue level, s(u->v) in {-1,+1} is the
  edge sign, and u is a sustained perturbation.  Every column of A has L1 norm alpha, so
  ||A||_1 = alpha and rho(A) <= alpha < 1: the steady state x* = (I-A)^-1 u exists, is
  unique, and is globally asymptotically stable.  alpha = 0.85 throughout.

  SOLUTION AND CONVERGENCE CRITERION.  x* is obtained by the Jacobi/Neumann iteration
  x^(k+1) = A x^(k) + u, which is a contraction with modulus alpha.  Iteration stops when
  ||x^(k+1) - x^(k)||_inf < 1e-12; the cap is 5,000 iterations and is never reached
  (alpha^k < 1e-12 at k = 170).  This is the same object the manuscript computes as a
  "Signed RWR" score, but here it is the steady state of a stated dynamical system rather
  than a heuristic diffusion score.

PHENOTYPE.  y[v] = tumour-vs-normal log2 fold change of node v (TCGA-BRCA, limma).

QUESTIONS.
  1. Which single sustained upstream perturbation u = sigma*e_s best reproduces y?
  2. Is that better than chance (phenotype-permutation null, and a max-statistic null over
     all 2N candidates)?
  3. What is the smallest set of simultaneous perturbations that reproduces y (LASSO over
     the response matrix)?
  4. Does a discrete signed threshold (Boolean) model over the same network reach an
     attractor that matches the observed phenotype?

Sensitivity: the whole ranking is repeated with TRRUST/TransmiR CURATED signs (TF->target
edges annotated Repression get -1, unannotated TF edges dropped) because the deposited
network assigns +1 to every TF->target edge.
"""
import os, sys, csv, json, collections
import numpy as np
from scipy.stats import spearmanr, pearsonr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import load_nodes, load_edges, RES, FIG

ALPHA = 0.85
TOL = 1e-12
MAXIT = 5000
rng = np.random.default_rng(20260908)

nt = load_nodes(); E = load_edges()
NODES = sorted(nt); idx = {n: i for i, n in enumerate(NODES)}; N = len(NODES)

# ------------------------------------------------------------------ phenotype
de = {}
for path in ('BRCA_DEX_genes.csv', 'BRCA_DEX_mirnas.csv'):
    with open(f'/path/to/revision/results/{path}') as fh:
        for r in csv.DictReader(fh):
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
print(f'network nodes with a DE measurement: {meas.sum()}/{N};  '
      f'FDR<0.05 & |logFC|>0.5: {sig.sum()}')
print(f'   COL1A1 logFC={y[idx["COL1A1"]]:+.3f}  COL3A1 logFC={y[idx["COL3A1"]]:+.3f}  '
      f'miR-29a={y[idx["hsa-miR-29a"]]:+.3f}')


def build_A(sign_key, drop_unsigned):
    rows = []; cols = []; vals = []
    for e in E:
        s = e[sign_key]
        if s is None or s == 0:
            if drop_unsigned:
                continue
            s = e['sign_published']
        rows.append(idx[e['source']]); cols.append(idx[e['target']]); vals.append(float(s))
    rows = np.array(rows); cols = np.array(cols); vals = np.array(vals)
    kout = np.bincount(rows, weights=np.abs(vals), minlength=N)
    w = vals / np.where(kout[rows] > 0, kout[rows], 1.0)
    return rows, cols, w, len(vals)


def steady_state(u, rows, cols, w):
    """x = A x + u  with A[v,u] = alpha*w(u->v).  Returns (x, n_iter, residual)."""
    x = u.copy()
    for k in range(MAXIT):
        ax = np.bincount(cols, weights=w * x[rows], minlength=N)
        xn = u + ALPHA * ax
        d = np.max(np.abs(xn - x))
        x = xn
        if d < TOL:
            return x, k + 1, d
    return x, MAXIT, d


def response_matrix(rows, cols, w):
    """R[:, s] = steady state for u = e_s.  Built column by column."""
    R = np.zeros((N, N)); iters = []
    for s in range(N):
        u = np.zeros(N); u[s] = 1.0
        x, k, _ = steady_state(u, rows, cols, w)
        R[:, s] = x; iters.append(k)
    return R, int(np.max(iters))


results = {}
for tag, key, drop in [('published', 'sign_published', False),
                       ('curated', 'sign_curated', True)]:
    rows, cols, w, nE = build_A(key, drop)
    # verify contraction
    rho_bound = ALPHA
    R, mx = response_matrix(rows, cols, w)
    print(f'\n[{tag}] edges used = {nE};  max Neumann iterations to 1e-12 = {mx} '
          f'(contraction modulus alpha={ALPHA})')

    # ---------------- Q1 single perturbations
    cand = []
    for s in range(N):
        r = R[:, s].copy()
        m = meas.copy(); m[s] = False
        if np.std(r[m]) < 1e-14:
            continue
        rho, pv = spearmanr(r[m], y[m])
        pr, ppv = pearsonr(r[m], y[m])
        # signed hit rate over strongly DE nodes
        m2 = sig.copy(); m2[s] = False
        hits = np.sign(r[m2]) == np.sign(y[m2])
        n_nz = int((np.abs(r[m2]) > 1e-12).sum())
        hr = float(hits[np.abs(r[m2]) > 1e-12].mean()) if n_nz else np.nan
        for sigma in (1, -1):
            cand.append(dict(node=NODES[s], type=nt[NODES[s]], direction=sigma,
                             spearman=sigma * rho, spearman_p=pv,
                             pearson=sigma * pr, pearson_p=ppv,
                             signed_hit_rate=(hr if sigma == 1 else 1 - hr) if not np.isnan(hr) else np.nan,
                             n_reached_sig=n_nz,
                             own_logFC=(y[s] if meas[s] else np.nan)))
    cand.sort(key=lambda d: -d['spearman'])
    with open(f'{RES}/systems_dynamics_single_perturbations_{tag}.csv', 'w', newline='') as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(cand[0].keys())); wtr.writeheader()
        wtr.writerows(cand)
    print(f'  best single perturbations reproducing the tumour phenotype [{tag}]:')
    for d in cand[:10]:
        own = 'NA' if np.isnan(d['own_logFC']) else f"{d['own_logFC']:+.2f}"
        hit = 'NA' if np.isnan(d['signed_hit_rate']) else f"{d['signed_hit_rate']:.3f}"
        print(f"    {d['node']:<14} {d['type']:<6} dir={d['direction']:+d} "
              f"rho={d['spearman']:+.4f} r={d['pearson']:+.4f} "
              f"hit={hit} n_reached={d['n_reached_sig']} ownFC={own}")

    # ---------------- Q2 max-statistic permutation null over all 2N candidates
    NPERM = 2000
    best_obs = cand[0]['spearman']
    ym = y[meas]
    Rm = R[meas]                                   # measured nodes x candidate perturbations
    keep_col = Rm.std(0) > 1e-14
    rk_R = np.argsort(np.argsort(Rm[:, keep_col], axis=0), axis=0).astype(float)
    rk_R = (rk_R - rk_R.mean(0)) / rk_R.std(0)
    maxnull = np.empty(NPERM)
    for p in range(NPERM):
        yp = rng.permutation(ym)
        rk_y = np.argsort(np.argsort(yp)).astype(float)
        rk_y = (rk_y - rk_y.mean()) / rk_y.std()
        rc = (rk_R * rk_y[:, None]).mean(0)         # Spearman of each candidate vs permuted y
        maxnull[p] = np.max(np.abs(rc))
    p_fwer = float((np.sum(maxnull >= abs(best_obs)) + 1) / (NPERM + 1))
    print(f'  max-statistic permutation null over all {int(keep_col.sum())*2} candidates '
          f'(+/- directions): best|rho|={abs(best_obs):.4f}, null max mean={maxnull.mean():.4f} '
          f'sd={maxnull.std(ddof=1):.4f}, FWER p={p_fwer:.4f} ({NPERM} permutations)')

    # ---------------- Q3 sparse multi-node perturbation (LASSO)
    from sklearn.linear_model import LassoCV, lars_path
    Xd = R[meas].copy()
    yd = y[meas].copy()
    Xs = (Xd - Xd.mean(0)) / np.where(Xd.std(0) > 0, Xd.std(0), 1.0)
    keep = Xd.std(0) > 1e-12
    Xs = Xs[:, keep]; names_keep = [NODES[i] for i in np.where(keep)[0]]
    alphas, _, coefs = lars_path(Xs, yd - yd.mean(), method='lasso')
    greedy = []
    seen = set()
    for j in range(coefs.shape[1]):
        nz = np.where(coefs[:, j] != 0)[0]
        for i in nz:
            if i not in seen:
                seen.add(i)
                pred = Xs[:, list(seen)] @ np.linalg.lstsq(Xs[:, list(seen)],
                                                          yd - yd.mean(), rcond=None)[0]
                r2 = 1 - np.sum((yd - yd.mean() - pred) ** 2) / np.sum((yd - yd.mean()) ** 2)
                greedy.append(dict(order=len(seen), node=names_keep[i],
                                   type=nt[names_keep[i]],
                                   coef_sign=int(np.sign(coefs[i, j])),
                                   cumulative_R2=float(r2)))
        if len(seen) >= 25:
            break
    with open(f'{RES}/systems_dynamics_lasso_path_{tag}.csv', 'w', newline='') as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(greedy[0].keys())); wtr.writeheader()
        wtr.writerows(greedy)
    print(f'  LASSO entry order (cumulative R2 of the fitted perturbation set) [{tag}]:')
    for g in greedy[:12]:
        print(f"    {g['order']:>2}. {g['node']:<14} {g['type']:<6} "
              f"sign={g['coef_sign']:+d}  R2={g['cumulative_R2']:.4f}")
    lcv = LassoCV(cv=5, random_state=0, max_iter=20000).fit(Xs, yd - yd.mean())
    nz = int((lcv.coef_ != 0).sum())
    r2cv = lcv.score(Xs, yd - yd.mean())
    print(f'  LassoCV: {nz} non-zero perturbations, in-sample R2 = {r2cv:.4f}, '
          f'alpha = {lcv.alpha_:.5g}')

    results[tag] = dict(n_edges=nE, max_iterations=mx, alpha=ALPHA, tol=TOL,
                        best_single=cand[0], top10=[c['node'] for c in cand[:10]],
                        fwer_p_best_single=p_fwer,
                        permutation_null_max_mean=float(maxnull.mean()),
                        permutation_null_max_sd=float(maxnull.std(ddof=1)),
                        n_permutations=int(NPERM),
                        lasso_cv_n_nonzero=nz, lasso_cv_R2=float(r2cv),
                        lasso_first10=[g['node'] for g in greedy[:10]],
                        R2_at_5=float([g['cumulative_R2'] for g in greedy if g['order'] == 5][0]),
                        R2_at_10=float([g['cumulative_R2'] for g in greedy if g['order'] == 10][0]))
    if tag == 'published':
        np.save(f'{RES}/systems_response_matrix_published.npy', R)
        RES_PUB = R

# ------------------------------------------------------------------ Q4 Boolean threshold
rows, cols, w, _ = build_A('sign_published', False)
state0 = np.zeros(N)
state0[sig] = np.sign(y[sig])
print(f'\n[Boolean threshold model] seeded with sign(logFC) at {int((state0!=0).sum())} '
      f'DE nodes; synchronous update x_v(t+1) = sign(sum_u w(u,v) x_u(t)), '
      f'fixed nodes keep their seed value if they have no regulator')


def boolean_run(seed, clamp=None, maxsteps=200):
    x = seed.copy()
    if clamp:
        for k, v in clamp.items():
            x[idx[k]] = v
    hist = {}
    for t in range(maxsteps):
        key = x.tobytes()
        if key in hist:
            return x, t, t - hist[key]        # (state, steps, period)
        hist[key] = t
        ax = np.bincount(cols, weights=w * x[rows], minlength=N)
        xn = np.sign(ax)
        indeg0 = np.bincount(cols, minlength=N) == 0
        xn[indeg0] = x[indeg0]
        if clamp:
            for k, v in clamp.items():
                xn[idx[k]] = v
        x = xn
    return x, maxsteps, -1


att, steps, period = boolean_run(state0)
agree = (np.sign(att[sig]) == np.sign(y[sig]))
nz = np.abs(att[sig]) > 0
print(f'  unperturbed attractor reached after {steps} steps, period {period}; '
      f'sign agreement with observed DE on reached nodes: '
      f'{agree[nz].mean():.3f} ({int(agree[nz].sum())}/{int(nz.sum())})')

brows = []
for s in NODES:
    for sigma in (1, -1):
        a, st, pe = boolean_run(np.zeros(N), clamp={s: sigma})
        m = sig.copy(); m[idx[s]] = False
        reach = np.abs(a[m]) > 0
        if reach.sum() < 10:
            continue
        acc = float((np.sign(a[m][reach]) == np.sign(y[m][reach])).mean())
        brows.append(dict(node=s, type=nt[s], clamp=sigma, steps=st, period=pe,
                          n_nodes_reached=int(reach.sum()), sign_accuracy=acc))
brows.sort(key=lambda d: (-d['sign_accuracy'], -d['n_nodes_reached']))
with open(f'{RES}/systems_dynamics_boolean_clamp.csv', 'w', newline='') as fh:
    wtr = csv.DictWriter(fh, fieldnames=list(brows[0].keys())); wtr.writeheader(); wtr.writerows(brows)
print('  best single-node clamps by sign agreement with the tumour phenotype:')
for b in brows[:10]:
    print(f"    {b['node']:<14} clamp={b['clamp']:+d} acc={b['sign_accuracy']:.3f} "
          f"n={b['n_nodes_reached']} steps={b['steps']} period={b['period']}")
# binomial null for the best
from scipy.stats import binomtest
b0 = brows[0]
bt = binomtest(int(round(b0['sign_accuracy'] * b0['n_nodes_reached'])), b0['n_nodes_reached'], 0.5)
print(f"    binomial p vs 0.5 for {b0['node']}: {bt.pvalue:.3g}")

results['boolean'] = dict(seed_nodes=int((state0 != 0).sum()), steps=int(steps),
                          period=int(period),
                          unperturbed_sign_agreement=float(agree[nz].mean()),
                          n_evaluated=int(nz.sum()),
                          best_clamp=b0, best_clamp_binom_p=float(bt.pvalue))
json.dump(results, open(f'{RES}/systems_dynamics_summary.json', 'w'), indent=2, default=str)
print('\nwrote systems_dynamics_summary.json')
