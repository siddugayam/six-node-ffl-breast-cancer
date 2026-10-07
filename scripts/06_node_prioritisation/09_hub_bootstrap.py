#!/usr/bin/env python3
"""
09 -- Bootstrap stability of the hub list.

Bootstrap resampling (n = 1,000) of the network edges, to assess the stability of the
hub list.

PROCEDURE (non-parametric edge bootstrap).
  For b = 1..1000:
     draw |E| = 6859 edges WITH REPLACEMENT from the observed edge list, collapse duplicates
     (this retains ~63.2% of the distinct edges, the standard bootstrap fraction), rebuild the
     directed graph on the full 587-node vertex set, and recompute every centrality.
  For each node type (TF / miRNA / Gene) and each metric, a node is "in the top 20" if its
  value is >= the 20th largest value among nodes of that type (ties admitted, so the top-20
  set can be slightly larger than 20 when values are tied; the realised size is reported).
  Stability = % of the 1000 bootstrap replicates in which the node is in the top 20.

METRICS
  degree            total degree (in + out) -- the "number of interactions" by which the
                    hubs are ranked (VEGFA 111, TP53 111, MYC 104, ...)
  betweenness       directed betweenness centrality
  ffl_participation number of 3-node FFL cores (miRNA-FFL + TF-FFL + Composite-FFL + TF-TF-FFL)
                    the node takes part in, in ANY role (regulator R, mediator M or target T)
                    -- i.e. how central the node is to the object the paper is actually about
"""
import csv, collections, os, random, time
import numpy as np
import igraph as ig

REV = '/path/to/revision'
LOG = f'{REV}/logs/hub_bootstrap.log'
N_BOOT = int(os.environ.get('N_BOOT', 1000))
SEED = 20260908
TOPN = 20

def say(*a):
    m = time.strftime('%H:%M:%S') + ' | ' + ''.join(str(x) for x in a)
    print(m, flush=True)
    with open(LOG, 'a') as fh: fh.write(m + '\n')

nodes = {r['name']: r['type'] for r in
         csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'), delimiter='\t')}
names = sorted(nodes); idx = {v: i for i, v in enumerate(names)}; n = len(names)
ntype = [nodes[v] for v in names]
E = []
for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'), delimiter='\t'):
    E.append((idx[r['source']], idx[r['target']]))
m = len(E)
NONMIR = np.array([t != 'miRNA' for t in ntype])

def ffl_participation(edges):
    outb = [0] * n
    for a, b in edges: outb[a] |= (1 << b)
    nm = 0
    for i in range(n):
        if ntype[i] != 'miRNA': nm |= (1 << i)
    eset = set(a * n + b for a, b in edges)
    part = np.zeros(n, dtype=np.int64)
    for a, b in edges:
        ta, tb = ntype[a], ntype[b]
        if tb == 'miRNA' and ta == 'miRNA': continue
        ok = ((ta == 'TF' and tb == 'miRNA') or (ta == 'miRNA' and tb == 'TF')
              or (ta == 'TF' and tb == 'TF'))
        if not ok: continue
        common = outb[a] & outb[b] & nm
        k = common.bit_count()
        if k:
            part[a] += k; part[b] += k
            t = common
            while t:                      # credit every target node of the core as well
                low = t & -t
                part[low.bit_length() - 1] += 1
                t ^= low
    return part

def metrics(edges):
    g = ig.Graph(n=n, edges=edges, directed=True)
    deg = np.array(g.degree(mode='all'), dtype=np.float64)
    btw = np.array(g.betweenness(directed=True), dtype=np.float64)
    fp = ffl_participation(edges).astype(np.float64)
    return {'degree': deg, 'betweenness': btw, 'ffl_participation': fp}

def top_sets(vals):
    """per node type -> (threshold, boolean mask of nodes in the top-20 of that type)"""
    out = {}
    for t in ('TF', 'miRNA', 'Gene'):
        sel = np.array([x == t for x in ntype])
        v = vals[sel]
        if len(v) == 0: continue
        thr = np.sort(v)[::-1][min(TOPN, len(v)) - 1]
        mask = sel & (vals >= thr) & (vals > 0)
        out[t] = (thr, mask)
    return out

def main():
    open(LOG, 'w').close()
    say('=== 09_hub_bootstrap :: START ===')
    say(f'network {n} nodes, {m} directed edges; N_BOOT={N_BOOT}; TOPN={TOPN}; SEED={SEED}')
    real = metrics(E)
    real_top = {k: top_sets(v) for k, v in real.items()}
    for k in real:
        for t, (thr, mask) in real_top[k].items():
            say(f'  real top-{TOPN} {k:<18} {t:<6}: threshold={thr:g}, set size={int(mask.sum())}')
    say('  real degree leaders: ' + ', '.join(
        f'{names[i]}={int(real["degree"][i])}' for i in np.argsort(-real['degree'])[:8]))

    rng = random.Random(SEED)
    hits = {k: np.zeros(n, dtype=np.int64) for k in real}
    ranksum = {k: np.zeros(n) for k in real}
    ranksq = {k: np.zeros(n) for k in real}
    t0 = time.time()
    for b in range(N_BOOT):
        samp = [E[rng.randrange(m)] for _ in range(m)]
        eb = sorted(set(samp))
        mm = metrics(eb)
        for k, v in mm.items():
            for t, (thr, mask) in top_sets(v).items():
                hits[k][mask] += 1
            # within-type rank (1 = highest)
            for t in ('TF', 'miRNA', 'Gene'):
                sel = np.where(np.array([x == t for x in ntype]))[0]
                order = sel[np.argsort(-v[sel], kind='stable')]
                rk = np.empty(n); rk[:] = np.nan
                for pos, node_i in enumerate(order, 1): rk[node_i] = pos
                ranksum[k][sel] += rk[sel]; ranksq[k][sel] += rk[sel] ** 2
        if (b + 1) % 100 == 0:
            el = time.time() - t0
            say(f'   bootstrap {b+1}/{N_BOOT}  elapsed {el/60:.1f} min  '
                f'eta {el/(b+1)*(N_BOOT-b-1)/60:.1f} min  (|E_uniq|={len(eb)})')
    say(f'bootstrap done in {(time.time()-t0)/60:.1f} min')

    rows = []
    for k in ('degree', 'betweenness', 'ffl_participation'):
        for t in ('TF', 'miRNA', 'Gene'):
            thr, mask = real_top[k][t]
            sel = np.where(np.array([x == t for x in ntype]))[0]
            order = sel[np.argsort(-real[k][sel], kind='stable')]
            realrank = {int(node_i): pos for pos, node_i in enumerate(order, 1)}
            for i in np.where(mask)[0]:
                mr = ranksum[k][i] / N_BOOT
                sdr = max(0.0, ranksq[k][i] / N_BOOT - mr ** 2) ** 0.5
                rows.append(dict(
                    node=names[i], node_type=t, metric=k,
                    real_value=(round(float(real[k][i]), 4)),
                    real_rank_within_type=realrank[int(i)],
                    stability_pct=round(100.0 * hits[k][i] / N_BOOT, 1),
                    n_boot_in_top20=int(hits[k][i]), n_boot=N_BOOT,
                    mean_boot_rank=round(mr, 2), sd_boot_rank=round(sdr, 2),
                    top20_threshold=round(float(thr), 4),
                    top20_set_size=int(mask.sum())))
    rows.sort(key=lambda r: (r['metric'], r['node_type'], r['real_rank_within_type']))
    fields = ['node', 'node_type', 'metric', 'real_value', 'real_rank_within_type',
              'stability_pct', 'n_boot_in_top20', 'n_boot', 'mean_boot_rank', 'sd_boot_rank',
              'top20_threshold', 'top20_set_size']
    with open(f'{REV}/results/hub_bootstrap_stability.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
        for r in rows: w.writerow(r)
    say(f'wrote results/hub_bootstrap_stability.csv ({len(rows)} rows)')
    for k in ('degree', 'betweenness', 'ffl_participation'):
        for t in ('TF', 'miRNA', 'Gene'):
            sub = [r for r in rows if r['metric'] == k and r['node_type'] == t]
            st = [r['stability_pct'] for r in sub]
            if not st:
                say(f'  {k:<18} {t:<6} n=0  (no node of this type ever enters the top-20 for '
                    f'this metric)')
                continue
            say(f'  {k:<18} {t:<6} n={len(sub):<3} stability median={np.median(st):.1f}% '
                f'min={min(st):.1f}% max={max(st):.1f}%  '
                f'>=95%: {sum(1 for x in st if x >= 95)}  <50%: {sum(1 for x in st if x < 50)}')
    say('=== 09 DONE ===')

if __name__ == '__main__':
    main()
