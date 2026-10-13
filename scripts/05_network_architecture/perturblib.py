#!/usr/bin/env python3
"""Fast knockout-damage kernel shared by the single- and double-knockout analyses.

Three damage components, each expressed as a fraction of the wild-type value computed on the
SAME surviving node universe (so a score reflects loss of routing, not loss of the node):

  dFFL   3-node FFL cores destroyed / 1,649
  dReach ordered reachable pairs lost / reachable pairs among survivors in the WT graph
  dCOL   signed change in the total absolute influence arriving at COL1A1 + COL3A1.
         Influence column x_c solves  x = e_c + (1-lambda) W x  with W[u,v] = sign(u,v)/outdeg(u).
         Solved by Neumann iteration; the iteration matrix has spectral radius <= 1-lambda = 0.85
         (measured 0.468 for this network) so the fixed point is unique and the iteration is a
         contraction.  CONVERGENCE CRITERION: ||x^(k+1) - x^(k)||_inf < 1e-12, cap 2,000 iterations.

composite_damage = mean(dFFL, dReach, |dCOL|)
"""
import collections
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

TOL = 1e-12
MAXIT = 2000
LAM = 0.15


class Kernel:
    def __init__(self, nodes, pairs, signs, node_type, cores, targets=('COL1A1', 'COL3A1'),
                 lam=LAM):
        self.nodes = list(nodes)
        self.idx = {n: i for i, n in enumerate(self.nodes)}
        self.N = len(self.nodes)
        self.lam = lam
        self.targets = [t for t in targets if t in self.idx]
        self.pairs = list(pairs)
        rows = np.array([self.idx[u] for u, v in self.pairs])
        cols = np.array([self.idx[v] for u, v in self.pairs])
        self.rows, self.cols = rows, cols
        self.sgn = np.array([signs[(u, v)] for u, v in self.pairs], float)
        self.node_type = node_type
        # FFL cores as index tuples
        self.cores = cores
        self.core_nodes = [tuple(self.idx[x] for x in c[:3]) for c in cores]
        self.core_is_comp = np.array([c[3] == 'Composite-FFL' for c in cores])
        self.core_by_node = collections.defaultdict(list)
        for k, cn in enumerate(self.core_nodes):
            for x in set(cn):
                self.core_by_node[x].append(k)
        self.FFL_WT = self._ffl_unique(np.ones(len(cores), bool))
        self.alive_all = np.ones(self.N, bool)
        self.R_WT_all = self.reach_pairs(self.alive_all, self.alive_all)
        self.INF_WT_all = self.influence(self.alive_all)

    # ---------------------------------------------------------------- FFL
    def _ffl_unique(self, keep):
        nc = int((keep & self.core_is_comp).sum())
        no = int((keep & ~self.core_is_comp).sum())
        return nc // 2 + no

    def ffl_after(self, ko_idx):
        dead = set()
        for x in ko_idx:
            dead.update(self.core_by_node.get(x, ()))
        keep = np.ones(len(self.core_nodes), bool)
        if dead:
            keep[list(dead)] = False
        return self._ffl_unique(keep)

    # ---------------------------------------------------------------- reachability
    def _reach_matrix(self, alive):
        """Boolean N x N reachability (path length >= 1) restricted to alive nodes."""
        m = alive[self.rows] & alive[self.cols]
        r, c = self.rows[m], self.cols[m]
        A = csr_matrix((np.ones(len(r), np.int8), (r, c)), shape=(self.N, self.N))
        ncomp, lab = connected_components(A, directed=True, connection='strong')
        # condensation DAG
        cr = lab[r]; cc = lab[c]
        keep = cr != cc
        # component -> node bitset as boolean rows
        Cnode = np.zeros((ncomp, self.N), bool)
        alive_idx = np.where(alive)[0]
        Cnode[lab[alive_idx], alive_idx] = True
        # topological order of condensation: scipy labels are in reverse topological order
        # for connection='strong' -- verify by construction instead of trusting it
        succ = collections.defaultdict(set)
        indeg = collections.Counter()
        for a, b in zip(cr[keep], cc[keep]):
            if b not in succ[a]:
                succ[a].add(b); indeg[b] += 1
        order = [c0 for c0 in range(ncomp) if indeg[c0] == 0]
        topo = []
        q = list(order)
        indeg2 = dict(indeg)
        while q:
            u = q.pop()
            topo.append(u)
            for v in succ.get(u, ()):
                indeg2[v] -= 1
                if indeg2[v] == 0:
                    q.append(v)
        Creach = np.zeros((ncomp, self.N), bool)
        selfloop = set(zip(cr[~keep], cc[~keep]))
        sizes = Cnode.sum(1)
        for u in reversed(topo):
            acc = np.zeros(self.N, bool)
            for v in succ.get(u, ()):
                acc |= Creach[v] | Cnode[v]
            if sizes[u] > 1 or (u, u) in selfloop:
                acc |= Cnode[u]
            Creach[u] = acc
        R = np.zeros((self.N, self.N), bool)
        R[alive_idx] = Creach[lab[alive_idx]]
        return R

    def reach_pairs(self, alive, restrict):
        R = self._reach_matrix(alive)
        sub = R[np.ix_(restrict, restrict)]
        np.fill_diagonal(sub, False)
        return int(sub.sum())

    # ---------------------------------------------------------------- influence
    def influence(self, alive, return_iters=False):
        m = alive[self.rows] & alive[self.cols]
        r, c, s = self.rows[m], self.cols[m], self.sgn[m]
        od = np.bincount(r, weights=np.abs(s), minlength=self.N)
        w = s / np.where(od[r] > 0, od[r], 1.0)
        # x = e_c + (1-lam) W x   -->   (W x)[u] = sum_v w[u,v] x[v]
        tot = 0.0; iters = 0
        for t in self.targets:
            if not alive[self.idx[t]]:
                continue
            e = np.zeros(self.N); e[self.idx[t]] = 1.0
            x = e.copy()
            for k in range(MAXIT):
                wx = np.bincount(r, weights=w * x[c], minlength=self.N)
                xn = e + (1 - self.lam) * wx
                xn[~alive] = 0.0
                d = np.max(np.abs(xn - x))
                x = xn
                if d < TOL:
                    break
            iters = max(iters, k + 1)
            x[self.idx[t]] = 0.0
            tot += float(np.abs(x[alive]).sum())
        return (tot, iters) if return_iters else tot

    # ---------------------------------------------------------------- damage
    def damage(self, ko):
        ki = [self.idx[n] for n in ko]
        alive = self.alive_all.copy()
        alive[ki] = False
        f = self.ffl_after(ki)
        r_wt = self.reach_pairs(self.alive_all, alive)
        r_ko = self.reach_pairs(alive, alive)
        inf_ko = self.influence(alive)
        if set(ko) & set(self.targets):
            dC = float('nan')
        else:
            dC = (self.INF_WT_all - inf_ko) / self.INF_WT_all
        dF = (self.FFL_WT - f) / self.FFL_WT
        dR = (r_wt - r_ko) / r_wt if r_wt else 0.0
        comp = float(np.mean([dF, dR, 0.0 if np.isnan(dC) else abs(dC)]))
        return dict(dFFL=dF, dReach=dR, dCOL=dC, composite=comp,
                    ffl_after=f, reach_after=r_ko, reach_wt=r_wt, influence_after=inf_ko)
