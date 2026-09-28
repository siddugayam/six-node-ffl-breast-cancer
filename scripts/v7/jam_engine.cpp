// jam_engine.cpp -- simulated-annealing active-subnetwork search
// Faithful re-implementation of the jActiveModules search of
// Ideker, Ozier, Schwikowski & Siegel (2002) Bioinformatics 18:S233-S240.
//   * state            = active/inactive assignment over all network nodes
//   * subnetwork score = z_A = (1/sqrt(k)) * sum_{i in A} z_i
//   * corrected score  = s_A = (z_A - mu_k) / sigma_k, mu_k/sigma_k from
//                        randomly sampled node sets of the same size k
//   * state score      = corrected score of the HIGHEST-SCORING CONNECTED
//                        COMPONENT of the subgraph induced by the active nodes
//   * move             = toggle one uniformly chosen node
//   * acceptance       = accept improvement; else accept with prob exp(delta/T)
//   * schedule         = geometric cooling from T0 to T1 over `iters` steps
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

struct Graph {
  const int *flat; const int *start; const int *len; int n;
};

// score of a state: max corrected score over connected components of active set
static double state_score(const Graph &G, const std::vector<char> &st,
                          const double *z, const double *mu, const double *sd,
                          std::vector<int> &comp, std::vector<int> &stack,
                          int *best_comp_id) {
  int n = G.n;
  std::fill(comp.begin(), comp.end(), -1);
  double best = R_NegInf; int bid = -1; int cid = 0;
  for (int s = 0; s < n; ++s) {
    if (!st[s] || comp[s] != -1) continue;
    // BFS/DFS from s over active nodes
    int top = 0; stack[top++] = s; comp[s] = cid;
    double sumz = 0.0; int k = 0;
    while (top > 0) {
      int v = stack[--top]; sumz += z[v]; ++k;
      int b = G.start[v], e = b + G.len[v];
      for (int j = b; j < e; ++j) {
        int u = G.flat[j];
        if (st[u] && comp[u] == -1) { comp[u] = cid; stack[top++] = u; }
      }
    }
    double zagg = sumz / std::sqrt((double)k);
    double sc = (zagg - mu[k]) / sd[k];
    if (sc > best) { best = sc; bid = cid; }
    ++cid;
  }
  if (best_comp_id) *best_comp_id = bid;
  return best;                       // -Inf if no node is active
}

// [[Rcpp::export]]
List sa_search_cpp(IntegerVector adj_flat, IntegerVector adj_start, IntegerVector adj_len,
                   NumericVector z, NumericVector mu_k, NumericVector sd_k,
                   LogicalVector init_state, LogicalVector allowed,
                   int iters, double T0, double T1, int trace_every) {
  int n = z.size();
  Graph G; G.flat = adj_flat.begin(); G.start = adj_start.begin();
  G.len = adj_len.begin(); G.n = n;
  const double *zp = z.begin(); const double *mu = mu_k.begin(); const double *sd = sd_k.begin();

  std::vector<char> st(n);
  for (int i = 0; i < n; ++i) st[i] = (init_state[i] && allowed[i]) ? 1 : 0;
  std::vector<int> comp(n), stack(n);

  double cur = state_score(G, st, zp, mu, sd, comp, stack, NULL);
  if (!R_finite(cur)) cur = 0.0;
  std::vector<char> best_state(st);
  double best = cur;

  std::vector<int> movable;
  for (int i = 0; i < n; ++i) if (allowed[i]) movable.push_back(i);
  int nm = movable.size();

  int ntr = trace_every > 0 ? (iters / trace_every + 1) : 0;
  NumericVector trace(ntr), temps(ntr);
  int ti = 0;
  double ratio = (iters > 1) ? std::pow(T1 / T0, 1.0 / (double)(iters - 1)) : 1.0;
  double T = T0;
  long naccept = 0, nworse = 0;

  for (int it = 0; it < iters; ++it) {
    int v = movable[(int)(unif_rand() * nm)];
    st[v] = !st[v];
    double cand = state_score(G, st, zp, mu, sd, comp, stack, NULL);
    if (!R_finite(cand)) cand = R_NegInf;
    double d = cand - cur;
    bool acc;
    if (d >= 0) acc = true;
    else { nworse++; acc = (unif_rand() < std::exp(d / T)); }
    if (acc) { cur = cand; naccept++; if (cur > best) { best = cur; best_state = st; } }
    else st[v] = !st[v];
    if (trace_every > 0 && it % trace_every == 0 && ti < ntr) {
      trace[ti] = cur; temps[ti] = T; ti++;
    }
    T *= ratio;
  }
  LogicalVector fs(n), bs(n);
  for (int i = 0; i < n; ++i) { fs[i] = (bool)st[i]; bs[i] = (bool)best_state[i]; }
  return List::create(_["final_state"] = fs, _["best_state"] = bs,
                      _["best_score"] = best, _["final_score"] = cur,
                      _["trace"] = trace, _["temps"] = temps,
                      _["n_accept"] = (double)naccept, _["n_worse"] = (double)nworse);
}

// Null model 2: randomly sampled CONNECTED subnetworks.
// Each replicate grows a connected subgraph from a random seed node, one node
// at a time, recording the aggregate z after every addition.  Returns a
// B x nmax matrix of aggregate z scores (row = replicate, column = size k).
// [[Rcpp::export]]
NumericMatrix grow_null_cpp(IntegerVector adj_flat, IntegerVector adj_start,
                            IntegerVector adj_len, NumericVector z, int B, int nmax) {
  int n = z.size();
  NumericMatrix out(B, nmax);
  std::vector<char> inset(n);
  std::vector<int> cand;
  for (int b = 0; b < B; ++b) {
    std::fill(inset.begin(), inset.end(), 0);
    cand.clear();
    int s = (int)(unif_rand() * n);
    inset[s] = 1; double sumz = z[s]; int k = 1;
    out(b, 0) = sumz;
    for (int j = adj_start[s]; j < adj_start[s] + adj_len[s]; ++j) cand.push_back(adj_flat[j]);
    while (k < nmax) {
      if (cand.empty()) {                       // component exhausted: restart elsewhere
        int t = -1; int guard = 0;
        while (guard++ < 10 * n) { int q = (int)(unif_rand() * n); if (!inset[q]) { t = q; break; } }
        if (t < 0) break;
        inset[t] = 1; sumz += z[t]; ++k; out(b, k - 1) = sumz;
        for (int j = adj_start[t]; j < adj_start[t] + adj_len[t]; ++j)
          if (!inset[adj_flat[j]]) cand.push_back(adj_flat[j]);
        continue;
      }
      int idx = (int)(unif_rand() * cand.size());
      int v = cand[idx];
      cand[idx] = cand.back(); cand.pop_back();
      if (inset[v]) continue;
      inset[v] = 1; sumz += z[v]; ++k; out(b, k - 1) = sumz;
      for (int j = adj_start[v]; j < adj_start[v] + adj_len[v]; ++j)
        if (!inset[adj_flat[j]]) cand.push_back(adj_flat[j]);
    }
  }
  // convert running sums to aggregate z = sum / sqrt(k)
  for (int b = 0; b < B; ++b)
    for (int k = 1; k <= nmax; ++k) out(b, k - 1) /= std::sqrt((double)k);
  return out;
}
