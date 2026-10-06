#!/usr/bin/env python3
"""A4 of SETTINGS.md: KnockTF 2.0 support for the TF-TF arcs of the six-node Bhat patterns, with the rule of P3b
(analyses/perturbation_tests), unchanged.

Arcs: the TF -> TF arcs the six-node networks use as their compulsory T1-T2 link (T1 -> T2 in miRNA-FFL and composite
instances, T2 -> T1 in TF-FFL and composite instances). Sign = the network edge files' sign (TRRUST Activation +1,
Repression -1, else none).
Per KnockTF dataset of a source TF (>= 2 control and >= 2 treated samples; TF Log2FC <= -0.74; duplicates of the same TF,
biosample, sample counts and TF Log2FC kept once):
  |log2FC| test: median |Log2FC| of the TF's targets in the arc set minus that of expressed genes that are not TRRUST
  targets of the TF in any mode (expressed = Mean_Control >= the dataset median); two-sided Wilcoxon; bootstrap SE
  (1,000; seeded from 20250908 and the dataset key, as P3b); fewer than 3 targets -> no test;
  sign agreement: Activation arcs expected down, Repression arcs up (Log2FC = 0 counts as disagreement).
Pooled per arc set: REML for the |log2FC| difference; agreement = agreements / signed arcs over the tested datasets,
two-sided binomial test. SUPPORTED if the pooled difference > 0 with 95 % CI excluding 0 AND agreement > 50 % with
binomial p < 0.05; otherwise NOT SUPPORTED.
Also, per arc: the target's Log2FC in every included dataset of its source TF, or "no perturbation data".

usage: python3 knocktf_tf_arcs.py <analysis root> <network deposit folder> <original SIF folder> <P3b folder>
                                  <KnockTF cache folder> <output folder>
"""
import os, sys, hashlib, subprocess, collections
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, binomtest
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bhat_common as bc


def rng_for(key):                                      # as analyses/perturbation_tests _seed.py
    return np.random.default_rng([20250908, int(hashlib.md5(key.encode()).hexdigest()[:8], 16)])


def main(root, netdir, sifdir, p3b, kcache, out):
    os.makedirs(out, exist_ok=True)
    L = []; P = lambda *a: (L.append(' '.join(map(str, a))), print(*a, flush=True))
    N = bc.networks(root, netdir, sifdir)
    use = collections.defaultdict(collections.Counter)
    for k in bc.CLASSES:
        for (m1, m2, t1, t2, g1, g2) in N['inst']['6node_' + k]:
            if k in ('miRNA_FFL', 'composite_FFL'): use[(t1, t2)]['T1->T2'] += 1
            if k in ('TF_FFL', 'composite_FFL'): use[(t2, t1)]['T2->T1'] += 1
    rows = []
    for (a, b), c in sorted(use.items()):
        et, sign, ssrc, tier, layer = N['mfn'].attributes((a, b), N['census'], N['table'])
        rows.append(dict(source=a, target=b, layer=layer, sign=sign, sign_source=ssrc, used_as_T2_T1=c['T2->T1'], used_as_T1_T2=c['T1->T2']))
    A = pd.DataFrame(rows)
    P(f'TF-TF arcs used by the six-node networks: {len(A)}; by layer {A.layer.value_counts().to_dict()}; used as T2->T1 '
      f'{int((A.used_as_T2_T1 > 0).sum())}, as T1->T2 {int((A.used_as_T1_T2 > 0).sum())}; signed {int(A.sign.isin(["1", "-1"]).sum())} '
      f'(+1 {int((A.sign == "1").sum())}, -1 {int((A.sign == "-1").sum())})')
    SIGN = {(r.source, r.target): int(r.sign) for r in A.itertuples() if r.sign in ('1', '-1')}
    SETS = {'all used TF-TF arcs': A, 'used as T2->T1': A[A.used_as_T2_T1 > 0], 'used as T1->T2': A[A.used_as_T1_T2 > 0]}
    for s in list(SETS):
        for lay in ('analysed network', 'census layer: TRRUST TF-target'):
            SETS[f'{s}, {lay.replace("census layer: ", "")}'] = SETS[s][SETS[s].layer == lay]
    srcs = set(A.source)
    DS = pd.read_csv(os.path.join(p3b, 'p3b_datasets.tsv'), sep='\t')
    assert srcs <= set(pd.read_csv(os.path.join(root, 'data/canonical_edges.tsv'), sep='\t').query('edge_type == "TF_target"').source)
    keep = DS[(DS.n_control >= 2) & (DS.n_treat >= 2)]
    keep = keep[keep.tf.isin(srcs)]
    TR = pd.read_csv(os.path.join(root, 'data/db/trrust_human.tsv'), sep='\t', header=None, names=['tf', 'target', 'mode', 'pmid'])
    K = pd.concat([pd.read_csv(os.path.join(kcache, f), sep='\t', dtype={'Gene': str})
                   for f in ('knocktf_network_tf_subset.tsv.gz', 'knocktf_network_tf_subset_extra.tsv.gz') if os.path.exists(os.path.join(kcache, f))])
    for c in ('Mean_Control', 'Log2FC'): K[c] = pd.to_numeric(K[c], errors='coerce')
    K = K.dropna(subset=['Log2FC'])
    R = []; sig_seen = set(); per_arc = []
    for r in keep.itertuples():
        k = K[K.Sample_ID == r.sample_id]
        base = dict(sample_id=r.sample_id, tf=r.tf, knock_method=r.knock_method, biosample=r.biosample, profile_id=r.profile_id,
                    n_control=r.n_control, n_treat=r.n_treat)
        if not len(k): R.append(dict(base, status='no rows in KnockTF differential expression')); continue
        k = k.groupby('Gene').agg(Mean_Control=('Mean_Control', 'mean'), Log2FC=('Log2FC', 'mean'))
        tfl = float(k.loc[r.tf, 'Log2FC']) if r.tf in k.index else np.nan
        base['tf_log2fc'] = tfl
        if not (tfl <= -0.74): R.append(dict(base, status='knockdown not verified (TF Log2FC > -0.74 or TF absent)')); continue
        sig = (r.tf, r.biosample, r.n_control, r.n_treat, round(tfl, 4))
        if sig in sig_seen: R.append(dict(base, status='duplicate of an earlier dataset (same TF, biosample, sample counts and TF Log2FC)')); continue
        sig_seen.add(sig)
        med = k.Mean_Control.median(); ex = k[k.Mean_Control >= med]
        for t in A[A.source == r.tf].target:
            per_arc.append(dict(source=r.tf, target=t, sample_id=r.sample_id, target_log2fc=float(k.loc[t, 'Log2FC']) if t in k.index else np.nan,
                                target_expressed=bool(t in ex.index), sign=SIGN.get((r.tf, t))))
        alltr = set(TR[TR.tf == r.tf].target)
        for sname, es in SETS.items():
            tg = set(es[es.source == r.tf].target)
            if not tg: continue
            a = ex.Log2FC[ex.index.isin(tg)].abs(); b = ex.Log2FC[~ex.index.isin(alltr) & ~ex.index.isin(tg)].abs()
            if len(a) < 3:
                R.append(dict(base, edge_set=sname, status='included', n_targets=len(a), n_nontargets=len(b))); continue
            d = float(a.median() - b.median())
            rng = rng_for(f'{r.sample_id}|{sname}')
            se = float(np.std([np.median(rng.choice(a.values, len(a))) - np.median(rng.choice(b.values, len(b))) for _ in range(1000)], ddof=1))
            p = float(mannwhitneyu(a, b, alternative='two-sided').pvalue)
            sg = [(g, SIGN[(r.tf, g)]) for g in ex.index[ex.index.isin(tg)] if (r.tf, g) in SIGN]
            agree = sum(1 for g, s_ in sg if (s_ == 1 and ex.Log2FC[g] < 0) or (s_ == -1 and ex.Log2FC[g] > 0))
            R.append(dict(base, edge_set=sname, status='included', n_targets=len(a), n_nontargets=len(b), abs_diff=d, abs_diff_se=se,
                          wilcoxon_p=p, n_signed=len(sg), n_agree=agree))
    R = pd.DataFrame(R); R.to_csv(os.path.join(out, 'A4_per_dataset.tsv'), sep='\t', index=False)
    inc = R[(R.status == 'included') & R.abs_diff.notna()] if 'abs_diff' in R else R.iloc[0:0]
    P(f'KnockTF 2.0 datasets of the {len(srcs)} source TFs with >= 2 control and >= 2 treated samples: {len(keep)} '
      f'({keep.tf.nunique()} TFs); knockdown verified and not duplicated: {R[R.status == "included"].sample_id.nunique()} '
      f'({R[R.status == "included"].tf.nunique()} TFs); not verified: {int(R.status.str.startswith("knockdown not verified").sum())}; '
      f'duplicates: {int(R.status.str.startswith("duplicate").sum())}; no rows: {int(R.status.str.startswith("no rows").sum())}')
    pin = pd.DataFrame([dict(group=s, yi=x.abs_diff, sei=x.abs_diff_se) for s in SETS for x in inc[inc.edge_set == s].itertuples()])
    PT = pd.DataFrame()
    if len(pin):
        pin.to_csv(os.path.join(out, 'A4_pool_in.tsv'), sep='\t', index=False)
        subprocess.run(['Rscript', os.path.join(HERE, '..', '..', 'perturbation_tests', '_rma.R'), os.path.join(out, 'A4_pool_in.tsv'), os.path.join(out, 'A4_pooled.tsv')], check=True)
        PT = pd.read_csv(os.path.join(out, 'A4_pooled.tsv'), sep='\t')
    dec = []
    for s, es in SETS.items():
        x = PT[PT.group == s] if len(PT) else PT; y = inc[inc.edge_set == s]
        if not len(x) or x.iloc[0].k == 0:
            P(f'{s}: {len(es)} arcs; no dataset with >= 3 of its targets expressed -> NOT FOUND'); dec.append(dict(arc_set=s, arcs=len(es), decision='NOT FOUND')); continue
        x = x.iloc[0]; na, ns = int(y.n_agree.sum()), int(y.n_signed.sum())
        bt = binomtest(na, ns, 0.5, alternative='two-sided') if ns else None
        ok = bool((x.ci_lb > 0) and ns and (na / ns > 0.5) and (bt.pvalue < 0.05))
        P(f'{s}: {len(es)} arcs; k = {int(x.k)} datasets ({y.tf.nunique()} TFs); pooled |log2FC| difference {x.est:+.4f} '
          f'(95 % CI {x.ci_lb:+.4f} to {x.ci_ub:+.4f}), p {x.p:.3g}; sign agreement {na}/{ns}'
          + (f' = {100 * na / ns:.1f} % (binomial p {bt.pvalue:.3g})' if ns else '') + f' -> {"SUPPORTED" if ok else "NOT SUPPORTED"}')
        dec.append(dict(arc_set=s, arcs=len(es), k=int(x.k), tfs=y.tf.nunique(), pooled_abs_diff=x.est, ci_lb=x.ci_lb, ci_ub=x.ci_ub, p=x.p,
                        n_agree=na, n_signed=ns, agreement=na / ns if ns else np.nan, binom_p=bt.pvalue if ns else np.nan,
                        decision='SUPPORTED' if ok else 'NOT SUPPORTED'))
    pd.DataFrame(dec).to_csv(os.path.join(out, 'A4_decision.csv'), index=False)
    PA = pd.DataFrame(per_arc)
    agg = []
    for r in A.itertuples():
        x = PA[(PA.source == r.source) & (PA.target == r.target)] if len(PA) else PA
        xv = x[x.target_log2fc.notna()] if len(x) else x
        s_ = SIGN.get((r.source, r.target))
        nag = int(sum((s_ == 1 and v < 0) or (s_ == -1 and v > 0) for v in xv.target_log2fc)) if (s_ and len(xv)) else None
        agg.append(dict(source=r.source, target=r.target, layer=r.layer, sign=r.sign, used_as_T2_T1=r.used_as_T2_T1, used_as_T1_T2=r.used_as_T1_T2,
                        datasets_of_source_tf=len(x), datasets_with_target=len(xv),
                        median_target_log2fc=float(xv.target_log2fc.median()) if len(xv) else np.nan,
                        n_agree=nag, status='no perturbation data' if not len(xv) else 'measured'))
    AG = pd.DataFrame(agg); AG.to_csv(os.path.join(out, 'A4_arcs.csv'), index=False)
    PA.to_csv(os.path.join(out, 'A4_per_arc_dataset.csv'), index=False)
    for s in ('all used TF-TF arcs', 'used as T2->T1', 'used as T1->T2'):
        es = SETS[s]; m = AG.merge(es[['source', 'target']], on=['source', 'target'])
        sg = m[m.sign.isin(['1', '-1']) & (m.status == 'measured')]
        P(f'per arc, {s}: {len(m)} arcs; with perturbation data {int((m.status == "measured").sum())}; no perturbation data '
          f'{int((m.status != "measured").sum())}; signed and measured {len(sg)}, of which the median Log2FC agrees with the sign for '
          f'{int(sum((x.sign == "1" and x.median_target_log2fc < 0) or (x.sign == "-1" and x.median_target_log2fc > 0) for x in sg.itertuples()))}')
    open(os.path.join(out, 'A4_summary.txt'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main(*sys.argv[1:7])
