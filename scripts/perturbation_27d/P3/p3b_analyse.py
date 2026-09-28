#!/usr/bin/env python3
"""P3b: validate the analysed network's TF->target edges with KnockTF 2.0 human knockdown datasets (SETTINGS.md, P3).
Inputs: p3b_datasets.tsv (datasets of network TFs with their control/treated counts); KnockTF's differential expression of
all human datasets (knocktf_v2_main_human.txt, processed values).
Per dataset (>= 2 control and >= 2 treated; TF Log2FC <= -0.74; duplicates of the same TF, biosample, sample counts and
TF Log2FC kept once):
  |log2FC| test: median |Log2FC| of the TF's network targets minus that of expressed genes that are not TRRUST targets of
  the TF in any mode (expressed = Mean_Control >= the dataset median); two-sided Wilcoxon; bootstrap SE (1,000, seed 20250908);
  sign agreement: Activation targets expected down, Repression targets up (Log2FC = 0 counts as disagreement).
Pooled: REML for the |log2FC| difference; agreement = sum of agreements / sum of signed targets, two-sided binomial test.
Decision (per edge set): SUPPORTED if the pooled difference > 0 with 95 % CI excluding 0 AND agreement > 50 % with
binomial p < 0.05; NOT SUPPORTED otherwise.  Edge sets: TRRUST, hTFtarget-only, all edges (the network's 770 TF->target
edges are all in TRRUST, so hTFtarget-only is empty).
Writes p3b_per_dataset.tsv, p3b_pooled.tsv and p3b_summary.txt (numbered lines)."""
import os, sys, csv, subprocess
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, binomtest
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); REV = os.path.dirname(INB); sys.path.insert(0, INB)
import _dl
from _seed import rng_for
D = f'{INB}/P3'
DS = pd.read_csv(f'{D}/p3b_datasets.tsv', sep='\t')
E = pd.read_csv(f'{REV}/data/canonical_edges.tsv', sep='\t'); E = E[E.edge_type == 'TF_target']
TR = pd.read_csv(f'{REV}/data/db/trrust_human.tsv', sep='\t', header=None, names=['tf', 'target', 'mode', 'pmid'])
MODE = {(a, b): m for a, b, m in zip(TR.tf, TR.target, TR['mode'])}
trs = set(zip(TR.tf, TR.target))
E['in_trrust'] = [(a, b) in trs for a, b in zip(E.source, E.target)]
SETS = {'TRRUST': E[E.in_trrust], 'hTFtarget-only': E[~E.in_trrust], 'all edges': E}
keep = DS[(DS.n_control >= 2) & (DS.n_treat >= 2)]
want = set(keep.sample_id)
cols = ['Sample_ID', 'TF', 'Gene', 'Mean_Control', 'Log2FC']
sub = os.path.join(_dl.CACHE, 'P3/knocktf/knocktf_network_tf_subset.tsv.gz')      # filtered rows, cached for re-runs
if os.path.exists(sub):
    K = pd.read_csv(sub, sep='\t', dtype={'Gene': str})
else:
    chunks = []
    for ch in pd.read_csv(os.path.join(_dl.CACHE, 'P3/knocktf/knocktf_v2_main_human.txt'), sep='\t', usecols=cols, chunksize=2_000_000, dtype=str):
        chunks.append(ch[ch.Sample_ID.isin(want)])
    K = pd.concat(chunks); K.to_csv(sub, sep='\t', index=False, compression='gzip')
extra = os.path.join(_dl.CACHE, 'P3/knocktf/knocktf_network_tf_subset_extra.tsv.gz')  # datasets absent from the bulk file (p3b_fetch_missing.py)
if os.path.exists(extra): K = pd.concat([K, pd.read_csv(extra, sep='\t', dtype={'Gene': str})])
for c in ('Mean_Control', 'Log2FC'): K[c] = pd.to_numeric(K[c], errors='coerce')     # a few non-numeric entries -> NaN
K = K.dropna(subset=['Log2FC']); print('rows kept', len(K), flush=True)
rows = []; sig_seen = set()
for r in keep.itertuples():
    k = K[K.Sample_ID == r.sample_id]
    if not len(k): rows.append(dict(sample_id=r.sample_id, tf=r.tf, status='no rows in KnockTF differential expression')); continue
    k = k.groupby('Gene').agg(Mean_Control=('Mean_Control', 'mean'), Log2FC=('Log2FC', 'mean'))
    tfl = float(k.loc[r.tf, 'Log2FC']) if r.tf in k.index else np.nan
    base = dict(sample_id=r.sample_id, tf=r.tf, knock_method=r.knock_method, biosample=r.biosample, profile_id=r.profile_id, n_control=r.n_control,
                n_treat=r.n_treat, tf_log2fc=tfl)
    if not (tfl <= -0.74):
        rows.append(dict(base, status='knockdown not verified (TF Log2FC > -0.74 or TF absent)')); continue
    sig = (r.tf, r.biosample, r.n_control, r.n_treat, round(tfl, 4))
    if sig in sig_seen: rows.append(dict(base, status='duplicate of an earlier dataset (same TF, biosample, sample counts and TF Log2FC)')); continue
    sig_seen.add(sig)
    ex = k[k.Mean_Control >= k.Mean_Control.median()]
    alltr = set(TR[TR.tf == r.tf].target)
    for sname, es in SETS.items():
        tg = set(es[es.source == r.tf].target)
        a = ex.Log2FC[ex.index.isin(tg)].abs(); b = ex.Log2FC[~ex.index.isin(alltr) & ~ex.index.isin(tg)].abs()
        if len(a) < 3:
            rows.append(dict(base, edge_set=sname, status='included', n_targets=len(a), n_nontargets=len(b))); continue
        d = float(a.median() - b.median())
        rng = rng_for(f'{r.sample_id}|{sname}')
        se = float(np.std([np.median(rng.choice(a.values, len(a))) - np.median(rng.choice(b.values, len(b))) for _ in range(1000)], ddof=1))
        p = float(mannwhitneyu(a, b, alternative='two-sided').pvalue)
        sg = [(g, MODE.get((r.tf, g))) for g in ex.index[ex.index.isin(tg)] if MODE.get((r.tf, g)) in ('Activation', 'Repression')]
        agree = sum(1 for g, mo in sg if (mo == 'Activation' and ex.Log2FC[g] < 0) or (mo == 'Repression' and ex.Log2FC[g] > 0))
        rows.append(dict(base, edge_set=sname, status='included', n_targets=len(a), n_nontargets=len(b), abs_diff=d, abs_diff_se=se, wilcoxon_p=p,
                         n_signed=len(sg), n_agree=agree))
R = pd.DataFrame(rows); R.to_csv(f'{D}/p3b_per_dataset.tsv', sep='\t', index=False)
inc = R[(R.status == 'included') & R.abs_diff.notna()] if 'abs_diff' in R else R.iloc[0:0]
L = []; P = L.append
P(f'KnockTF 2.0 datasets of network TFs: {len(DS)} ({DS.tf.nunique()} TFs); with >= 2 control and >= 2 treated samples: {len(keep)}; '
  f'knockdown verified and not duplicated: {R[R.status == "included"].sample_id.nunique()} ({R[R.status == "included"].tf.nunique()} TFs); '
  f'not verified: {int((R.status.str.startswith("knockdown not verified")).sum())}; duplicates: {int((R.status.str.startswith("duplicate")).sum())}')
pin = [dict(group=s, yi=x.abs_diff, sei=x.abs_diff_se) for s in SETS for x in inc[inc.edge_set == s].itertuples()]
PT = pd.DataFrame()
if pin:
    pd.DataFrame(pin).to_csv(f'{D}/p3b_pool_in.tsv', sep='\t', index=False)
    subprocess.run(['Rscript', f'{INB}/_rma.R', f'{D}/p3b_pool_in.tsv', f'{D}/p3b_pooled.tsv'], check=True)
    PT = pd.read_csv(f'{D}/p3b_pooled.tsv', sep='\t')
for s, es in SETS.items():
    x = PT[PT.group == s] if len(PT) else PT
    if len(es) == 0: P(f'{s}: no edge of this kind in the analysed network (all 770 TF->target edges are in TRRUST v2) -> NOT FOUND'); continue
    if not len(x) or x.iloc[0].k == 0: P(f'{s}: no dataset'); continue
    x = x.iloc[0]; y = inc[inc.edge_set == s]
    na, ns = int(y.n_agree.sum()), int(y.n_signed.sum())
    bt = binomtest(na, ns, 0.5, alternative='two-sided') if ns else None
    ok = (x.ci_lb > 0) and ns and (na / ns > 0.5) and (bt.pvalue < 0.05)
    P(f'{s}: k = {int(x.k)} datasets ({y.tf.nunique()} TFs); pooled |log2FC| difference {x.est:+.4f} (95 % CI {x.ci_lb:+.4f} to {x.ci_ub:+.4f}), p {x.p:.3g}; '
      f'sign agreement {na}/{ns} = {100 * na / ns:.1f} % (binomial p {bt.pvalue:.3g}) -> {"SUPPORTED" if ok else "NOT SUPPORTED"}')
open(f'{D}/p3b_summary.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
