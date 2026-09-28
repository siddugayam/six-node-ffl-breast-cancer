#!/usr/bin/env python3
"""P4 summary: applies the revised decision rule (SETTINGS.md, change 1) to the deposited Enformer values (sanity check)
and, when the Borzoi scores exist, to Borzoi.  Writes p4_calls.csv (one row per call), p4_controls.csv and
p4_summary.txt (numbered lines).
Effect of a call (Borzoi) = mean over its 5 shuffles of (score - reference), per replicate; then the mean over the
replicates present.  Primary score: fibroblast RNA-seq exon sum; CAGE (TSS, 384 bp) reported alongside.
Most influential = most negative effect; rank 1 = most negative (ties: minimum rank), as in 44_seqreg_ext_synthesis.py.
Families (44_ lines 14-16): NF-kB = NFKB1, NFKB2, RELA, RELB, REL; SP = SP1, SP2, SP3; ETS = ETS1, ETS2, ELK1, ELK4, GABPA."""
import os, glob
import numpy as np, pandas as pd
from scipy.stats import spearmanr
OUT = os.path.dirname(os.path.abspath(__file__)); REV = os.path.dirname(os.path.dirname(os.path.dirname(OUT))); RES = f'{REV}/results/v3'
FAM = {"NF-kB": ["NFKB1", "NFKB2", "RELA", "RELB", "REL"], "SP": ["SP1", "SP2", "SP3"], "ETS": ["ETS1", "ETS2", "ELK1", "ELK4", "GABPA"]}
DONOR = (184, 210)
L = []; P = L.append

def rule(calls, eff, ctrl_donor, ctrl_ggaa, label):
    """calls: DataFrame with region, motif, tf_name, rel_to_TSS, site_len and the effect column eff (negative = drop)."""
    out = {}
    g = calls[calls.region == 'COL1A1'].copy(); g['rank'] = g[eff].rank(method='min').astype(int)
    top = g.nsmallest(1, eff).iloc[0]
    top_fam = next((f for f, t in FAM.items() if top.tf_name in t), None)
    ok1 = top_fam is None; parts = []
    for f, tfs in FAM.items():
        s = g[g.tf_name.isin(tfs)]
        b = s.nsmallest(1, eff).iloc[0]; ratio = b[eff] / top[eff]
        ok1 = ok1 and ratio <= 1 / 3
        parts.append(f'{f}: {b.motif} at TSS{int(b.rel_to_TSS):+d}, rank {int(b["rank"])}/{len(g)}, effect {b[eff]:.4g} ({b["pct"]:+.2f} %), ratio to top {ratio:.3f}')
    P(f'{label} criterion 1 (COL1A1): top call {top.motif} ({top.tf_name}) at TSS{int(top.rel_to_TSS):+d}, effect {top[eff]:.4g} ({top["pct"]:+.2f} %), '
      f'family {top_fam or "none of NF-kB/SP/ETS"}; ' + '; '.join(parts) + f' -> {"PASS" if ok1 else "FAIL"}')
    h = calls[calls.region == 'COL3A1'].copy(); h['rank'] = h[eff].rank(method='min').astype(int)
    t3 = h.nsmallest(1, eff).iloc[0]
    lo, hi = int(t3.rel_to_TSS), int(t3.rel_to_TSS) + int(t3.site_len) - 1
    ok2 = not (hi < DONOR[0] or lo > DONOR[1])
    P(f'{label} criterion 2 (COL3A1): top call {t3.motif} ({t3.tf_name}) at TSS+{lo}..+{hi}, effect {t3[eff]:.4g} ({t3["pct"]:+.2f} %); '
      f'donor window +{DONOR[0]}..+{DONOR[1]} -> {"PASS" if ok2 else "FAIL"}')
    ok3 = ctrl_donor < ctrl_ggaa          # both negative: the donor kill lowers the prediction more
    P(f'{label} criterion 3 (COL3A1): donor GT->AA effect {ctrl_donor:.4g}; GGAA->CCTT effect {ctrl_ggaa:.4g} -> {"PASS" if ok3 else "FAIL"}')
    P(f'{label} decision: {"REPLICATED" if (ok1 and ok2 and ok3) else "NOT REPLICATED"} (criteria {int(ok1)}{int(ok2)}{int(ok3)})')
    return ok1 and ok2 and ok3

# ------------------------------------------------------------------ Enformer (deposited)
E = pd.read_csv(f'{RES}/seqreg_ext_occlusion_allsites.csv')
E['pct'] = E.pct_change_fib_cage
EC = pd.read_csv(f'{RES}/seqreg_ext_splicedonor_control.csv')
ec = lambda v: float(EC[(EC.region == 'COL3A1') & (EC.variant == v)].d_fib_cage.iloc[0])
rule(E, 'd_fib_cage', ec('donor_kill_GT_to_AA'), ec('ets_core_kill_GGAA_to_CCTT'), 'Enformer (deposited, fibroblast CAGE)')

# ------------------------------------------------------------------ Borzoi
reps = sorted(glob.glob(f'{OUT}/p4_scores_rep*.tsv'))
full = []
def read_scores(p):
    # p4_borzoi_occlusion.py wrote each score as a NumPy repr, 'np.float64(<value>)'; the value inside is the full number
    d = pd.read_csv(p, sep='\t', dtype=str)
    for c in ('rna_exon_sum', 'cage_tss_sum'):
        d[c] = d[c].str.replace(r'^np\.float(?:32|64)\((.*)\)$', r'\1', regex=True).astype(float)
    return d
for p in reps:
    d = read_scores(p)
    if len(d) >= 8140: full.append((int(p.split('rep')[-1].split('.')[0]), d))
P(f'Borzoi replicates complete: {[r for r, _ in full]}')
if full:
    sites = E[['region', 'motif', 'tf_name', 'rel_to_TSS', 'site_len', 'd_fib_cage', 'pct_change_fib_cage']].copy()
    sites['n'] = sites.groupby('region').cumcount()
    per = []
    for r, d in full:
        d = d.set_index('key')
        for reg in ('COL1A1', 'COL3A1'):
            ref = d.loc[f'{reg}|ref']
            c = d[d.index.str.startswith(f'{reg}|call|')].copy()
            c['n'] = c.index.str.split('|').str[2].astype(int)
            m = c.groupby('n')[['rna_exon_sum', 'cage_tss_sum']].mean()
            m['d_rna'] = m.rna_exon_sum - ref.rna_exon_sum; m['d_cage'] = m.cage_tss_sum - ref.cage_tss_sum
            m['ref_rna'] = ref.rna_exon_sum; m['ref_cage'] = ref.cage_tss_sum; m['rep'] = r; m['region'] = reg
            per.append(m.reset_index())
    per = pd.concat(per)
    agg = per.groupby(['region', 'n']).agg(d_rna=('d_rna', 'mean'), d_cage=('d_cage', 'mean'), ref_rna=('ref_rna', 'mean'), ref_cage=('ref_cage', 'mean'),
                                           d_rna_sd_reps=('d_rna', 'std'), n_reps=('rep', 'nunique')).reset_index()
    B = sites.merge(agg, on=['region', 'n'])
    B['pct'] = 100 * B.d_rna / B.ref_rna; B['pct_cage'] = 100 * B.d_cage / B.ref_cage
    for reg in ('COL1A1', 'COL3A1'):
        m = B.region == reg
        B.loc[m, 'rank_rna'] = B.loc[m, 'd_rna'].rank(method='min'); B.loc[m, 'rank_cage'] = B.loc[m, 'd_cage'].rank(method='min')
        B.loc[m, 'rank_enformer'] = B.loc[m, 'd_fib_cage'].rank(method='min')
    B.to_csv(f'{OUT}/p4_calls.csv', index=False)
    # controls
    rows = []
    for r, d in full:
        d = d.set_index('key')
        for reg in ('COL1A1', 'COL3A1'):
            ref = d.loc[f'{reg}|ref']
            for k in d.index[d.index.str.startswith(f'{reg}|ctrl|')]:
                rows.append(dict(region=reg, variant=k.split('|')[2], description=k.split('|')[3], rep=r,
                                 d_rna=d.loc[k].rna_exon_sum - ref.rna_exon_sum, pct_rna=100 * (d.loc[k].rna_exon_sum - ref.rna_exon_sum) / ref.rna_exon_sum,
                                 d_cage=d.loc[k].cage_tss_sum - ref.cage_tss_sum, pct_cage=100 * (d.loc[k].cage_tss_sum - ref.cage_tss_sum) / ref.cage_tss_sum))
    C = pd.DataFrame(rows).groupby(['region', 'variant', 'description']).agg(d_rna=('d_rna', 'mean'), pct_rna=('pct_rna', 'mean'), d_cage=('d_cage', 'mean'),
                                                                             pct_cage=('pct_cage', 'mean'), n_reps=('rep', 'nunique')).reset_index()
    C.to_csv(f'{OUT}/p4_controls.csv', index=False)
    bc = lambda v, col: float(C[(C.region == 'COL3A1') & (C.variant == v)][col].iloc[0])
    lab = f'Borzoi ({len(full)} replicates, fibroblast RNA-seq exon sum)'
    rule(B, 'd_rna', bc('donor_kill_GT_to_AA', 'd_rna'), bc('ets_core_kill_GGAA_to_CCTT', 'd_rna'), lab)
    Bc = B.copy(); Bc['pct'] = Bc.pct_cage
    rule(Bc, 'd_cage', bc('donor_kill_GT_to_AA', 'd_cage'), bc('ets_core_kill_GGAA_to_CCTT', 'd_cage'), f'Borzoi ({len(full)} replicates, fibroblast CAGE at the TSS; secondary)')
    for reg in ('COL1A1', 'COL3A1'):
        g = B[B.region == reg]
        r1 = spearmanr(g.d_fib_cage, g.d_rna); r2 = spearmanr(g.d_fib_cage, g.d_cage)
        P(f'{reg}: Spearman rho, Enformer fibroblast CAGE vs Borzoi RNA exon sum over {len(g)} calls = {r1.statistic:.3f} (p {r1.pvalue:.2g}); '
          f'vs Borzoi CAGE = {r2.statistic:.3f} (p {r2.pvalue:.2g})')
        for col, nm in (('d_rna', 'RNA'), ('d_cage', 'CAGE')):
            t = g.nsmallest(3, col)
            P(f'{reg}: top 3 calls by Borzoi {nm}: ' + '; '.join(f'{x.motif} at TSS{int(x.rel_to_TSS):+d} ({x.pct if nm == "RNA" else x.pct_cage:+.2f} %)' for x in t.itertuples()))
    for x in C.itertuples():
        P(f'control {x.region} {x.variant} ({x.description}): RNA {x.pct_rna:+.2f} %, CAGE {x.pct_cage:+.2f} % ({x.n_reps} replicates)')
open(f'{OUT}/p4_summary.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
