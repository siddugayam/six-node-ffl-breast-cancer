#!/usr/bin/env python3
"""
10 -- SPECIFICITY control: is the FFL structure specific to breast cancer, or does it merely
      reflect the general density of known molecular interactions?


DESIGN.  Every network compared here -- breast cancer and control alike -- is assembled by ONE
identical pipeline from the same three primary sources, so that nothing differs except which
nodes were chosen:
    TF -> target   : TRRUST v2 curated pairs
    TF -> miRNA    : TransmiR v2 curated pairs
    miRNA -> target: multiMiR *validated* interactions (miRTarBase / TarBase / miRecords),
                     retrieved BY TARGET GENE so that each gene's miRNA in-degree is complete
                     and genome-wide, hence directly comparable across gene sets.
The authors' miRNA-miRNA and gene-gene layers are excluded from BOTH sides: the audit showed
28 of their 30 miRNA-miRNA edges have no co-transcription support and their single gene-gene
edge (COL1A1-COL3A1) is an undirected STRING association, so including them would advantage
the BRCA network by construction.

FOUR NETWORK CLASSES
  BRCA_published   the deposited canonical network exactly as it stands (6859 edges; includes
                   miRWalk-predicted miRNA-target edges).  Reported for reference only -- it
                   is NOT comparable to the controls because its miRNA-target layer is 37.4%
                   pure sequence prediction.
  BRCA_rebuilt     the 587 breast-cancer nodes run through the control pipeline (validated
                   miRNA-target edges only).  THIS is the correct comparator.
  CTRL_random      157 TFs drawn uniformly from TRRUST regulators not in the BRCA network,
                   207 genes drawn uniformly from protein-coding genes measured in TCGA-BRCA
                   and not in the BRCA network, 223 miRNAs drawn uniformly from miRBase mature
                   stems not in the BRCA network.  Node-type composition matched exactly.
  CTRL_topann      same sizes, drawn from the MOST HEAVILY ANNOTATED non-BRCA nodes of each
                   type (top 2x oversampled pool).  This is the strongest control that the
                   annotation databases permit, and it exists because CTRL_matched cannot
                   actually reach the BRCA annotation level -- see the annotation-rank table
                   written to results/motif_specificity_annotation_rank.csv.
  CTRL_matched     same sizes, but each control node is drawn from the 15 nearest neighbours
                   of the corresponding BRCA node in ANNOTATION-DEGREE space (log1p of the
                   node's genome-wide TRRUST / TransmiR / multiMiR-validated degrees).  This
                   control removes the "well-studied genes are better annotated" confound,
                   which is precisely the "general density of known interactions" the concern
                   is worried about.

OUTPUT: results/motif_specificity.csv -- one row per network.
"""
import csv, collections, json, math, os, random, time
import numpy as np

REV = '/path/to/revision'
D = f'{REV}/data/control'
LOG = f'{REV}/logs/motif_specificity.log'
N_CTRL = int(os.environ.get('N_CTRL', 30))
KNN = int(os.environ.get('KNN', 15))
SEED = 20260908

def say(*a):
    m = time.strftime('%H:%M:%S') + ' | ' + ''.join(str(x) for x in a)
    print(m, flush=True)
    with open(LOG, 'a') as fh: fh.write(m + '\n')

def rd(path, sep='\t'):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter=sep))

# ------------------------------------------------------------------ FFL counting
def count_ffl(edges, ntype):
    """edges: iterable of (src, tgt) name pairs.  ntype: name -> 'TF'|'miRNA'|'Gene'.
    Same 3-node FFL definition as scripts/08 and 02x."""
    ns = sorted({x for e in edges for x in e})
    idx = {v: i for i, v in enumerate(ns)}; n = len(ns)
    outb = [0] * n
    eset = set()
    for a, b in edges:
        ia, ib = idx[a], idx[b]
        outb[ia] |= (1 << ib); eset.add(ia * n + ib)
    nm = 0
    for v, i in idx.items():
        if ntype[v] != 'miRNA': nm |= (1 << i)
    c = collections.Counter()
    for a, b in edges:
        ia, ib = idx[a], idx[b]
        com = outb[ia] & outb[ib]
        if not com: continue
        c['Any-3node-FFL'] += com.bit_count()
        k = (com & nm).bit_count()
        if not k: continue
        ta, tb = ntype[a], ntype[b]
        cls = None
        if ta == 'TF' and tb == 'miRNA':
            cls = 'Composite-FFL' if (ib * n + ia) in eset else 'TF-FFL'
        elif ta == 'miRNA' and tb == 'TF':
            cls = 'Composite-FFL' if (ib * n + ia) in eset else 'miRNA-FFL'
        elif ta == 'TF' and tb == 'TF':
            cls = 'TF-TF-FFL'
        if cls: c[cls] += k
    c['n_connected_nodes'] = n
    return c

# ------------------------------------------------------------------ pipeline
def build_network(tfs, genes, mirs, TRR, TMR, VAL):
    """The one pipeline used for every network in this analysis."""
    tfs = set(tfs); genes = set(genes); mirs = set(mirs)
    prot = tfs | genes
    edges = []                       # iteration is over SORTED sets so the edge ORDER, and
    for t in sorted(tfs):            # hence any downstream subsampling, is reproducible
        for g in sorted(TRR.get(t, ())):
            if g in prot and g != t: edges.append((t, g))
    for t in sorted(tfs):
        for mi in sorted(TMR.get(t, ())):
            if mi in mirs: edges.append((t, mi))
    for mi in sorted(mirs):
        for g in sorted(VAL.get(mi, ())):
            if g in prot: edges.append((mi, g))
    return list(dict.fromkeys(edges))

def main():
    say('=== 10_motif_specificity :: START ===')
    say(f'N_CTRL={N_CTRL} per control class, KNN={KNN}, SEED={SEED}')
    nodes = {r['name']: r['type'] for r in rd(f'{REV}/data/canonical_nodes.tsv')}
    brca_tf = sorted(v for v, t in nodes.items() if t == 'TF')
    brca_gn = sorted(v for v, t in nodes.items() if t == 'Gene')
    brca_mi = sorted(v for v, t in nodes.items() if t == 'miRNA')
    say(f'BRCA node sizes: TF={len(brca_tf)} Gene={len(brca_gn)} miRNA={len(brca_mi)}')

    TRR = collections.defaultdict(set); TRR_in = collections.Counter()
    for r in rd(f'{D}/trrust_pairs.tsv'):
        TRR[r['TF']].add(r['target']); TRR_in[r['target']] += 1
    TMR = collections.defaultdict(set); TMR_in = collections.Counter()
    for r in rd(f'{D}/transmir_pairs.tsv'):
        TMR[r['TF']].add(r['mirna']); TMR_in[r['mirna']] += 1
    VAL = collections.defaultdict(set); VAL_in = collections.Counter()
    for r in rd(f'{D}/validated_mirna_target_pairs.tsv'):
        VAL[r['mirna']].add(r['gene']); VAL_in[r['gene']] += 1
    queried = {r['gene'] for r in rd(f'{D}/queried_gene_universe.tsv')}
    say(f'TRRUST regulators={len(TRR)}  TransmiR TFs={len(TMR)}  '
        f'validated miRNAs={len(VAL)}  queried gene universe={len(queried)}')

    # ---------- universes for control sampling (all disjoint from the BRCA node set) ----------
    brca_all = set(nodes)
    tf_uni = sorted(set(readlines(f'{REV}/cache/control_multimir/pool_tfs.txt')) - brca_all)
    gn_uni = sorted(set(readlines(f'{REV}/cache/control_multimir/pool_genes.txt')) - brca_all)
    mi_uni = sorted(set(readlines(f'{D}/mirbase_mature_stems.txt')) - brca_all)
    say(f'control universes: TF={len(tf_uni)}  Gene={len(gn_uni)}  miRNA={len(mi_uni)}')
    assert len(tf_uni) >= len(brca_tf) and len(gn_uni) >= len(brca_gn) and len(mi_uni) >= len(brca_mi)
    for g in gn_uni[:50] + tf_uni[:50]:
        assert g in queried, f'{g} was never queried against multiMiR'

    # ---------- annotation-degree feature vectors ----------
    def feat_tf(v):  return (len(TRR.get(v, ())), len(TMR.get(v, ())), TRR_in[v])
    def feat_gn(v):  return (TRR_in[v], VAL_in[v])
    def feat_mi(v):  return (TMR_in[v], len(VAL.get(v, ())))
    say('annotation degree, BRCA vs control universe (median [IQR]):')
    for lbl, fn, a, b in (('TF  (TRRUST-out, TransmiR-out, TRRUST-in)', feat_tf, brca_tf, tf_uni),
                          ('Gene(TRRUST-in, validated-miRNA-in)',       feat_gn, brca_gn, gn_uni),
                          ('miR (TransmiR-in, validated-target-out)',   feat_mi, brca_mi, mi_uni)):
        A = np.array([fn(x) for x in a], float); B = np.array([fn(x) for x in b], float)
        say(f'   {lbl}')
        say(f'      BRCA    n={len(a):<5} ' + ' | '.join(
            f'{np.median(A[:,j]):.0f} [{np.percentile(A[:,j],25):.0f}-{np.percentile(A[:,j],75):.0f}]'
            for j in range(A.shape[1])))
        say(f'      control n={len(b):<5} ' + ' | '.join(
            f'{np.median(B[:,j]):.0f} [{np.percentile(B[:,j],25):.0f}-{np.percentile(B[:,j],75):.0f}]'
            for j in range(B.shape[1])))

    def matched_draw(targets, universe, fn, rng):
        T = np.log1p(np.array([fn(x) for x in targets], float))
        U = np.log1p(np.array([fn(x) for x in universe], float))
        sd = U.std(axis=0); sd[sd == 0] = 1.0
        Tn = T / sd; Un = U / sd
        used = set(); out = []
        order = rng.sample(range(len(targets)), len(targets))
        for ti in order:
            d = np.sqrt(((Un - Tn[ti]) ** 2).sum(axis=1))
            cand = np.argsort(d, kind='stable')
            pool = [int(ci) for ci in cand if ci not in used][:KNN]
            picked = rng.choice(pool)
            used.add(picked); out.append(universe[picked])
        return out

    # ---------- assemble every network ----------
    rows = []
    nodesets = {}
    def add(label, cls, tfs, genes, mirs, edges=None):
        nodesets[label] = dict(cls=cls, TF=list(tfs), Gene=list(genes), miRNA=list(mirs))
        if edges is not None:                 # networks whose edge list is not reproducible
            nodesets[label]['edges'] = [list(e) for e in edges]   # from the node set alone
        ntl = {}
        for v in tfs: ntl[v] = 'TF'
        for v in genes: ntl[v] = 'Gene'
        for v in mirs: ntl[v] = 'miRNA'
        if edges is None: edges = build_network(tfs, genes, mirs, TRR, TMR, VAL)
        c = count_ffl(edges, ntl) if edges else collections.Counter({'n_connected_nodes': 0})
        ne = len(edges); nn = len(tfs) + len(genes) + len(mirs)
        nf = c['Any-3node-FFL']
        rows.append(dict(network=label, network_class=cls, n_nodes=nn,
                         n_nodes_connected=c['n_connected_nodes'], n_edges=ne, n_FFL=nf,
                         FFL_per_node=round(nf / nn, 6) if nn else 0.0,
                         FFL_per_edge=round(nf / ne, 6) if ne else 0.0,
                         n_miRNA_FFL=c['miRNA-FFL'], n_TF_FFL=c['TF-FFL'],
                         n_Composite_FFL=c['Composite-FFL'], n_TF_TF_FFL=c['TF-TF-FFL'],
                         n_TF=len(tfs), n_Gene=len(genes), n_miRNA=len(mirs)))
        return rows[-1]

    pub = [(r['source'], r['target']) for r in rd(f'{REV}/data/canonical_edges.tsv')
           if r['edge_type'] in ('miRNA_target', 'TF_miRNA', 'TF_target')]
    r0 = add('BRCA_published', 'BRCA_published', brca_tf, brca_gn, brca_mi, edges=pub)
    say(f"BRCA_published (miRNA_target+TF_miRNA+TF_target only): {r0['n_edges']} edges, "
        f"{r0['n_FFL']} FFL cores, {r0['FFL_per_node']:.4f}/node, {r0['FFL_per_edge']:.4f}/edge")
    r1 = add('BRCA_rebuilt', 'BRCA_rebuilt', brca_tf, brca_gn, brca_mi)
    say(f"BRCA_rebuilt (control pipeline): {r1['n_edges']} edges, {r1['n_FFL']} FFL cores, "
        f"{r1['FFL_per_node']:.4f}/node, {r1['FFL_per_edge']:.4f}/edge")

    rng = random.Random(SEED)
    for i in range(1, N_CTRL + 1):
        add(f'CTRL_random_{i:02d}', 'CTRL_random',
            rng.sample(tf_uni, len(brca_tf)), rng.sample(gn_uni, len(brca_gn)),
            rng.sample(mi_uni, len(brca_mi)))
    say(f'built {N_CTRL} CTRL_random networks')
    # ---- how much of the genome-wide annotation top the BRCA node set already occupies ----
    annrows = []
    def annrank(label, universe, degfn, brcaset, N):
        ranked = sorted(universe, key=lambda x: -degfn(x))
        ov = sum(1 for x in ranked[:N] if x in brcaset)
        annrows.append(dict(stratum=label, top_N=N, n_in_BRCA=ov,
                            pct_in_BRCA=round(100.0 * ov / N, 1), universe=len(universe),
                            BRCA_set_size=len(brcaset),
                            BRCA_median_degree=float(np.median([degfn(x) for x in brcaset])),
                            nonBRCA_median_degree=float(np.median(
                                [degfn(x) for x in universe if x not in brcaset]))))
        say(f'   {label}: {ov}/{N} ({100.0*ov/N:.1f}%) of the most-annotated are BRCA nodes; '
            f'median degree BRCA={annrows[-1]["BRCA_median_degree"]:.0f} vs '
            f'non-BRCA={annrows[-1]["nonBRCA_median_degree"]:.0f}')
    say('ANNOTATION-RANK OCCUPANCY of the BRCA node set:')
    annrank('TF by TRRUST out-degree', list(TRR), lambda x: len(TRR.get(x, ())), set(brca_tf), 157)
    annrank('TF by TransmiR out-degree', list(TMR), lambda x: len(TMR.get(x, ())), set(brca_tf), 157)
    annrank('miRNA by validated-target count', list(VAL), lambda x: len(VAL.get(x, ())), set(brca_mi), 223)
    annrank('miRNA by TransmiR in-degree', list(TMR_in), lambda x: TMR_in[x], set(brca_mi), 223)
    annrank('protein-coding by validated-miRNA in-degree', sorted(queried),
            lambda x: VAL_in[x], set(brca_gn) | set(brca_tf), 364)
    with open(f'{REV}/results/motif_specificity_annotation_rank.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(annrows[0].keys())); w.writeheader()
        for r in annrows: w.writerow(r)
    say('wrote results/motif_specificity_annotation_rank.csv')

    # ---- CTRL_topann : the most heavily annotated non-BRCA nodes available ----
    def toppool(universe, degfn, k, mult=2):
        return sorted(universe, key=lambda x: (-degfn(x), x))[:min(len(universe), mult * k)]
    tf_top = toppool(tf_uni, lambda x: len(TRR.get(x, ())) + len(TMR.get(x, ())), len(brca_tf))
    gn_top = toppool(gn_uni, lambda x: VAL_in[x] + TRR_in[x], len(brca_gn))
    mi_top = toppool(mi_uni, lambda x: len(VAL.get(x, ())) + TMR_in[x], len(brca_mi))
    say(f'CTRL_topann pools: TF={len(tf_top)} Gene={len(gn_top)} miRNA={len(mi_top)}')
    rng3 = random.Random(SEED + 13)
    for i in range(1, N_CTRL + 1):
        add(f'CTRL_topann_{i:02d}', 'CTRL_topann',
            rng3.sample(tf_top, len(brca_tf)), rng3.sample(gn_top, len(brca_gn)),
            rng3.sample(mi_top, len(brca_mi)))
    say(f'built {N_CTRL} CTRL_topann networks')

    rng2 = random.Random(SEED + 7)
    for i in range(1, N_CTRL + 1):
        add(f'CTRL_matched_{i:02d}', 'CTRL_matched',
            matched_draw(brca_tf, tf_uni, feat_tf, rng2),
            matched_draw(brca_gn, gn_uni, feat_gn, rng2),
            matched_draw(brca_mi, mi_uni, feat_mi, rng2))
    say(f'built {N_CTRL} CTRL_matched networks')

    # ---- edge-count-matched BRCA subsamples: the same BRCA nodes and the same edge
    #      sources, randomly thinned to the mean edge count of the matched controls, so that
    #      "more FFLs" cannot be explained by "more edges" alone.
    tgt = int(round(np.mean([r['n_edges'] for r in rows
                             if r['network_class'] == 'CTRL_matched'])))
    brca_edges = build_network(brca_tf, brca_gn, brca_mi, TRR, TMR, VAL)
    say(f'edge-count-matched BRCA subsamples: thinning {len(brca_edges)} -> {tgt} edges')
    rng4 = random.Random(SEED + 29)
    for i in range(1, N_CTRL + 1):
        add(f'BRCA_edgematched_{i:02d}', 'BRCA_edgematched', brca_tf, brca_gn, brca_mi,
            edges=rng4.sample(brca_edges, tgt))
    say(f'built {N_CTRL} BRCA_edgematched networks')

    fields = ['network', 'n_nodes', 'n_edges', 'n_FFL', 'FFL_per_node', 'FFL_per_edge',
              'network_class', 'n_nodes_connected', 'n_miRNA_FFL', 'n_TF_FFL',
              'n_Composite_FFL', 'n_TF_TF_FFL', 'n_TF', 'n_Gene', 'n_miRNA']
    with open(f'{REV}/results/motif_specificity.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
        for r in rows: w.writerow(r)
    say(f'wrote results/motif_specificity.csv ({len(rows)} rows)')
    json.dump(nodesets, open(f'{REV}/cache/specificity_nodesets.json', 'w'))
    say(f'wrote cache/specificity_nodesets.json ({len(nodesets)} networks)')

    # ---------- statistics ----------
    stats = []
    for cls in ('CTRL_random', 'CTRL_matched', 'CTRL_topann'):
        sub = [r for r in rows if r['network_class'] == cls]
        for metric in ('n_edges', 'n_FFL', 'FFL_per_node', 'FFL_per_edge'):
            v = np.array([r[metric] for r in sub], float)
            refs = [r1, r0]
            for ref in refs:
                x = float(ref[metric])
                p_over = (1 + int((v >= x).sum())) / (1 + len(v))
                z = (x - v.mean()) / v.std(ddof=1) if v.std(ddof=1) > 0 else float('nan')
                stats.append(dict(reference=ref['network'], control_class=cls, metric=metric,
                                  ref_value=round(x, 6), ctrl_mean=round(float(v.mean()), 6),
                                  ctrl_sd=round(float(v.std(ddof=1)), 6),
                                  ctrl_min=round(float(v.min()), 6),
                                  ctrl_max=round(float(v.max()), 6),
                                  Z=round(z, 3), p_emp_over=p_over, n_ctrl=len(v)))
    # ---- the density-controlled comparison: BRCA thinned to the control edge count ----
    from scipy.stats import mannwhitneyu
    em = [r for r in rows if r['network_class'] == 'BRCA_edgematched']
    for cls in ('CTRL_matched', 'CTRL_topann', 'CTRL_random'):
        sub = [r for r in rows if r['network_class'] == cls]
        for metric in ('n_edges', 'n_FFL', 'FFL_per_node', 'FFL_per_edge'):
            a = np.array([r[metric] for r in em], float)
            b = np.array([r[metric] for r in sub], float)
            try:
                u, pmw = mannwhitneyu(a, b, alternative='greater')
            except ValueError:
                pmw = float('nan')
            stats.append(dict(reference='BRCA_edgematched', control_class=cls, metric=metric,
                              ref_value=round(float(a.mean()), 6),
                              ctrl_mean=round(float(b.mean()), 6),
                              ctrl_sd=round(float(b.std(ddof=1)), 6),
                              ctrl_min=round(float(b.min()), 6),
                              ctrl_max=round(float(b.max()), 6),
                              Z=(round(float((a.mean() - b.mean()) / b.std(ddof=1)), 3)
                                 if b.std(ddof=1) > 0 else 'NA'),
                              p_emp_over=round(float(pmw), 6), n_ctrl=len(b)))
            if metric == 'n_FFL':
                say(f'   DENSITY-CONTROLLED  BRCA thinned to {int(a.mean()) and int(np.mean([r["n_edges"] for r in em]))} '
                    f'edges: {a.mean():.0f} +-{a.std(ddof=1):.0f} FFL cores vs {cls} '
                    f'{b.mean():.0f} +-{b.std(ddof=1):.0f}  ratio={a.mean()/b.mean():.2f}x  '
                    f'Mann-Whitney p={pmw:.3g}')

    with open(f'{REV}/results/motif_specificity_stats.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(stats[0].keys())); w.writeheader()
        for r in stats: w.writerow(r)
    say('wrote results/motif_specificity_stats.csv')
    for s in stats:
        if s['reference'] != 'BRCA_rebuilt': continue
        say(f"   {s['control_class']:<13} {s['metric']:<14} BRCA_rebuilt={s['ref_value']:<12} "
            f"ctrl={s['ctrl_mean']:.4g} +-{s['ctrl_sd']:.4g} "
            f"[{s['ctrl_min']:.4g},{s['ctrl_max']:.4g}]  Z={s['Z']}  p={s['p_emp_over']:.4g}")
    say('=== 10 DONE ===')

def readlines(p):
    with open(p) as fh:
        return [x.strip() for x in fh if x.strip()]

if __name__ == '__main__':
    main()
