"""Assemble the portal cross-check tables.
   E1: 16 protein-coding genes  - our TCGA value vs KM Plotter / UALCAN / bc-GenExMiner / GEPIA2.
   E2: miRNAs                   - our TCGA value vs KM Plotter breast-miRNA module.
   Every cell is copied from a portals_*.csv produced earlier in this run; nothing is invented."""
import csv, re, os, math, time

RES = '/path/to/revision/results/v4/'
def load(f):
    p = RES + f
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []

GENES = ["PTEN", "SMAD4", "TGFBR2", "DICER1", "KLF4", "XIAP", "COL1A1", "COL3A1",
         "NFKB1", "RELA", "SP1", "ETS1", "VEGFA", "CCND2", "MYC", "EZH2"]

ours   = {r['gene']: r for r in load('portals_our_tcga_reference_16genes.csv')}
kmg    = load('portals_kmplot_mrna_breast.csv')
ualx   = load('portals_ualcan_expression.csv')
ualm   = load('portals_ualcan_methylation.csv')
ualp   = load('portals_ualcan_cptac_protein.csv')
uals   = load('portals_ualcan_stats.csv')
bcde   = load('portals_bcgenexminer_diffexpr.csv')
bcpr   = load('portals_bcgenexminer_prognosis.csv')
g2e    = load('portals_gepia2_expression.csv')
g2s    = load('portals_gepia2_survival.csv')

def ual_pair(rows, gene):
    d = {r['category']: r for r in rows if r['gene'] == gene and r['grouping'] == 'Sample types'}
    if 'Normal' in d and 'Primary tumor' in d:
        return float(d['Normal']['median']), float(d['Primary tumor']['median']), \
               d['Normal']['n'], d['Primary tumor']['n']
    return None, None, None, None

def ual_p(gene, assay):
    for s in uals:
        if s['gene'] == gene and s['assay'] == assay and s['comparison'] == 'Normal-vs-Primary':
            return s['p_value']
    return ''

def km(gene, endpoint):
    for r in kmg:
        if (r['gene'] == gene and r['endpoint'] == endpoint and r['restriction'] == 'all'
                and r['cutoff_method'] == 'median'):
            return r
    return {}

def bcgem_pr(gene, ep):
    for r in bcpr:
        if (r['gene'] == gene and r['endpoint'] == ep and r['population_ER'] == 'ER all'
                and r['population_PR'] == 'PR all' and r['population_HER2'] == 'HER2 all'):
            return r
    return {}

def bcgem_de(gene):
    d = {r['comparison_group']: r for r in bcde if r['gene'] == gene}
    if 'Healthy' in d and 'Tumour' in d:
        h, t = float(d['Healthy']['median']), float(d['Tumour']['median'])
        a = float(d['Tumour-adjacent']['median']) if 'Tumour-adjacent' in d else None
        return h, t, d['Tumour'].get('DTK_p', ''), a
    return None, None, '', None

def g2expr(gene):
    for r in g2e:
        if r['gene'] == gene and r['cancer'] == 'BRCA':
            return r
    return {}

def g2surv(gene, ep):
    for r in g2s:
        if r['gene'] == gene and r['endpoint'] == ep and r['split'] == 'median':
            return r
    return {}

def dirn(a, b):                       # direction of b (tumour) relative to a (normal)
    if a is None or b is None: return ''
    return "UP in tumour" if b > a else ("DOWN in tumour" if b < a else "no change")

def pnum(p):
    """kmplot/bc-GenExMiner print '<1e-16', '< 0.0001', '> 0.10'; return a usable float."""
    if p is None: return None
    p = str(p).strip().replace(' ', '')
    if not p: return None
    m = re.match(r'^([<>])?(=)?([0-9.]+(?:[eE][+\-]?[0-9]+)?)$', p)
    if not m: return None
    v = float(m.group(3))
    if m.group(1) == '>': return v * 1.0000001      # ">0.10" -> treat as 0.10, i.e. not sig
    return v

def hrdir(hr, p=None, alpha=0.05):
    """Direction only counts when the association is significant; otherwise it is 'null'."""
    try: hr = float(hr)
    except Exception: return ''
    pv = pnum(p)
    if pv is not None and pv >= alpha:
        return "null (p=%s)" % str(p).strip()
    return "high = worse" if hr > 1 else ("high = better" if hr < 1 else "null")

F = ["feature", "class",
     # expression direction, four sources
     "our_TCGA_direction", "our_TCGA_log2FC", "our_TCGA_FDR",
     "our_TCGA_median_tumour", "our_TCGA_median_normal",
     "UALCAN_direction", "UALCAN_median_normal_TPM", "UALCAN_median_tumour_TPM", "UALCAN_p",
     "bcGenExMiner_direction", "bcGenExMiner_median_healthy",
     "bcGenExMiner_median_tumour_adjacent", "bcGenExMiner_median_tumour",
     "bcGenExMiner_direction_vs_tumour_adjacent", "bcGenExMiner_DTK_p",
     "GEPIA2_direction", "GEPIA2_median_normal_TPM", "GEPIA2_median_tumour_TPM", "GEPIA2_log2FC",
     "EXPRESSION_AGREEMENT",
     # survival, four sources
     "our_TCGA_HRperSD_OS", "our_TCGA_p_OS", "our_TCGA_HRperSD_PFI", "our_TCGA_p_PFI",
     "KMplotter_probe", "KMplotter_RFS_HR", "KMplotter_RFS_CI", "KMplotter_RFS_p", "KMplotter_RFS_n",
     "KMplotter_OS_HR", "KMplotter_OS_CI", "KMplotter_OS_p", "KMplotter_OS_n",
     "bcGenExMiner_DFS_HR", "bcGenExMiner_DFS_CI", "bcGenExMiner_DFS_p", "bcGenExMiner_DFS_n",
     "bcGenExMiner_OS_HR", "bcGenExMiner_OS_CI", "bcGenExMiner_OS_p", "bcGenExMiner_OS_n",
     "GEPIA2_OS_HRhigh", "GEPIA2_OS_logrank_p", "GEPIA2_DFS_HRhigh", "GEPIA2_DFS_logrank_p",
     "SURVIVAL_DIRECTION_AGREEMENT", "normal_reference_note",
     # extras
     "UALCAN_promoter_meth_normal", "UALCAN_promoter_meth_tumour", "UALCAN_promoter_meth_p",
     "UALCAN_CPTAC_protein_normal_Z", "UALCAN_CPTAC_protein_tumour_Z", "UALCAN_CPTAC_p",
     "FLAG"]

rows = []
for g in GENES:
    o = ours.get(g, {})
    un, ut, unn, utn = ual_pair(ualx, g)
    mn, mt, _, _ = ual_pair(ualm, g)
    pn, pt, pnn, ptn = ual_pair(ualp, g)
    bh, bt, bp, ba = bcgem_de(g)
    ge = g2expr(g)
    kR, kO = km(g, 'RFS'), km(g, 'OS')
    bD, bO = bcgem_pr(g, 'DFS'), bcgem_pr(g, 'OS')
    sO, sD = g2surv(g, 'OS'), g2surv(g, 'DFS/RFS')

    exd = [o.get('our_direction', ''), dirn(un, ut), dirn(bh, bt), ge.get('direction', '')]
    exd = [x for x in exd if x]
    agree_ex = "ALL AGREE (%s)" % exd[0] if len(set(exd)) == 1 else \
               "DISAGREE: " + " | ".join("%s=%s" % (s, d) for s, d in
                    zip(["ours", "UALCAN", "bc-GenExMiner", "GEPIA2"], exd))
    srcs = [("ours TCGA OS", o.get('our_HRperSD_OS'), o.get('our_p_OS')),
            ("ours TCGA PFI", o.get('our_HRperSD_PFI'), o.get('our_p_PFI')),
            ("KM Plotter RFS", kR.get('HR'), kR.get('logrank_p')),
            ("KM Plotter OS", kO.get('HR'), kO.get('logrank_p')),
            ("bc-GenExMiner DFS", bD.get('HR'), bD.get('p_value')),
            ("bc-GenExMiner OS", bO.get('HR'), bO.get('p_value')),
            ("GEPIA2 OS", sO.get('HR_high'), sO.get('logrank_p')),
            ("GEPIA2 DFS", sD.get('HR_high'), sD.get('logrank_p'))]
    calls = [(nm, hrdir(h, pp)) for nm, h, pp in srcs if h]
    sig = [(nm, d) for nm, d in calls if d.startswith("high")]
    nulls = [nm for nm, d in calls if not d.startswith("high")]
    if not sig:
        agree_sv = "no significant association in any source (%d analyses, all p>=0.05)" % len(calls)
    elif len(set(d for _, d in sig)) == 1:
        agree_sv = ("CONCORDANT among significant: %s in %s"
                    % (sig[0][1], ", ".join(nm for nm, _ in sig))
                    + ("; null in " + ", ".join(nulls) if nulls else ""))
    else:
        agree_sv = ("CONFLICT among significant: "
                    + " | ".join("%s=%s" % (nm, d) for nm, d in sig)
                    + ("; null in " + ", ".join(nulls) if nulls else ""))
    flags = []
    if "DISAGREE" in agree_ex:
        tcga_based = {o.get('our_direction', ''), dirn(un, ut)} - {''}
        gtex_based = {dirn(bh, bt), ge.get('direction', '')} - {''}
        if (len(tcga_based) == 1 and len(gtex_based) == 1 and tcga_based != gtex_based):
            msg = ("expression direction splits exactly on the NORMAL REFERENCE: "
                   "TCGA solid-tissue normal (ours + UALCAN) says %s, GTEx healthy breast "
                   "(bc-GenExMiner GTEx&TCGA + GEPIA2) says %s. GTEx breast is largely "
                   "adipose/stroma, so this is a tissue-composition artefact of the reference, "
                   "not a contradiction about tumours"
                   % (list(tcga_based)[0], list(gtex_based)[0]))
            if ba is not None and dirn(ba, bt) == list(tcga_based)[0]:
                msg += ("; RESOLVED - bc-GenExMiner's own tumour-ADJACENT arm (median %s vs "
                        "tumour %s) reproduces the TCGA-normal answer (%s), so all four sources "
                        "agree once the reference tissue is matched" % (ba, bt, dirn(ba, bt)))
            flags.append(msg)
        else:
            small = []
            try:
                if abs(float(o.get('our_logFC_TvsN', 0))) < 0.15 or \
                        float(o.get('our_FDR_TvsN', 1)) > 0.05:
                    small.append("ours (log2FC %s, FDR %s)"
                                 % (o.get('our_logFC_TvsN'), o.get('our_FDR_TvsN')))
            except Exception:
                pass
            try:
                if abs(float(ge.get('log2FC_tumour_vs_normal', 0))) < 0.15:
                    small.append("GEPIA2 (log2FC %s)" % ge.get('log2FC_tumour_vs_normal'))
            except Exception:
                pass
            if small:
                flags.append("expression direction differs between portals, but the dissenting "
                             "source(s) report a negligible / non-significant fold change: "
                             + "; ".join(small) + " - treat as 'no reliable change', not as a "
                             "contradiction")
            else:
                flags.append("expression direction differs between portals")
    if agree_sv.startswith("CONFLICT"):
        flags.append("survival direction conflicts between significant sources")
    if pt is not None and abs(pt) < 0.06:
        flags.append("UALCAN CPTAC tumour median Z is ~0 by construction (Z is centred on the "
                     "full cancer-type cohort, which is 125 tumours vs 18 normals) - the contrast "
                     "is really 'where the 18 normals sit', not a fold change")
    if not ge:
        flags.append("no GEPIA2 value retrieved")
    if pn is None:
        flags.append("gene not quantified in the CPTAC breast proteome")

    rows.append(dict(feature=g, **{"class": "protein-coding"},
        our_TCGA_direction=o.get('our_direction', ''), our_TCGA_log2FC=o.get('our_logFC_TvsN', ''),
        our_TCGA_FDR=o.get('our_FDR_TvsN', ''),
        our_TCGA_median_tumour=o.get('our_median_tumour', ''),
        our_TCGA_median_normal=o.get('our_median_normal', ''),
        UALCAN_direction=dirn(un, ut), UALCAN_median_normal_TPM=un,
        UALCAN_median_tumour_TPM=ut, UALCAN_p=ual_p(g, 'expression'),
        bcGenExMiner_direction=dirn(bh, bt), bcGenExMiner_median_healthy=bh,
        bcGenExMiner_median_tumour=bt, bcGenExMiner_median_tumour_adjacent=ba,
        bcGenExMiner_direction_vs_tumour_adjacent=dirn(ba, bt), bcGenExMiner_DTK_p=bp,
        GEPIA2_direction=ge.get('direction', ''),
        GEPIA2_median_normal_TPM=ge.get('normal_median_TPM', ''),
        GEPIA2_median_tumour_TPM=ge.get('tumour_median_TPM', ''),
        GEPIA2_log2FC=ge.get('log2FC_tumour_vs_normal', ''),
        EXPRESSION_AGREEMENT=agree_ex,
        our_TCGA_HRperSD_OS=o.get('our_HRperSD_OS', ''), our_TCGA_p_OS=o.get('our_p_OS', ''),
        our_TCGA_HRperSD_PFI=o.get('our_HRperSD_PFI', ''), our_TCGA_p_PFI=o.get('our_p_PFI', ''),
        KMplotter_probe=kR.get('probe', ''), KMplotter_RFS_HR=kR.get('HR', ''),
        KMplotter_RFS_CI="%s-%s" % (kR.get('CI_low', ''), kR.get('CI_high', '')),
        KMplotter_RFS_p=kR.get('logrank_p', ''), KMplotter_RFS_n=kR.get('n_total', ''),
        KMplotter_OS_HR=kO.get('HR', ''),
        KMplotter_OS_CI="%s-%s" % (kO.get('CI_low', ''), kO.get('CI_high', '')),
        KMplotter_OS_p=kO.get('logrank_p', ''), KMplotter_OS_n=kO.get('n_total', ''),
        bcGenExMiner_DFS_HR=bD.get('HR', ''), bcGenExMiner_DFS_CI=bD.get('CI_95', ''),
        bcGenExMiner_DFS_p=bD.get('p_value', ''), bcGenExMiner_DFS_n=bD.get('n_patients', ''),
        bcGenExMiner_OS_HR=bO.get('HR', ''), bcGenExMiner_OS_CI=bO.get('CI_95', ''),
        bcGenExMiner_OS_p=bO.get('p_value', ''), bcGenExMiner_OS_n=bO.get('n_patients', ''),
        GEPIA2_OS_HRhigh=sO.get('HR_high', ''), GEPIA2_OS_logrank_p=sO.get('logrank_p', ''),
        GEPIA2_DFS_HRhigh=sD.get('HR_high', ''), GEPIA2_DFS_logrank_p=sD.get('logrank_p', ''),
        SURVIVAL_DIRECTION_AGREEMENT=agree_sv,
        normal_reference_note=("ours + UALCAN compare against 114 TCGA-BRCA solid-tissue "
                               "normals; bc-GenExMiner (GTEx & TCGA source) and GEPIA2 compare "
                               "against GTEx healthy breast, which is adipose-rich"),
        UALCAN_promoter_meth_normal=mn, UALCAN_promoter_meth_tumour=mt,
        UALCAN_promoter_meth_p=ual_p(g, 'promoter_methylation'),
        UALCAN_CPTAC_protein_normal_Z=pn, UALCAN_CPTAC_protein_tumour_Z=pt,
        UALCAN_CPTAC_p=ual_p(g, 'protein_CPTAC'),
        FLAG="; ".join(flags)))

with open(RES + 'portals_CROSSCHECK_genes.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=F); w.writeheader(); w.writerows(rows)
print("wrote portals_CROSSCHECK_genes.csv", len(rows), "rows")

# ---------------- miRNA cross-check ----------------
mirref = load('portals_our_tcga_mirna_reference.csv')
kmm = load('portals_kmplot_mirna_breast.csv')
MIRS = [("hsa-miR-130a", "hsa-miR-130a-3p"), ("hsa-miR-130a_5p", "hsa-miR-130a-5p"), ("hsa-miR-29a", "hsa-miR-29a-3p"),
        ("hsa-miR-29b", "hsa-miR-29b-3p"), ("hsa-miR-29c", "hsa-miR-29c-3p"),
        ("hsa-let-7b", "hsa-let-7b-5p"), ("hsa-let-7e", "hsa-let-7e-5p"),
        ("hsa-miR-130b", "hsa-miR-130b-3p"), ("hsa-miR-301a", "hsa-miR-301a-3p")]
# arm-level DE (recomputed: results/BRCA_DEX_mirnas.csv is keyed on precursor names)
armde = {r['feature']: r for r in load('portals_our_tcga_mirna_arm_DE.csv')}
# precursor-level DE as shipped in the manuscript's own table, for reference
dexm = {r['feature']: r for r in csv.DictReader(
    open('/path/to/revision/results/BRCA_DEX_mirnas.csv'))}
MF = ["feature", "kmplot_id", "our_arm_id",
      "our_TCGA_logFC_TvsN", "our_TCGA_FDR_TvsN", "our_TCGA_direction",
      "our_TCGA_logFC_precursor_level", "our_TCGA_FDR_precursor_level",
      "our_TCGA_HRperSD_OS", "our_TCGA_p_OS", "our_TCGA_HRperSD_PFI", "our_TCGA_p_PFI",
      "our_TCGA_HRperSD_DSS", "our_TCGA_p_DSS"]
for ds in ["METABRIC", "TCGA", "GSE40267", "GSE19783"]:
    for cut in ["median", "auto_percentile"]:
        MF += ["KM_%s_%s_HR" % (ds, cut), "KM_%s_%s_CI" % (ds, cut),
               "KM_%s_%s_p" % (ds, cut), "KM_%s_%s_n" % (ds, cut),
               "KM_%s_%s_FDR" % (ds, cut)]
MF += ["KM_direction_consistency", "NOTE"]

mrows = []
for kid, arm in MIRS:
    r = dict(feature=kid, kmplot_id=kid, our_arm_id=arm or "")
    armkey = arm or (kid + "-3p")
    a = armde.get(armkey, {})
    d = dexm.get(kid, {})
    r["our_arm_id"] = armkey
    r["our_TCGA_logFC_TvsN"] = round(float(a['logFC']), 4) if a else ''
    r["our_TCGA_FDR_TvsN"] = float(a['adj.P.Val']) if a else ''
    r["our_TCGA_direction"] = a.get('direction', '')
    r["our_TCGA_logFC_precursor_level"] = round(float(d['logFC']), 4) if d else ''
    r["our_TCGA_FDR_precursor_level"] = float(d['adj.P.Val']) if d else ''
    for ep in ("OS", "PFI", "DSS"):
        m = [x for x in mirref if x['feature'] == armkey and x['endpoint'] == ep]
        r["our_TCGA_HRperSD_%s" % ep] = round(float(m[0]['HR_per_SD']), 4) if m else ''
        r["our_TCGA_p_%s" % ep] = round(float(m[0]['p_value']), 5) if m else ''
    hrs = []
    for ds in ["METABRIC", "TCGA", "GSE40267", "GSE19783"]:
        for cut in ["median", "auto_percentile"]:
            k = [x for x in kmm if x['mirna'] == kid.replace("_5p", "") and x['dataset'] == ds
                 and x['cutoff_method'] == cut and x['restriction'] == 'all']
            k = k[0] if k else {}
            r["KM_%s_%s_HR" % (ds, cut)] = k.get('HR', '')
            r["KM_%s_%s_CI" % (ds, cut)] = "%s-%s" % (k.get('CI_low', ''), k.get('CI_high', ''))
            r["KM_%s_%s_p" % (ds, cut)] = k.get('logrank_p', '') or k.get('error', '')
            r["KM_%s_%s_n" % (ds, cut)] = k.get('n_total', '')
            r["KM_%s_%s_FDR" % (ds, cut)] = k.get('FDR', '')
            if k.get('HR'):
                hrs.append(float(k['HR']))
    if hrs:
        lo = sum(1 for h in hrs if h < 1)
        r["KM_direction_consistency"] = ("%d/%d KM Plotter analyses give HR<1 "
                                         "(higher miRNA = better OS)" % (lo, len(hrs)))
    else:
        r["KM_direction_consistency"] = ""
    r["NOTE"] = ("KM Plotter's breast miRNA module offers OS only (no RFS/DMFS). Its TCGA "
                 "dataset is the same TCGA-BRCA cohort as our own analysis, so only METABRIC, "
                 "GSE40267 and GSE19783 are independent of us.")
    mrows.append(r)

with open(RES + 'portals_CROSSCHECK_mirnas.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=MF); w.writeheader(); w.writerows(mrows)
print("wrote portals_CROSSCHECK_mirnas.csv", len(mrows), "rows")
