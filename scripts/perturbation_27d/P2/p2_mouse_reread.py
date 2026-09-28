#!/usr/bin/env python3
"""P2 mouse, second reading (SETTINGS change 13): the mouse cluster perturbations found on re-reading the P2 hits.
RNA-seq series (no NCBI-generated counts for any of them): the submitters' count or normalised tables, built here into
GSM-labelled matrices and run through _de.R (raw counts: filterByExpr, TMM, voom; normalised values: log2(x + 1), limma,
flagged 'processed').  Array series: written to P2/p2_mouse_arms.tsv and run by P1/p1_run_dataset.py --mouse
(CEL: RMA; Agilent: normexp + quantile from the Feature Extraction files; 3D-Gene: series matrix, processed).
Mouse genes -> human symbols by one-to-one orthologues (_mouse.py).  One contrast per series and tissue / cell type:
the whole-cluster perturbation where a series also has partial deletions.  Writes P2/de/<name>.tsv + .meta.json (cls
'mouse'), P2/p2_mouse_arms.tsv and P2/batch_mouse_arrays.tsv."""
import os, sys, io, re, gzip, json, tarfile, subprocess
import numpy as np, pandas as pd
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); REV = os.path.dirname(INB); sys.path.insert(0, INB)
import _dl, _geo, _mouse
OUT = f'{INB}/P2/de'; W = os.path.join(_dl.CACHE, 'P2/work'); os.makedirs(W, exist_ok=True)
L = pd.read_csv(f'{REV}/data/layer_miRNA_miRNA.tsv', sep='\t')
C420 = sorted(set(L[L.cluster_id == 'cl10kb_420'].miRNA_1) | set(L[L.cluster_id == 'cl10kb_420'].miRNA_2))
MEM = {'miR-17~92': 'hsa-miR-17;hsa-miR-18a;hsa-miR-19a;hsa-miR-19b;hsa-miR-20a;hsa-miR-92a',
       'miR-183~96~182': 'hsa-miR-183;hsa-miR-96;hsa-miR-182', 'DLK1-DIO3 (miR-379/410)': ';'.join(C420),
       'miR-19a/b': 'hsa-miR-19a;hsa-miR-19b',
       'miR-17~92 + miR-106a~363 + miR-106b~25': 'hsa-miR-17;hsa-miR-18a;hsa-miR-19a;hsa-miR-19b;hsa-miR-20a;hsa-miR-92a;hsa-miR-106a;hsa-miR-18b;'
                                                'hsa-miR-20b;hsa-miR-363;hsa-miR-106b;hsa-miR-93;hsa-miR-25',
       'miR-200b~200a~429': 'hsa-miR-200b;hsa-miR-200a;hsa-miR-429', 'miR-200c~141': 'hsa-miR-200c;hsa-miR-141',
       'both miR-200 clusters': 'hsa-miR-200c;hsa-miR-141;hsa-miR-200b;hsa-miR-200a;hsa-miR-429',
       'miR-23~27~24 (both clusters)': 'hsa-miR-23a;hsa-miR-27a;hsa-miR-24;hsa-miR-23b;hsa-miR-27b'}
def run_de(name, M, samples, kind, meta):
    """M: DataFrame indexed by mouse gene id (Ensembl, Entrez or symbol), columns = GSMs; samples: gsm -> arm."""
    M = M[[g for g in samples]]
    m = f'{W}/{name}.matrix.tsv'; M.rename_axis('id').reset_index().to_csv(m, sep='\t', index=False)
    ids = M.index.astype(str); hs = _mouse.to_human(ids)
    gm = f'{W}/{name}.genemap.tsv'; pd.DataFrame({'id': ids, 'symbol': hs})[lambda d: d.symbol != ''].to_csv(gm, sep='\t', index=False)
    sp = f'{W}/{name}.spec.json'; json.dump(dict(kind=kind, matrix=m, gene_map=gm, samples=samples), open(sp, 'w'))
    r = subprocess.run(['Rscript', f'{INB}/_de.R', sp, f'{OUT}/{name}.tsv'], capture_output=True, text=True)
    if r.returncode: print('FAILED', name, r.stderr[-400:]); return
    meta.update(de=r.stdout.strip(), n_pert=sum(v == 'pert' for v in samples.values()), n_ctrl=sum(v == 'ctrl' for v in samples.values()),
                cls='mouse', tag='primary', note=meta.get('note', '') + '; mouse -> human one-to-one orthologues (MGI)')
    json.dump(meta, open(f'{OUT}/{name}.meta.json', 'w')); print(name, meta['de'], meta['note'][:60], flush=True)
def arms(S, pert, ctrl):
    out = {}
    for s in S:
        if re.search(pert, s['title']): out[s['gsm']] = 'pert'
        elif re.search(ctrl, s['title']): out[s['gsm']] = 'ctrl'
    return out
def per_sample_tar(gse, reader):
    tf = tarfile.open(_geo.suppl_fetch(gse, f'{gse}_RAW.tar', 'P2')); cols = {}
    for mb in tf.getmembers():
        raw = tf.extractfile(mb).read(); raw = gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw
        cols[mb.name.split('_')[0]] = reader(raw)
    return pd.DataFrame(cols)
base = lambda gse, tf, cell, direction, route, note: dict(gse=gse, tf=tf, cell=cell, direction=direction, mirna=MEM[tf], route=route, platform='', note=note)

# GSE63813: miR-17~92 full knockout vs wild type, E9.5, per tissue (partial-deletion arms of the series not used)
S = _geo.samples('GSE63813', 'P2'); src = {s['description'].split('|')[-1].strip(): s['gsm'] for s in S}   # library code (e.g. MA1) = count-table column
T = pd.read_csv(_geo.suppl_fetch('GSE63813', 'GSE63813_countTable.txt.gz', 'P2'), sep='\t', index_col=0).rename(columns=src)
for code, tis in (('TB', 'tail bud'), ('H', 'heart'), ('E', 'whole embryo without heart and tail bud')):
    run_de(f'GSE63813__miR-17_92__mouse_E9.5_{code}__primary', T, arms(S, rf'^Full-KO_{code}_\d', rf'^Wild-type_{code}_\d'), 'counts',
           base('GSE63813', 'miR-17~92', f'mouse E9.5 {tis} (miR-17~92 KO)', 'loss', 'submitter counts', 'submitter count table (Ensembl ids)'))
# GSE113162: RJ423 miR-200b/200a/429 vs empty vector
S = _geo.samples('GSE113162', 'P2')
T = pd.read_csv(_geo.suppl_fetch('GSE113162', 'GSE113162_RJ423EV_RJ423_miR200ab429_rawCountMatrix.csv.gz', 'P2'), index_col=0).drop(columns=['Symbol'])
ev = [s['gsm'] for s in S if 'EV' in s['title']]; mi = [s['gsm'] for s in S if 'miR200' in s['title']]
T = T.rename(columns=dict(zip(['RJ423EV1', 'RJ423EV2', 'RJ423EV3'], ev)) | dict(zip(['RJ423200ab429A', 'RJ423200ab429B', 'RJ423200ab429C'], mi)))
run_de('GSE113162__miR-200b_200a_429__mouse_RJ423__primary', T, {**{g: 'ctrl' for g in ev}, **{g: 'pert' for g in mi}}, 'counts',
       base('GSE113162', 'miR-200b~200a~429', 'mouse RJ423 mammary tumour cells', 'gain', 'submitter counts', 'submitter raw count matrix; columns assigned to GSMs by arm'))
# GSE150101: RJ423 miR-200c/141 vs empty vector (per-sample raw counts)
S = _geo.samples('GSE150101', 'P2')
T = per_sample_tar('GSE150101', lambda raw: pd.read_csv(io.BytesIO(raw), sep='\t', index_col=0).iloc[:, 0])
run_de('GSE150101__miR-200c_141__mouse_RJ423__primary', T, arms(S, r'c141', r'EV'), 'counts',
       base('GSE150101', 'miR-200c~141', 'mouse RJ423 mammary tumour cells', 'gain', 'submitter counts', 'submitter per-sample raw counts'))
# GSE180264: mammary glands, MTB-200ba429 vs control; and on the MTB-IGFIR background
S = _geo.samples('GSE180264', 'P2')
fc = lambda raw: pd.read_csv(io.BytesIO(raw), sep='\t', index_col=0).RAW_COUNT
T = per_sample_tar('GSE180264', fc)
run_de('GSE180264__miR-200b_200a_429__mouse_mammary__primary', T, arms(S, r'^MTB-200ba429', r'^Control\d'), 'counts',
       base('GSE180264', 'miR-200b~200a~429', 'mouse mammary gland (MTB-200ba429)', 'gain', 'submitter counts', 'submitter per-sample featureCounts'))
run_de('GSE180264__miR-200b_200a_429__mouse_mammary_IGFIR__primary', T, arms(S, r'^MTB-IGFIRba429', r'^MTB-IGFIR\d'), 'counts',
       base('GSE180264', 'miR-200b~200a~429', 'mouse mammary gland, MTB-IGFIR background', 'gain', 'submitter counts', 'submitter per-sample featureCounts'))
# GSE270979: mammary glands (55 d), MTB-TANba429 vs MTB-TAN
S = _geo.samples('GSE270979', 'P2'); T = per_sample_tar('GSE270979', fc)
run_de('GSE270979__miR-200b_200a_429__mouse_mammary_Neu__primary', T, arms(S, r'^MTB-TANba429 55D', r'^MTB-TAN 55D'), 'counts',
       base('GSE270979', 'miR-200b~200a~429', 'mouse mammary gland, Neu background (55 d)', 'gain', 'submitter counts', 'submitter per-sample featureCounts'))
# GSE250195: ID8 and 28-2 ovarian cancer cells, both miR-200 clusters vs control lentivirus
S = _geo.samples('GSE250195', 'P2')
T = per_sample_tar('GSE250195', lambda raw: pd.read_excel(io.BytesIO(raw)).set_index('FEATURE_ID').RAW_COUNT)
for line in ('ID8', '28-2'):
    run_de(f'GSE250195__miR-200__mouse_{line}__primary', T, arms(S, rf'^{line}_lentivirus containing', rf'^{line}_control'), 'counts',
           base('GSE250195', 'both miR-200 clusters', f'mouse {line} ovarian cancer cells', 'gain', 'submitter counts', 'submitter per-sample featureCounts'))
# GSE268482: C2C12 myoblasts (MB) and myotubes (MT), Dlk1-Dio3 over-expressing lines (2-4, 2-7) vs control (px)
S = _geo.samples('GSE268482', 'P2'); t2g = {s['title']: s['gsm'] for s in S}
T = pd.read_csv(_geo.suppl_fetch('GSE268482', 'GSE268482_counts.txt.gz', 'P2'), sep='\t', index_col=0).rename(columns=t2g)
for st, lab in (('MB', 'myoblasts'), ('MT', 'myotubes')):
    run_de(f'GSE268482__DLK1-DIO3__mouse_C2C12_{st}__primary', T, arms(S, rf'^{st} 2-[47] ', rf'^{st} px '), 'counts',
           base('GSE268482', 'DLK1-DIO3 (miR-379/410)', f'mouse C2C12 {lab} (Dlk1-Dio3 over-expression)', 'gain', 'submitter counts', 'submitter count table'))
# GSE118698: T-cell-specific miR-23~27~24 family KO vs WT, per subset (RPKM, processed; columns named by sample title)
S = _geo.samples('GSE118698', 'P2'); t2g = {s['title']: s['gsm'] for s in S}
T = pd.read_csv(_geo.suppl_fetch('GSE118698', 'GSE118698_miR-23_Tfh_RPKM.txt.gz', 'P2'), sep='\t')
T.index = T['Annotation/Divergence'].astype(str).str.split('|').str[0]
for sub in ('Tn', 'Th1', 'Tfh', 'GCTfh'):
    cols = [c for c in T.columns if re.fullmatch(rf'{sub}_\d_(WT|KO)', c)]
    M = T[cols].rename(columns=t2g).groupby(level=0).max()
    run_de(f'GSE118698__miR-23_27_24__mouse_{sub}__primary', M, {t2g[c]: ('pert' if c.endswith('_KO') else 'ctrl') for c in cols}, 'linear_matrix',
           base('GSE118698', 'miR-23~27~24 (both clusters)', f'mouse CD4+ T cells, {sub}', 'loss', 'processed', 'submitter RPKM (processed)'))
# GSE222351: adipose-tissue macrophages, Mirc11/Mirc22 (miR-23~27~24 clusters) KO vs fl/fl (normalised counts, processed)
S = _geo.samples('GSE222351', 'P2')
T = pd.read_csv(_geo.suppl_fetch('GSE222351', 'GSE222351_Table_S1-Normalized_gene_counts.txt.gz', 'P2'), sep='\t', index_col=0)
g_fl = [s['gsm'] for s in S if s['title'].startswith('ATMs, fl/fl')]; g_ko = [s['gsm'] for s in S if s['title'].startswith('ATMs, KO')]
M = T[[f'fl/fl {i}' for i in range(1, 5)] + [f'KO {i}' for i in range(1, 5)]].copy(); M.columns = g_fl + g_ko
run_de('GSE222351__miR-23_27_24__mouse_ATM__primary', M, {**{g: 'ctrl' for g in g_fl}, **{g: 'pert' for g in g_ko}}, 'linear_matrix',
       base('GSE222351', 'miR-23~27~24 (both clusters)', 'mouse adipose-tissue macrophages (obese)', 'loss', 'processed', 'submitter normalised counts (processed)'))
# GSE126888: heart after myocardial infarction, miR-19a/b mimic vs control mimic (normalised expression, processed; stimulated)
S = _geo.samples('GSE126888', 'P2')
T = pd.read_csv(_geo.suppl_fetch('GSE126888', 'GSE126888_miR19_4dvsCtrl_4d_Expression_Profile.txt.gz', 'P2'), sep='\t', index_col=0)
g_c = [s['gsm'] for s in S if s['title'].startswith('Ctrl')]; g_m = [s['gsm'] for s in S if s['title'].startswith('miR-19')]
M = T[['Ctrl_1', 'Ctrl_2', 'Ctrl_3', 'miR19_1', 'miR19_2', 'miR19_3']].copy(); M.columns = g_c + g_m
run_de('GSE126888__miR-19a_b__mouse_heart_MI__primary', M, {**{g: 'ctrl' for g in g_c}, **{g: 'pert' for g in g_m}}, 'linear_matrix',
       base('GSE126888', 'miR-19a/b', 'mouse heart after myocardial infarction (both arms)', 'gain', 'processed', 'submitter normalised expression (processed); stimulated'))
# GSE110714: E12.5 dorsal root ganglia, miR-183 cluster conditional KO vs controls (processed; decimal commas)
S = _geo.samples('GSE110714', 'P2')
T = pd.read_csv(_geo.suppl_fetch('GSE110714', 'GSE110714_Peng_et_al_expression.csv.gz', 'P2'), dtype=str, index_col=0)
T = T.apply(lambda c: pd.to_numeric(c.str.replace(',', '.'), errors='coerce'))
g_c = [s['gsm'] for s in S if s['title'].startswith('Control')]; g_k = [s['gsm'] for s in S if s['title'].startswith('miR-CKO')]
M = T[['Control 1', 'Control 2', 'Control 3', 'miR CKO1', 'miR CKO2', 'miR CKO3']].copy(); M.columns = g_c + g_k; M = M.groupby(level=0).max()
run_de('GSE110714__miR-183_96_182__mouse_DRG__primary', M, {**{g: 'ctrl' for g in g_c}, **{g: 'pert' for g in g_k}}, 'linear_matrix',
       base('GSE110714', 'miR-183~96~182', 'mouse E12.5 dorsal root ganglia (miR-183 cluster cKO)', 'loss', 'processed', 'submitter expression values (processed)'))
# GSE114495: calvarial cells, mir17 fl/fl;106b-/- with Cre vs YFP (FPKM, processed; the sample number is in the GSM description)
S = _geo.samples('GSE114495', 'P2')
T = pd.read_excel(_geo.suppl_fetch('GSE114495', 'GSE114495_akozlova_GeneExpression.xlsx', 'P2')).set_index('gene_id')
T.columns = [str(c).strip() for c in T.columns]; T.index = T.index.astype(str)
num = {s['gsm']: s['description'].strip() for s in S}
M = T[[f'{num[s["gsm"]]}_FPKM' for s in S]].copy(); M.columns = [s['gsm'] for s in S]
run_de('GSE114495__miR-17_92__mouse_calvaria__primary', M, {s['gsm']: ('pert' if s['title'].startswith('Cre') else 'ctrl') for s in S}, 'linear_matrix',
       base('GSE114495', 'miR-17~92', 'mouse calvarial cells (miR-17~92 deletion on a miR-106b~25-null background)', 'loss', 'processed', 'submitter FPKM (processed)'))

# GSE116682: arcuate nucleus, Pomc-Cre;miR-17~92 fl/fl vs wild type (FPKM sheet of the submitter's table, processed)
S = _geo.samples('GSE116682', 'P2'); t2g = {s['title']: s['gsm'] for s in S}
X = pd.read_excel(_geo.suppl_fetch('GSE116682', 'GSE116682_All_Differentially_Expressed_Genes.xlsx', 'P2'), sheet_name='FPKM(17-92 and 17-92-C)', header=1)
X.index = X['Track_id'].astype(str).str.split('.').str[0]
M = X[['17-92-1', '17-92-2', '17-92-C1', '17-92-C2']].rename(columns=t2g)
run_de('GSE116682__miR-17_92__mouse_arcuate__primary', M, {t2g['17-92-1']: 'pert', t2g['17-92-2']: 'pert', t2g['17-92-C1']: 'ctrl', t2g['17-92-C2']: 'ctrl'},
       'linear_matrix', base('GSE116682', 'miR-17~92', 'mouse arcuate nucleus (POMC-neuron miR-17~92 KO)', 'loss', 'processed', 'submitter FPKM (processed)'))

# ---------------- array series: arms for P1/p1_run_dataset.py --mouse
rows = []
def arm_rows(gse, tf, cell, direction, pert, ctrl, platform=None, flt=None):
    S = [s for s in _geo.samples(gse, 'P2') if (flt is None or flt(s))]
    for s in S:
        a = 'pert' if re.search(pert, s['title']) else ('ctrl' if re.search(ctrl, s['title']) else '')
        if a: rows.append(dict(gse=gse, tf=tf, gsm=s['gsm'], arm=a, title=s['title'], platform=platform or s['platform'], strategy=s['strategy'],
                               cell=cell, cls='mouse', direction=direction, mirna=MEM[tf], organism='mouse'))
arm_rows('GSE30877', 'miR-17~92', 'mouse follicular B cells (miR-17~92 transgene)', 'gain', r'^trans gene', r'^wild type')
arm_rows('GSE56379', 'miR-17~92 + miR-106a~363 + miR-106b~25', 'mouse follicular B cells (triple cluster KO)', 'loss', r'^miR17_92_tKO\d$', r'^Control_\d$',
         flt=lambda s: s['platform'] == 'GPL6246' and 'hr' not in s['title'])
arm_rows('GSE56377', 'miR-17~92', 'mouse follicular B cells, LPS/IL-4 25.5 h (stimulated)', 'gain', r'^miR-17.92_TG\d_25\.5hr', r'^WT\d_25\.5hr')
arm_rows('GSE50559', 'miR-17~92', 'mouse spleen B cells (miR-17~92 transgenes)', 'gain', r'^Tg', r'^WT', flt=lambda s: s['platform'] == 'GPL5642')
arm_rows('GSE42760', 'miR-17~92', 'mouse TFH cells (miR-17~92 deletion)', 'loss', r'^KO TFH', r'^Control TFH')
arm_rows('GSE47155', 'DLK1-DIO3 (miR-379/410)', 'mouse liver (miR-379/410 KO)', 'loss', r'^Liver_KO', r'^Liver_WildType')
arm_rows('GSE47156', 'DLK1-DIO3 (miR-379/410)', 'mouse liver, P0 (miR-379/410 KO)', 'loss', r'^Liver_KO_P0', r'^Liver_WildType_P0')
arm_rows('GSE57112', 'DLK1-DIO3 (miR-379/410)', 'mouse liver, P1 (miR-379/410 KO)', 'loss', r'^Liver_KO_', r'^Liver_WildType_P1')
arm_rows('GSE32533', 'miR-17~92', 'mouse CD4+CD25- T cells (miR-17~92 deletion)', 'loss', r'^KO-mock', r'^WT-mock')
arm_rows('GSE15741', 'miR-200b~200a~429', 'mouse 344SQ lung adenocarcinoma cells', 'gain', r'^429 ', r'^Control ')
arm_rows('GSE19631', 'miR-200b~200a~429', 'mouse 4TO7 cells (cluster 1)', 'gain', r'^4TO7 Cluster 1,', r'^4TO7 Puro vec control', platform='GPL7202')
arm_rows('GSE19631', 'miR-200c~141', 'mouse 4TO7 cells (cluster 2)', 'gain', r'^4TO7 Cluster 2,', r'^4TO7 Hygro vec, ', platform='GPL7202')
arm_rows('GSE19631', 'both miR-200 clusters', 'mouse 4TO7 cells (clusters 1+2)', 'gain', r'^4TO7 Clusters 1\+2', r'^4TO7 Hygro vec \+ Puro vec', platform='GPL7202')
A = pd.DataFrame(rows); A.to_csv(f'{INB}/P2/p2_mouse_arms.tsv', sep='\t', index=False)
print(A.groupby(['gse', 'cell', 'arm']).size().unstack().to_string())
route = {'GPL6246': 'cel', 'GPL11533': 'cel', 'GPL8321': 'cel', 'GPL1261': 'cel', 'GPL5642': 'matrix', 'GPL7202': 'agilent1', 'GPL10787': 'agilent1'}
with open(f'{INB}/P2/batch_mouse_arrays.tsv', 'w') as f:
    for (g, tf, cell), d in A.groupby(['gse', 'tf', 'cell'], sort=False):
        f.write(f"P2\tP2/p2_mouse_arms.tsv\t{g}\t{tf}\t{cell}\tprimary\t--route {route[d.platform.iloc[0]]} --mouse\n")
