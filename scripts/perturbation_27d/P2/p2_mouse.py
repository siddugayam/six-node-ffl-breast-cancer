#!/usr/bin/env python3
"""P2 mouse (secondary): differential expression for the mouse cluster-perturbation datasets, with mouse genes mapped to
human genes by one-to-one orthologues (MGI HOM_MouseHumanSequence.rpt; a homology class with exactly one mouse and one
human gene).  No NCBI-generated counts exist for these series, so the submitters' files are used:
  raw counts -> _de.R counts (filterByExpr, TMM, voom); processed values (normalised counts, FPKM) -> log2(x+1), limma
  (flagged 'processed'); a submitter DESeq2 table -> its log2FoldChange taken as is (flagged 'processed DE'; expressed =
  MeanNormCts >= median).
Members are the human network nodes of the perturbed cluster (miR-379/410 = the network's 10 kb cluster cl10kb_420, 28
nodes).  Writes P2/de/<name>.tsv and .meta.json in the P2 format (cls = 'mouse')."""
import os, sys, io, re, gzip, json, tarfile, subprocess
import numpy as np, pandas as pd
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); REV = os.path.dirname(INB); sys.path.insert(0, INB)
import _dl, _geo
OUT = f'{INB}/P2/de'; W = os.path.join(_dl.CACHE, 'P2/work'); os.makedirs(W, exist_ok=True)
H = pd.read_csv(os.path.join(_dl.CACHE, 'common/HOM_MouseHumanSequence.rpt'), sep='\t', dtype=str)
mm, hs = H[H['NCBI Taxon ID'] == '10090'], H[H['NCBI Taxon ID'] == '9606']
one = set(mm['DB Class Key'].value_counts()[lambda x: x == 1].index) & set(hs['DB Class Key'].value_counts()[lambda x: x == 1].index)
hmap = dict(zip(hs['DB Class Key'], hs.Symbol))
M2H = {s: hmap[k] for s, k in zip(mm.Symbol, mm['DB Class Key']) if k in one}
E2H = {e: hmap[k] for e, k in zip(mm['EntrezGene ID'], mm['DB Class Key']) if k in one}
gi = pd.read_csv(_dl.fetch('https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz', 'common', 'Mus_musculus.gene_info.gz', 'P2',
                            'NCBI Gene (gene_info, mouse)', 'Mus_musculus.gene_info', 'NCBI: no restrictions on use or distribution (NCBI website and data usage policies)'),
                 sep='\t', usecols=['Symbol', 'dbXrefs'], dtype=str)
ens = gi.dbXrefs.str.extract(r'Ensembl:(ENSMUSG\d+)')[0]
ENS2H = {e: M2H[sym] for e, sym in zip(ens, gi.Symbol) if isinstance(e, str) and sym in M2H}
L = pd.read_csv(f'{REV}/data/layer_miRNA_miRNA.tsv', sep='\t')
C420 = sorted(set(L[L.cluster_id == 'cl10kb_420'].miRNA_1) | set(L[L.cluster_id == 'cl10kb_420'].miRNA_2))
MEM = {'miR-17~92': 'hsa-miR-17;hsa-miR-18a;hsa-miR-19a;hsa-miR-19b;hsa-miR-20a;hsa-miR-92a', 'miR-29a/b-1': 'hsa-miR-29a;hsa-miR-29b',
       'miR-183~96~182': 'hsa-miR-183;hsa-miR-96;hsa-miR-182', 'DLK1-DIO3 (miR-379/410)': ';'.join(C420)}
def write(name, tab, meta):
    tab.to_csv(f'{OUT}/{name}.tsv', sep='\t', index=False); json.dump(meta, open(f'{OUT}/{name}.meta.json', 'w')); print(name, len(tab), meta['note'][:80])
def run_de(name, counts, samples, kind, meta, idmap):
    """counts: DataFrame (first column id); idmap: id -> human symbol."""
    m = f'{W}/{name}.matrix.tsv'; counts.to_csv(m, sep='\t', index=False)
    gm = f'{W}/{name}.genemap.tsv'; pd.DataFrame({'id': list(idmap), 'symbol': list(idmap.values())}).to_csv(gm, sep='\t', index=False)
    sp = f'{W}/{name}.spec.json'; json.dump(dict(kind=kind, matrix=m, gene_map=gm, samples=samples), open(sp, 'w'))
    r = subprocess.run(['Rscript', f'{INB}/_de.R', sp, f'{OUT}/{name}.tsv'], capture_output=True, text=True)
    if r.returncode: print('FAILED', name, r.stderr[-400:]); return
    meta['de'] = r.stdout.strip(); json.dump(meta, open(f'{OUT}/{name}.meta.json', 'w')); print(name, meta['de'], meta['note'][:70])
base = dict(cls='mouse', tag='primary', platform='', n_pert=0, n_ctrl=0)
# GSE83636: miR-17~92 KO vs WT hippocampus, per-sample HTSeq counts (mouse symbols)
tf = tarfile.open(_geo.suppl_fetch('GSE83636', 'GSE83636_RAW.tar', 'P2')); cols = {}
for mbr in tf.getmembers():
    x = pd.read_csv(io.BytesIO(gzip.decompress(tf.extractfile(mbr).read())), sep='\t', header=None); x = x[~x[0].astype(str).str.startswith('__')]
    cols[mbr.name.split('_')[0]] = x.set_index(0)[1]
C = pd.DataFrame(cols).reset_index().rename(columns={'index': 'id'})
S = {'GSM2211615': 'ctrl', 'GSM2211616': 'ctrl', 'GSM2211617': 'ctrl', 'GSM2211618': 'pert', 'GSM2211619': 'pert', 'GSM2211620': 'pert'}
run_de('GSE83636__miR-17_92__mouse_hippocampus__primary', C, S, 'counts', dict(base, gse='GSE83636', tf='miR-17~92', cell='mouse hippocampus (miR-17~92 KO)',
       direction='loss', mirna=MEM['miR-17~92'], route='submitter counts', n_pert=3, n_ctrl=3, note='submitter HTSeq counts (Ensembl ids); mouse -> human one-to-one orthologues'), ENS2H)
# GSE77010: del17-92 vs floxed lymphoma cells, counts with mouse Entrez ids
C = pd.read_csv(_geo.suppl_fetch('GSE77010', 'GSE77010_counts.tsv.gz', 'P2'), sep='\t', dtype={0: str}); C = C.rename(columns={C.columns[0]: 'id'})
S = {'fl_rep1': 'ctrl', 'fl_rep2': 'ctrl', 'fl_rep3': 'ctrl', 'del_rep1': 'pert', 'del_rep2': 'pert', 'del_rep3': 'pert'}
run_de('GSE77010__miR-17_92__mouse_lymphoma__primary', C, S, 'counts', dict(base, gse='GSE77010', tf='miR-17~92', cell='mouse Eu-Myc lymphoma cells (miR-17~92 deletion)',
       direction='loss', mirna=MEM['miR-17~92'], route='submitter counts', n_pert=3, n_ctrl=3, note='submitter counts (Entrez ids); mouse -> human one-to-one orthologues'), E2H)
# GSE253734: miR-29a/b-1 KO vs WT crypts (untreated comparison), counts with Ensembl ids and gene names
t = pd.read_csv(_geo.suppl_fetch('GSE253734', 'GSE253734_all_compare.txt.gz', 'P2'), sep='\t')
C = t[['gene_name'] + [c for c in t.columns if c.endswith('_count')]].rename(columns={'gene_name': 'id'}).groupby('id', as_index=False).sum()
S = {f'KO{i}_count': 'pert' for i in range(1, 5)} | {f'WT{i}_count': 'ctrl' for i in range(1, 5)}
run_de('GSE253734__miR-29a_b-1__mouse_crypts__primary', C, S, 'counts', dict(base, gse='GSE253734', tf='miR-29a/b-1', cell='mouse small-intestinal crypts (miR-29a/b-1 KO)',
       direction='loss', mirna=MEM['miR-29a/b-1'], route='submitter counts', n_pert=4, n_ctrl=4, note='submitter counts (untreated comparison); mouse -> human one-to-one orthologues'), M2H)
# GSE151844: miR-29a/b-1 transgene (gain) vs WT islets, FPKM (processed)
t = pd.read_csv(_geo.suppl_fetch('GSE151844', 'GSE151844_gene_fpkm.txt.gz', 'P2'), sep='\t')
C = t[['gene_name', 'wt1', 'wt2', 'miR29ab1', 'miR29ab2']].rename(columns={'gene_name': 'id'}).groupby('id', as_index=False).sum()
run_de('GSE151844__miR-29a_b-1__mouse_islets__primary', C, {'wt1': 'ctrl', 'wt2': 'ctrl', 'miR29ab1': 'pert', 'miR29ab2': 'pert'}, 'linear_matrix',
       dict(base, gse='GSE151844', tf='miR-29a/b-1', cell='mouse pancreatic islets (miR-29a/b-1 transgene)', direction='gain', mirna=MEM['miR-29a/b-1'],
            route='processed (FPKM)', n_pert=2, n_ctrl=2, note='processed FPKM, log2(x+1), limma; mouse -> human one-to-one orthologues'), M2H)
# GSE247439: miR-379/410 cKO vs WT hippocampus, normalised read counts (processed)
t = pd.read_csv(_geo.suppl_fetch('GSE247439', 'GSE247439_readcount_matrix.tab.gz', 'P2'), sep='\t', index_col=0).reset_index().rename(columns={'index': 'id'})
S = {f'ko{i}': 'pert' for i in range(1, 7)} | {f'Wt{i}': 'ctrl' for i in range(1, 7)}
run_de('GSE247439__DLK1-DIO3__mouse_hippocampus__primary', t.rename(columns={t.columns[0]: 'id'}), S, 'linear_matrix',
       dict(base, gse='GSE247439', tf='DLK1-DIO3 (miR-379/410)', cell='mouse hippocampus (miR-379/410 cKO)', direction='loss', mirna=MEM['DLK1-DIO3 (miR-379/410)'],
            route='processed (normalised read counts)', n_pert=6, n_ctrl=6, note='processed (non-integer read counts), log2(x+1), limma; mouse -> human orthologues'), M2H)
# GSE140838: Mirg deletion vs control ESC, submitter DESeq2 table (processed DE)
t = pd.read_csv(_geo.suppl_fetch('GSE140838', 'GSE140838_RNAseq_DESeq2_Del.ESC_vs_Ctl.ESC.txt.gz', 'P2'), sep='\t')
t['human'] = t.gene_name.map(M2H); t = t.dropna(subset=['human', 'log2FoldChange']).sort_values('MeanNormCts', ascending=False).drop_duplicates('human')
tab = pd.DataFrame({'gene': t.human, 'logFC': t.log2FoldChange, 'SE': np.nan, 't': np.nan, 'P': np.nan, 'adjP': t.padj, 'AveExpr': np.log2(t.MeanNormCts + 1),
                    'ctrl_mean_log2': np.log2(t.MeanNormCts + 1), 'expressed': t.MeanNormCts >= t.MeanNormCts.median(), 'n_pert': np.nan, 'n_ctrl': np.nan})
write('GSE140838__DLK1-DIO3__mouse_ESC__primary', tab, dict(base, gse='GSE140838', tf='DLK1-DIO3 (miR-379/410)', cell='mouse ESC (Mirg deletion)', direction='loss',
      mirna=MEM['DLK1-DIO3 (miR-379/410)'], route='processed DE', note='submitter DESeq2 log2FoldChange (Del vs Ctl ESC); expressed = MeanNormCts >= median; mouse -> human orthologues'))
# GSE78505: miR-183C KO vs WT Th17 cells, gene counts (xlsx)
x = pd.read_excel(_geo.suppl_fetch('GSE78505', 'GSE78505_Xiao_NS107_gene_cts.xlsx', 'P2')).rename(columns={'Gene': 'id'})
x = x[['id', 'AGM-E', 'AGM-G', 'AGM-F', 'AGM-H']].groupby('id', as_index=False).sum()      # E, G = WT; F, H = miR-183C KO (GSM descriptions)
idmap = ENS2H if x.id.astype(str).str.startswith('ENSMUSG').mean() > 0.5 else M2H
run_de('GSE78505__miR-183_96_182__mouse_Th17__primary', x, {'AGM-E': 'ctrl', 'AGM-G': 'ctrl', 'AGM-F': 'pert', 'AGM-H': 'pert'}, 'counts',
       dict(base, gse='GSE78505', tf='miR-183~96~182', cell='mouse Th17 cells (miR-183C KO)', direction='loss', mirna=MEM['miR-183~96~182'],
            route='submitter counts', n_pert=2, n_ctrl=2, note='submitter gene counts (xlsx); mouse -> human one-to-one orthologues'), idmap)
