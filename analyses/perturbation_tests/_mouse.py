#!/usr/bin/env python3
"""Mouse -> human gene mapping by one-to-one orthologues (MGI HOM_MouseHumanSequence.rpt: a homology class with exactly
one mouse and one human gene), as in P2/p2_mouse.py.  M2H: mouse symbol -> human symbol; E2H: mouse Entrez id -> human
symbol; ENS2H: mouse Ensembl gene id -> human symbol (via NCBI mouse gene_info dbXrefs)."""
import os
import pandas as pd
import _dl
H = pd.read_csv(os.path.join(_dl.CACHE, 'common/HOM_MouseHumanSequence.rpt'), sep='\t', dtype=str)
mm, hs = H[H['NCBI Taxon ID'] == '10090'], H[H['NCBI Taxon ID'] == '9606']
one = set(mm['DB Class Key'].value_counts()[lambda x: x == 1].index) & set(hs['DB Class Key'].value_counts()[lambda x: x == 1].index)
hmap = dict(zip(hs['DB Class Key'], hs.Symbol))
M2H = {s: hmap[k] for s, k in zip(mm.Symbol, mm['DB Class Key']) if k in one}
E2H = {e: hmap[k] for e, k in zip(mm['EntrezGene ID'], mm['DB Class Key']) if k in one}
gi = pd.read_csv(os.path.join(_dl.CACHE, 'common', 'Mus_musculus.gene_info.gz'), sep='\t', usecols=['Symbol', 'dbXrefs'], dtype=str)
ens = gi.dbXrefs.str.extract(r'Ensembl:(ENSMUSG\d+)')[0]
ENS2H = {e: M2H[sym] for e, sym in zip(ens, gi.Symbol) if isinstance(e, str) and sym in M2H}
def to_human(ids, symbols=None):
    """Human symbol for each mouse id (Ensembl, Entrez or symbol) or, if given, mouse symbol; '' when not one-to-one."""
    out = []
    for i, s in zip(ids, symbols if symbols is not None else [None] * len(ids)):
        i = str(i).split('.')[0]
        h = ENS2H.get(i) if i.startswith('ENSMUSG') else (E2H.get(i) if i.isdigit() else M2H.get(i))
        if not h and isinstance(s, str): h = M2H.get(s, '')
        out.append(h or '')
    return out
