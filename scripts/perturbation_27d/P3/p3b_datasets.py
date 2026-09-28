#!/usr/bin/env python3
"""P3b: KnockTF 2.0 human datasets for the network TFs (sources of the analysed network's TF->target edges), with the
number of control and treated samples read from each dataset's KnockTF page.  Writes p3b_datasets.tsv."""
import os, sys, re, csv, time
import pandas as pd
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); REV = os.path.dirname(INB); sys.path.insert(0, INB)
import _dl
E = pd.read_csv(f'{REV}/data/canonical_edges.tsv', sep='\t'); TFS = set(E[E.edge_type == 'TF_target'].source)
K = pd.read_csv(f'{INB}/P1/knocktf_datasets.tsv', sep='\t')
K['sample_id'] = K.profile_file.str.extract(r'(DataSet_[0-9_]+)')[0]
K = K[K.tf.str.upper().isin({t.upper() for t in TFS})].copy()
rows = []
for r in K.itertuples():
    url = f'http://www.licpathway.net/KnockTFv2/search/search_sample_result.php?sample_id={r.sample_id}'
    try:
        p = _dl.fetch(url, 'P3/knocktf_pages', f'{r.sample_id}.html', 'P3b', 'KnockTF 2.0', r.sample_id, 'not stated on the download page')
        t = re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', ' ', open(p, encoding='utf-8', errors='replace').read()))
        nc = re.search(r'# of Control:\s*(\d+)', t); nt = re.search(r'# of Treat:\s*(\d+)', t)
        rows.append(dict(sample_id=r.sample_id, tf=r.tf, knock_method=r.knock_method, biosample=r.biosample, tissue=r.tissue, biosample_type=r.biosample_type,
                         profile_id=r.profile_id, platform=r.platform, n_control=int(nc.group(1)) if nc else None, n_treat=int(nt.group(1)) if nt else None))
    except Exception as e:
        rows.append(dict(sample_id=r.sample_id, tf=r.tf, knock_method=r.knock_method, biosample=r.biosample, tissue=r.tissue, biosample_type=r.biosample_type,
                         profile_id=r.profile_id, platform=r.platform, n_control=None, n_treat=None, error=str(e)[:100]))
    time.sleep(0.3)
pd.DataFrame(rows).to_csv(f'{INB}/P3/p3b_datasets.tsv', sep='\t', index=False)
D = pd.DataFrame(rows)
print(len(D), 'KnockTF datasets for', D.tf.nunique(), 'network TFs; with >= 2 control and >= 2 treated:', int(((D.n_control >= 2) & (D.n_treat >= 2)).sum()))
