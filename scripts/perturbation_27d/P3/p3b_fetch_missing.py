#!/usr/bin/env python3
"""P3b: KnockTF's bulk file (knocktf_v2_main_human.txt) holds only the DataSet_01/02 series and 4 DataSet_03 series.
The differential expression of every other network-TF dataset is fetched from KnockTF's own per-dataset endpoint
(search_sample_server.php, the table behind each dataset page), and written to knocktf_network_tf_subset_extra.tsv.gz."""
import os, sys, json, urllib.parse, time
import pandas as pd
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, INB)
import _dl
D = f'{INB}/P3'
R = pd.read_csv(f'{D}/p3b_per_dataset.tsv', sep='\t')
miss = R[R.status.str.startswith('no rows')].drop_duplicates('sample_id').sample_id.tolist()
rows = []
for sid in miss:
    q = urllib.parse.urlencode({'sample_id': sid, 'species': 'Homo sapiens', 'draw': 1, 'start': 0, 'length': 100000})
    try:
        p = _dl.fetch(f'http://www.licpathway.net/KnockTFv2/search/search_sample_server.php?{q}', 'P3/knocktf_json', f'{sid}.json', 'P3b', 'KnockTF 2.0',
                      f'{sid} differential expression (dataset endpoint)', 'not stated on the download page', timeout=300)
        d = json.load(open(p))['data']
        for x in d: rows.append((sid, x['tf_name'], x['gene_name'], x['mean_control'], x['log2fc']))
    except Exception as e:
        print('failed', sid, str(e)[:120], flush=True)
    time.sleep(0.3)
pd.DataFrame(rows, columns=['Sample_ID', 'TF', 'Gene', 'Mean_Control', 'Log2FC']).to_csv(os.path.join(_dl.CACHE, 'P3/knocktf/knocktf_network_tf_subset_extra.tsv.gz'),
                                                                                       sep='\t', index=False, compression='gzip')
print('datasets', len(miss), 'rows', len(rows))
