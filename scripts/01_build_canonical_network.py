#!/usr/bin/env python3
"""
01_build_canonical_network.py
Rebuild a single, canonical, DIRECTED, SIGNED, type-annotated regulatory network
from the six Cytoscape exports in miRNA_Github_GPR/.

Repairs applied (each logged):
  R1  miRBase case collapse:  hsa-mir-X (precursor) and hsa-miR-X (mature) were
      deposited as two disconnected nodes, which severs the TF->miRNA and
      miRNA->target arms of every FFL. Collapsed to one canonical mature node.
  R2  Duplicate edges removed (SIF files contain repeated rows).
  R3  Node-type conflicts resolved (10 bona fide TFs were labelled 'Gene' in
      some exports): TF assignment is the union across all six files.
  R4  Edges made DIRECTED (the published pipeline built an undirected igraph).
  R5  Edge signs assigned by regulator class, later refined with TRRUST mode.
"""
import csv, re, json, collections, os, sys

REPO = "/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR"
OUT  = "/path/to/revision/data"
PREFIXES = ['3-miR', '3-TF', '3-Comp', '4-TF', '5-TF', '6-TF']
os.makedirs(OUT, exist_ok=True)

log = []
def L(msg):
    log.append(msg); print(msg)

# ---------------------------------------------------------------- node types
raw_types = collections.defaultdict(set)
for p in PREFIXES:
    for r in csv.DictReader(open(f'{REPO}/node_attributes/{p}.csv')):
        raw_types[r['name'].strip()].add(r['Type'].strip())

TF_set = {n for n, t in raw_types.items() if 'TF' in t}
conflicts = sorted(n for n, t in raw_types.items() if len(t) > 1)
L(f"[R3] node-type conflicts resolved as TF: {len(conflicts)} -> {conflicts}")

def is_mir(n):
    return bool(re.match(r'^hsa-(mir|miR|let)-', n, flags=re.I))

# ------------------------------------------------- canonical miRNA naming
def canon(n):
    """hsa-mir-130a -> hsa-miR-130a ; hsa-let-7a stays ; genes untouched."""
    if not is_mir(n):
        return n
    m = re.match(r'^hsa-(mir|miR)-(.+)$', n)
    if m:
        return 'hsa-miR-' + m.group(2)
    return n.replace('hsa-LET-', 'hsa-let-')

# collapse precursor copy-number suffixes (hsa-miR-129-2 -> hsa-miR-129) only
# when the un-suffixed node also exists in the data set
all_names = set(raw_types)
canon_names = {canon(n) for n in all_names}
def canon2(n):
    c = canon(n)
    m = re.match(r'^(hsa-miR-\d+[a-z]?)-([12])$', c)
    if m and m.group(1) in canon_names:
        return m.group(1)
    return c

merge_map = {n: canon2(n) for n in all_names}
merged = collections.defaultdict(list)
for k, v in merge_map.items():
    merged[v].append(k)
n_collapsed = sum(1 for v in merged.values() if len(v) > 1)
L(f"[R1] miRBase case/precursor collapse: {len(all_names)} raw node labels -> "
  f"{len(merged)} canonical nodes ({n_collapsed} labels merged from >1 alias)")

def ntype(c):
    if is_mir(c): return 'miRNA'
    if any(a in TF_set for a in merged[c]): return 'TF'
    return 'Gene'

node_type = {c: ntype(c) for c in merged}
L("[nodes] " + str(collections.Counter(node_type.values())))

# ------------------------------------------------------------------- edges
# source column of the SIF is the regulator (verified: mature miRNA nodes are
# never targets, TransmiR TF->miRNA rows always have the precursor as target)
edges = {}                       # (src,dst) -> dict
prov  = collections.defaultdict(set)
raw_rows = 0
for p in PREFIXES:
    for line in open(f'{REPO}/SIF_files/{p}.sif'):
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 3:
            parts = line.split()
        if len(parts) < 3:
            continue
        raw_rows += 1
        s, d = merge_map.get(parts[0], canon2(parts[0])), merge_map.get(parts[-1], canon2(parts[-1]))
        if s == d:
            continue
        ts, td = node_type.get(s, 'Gene'), node_type.get(d, 'Gene')
        if ts == 'miRNA' and td in ('Gene', 'TF'):
            et, sign = 'miRNA_target', -1
        elif ts in ('TF',) and td == 'miRNA':
            et, sign = 'TF_miRNA', +1
        elif ts == 'Gene' and td == 'miRNA':
            et, sign = 'TF_miRNA', +1          # 10 mis-typed TFs, resolved above
        elif ts == 'TF' and td in ('Gene', 'TF'):
            et, sign = 'TF_target', +1
        elif ts == 'miRNA' and td == 'miRNA':
            et, sign = 'miRNA_miRNA', +1
        elif ts == 'Gene' and td in ('Gene', 'TF'):
            et, sign = 'gene_gene', +1
        else:
            et, sign = 'other', +1
        edges.setdefault((s, d), {'src': s, 'dst': d, 'etype': et, 'sign': sign})
        prov[(s, d)].add(p)

L(f"[R2] SIF rows read: {raw_rows}  ->  unique directed canonical edges: {len(edges)}")
et_counts = collections.Counter(e['etype'] for e in edges.values())
L("[edges] " + str(dict(et_counts)))

# reciprocal TF<->miRNA pairs (composite arm)
recip = sum(1 for (s, d) in edges
            if node_type[s] == 'TF' and node_type[d] == 'miRNA' and (d, s) in edges)
L(f"[edges] reciprocal TF<->miRNA pairs (composite arm): {recip}")

# ------------------------------------------------------------------ write out
with open(f'{OUT}/canonical_edges.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['source', 'target', 'edge_type', 'sign', 'source_type',
                'target_type', 'in_networks'])
    for (s, d), e in sorted(edges.items()):
        w.writerow([s, d, e['etype'], e['sign'], node_type[s], node_type[d],
                    ';'.join(sorted(prov[(s, d)]))])

with open(f'{OUT}/canonical_nodes.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['name', 'type', 'aliases'])
    for c in sorted(merged):
        w.writerow([c, node_type[c], ';'.join(sorted(merged[c]))])

json.dump({'merge_map': merge_map, 'node_type': node_type},
          open(f'{OUT}/name_map.json', 'w'), indent=1)
open(f'{OUT}/01_build_log.txt', 'w').write('\n'.join(log) + '\n')
print("\nWritten to", OUT)
