#!/usr/bin/env python3
"""
INBOX_2026-09-27c, Q7 (descriptive): what drives the co-transcription effect, and where is it concentrated?
Inputs (read-only): the S1 census graph of INBOX_2026-09-27b (S1/graph_nolegacy_null.txt, graph_nolegacy_labels.txt,
node_names.txt; label 0 = TF->miRNA, 1 = miRNA->target, 6 = 10 kb co-transcription pair, exactly as
S1/s1_cotranscription_overlap.py reads them); REV/data/layer_TF_miRNA.tsv (TransmiR PMIDs per TF->miRNA edge);
REV/cache/direct/ts/miR_Family_Info.txt (TargetScan: human mature miRNAs, seed = nt 2-8 'Seed+m8', 'miR family');
REV/data/mirna_id_map.tsv (the project's mature -> canonical map); S1/obs/nolegacy_bhat6_composite_instances.tsv
(the 3,898 BHAT6 composite instances: M1, M2, T1, T2, G1, G2 node indices).
(a) co-transcribed pairs sharing >= 1 TF regulator; for each shared TF, whether the two TransmiR records carry the same
    PMID set (identical), share a PMID, or share none.
(b) seed sharing: a canonical node's seeds are those of every human mature form mapped to it (mirna_id_map, else the
    mature ID without its -3p/-5p suffix); a pair shares a seed if the two sets intersect; the same for TargetScan family
    names.  A miRBase family file (miFam.dat) is not in the project.  Target Jaccard as in S1/s1_cotranscription_overlap.py.
(c) BHAT6 composite instances per M1-M2 pair and per co-transcription cluster: primary = the genomic cluster_id of the
    pair in REV/data/layer_miRNA_miRNA.tsv (10 kb); secondary = connected components of the 10 kb pairs, which merge
    genomic clusters through canonical nodes that collapse paralog copies (e.g. miR-19b, miR-92a on chr13 and chrX).
"""
import os, csv, re, statistics, collections
H = os.path.dirname(os.path.abspath(__file__)); REV = '/path/to/revision'; S1 = f'{REV}/INBOX_2026-09-27b/S1'
names = open(f'{S1}/node_names.txt').read().split('\n')[:587]
L = open(f'{S1}/graph_nolegacy_null.txt').read().rstrip('\n').split('\n'); ty = list(map(int, L[1].split()))
labs = list(map(int, open(f'{S1}/graph_nolegacy_labels.txt').read().split()))
arcs = [tuple(map(int, l.split()[:2])) for l in L[2:]]
tg = collections.defaultdict(set); rg = collections.defaultdict(set); pairs = set()
for (u, v), lb in zip(arcs, labs):
    if lb == 1: tg[u].add(v)
    elif lb == 0: rg[v].add(u)
    elif lb == 6: pairs.add(frozenset((u, v)))
pairs = sorted(tuple(sorted(p)) for p in pairs)
jac = lambda a, b: len(a & b) / len(a | b) if (a | b) else 0.0
out = []
P = lambda *a: (out.append(' '.join(str(x) for x in a)), print(*a))
P(f'Q7 co-transcribed (10 kb) pairs {len(pairs)} among {len({x for p in pairs for x in p})} miRNA nodes')
# ------------------------------------------------------------ (a)
pm = collections.defaultdict(set)
for r in csv.DictReader(open(f'{REV}/data/layer_TF_miRNA.tsv'), delimiter='\t'):
    pm[(r['source'], r['target'])] |= {x.strip() for x in r['pmid'].split(';') if x.strip() not in ('', 'NA')}
shared_pairs = [(a, b) for a, b in pairs if rg[a] & rg[b]]
cnt = collections.Counter(); per_pair_all_identical = 0; rows = []
for a, b in pairs:
    sh = sorted(rg[a] & rg[b]); kinds = []
    for t in sh:
        A, B = pm.get((names[t], names[a]), set()), pm.get((names[t], names[b]), set())
        k = 'no record' if not A or not B else 'identical PMID set' if A == B else 'share a PMID' if A & B else 'no common PMID'
        cnt[k] += 1; kinds.append(k)
    if sh and all(k == 'identical PMID set' for k in kinds): per_pair_all_identical += 1
    rows.append(dict(miRNA_1=names[a], miRNA_2=names[b], shared_TF_regulators=len(sh), shared_TFs=';'.join(names[t] for t in sh),
                     shared_TF_identical_PMIDs=sum(k == 'identical PMID set' for k in kinds), target_jaccard=round(jac(tg[a], tg[b]), 4)))
P(f'Q7(a) co-transcribed pairs sharing >= 1 TF regulator: {len(shared_pairs)} of {len(pairs)}; shared TF->miRNA edge pairs {sum(cnt.values())}: '
  + '; '.join(f'{k} {v}' for k, v in cnt.most_common()) + f'; pairs whose every shared TF has an identical PMID set for both miRNAs: {per_pair_all_identical}')
# ------------------------------------------------------------ (b)
idmap = {r['raw_id']: r['canonical_id'] for r in csv.DictReader(open(f'{REV}/data/mirna_id_map.tsv'), delimiter='\t')}
seeds = collections.defaultdict(set); fams = collections.defaultdict(set)
node_ix = {n: i for i, n in enumerate(names)}
for r in csv.DictReader(open(f'{REV}/cache/direct/ts/miR_Family_Info.txt'), delimiter='\t'):
    if r['Species ID'] != '9606': continue
    mid = r['MiRBase ID']; c = idmap.get(mid) or re.sub(r'-[35]p$', '', mid)
    if c in node_ix and ty[node_ix[c]] == 0:
        seeds[node_ix[c]].add(r['Seed+m8']); fams[node_ix[c]].add(r['miR family'])
cot = {x for p in pairs for x in p}
P(f'Q7(b) co-transcribed miRNA nodes with a TargetScan human seed: {sum(1 for x in cot if seeds[x])} of {len(cot)}')
known = [(a, b) for a, b in pairs if seeds[a] and seeds[b]]
ss = [(a, b) for a, b in known if seeds[a] & seeds[b]]; ns = [(a, b) for a, b in known if not seeds[a] & seeds[b]]
fs = [(a, b) for a, b in known if fams[a] & fams[b]]
jm = lambda L_: statistics.mean(jac(tg[a], tg[b]) for a, b in L_) if L_ else float('nan')
P(f'Q7(b) pairs with seeds for both miRNAs {len(known)} of {len(pairs)}; sharing an identical nt 2-8 seed {len(ss)}; sharing a TargetScan family {len(fs)}; '
  f'mean target Jaccard: seed-sharing {jm(ss):.3f} (n = {len(ss)}), not seed-sharing {jm(ns):.3f} (n = {len(ns)}), pairs without seed data '
  f'{jm([p for p in pairs if p not in known]):.3f} (n = {len(pairs) - len(known)}); miRBase family (miFam.dat): NOT FOUND in the project')
for r, (a, b) in zip(rows, pairs):
    r['seed_known_both'] = bool(seeds[a] and seeds[b]); r['share_seed_nt2_8'] = bool(seeds[a] & seeds[b]); r['share_targetscan_family'] = bool(fams[a] & fams[b])
with open(f'{H}/q7_pairs.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
# ------------------------------------------------------------ (c)
adj = collections.defaultdict(set)
for a, b in pairs: adj[a].add(b); adj[b].add(a)
comp = {}
for s in sorted(adj):
    if s in comp: continue
    st = [s]; comp[s] = s
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in comp: comp[y] = s; st.append(y)
members = collections.defaultdict(list)
for x, c in comp.items(): members[c].append(names[x])
inst = [tuple(map(int, l.split())) for l in open(f'{S1}/obs/nolegacy_bhat6_composite_instances.tsv') if l.strip()]
bypair = collections.Counter(tuple(sorted((names[i[0]], names[i[1]]))) for i in inst)
bycl = collections.Counter(comp[i[0]] for i in inst)
assert all(comp[i[0]] == comp[i[1]] for i in inst)
# genomic clusters: cluster_id of REV/data/layer_miRNA_miRNA.tsv for the 10 kb pairs (primary definition); the connected
# components above merge genomic clusters through canonical nodes that collapse paralog copies (e.g. miR-19b, miR-92a)
LY = [r for r in csv.DictReader(open(f'{REV}/data/layer_miRNA_miRNA.tsv'), delimiter='\t') if r['threshold_10kb'] == 'TRUE']
cid = collections.defaultdict(set)
for r in LY: cid[tuple(sorted((r['miRNA_1'], r['miRNA_2'])))].add((r['cluster_id'], r['chrom']))
P(f'Q7(c) 10 kb pairs in layer_miRNA_miRNA.tsv {len(cid)}; equal to the census pairs: {set(cid) == {tuple(sorted((names[a], names[b]))) for a, b in pairs}}; '
  f'pairs with more than one genomic cluster id: {sum(1 for v in cid.values() if len(v) > 1)}')
gmem = collections.defaultdict(set)
for k, v in cid.items():
    for c in v: gmem[c] |= set(k)
bygc = collections.Counter()
for i in inst:
    for c in cid[tuple(sorted((names[i[0]], names[i[1]])))]: bygc[c] += 1
topg = bygc.most_common(10)
P(f'Q7(c) genomic clusters (layer cluster_id) used {len(bygc)} of {len(gmem)}; top 10 (instances, share of {len(inst)}; cluster id, chromosome; members): ' + ' | '.join(
  f'{v} ({100 * v / len(inst):.1f} %): {c[0]}, {c[1]}; {", ".join(sorted(gmem[c]))}' for c, v in topg))
P(f'Q7(c) share of all {len(inst)} instances in the top 5 genomic clusters: {sum(v for c, v in topg[:5])} ({100 * sum(v for c, v in topg[:5]) / len(inst):.1f} %)')
with open(f'{H}/q7c_bhat6_by_genomic_cluster.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['rank', 'cluster_id', 'chrom', 'instances', 'pct_of_all', 'members'])
    for k, (c, v) in enumerate(bygc.most_common(), 1): w.writerow([k, c[0], c[1], v, round(100 * v / len(inst), 2), ';'.join(sorted(gmem[c]))])
P(f'Q7(c) BHAT6 composite instances {len(inst)}; M1-M2 pairs used {len(bypair)}; connected components of the 10 kb pairs {len(members)}, used {len(bycl)} (secondary view; see the genomic clusters)')
top = bycl.most_common(10); tot = len(inst)
P('Q7(c) secondary view, top 10 connected components (instances; members): ' + ' | '.join(f'{v} ({100 * v / tot:.1f} %): {", ".join(sorted(members[c]))}' for c, v in top))
P(f'Q7(c) secondary view, share in the top 5 connected components: {sum(v for c, v in top[:5])} ({100 * sum(v for c, v in top[:5]) / tot:.1f} %); top 1: {top[0][1]} ({100 * top[0][1] / tot:.1f} %)')
P('Q7(c) top 10 M1-M2 pairs: ' + '; '.join(f'{a}/{b} {v}' for (a, b), v in bypair.most_common(10)))
with open(f'{H}/q7c_bhat6_by_cluster.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['rank', 'instances', 'pct_of_all', 'n_members', 'members'])
    for k, (c, v) in enumerate(bycl.most_common(), 1): w.writerow([k, v, round(100 * v / tot, 2), len(members[c]), ';'.join(sorted(members[c]))])
with open(f'{H}/q7c_bhat6_by_pair.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['miRNA_1', 'miRNA_2', 'instances'])
    for (a, b), v in bypair.most_common(): w.writerow([a, b, v])
open(f'{H}/q7_cotranscription.txt', 'w').write('\n'.join(out) + '\n')
