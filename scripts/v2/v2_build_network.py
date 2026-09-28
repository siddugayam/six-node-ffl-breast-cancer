#!/usr/bin/env python3
"""
Independent re-derivation of the canonical miRNA-TF-gene FFL network from the
RAW deposited GitHub files.  Written from scratch; does NOT read
revision/scripts/01_build_canonical_network.py nor revision/data/canonical_*.tsv
as an input (canonical_*.tsv is read ONLY at the very end for the diff).
"""
import csv, os, re, sys, json, itertools, collections

SIF_DIR  = "/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR/SIF_files"
ATTR_DIR = "/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR/node_attributes"
OUT      = "/path/to/revision/results/v2"
NETS     = ["3-miR", "3-TF", "3-Comp", "4-TF", "5-TF", "6-TF"]

os.makedirs(OUT, exist_ok=True)
log_lines = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); log_lines.append(s)

# ----------------------------------------------------------------- read SIFs
raw_edges = {}          # net -> list of (src, tgt) raw labels, in file order
for net in NETS:
    rows = []
    with open(os.path.join(SIF_DIR, net + ".sif")) as fh:
        for ln, line in enumerate(fh, 1):
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            f = line.split("\t")
            assert len(f) == 3, (net, ln, f)
            rows.append((f[0].strip(), f[2].strip()))
    raw_edges[net] = rows
    P(f"[read] {net}.sif  raw lines = {len(rows)}")

# ------------------------------------------------------- read node attributes
attr_type = collections.defaultdict(set)   # raw label -> set of Type strings seen
attr_nodes_per_net = {}
attr_degree = {}                           # net -> {label: Degree}
for net in NETS:
    names = []
    deg = {}
    with open(os.path.join(ATTR_DIR, net + ".csv"), newline="") as fh:
        for r in csv.DictReader(fh):
            nm = r["name"].strip()
            names.append(nm)
            attr_type[nm].add(r["Type"].strip())
            deg[nm] = int(r["Degree"])
    attr_nodes_per_net[net] = names
    attr_degree[net] = deg
    P(f"[read] {net}.csv  attribute rows = {len(names)}  unique names = {len(set(names))}")

# ------------------------------------------------------------ raw node labels
sif_labels = set()
for net in NETS:
    for s, t in raw_edges[net]:
        sif_labels.add(s); sif_labels.add(t)
attr_labels = set(attr_type)
P(f"\n[raw labels] union over SIF endpoints          = {len(sif_labels)}")
P(f"[raw labels] union over attribute-file names   = {len(attr_labels)}")
P(f"[raw labels] union of both                     = {len(sif_labels | attr_labels)}")
P(f"[raw labels] in attr but not in any SIF        = {len(attr_labels - sif_labels)} -> {sorted(attr_labels - sif_labels)[:10]}")
P(f"[raw labels] in SIF but not in any attr file   = {len(sif_labels - attr_labels)} -> {sorted(sif_labels - attr_labels)[:10]}")

RAW = sif_labels | attr_labels

# --------------------------------------------------------------- node typing
MIR_RE = re.compile(r"^hsa-(mir|miR|let)-")
tf_labels = {n for n, ts in attr_type.items() if "TF" in ts}

def raw_type(label):
    if label in tf_labels:
        return "TF"
    if MIR_RE.match(label):
        return "miRNA"
    return "Gene"

raw_types = {l: raw_type(l) for l in RAW}
P(f"\n[typing] labels called TF in >=1 attribute file = {len(tf_labels)}")
P(f"[typing] raw label type breakdown              = {dict(collections.Counter(raw_types.values()))}")
# sanity: any label with conflicting Type across files?
conflict = {n: sorted(ts) for n, ts in attr_type.items() if len(ts) > 1}
P(f"[typing] labels with conflicting Type across attribute files = {len(conflict)}")
if conflict:
    P("         e.g. " + json.dumps(dict(itertools.islice(conflict.items(), 8))))

# ------------------------------------------------ miRNA identifier harmonisation
# step 1: hsa-mir-X -> hsa-miR-X   (leave hsa-let-* alone)
def case_norm(label):
    if label.startswith("hsa-mir-"):
        return "hsa-miR-" + label[len("hsa-mir-"):]
    return label

case_normed = {l: case_norm(l) for l in RAW}
cn_set = set(case_normed.values())

# step 2: collapse a trailing precursor index -1 / -2 onto the base name,
#         ONLY when the base name also exists in the (case-normalised) data
PREC_RE = re.compile(r"^(hsa-(?:miR|let)-.+?)-([12])$")
def collapse(label):
    m = PREC_RE.match(label)
    if m and m.group(1) in cn_set:
        return m.group(1)
    return label

canon = {l: collapse(case_normed[l]) for l in RAW}
n_case_changed = sum(1 for l in RAW if case_normed[l] != l)
n_prec_changed = sum(1 for l in RAW if collapse(case_normed[l]) != case_normed[l])
P(f"\n[harmonise] labels changed by mir->miR case fix        = {n_case_changed}")
P(f"[harmonise] labels changed by -1/-2 precursor collapse = {n_prec_changed}")
prec_examples = sorted({(l, canon[l]) for l in RAW if collapse(case_normed[l]) != case_normed[l]})
P(f"[harmonise] precursor collapses: {prec_examples}")

CANON_NODES = sorted(set(canon.values()))
# canonical node type: from the raw types of all labels mapping to it
canon_types = {}
type_conflicts = []
for c in CANON_NODES:
    ts = {raw_types[l] for l in RAW if canon[l] == c}
    if len(ts) > 1:
        type_conflicts.append((c, sorted(ts)))
        # precedence TF > miRNA > Gene
        canon_types[c] = "TF" if "TF" in ts else ("miRNA" if "miRNA" in ts else "Gene")
    else:
        canon_types[c] = ts.pop()
P(f"\n[nodes] canonical nodes = {len(CANON_NODES)}")
P(f"[nodes] type breakdown  = {dict(collections.Counter(canon_types.values()))}")
P(f"[nodes] canonical nodes with conflicting collapsed types = {len(type_conflicts)} {type_conflicts[:5]}")

# --------------------------------------------------------------------- edges
CLASS = {
    ("miRNA", "Gene"): "miRNA_target",
    ("miRNA", "TF"):   "miRNA_target",
    ("TF", "miRNA"):   "TF_miRNA",
    ("TF", "Gene"):    "TF_target",
    ("TF", "TF"):      "TF_target",
    ("miRNA", "miRNA"): "miRNA_miRNA",
    ("Gene", "Gene"):  "gene_gene",
    ("Gene", "TF"):    "gene_gene",
    ("Gene", "miRNA"): "gene_miRNA",   # not in spec; flag if it occurs
}

edge_nets = collections.defaultdict(set)   # (s,t) canonical -> set of nets
for net in NETS:
    for s, t in raw_edges[net]:
        edge_nets[(canon[s], canon[t])].add(net)

edges = sorted(edge_nets)
edge_class = {e: CLASS[(canon_types[e[0]], canon_types[e[1]])] for e in edges}
P(f"\n[edges] unique directed canonical edges = {len(edges)}")
P(f"[edges] per class = {json.dumps(dict(collections.Counter(edge_class.values())), indent=None)}")
selfloops = [e for e in edges if e[0] == e[1]]
P(f"[edges] self-loops = {len(selfloops)} {selfloops[:10]}")

# reciprocal TF<->miRNA pairs
eset = set(edges)
recip = set()
for (s, t) in eset:
    if canon_types[s] == "TF" and canon_types[t] == "miRNA" and (t, s) in eset:
        recip.add((s, t))
P(f"[edges] reciprocal TF<->miRNA pairs = {len(recip)}")

# raw-line accounting per network
P("\n[per-network] raw lines / unique raw pairs / unique canonical pairs")
for net in NETS:
    rows = raw_edges[net]
    P(f"   {net:7s} raw={len(rows):5d}  uniq_raw={len(set(rows)):5d}  uniq_canon={len(set((canon[s],canon[t]) for s,t in rows)):5d}")

# ------------------------------------------------- (c) miRBase case duplicates
mir_raw = {l for l in RAW if raw_types[l] == "miRNA"}
P(f"\n[case] raw miRNA labels = {len(mir_raw)}")
lower_groups = collections.defaultdict(set)
for l in mir_raw:
    lower_groups[l.lower()].add(l)
P(f"[case] distinct case-insensitive miRNA names = {len(lower_groups)}")
both_cases = {k: v for k, v in lower_groups.items() if len(v) > 1}
P(f"[case] names present in BOTH miRBase cases   = {len(both_cases)}")

# Jaccard of target sets for the two copies
out_targets = collections.defaultdict(set)
for net in NETS:
    for s, t in raw_edges[net]:
        out_targets[s].add(t)          # RAW labels, un-harmonised

jac_rows = []
for base in ["130a", "21", "124", "29a", "34a"]:
    a, b = f"hsa-mir-{base}", f"hsa-miR-{base}"
    A, B = out_targets.get(a, set()), out_targets.get(b, set())
    j = len(A & B) / len(A | B) if (A | B) else float("nan")
    jac_rows.append((a, b, len(A), len(B), len(A & B), len(A | B), j))
    P(f"[case] {a:14s} n={len(A):4d} | {b:14s} n={len(B):4d} | inter={len(A&B):3d} union={len(A|B):4d} Jaccard={j:.4f}")

# ==================================================== (d) 3-node FFL cores
def build_graph(pairs, types):
    out = collections.defaultdict(set)
    for s, t in pairs:
        if s != t:
            out[s].add(t)
    return out

def count_ffl(pairs, types):
    """Formal 3-node FFL cores.
       TF-FFL       : TF -> miRNA, TF -> g,  miRNA -> g
       miRNA-FFL    : miRNA -> TF, miRNA -> g, TF -> g
       Composite-FFL: TF <-> miRNA (both directions) and both -> g
       Reported both mutually-exclusive (composite removed from the other two)
       and inclusive (overlapping).
    """
    E = {(s, t) for s, t in pairs if s != t}
    out = collections.defaultdict(set)
    for s, t in E:
        out[s].add(t)
    tfs   = [n for n in types if types.get(n) == "TF"]
    mirs  = [n for n in types if types.get(n) == "miRNA"]
    res = dict(tf_incl=0, mir_incl=0, comp=0, tf_excl=0, mir_excl=0)
    inst = []
    for tf in tfs:
        for m in out.get(tf, ()):            # tf -> m
            if types.get(m) != "miRNA":
                continue
            mutual = (m, tf) in E
            shared = (out.get(tf, set()) & out.get(m, set())) - {tf, m}
            k = len(shared)
            res["tf_incl"] += k
            if mutual:
                res["comp"] += k
            else:
                res["tf_excl"] += k
            for g in shared:
                inst.append(("composite" if mutual else "TF-FFL", tf, m, g))
    for m in mirs:
        for tf in out.get(m, ()):            # m -> tf
            if types.get(tf) != "TF":
                continue
            mutual = (tf, m) in E
            shared = (out.get(tf, set()) & out.get(m, set())) - {tf, m}
            k = len(shared)
            res["mir_incl"] += k
            if not mutual:
                res["mir_excl"] += k
                for g in shared:
                    inst.append(("miRNA-FFL", tf, m, g))
    return res, inst

P("\n================ (d) 3-node FFL cores ================")
P(f"{'network':9s} {'TFFFL_x':>8s} {'miRFFL_x':>9s} {'comp':>7s} {'TOTAL_excl':>11s} | {'TFFFL_i':>8s} {'miRFFL_i':>9s} {'TOTAL_incl':>11s}")
ffl_table = {}
for net in NETS:
    pairs = {(canon[s], canon[t]) for s, t in raw_edges[net]}
    ntypes = {n: canon_types[n] for n in {x for p in pairs for x in p}}
    r, _ = count_ffl(pairs, ntypes)
    tot_x = r["tf_excl"] + r["mir_excl"] + r["comp"]
    tot_i = r["tf_incl"] + r["mir_incl"]
    ffl_table[net] = dict(r, total_excl=tot_x, total_incl=tot_i)
    P(f"{net:9s} {r['tf_excl']:8d} {r['mir_excl']:9d} {r['comp']:7d} {tot_x:11d} | {r['tf_incl']:8d} {r['mir_incl']:9d} {tot_i:11d}")

r, _ = count_ffl(eset, canon_types)
tot_x = r["tf_excl"] + r["mir_excl"] + r["comp"]
tot_i = r["tf_incl"] + r["mir_incl"]
ffl_table["UNION"] = dict(r, total_excl=tot_x, total_incl=tot_i)
P(f"{'UNION':9s} {r['tf_excl']:8d} {r['mir_excl']:9d} {r['comp']:7d} {tot_x:11d} | {r['tf_incl']:8d} {r['mir_incl']:9d} {tot_i:11d}")

# also: un-harmonised (raw label) counts, for contrast
P("\n[d, contrast] same counts on RAW un-harmonised labels:")
for net in NETS:
    pairs = set(raw_edges[net])
    ntypes = {n: raw_types[n] for n in {x for p in pairs for x in p}}
    r2, _ = count_ffl(pairs, ntypes)
    P(f"   {net:9s} excl_total={r2['tf_excl']+r2['mir_excl']+r2['comp']:6d} (TF {r2['tf_excl']}, miR {r2['mir_excl']}, comp {r2['comp']})")

# ---------------------------------------------------------------- write files
with open(os.path.join(OUT, "v2_canonical_edges.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["source", "target", "edge_type", "source_type", "target_type", "in_networks"])
    for e in edges:
        w.writerow([e[0], e[1], edge_class[e], canon_types[e[0]], canon_types[e[1]],
                    ";".join(sorted(edge_nets[e]))])
with open(os.path.join(OUT, "v2_canonical_nodes.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["name", "type", "aliases"])
    alias = collections.defaultdict(set)
    for l in RAW:
        alias[canon[l]].add(l)
    for c in CANON_NODES:
        w.writerow([c, canon_types[c], ";".join(sorted(alias[c]))])

json.dump({"ffl_table": ffl_table,
           "n_raw_labels_sif": len(sif_labels),
           "n_raw_labels_attr": len(attr_labels),
           "n_raw_labels_union": len(RAW),
           "n_canon_nodes": len(CANON_NODES),
           "node_types": dict(collections.Counter(canon_types.values())),
           "n_edges": len(edges),
           "edge_classes": dict(collections.Counter(edge_class.values())),
           "n_recip": len(recip),
           "n_mir_raw": len(mir_raw),
           "n_mir_ci": len(lower_groups),
           "n_both_cases": len(both_cases),
           "jaccard": [list(x) for x in jac_rows]},
          open(os.path.join(OUT, "v2_network_summary.json"), "w"), indent=2)

with open("/path/to/revision/logs/v2/v2_build_network.log", "w") as fh:
    fh.write("\n".join(log_lines) + "\n")
