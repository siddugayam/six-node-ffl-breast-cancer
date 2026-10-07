#!/usr/bin/env python3
# Sub-type numbering of Mangan and Alon (2003), signs (R->M, M->T, R->T).
"""
02_ffl_census.py -- formal definition, exact census and Alon coherence typing of
n-node feed-forward loops in the canonical breast-cancer miRNA-TF-gene network.

FORMAL DEFINITION
-------------------------------------------------------------------
Let G be the canonical DIRECTED, SIGNED, TYPED graph built by scripts/03_ffl_census/ffl_graph.py
(canonical_edges.tsv augmented with layer_TF_target / layer_TF_miRNA /
layer_gene_gene / layer_miRNA_miRNA).  Undirected evidence (miRNA-miRNA
polycistronic co-transcription; optionally STRING association) is stored as a
reciprocal ARC PAIR, and no single directed path may traverse more than one
undirected edge.

An n-node FFL is a vertex-induced connected subgraph S on n nodes for which
there EXISTS an orientation of the reciprocal arc pairs ("collapsing" each
reciprocal TF<->miRNA pair into one composite edge whose effective direction is
fixed by the circuit) such that

  (i)   S contains at least one 3-node FFL core: R->M, R->T, M->T ;
  (ii)  S is a directed acyclic graph ;
  (iii) S has exactly one source (in-degree 0 in S) and exactly one sink ;
  (iv)  >= 2 internally vertex-disjoint directed paths run from that source to
        that sink.

Condition (iv) preserves the defining feed-forward property (parallel paths from
one source to one target) and is what excludes a cascade with a side branch.

Enumerating admissible orientations is equivalent to enumerating topological
orders with the source first and the sink last and orienting every edge forward;
this is complete and sound (any valid orientation is recovered by its own
topological order) and is how both implementations work.

IMPLEMENTATIONS AND CROSS-VALIDATION
------------------------------------
  scripts/03_ffl_census/ffl_enum.c   exact parallel enumerator, drives on ordered (source,sink)
                       pairs and (n-2)-subsets of the nodes that lie on a
                       source->sink path of length <= n-1.
  scripts/03_ffl_census/ffl_def.py   independent reference implementation used for the
                       hand-worked checks, the exemplar audit and cross-validation.
The two agree exactly at n=3 (6725) and n=4 (71307) although they enumerate over
different candidate spaces.

OUTPUTS (results/)
  ffl_census.csv          per architecture class: n_nodes, ffl_class, architecture,
                          count, n_distinct_nodes
  ffl_3node_cores.csv     every 3-node core with signs, C1-C4/I1-I4 type, evidence
  ffl_higher_order.csv    4/5/6-node modules: ALL 4-node modules, and up to 5000
                          instances per architecture class for n=5 and n=6 (the
                          exact totals are in ffl_census.csv; a complete listing
                          is impossible at n=6).
  ffl_definition_check.md how the exemplar circuits of the original submission score
"""
import csv, os, sys, subprocess, itertools, json
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(os.path.dirname(HERE))   # repository root
RES = os.path.join(ROOT, "results")
LOG = os.path.join(ROOT, "logs")
SCRATCH = os.environ.get("FFL_SCRATCH",
    "/path/to/scratch")

import ffl_graph as fg
import ffl_cores
import ffl_def

ICN = ["TF-TF", "TF-gene", "TF-miRNA", "miRNA-gene", "miRNA-TF",
       "miRNA-miRNA", "gene-TF", "gene-miRNA", "gene-gene"]
FCN = ["Composite-FFL", "TF-FFL", "miRNA-FFL", "gene-FFL"]
NAMED6 = {"TF-TF", "TF-miRNA", "miRNA-miRNA", "miRNA-gene", "TF-gene", "gene-gene"}

logfh = open(os.path.join(LOG, "ffl_census.log"), "a")
def log(*a):
    s = " ".join(str(x) for x in a)
    print(s); logfh.write(s + "\n"); logfh.flush()


# ------------------------------------------------------------------ graph ---
def export_graph(path, **kw):
    nodes, ntype, arcs, meta = fg.build(verbose=False, **kw)
    tmap = {"TF": 0, "Gene": 1, "miRNA": 2}
    idx = {n: i for i, n in enumerate(nodes)}
    with open(path, "w") as fh:
        fh.write("%d %d\n" % (len(nodes), len(arcs)))
        for n in nodes:
            fh.write("%s %d\n" % (n, tmap[ntype[n]]))
        for (u, v), a in arcs.items():
            fh.write("%d %d %d\n" % (idx[u], idx[v], 1 if a["undirected"] else 0))
    return nodes, ntype, arcs, meta


# ------------------------------------------------------- 3-node core table ---
def sign_of(u, v, arcs, ann):
    """(sign, source_of_sign, is_assumed)"""
    a = arcs[(u, v)]
    if a["layer"] == "miRNA_target":
        return -1, "miRNA_repression", False
    if a["layer"].startswith("miRNA_miRNA"):
        return 1, "assumed_cotranscription", True
    k = ann.get((u, v))
    if k and k["curated_sign"] is not None:
        return k["curated_sign"], ("TRRUST" if k["layer"] == "TF_target" else "TransmiR"), False
    return 1, "assumed_activation", True


def write_cores(nodes, ntype, arcs):
    ann = ffl_cores.load_edge_annotation()
    cores = ffl_cores.enumerate_cores(nodes, ntype, arcs)
    und = set(e for e, a in arcs.items() if a["undirected"])
    A = set(arcs)
    valid_sets = set()
    path = os.path.join(RES, "ffl_3node_cores.csv")
    cnt_all, cnt_cur = Counter(), Counter()
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["regulator", "intermediate", "target",
                    "regulator_type", "intermediate_type", "target_type",
                    "ffl_class", "coherence_type", "coherent",
                    "sign_R_M", "sign_M_T", "sign_R_T",
                    "sign_source_R_M", "sign_source_M_T", "sign_source_R_T",
                    "any_assumed_sign", "coherence_type_curated_only",
                    "evidence_R_M", "evidence_M_T", "evidence_R_T",
                    "layer_R_M", "layer_M_T", "layer_R_T",
                    "uses_undirected_edge", "satisfies_full_definition"])
        cache = {}
        for (R, M, T) in cores:
            s1, o1, a1 = sign_of(R, M, arcs, ann)
            s2, o2, a2 = sign_of(M, T, arcs, ann)
            s3, o3, a3 = sign_of(R, T, arcs, ann)
            ty = ffl_cores.TYPE_TABLE[(s1, s2, s3)]
            assumed = a1 or a2 or a3
            cnt_all[ty] += 1
            cnt_cur["unclassifiable" if assumed else ty] += 1
            key = frozenset((R, M, T))
            if key not in cache:
                cache[key] = ffl_def.test_set([R, M, T], A, und) is not None
            if cache[key]:
                valid_sets.add(key)
            uu = any(arcs[e]["undirected"] for e in ((R, M), (M, T), (R, T)))
            ev = lambda u, v: (ann.get((u, v), {}).get("ev", "NA"))
            w.writerow([R, M, T, ntype[R], ntype[M], ntype[T],
                        ffl_cores.core_class(R, M, T, ntype, arcs), ty,
                        "TRUE" if ty.startswith("C") else "FALSE",
                        s1, s2, s3, o1, o2, o3,
                        "TRUE" if assumed else "FALSE",
                        "unclassifiable" if assumed else ty,
                        ev(R, M), ev(M, T), ev(R, T),
                        arcs[(R, M)]["layer"], arcs[(M, T)]["layer"], arcs[(R, T)]["layer"],
                        "TRUE" if uu else "FALSE",
                        "TRUE" if cache[key] else "FALSE"])
    log("[cores] %d ordered 3-node cores spanning %d distinct node sets; %d node sets "
        "satisfy the full definition" % (len(cores), len(cache), len(valid_sets)))
    log("[cores] coherence (assumption-filled signs): " +
        " ".join("%s=%d" % (t, cnt_all[t]) for t in ["C1","C2","C3","C4","I1","I2","I3","I4"]))
    write_coherence_summary()
    log("[cores] coherence (curated signs only)     : " +
        " ".join("%s=%d" % (t, cnt_cur[t]) for t in ["C1","C2","C3","C4","I1","I2","I3","I4"]) +
        "  unclassifiable=%d" % cnt_cur["unclassifiable"])
    return cores, cnt_all, cnt_cur


# ------------------------------------------------------------- census file ---
SIGNS = {"C1":"(+,+,+)","C2":"(-,+,-)","C3":"(+,-,-)","C4":"(-,-,+)",
         "I1":"(+,-,+)","I2":"(-,-,-)","I3":"(+,+,-)","I4":"(-,+,+)"}


def write_coherence_summary():
    """results/ffl_coherence_summary.tsv -- C1-C4 / I1-I4 counts overall and per class."""
    order = ["C1","C2","C3","C4","I1","I2","I3","I4"]
    rows = list(csv.DictReader(open(os.path.join(RES, "ffl_3node_cores.csv"))))
    tot, cur = Counter(), Counter()
    byclass = defaultdict(Counter)
    canon = Counter()
    for r in rows:
        tot[r["coherence_type"]] += 1
        if r["any_assumed_sign"] == "FALSE":
            cur[r["coherence_type"]] += 1
        byclass[r["ffl_class"]][r["coherence_type"]] += 1
        if r["regulator_type"] == "TF" and r["intermediate_type"] == "miRNA" and r["target_type"] == "Gene":
            canon[r["coherence_type"]] += 1
    with open(os.path.join(RES, "ffl_coherence_summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["ffl_class", "coherence_type", "signs_R_M__M_T__R_T", "coherent",
                    "n_cores_all", "n_cores_fully_curated_signs"])
        for t in order:
            w.writerow(["ALL", t, SIGNS[t], "TRUE" if t[0] == "C" else "FALSE", tot[t], cur[t]])
        w.writerow(["ALL", "COHERENT_TOTAL", "", "TRUE",
                    sum(tot[t] for t in order[:4]), sum(cur[t] for t in order[:4])])
        w.writerow(["ALL", "INCOHERENT_TOTAL", "", "FALSE",
                    sum(tot[t] for t in order[4:]), sum(cur[t] for t in order[4:])])
        w.writerow(["ALL", "UNCLASSIFIABLE_assumed_sign", "", "NA",
                    sum(1 for r in rows if r["any_assumed_sign"] == "TRUE"), 0])
        for t in order:
            if canon[t]:
                w.writerow(["canonical_TF->miRNA->gene_with_TF->gene", t, SIGNS[t],
                            "TRUE" if t[0] == "C" else "FALSE", canon[t], ""])
        for k in sorted(byclass, key=lambda k: -sum(byclass[k].values())):
            for t in order:
                if byclass[k][t]:
                    w.writerow([k, t, SIGNS[t], "TRUE" if t[0] == "C" else "FALSE", byclass[k][t], ""])
    log("[cores] canonical TF->miRNA->gene + TF->gene cores: %d, of which I1 (textbook "
        "incoherent type-1 miRNA-mediated FFL) = %d" % (sum(canon.values()), canon["I1"]))


def decode(key):
    icl = key & 0x1FF
    nTF = (key >> 9) & 15
    nmiR = (key >> 13) & 15
    fc = (key >> 17) & 3
    return icl, nTF, nmiR, fc


def build_census(sizes, prefix):
    rows = []
    for n in sizes:
        f = "%s%d_census.tsv" % (prefix, n)
        if not os.path.exists(f):
            log("[census] MISSING %s -- skipped" % f); continue
        for r in csv.DictReader(open(f), delimiter="\t"):
            icl, nTF, nmiR, fc = decode(int(r["arch_key"]))
            classes = [ICN[i] for i in range(9) if icl >> i & 1]
            rows.append(dict(
                n_nodes=n, ffl_class=FCN[fc],
                architecture="+".join(classes),
                count=int(r["count"]), n_distinct_nodes=int(r["n_distinct_nodes"]),
                node_composition="TF=%d;Gene=%d;miRNA=%d" % (nTF, n - nTF - nmiR, nmiR),
                n_interaction_classes=len(classes),
                has_all_six_named_classes="TRUE" if NAMED6 <= set(classes) else "FALSE",
                arch_key=int(r["arch_key"])))
    rows.sort(key=lambda r: (r["n_nodes"], -r["count"]))
    p = os.path.join(RES, "ffl_census.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows: w.writerow(r)
    for n in sizes:
        sub = [r for r in rows if r["n_nodes"] == n]
        if not sub: continue
        log("[census] n=%d : %d architecture classes, %d modules, max interaction "
            "classes in one module = %d, modules with all six named classes = %d"
            % (n, len(sub), sum(r["count"] for r in sub),
               max(r["n_interaction_classes"] for r in sub),
               sum(r["count"] for r in sub if r["has_all_six_named_classes"] == "TRUE")))
    return rows


def build_higher_order(sizes, prefix, caps):
    p = os.path.join(RES, "ffl_higher_order.csv")
    seen = Counter()
    n_written = Counter()
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["n_nodes", "members", "source", "sink", "n_disjoint_paths",
                    "ffl_class", "architecture", "node_composition", "arch_key",
                    "listing"])
        for n in sizes:
            f = "%s%d_instances.tsv" % (prefix, n)
            if not os.path.exists(f):
                log("[higher] MISSING %s" % f); continue
            cap = caps.get(n, 500)
            for r in csv.DictReader(open(f), delimiter="\t"):
                key = int(r["arch_key"])
                if cap is not None and seen[(n, key)] >= cap:
                    seen[(n, key)] += 1
                    continue
                seen[(n, key)] += 1
                icl, nTF, nmiR, fc = decode(key)
                mem = [r["source"]] + r["members"].split(";") + [r["sink"]]
                w.writerow([n, ";".join(mem), r["source"], r["sink"],
                            r["n_disjoint_paths"], FCN[fc],
                            "+".join(ICN[i] for i in range(9) if icl >> i & 1),
                            "TF=%d;Gene=%d;miRNA=%d" % (nTF, n - nTF - nmiR, nmiR),
                            key,
                            "complete" if cap is None else "capped_%d_per_architecture" % cap])
                n_written[n] += 1
    for n in sizes:
        log("[higher] n=%d : %d instance rows written" % (n, n_written[n]))
    return n_written


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    gpath = os.path.join(SCRATCH, "graph_primary.txt")
    nodes, ntype, arcs, meta = export_graph(gpath)
    log("[graph] " + "; ".join(meta["log"]))
    log("[graph] arcs=%d nodes=%d" % (len(arcs), len(nodes)))
    if what in ("all", "cores"):
        write_cores(nodes, ntype, arcs)
    if what in ("all", "census"):
        sizes = [int(x) for x in os.environ.get("FFL_SIZES", "3,4,5,6").split(",")]
        prefix = os.path.join(SCRATCH, "n")
        build_census(sizes, prefix)
        build_higher_order([s for s in sizes if s >= 4], prefix,
                           {4: None, 5: 5000, 6: 5000})
