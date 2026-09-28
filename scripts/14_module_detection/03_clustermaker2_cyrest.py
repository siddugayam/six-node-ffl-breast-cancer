#!/usr/bin/env python3
"""
v7 / 02 : drive the REAL clusterMaker2 (v2.3.4) through Cytoscape 3.10.0 + CyREST.

clusterMaker2 is a Cytoscape Java app, not an R/Bioconductor package, so it
cannot be "installed" from CRAN. It was downloaded from the Cytoscape App Store
(https://apps.cytoscape.org/download/clustermaker2/2.3.4), dropped into
~/CytoscapeConfiguration/3/apps/installed/, and Cytoscape was started headless-ish
with `cytoscape.sh -R 1234`. Every clustering below is therefore produced by the
authors' own clusterMaker2 code, not by a reimplementation.

Note on the battery: clusterMaker2 2.3.4 exposes mcl, multilevel (= Louvain),
fastgreedy, infomap, labelpropagation, leadingeigenvector, leiden, glay and
others, but NOT walktrap. Walktrap is therefore run in igraph (script 03).
"""
import json, sys, time, urllib.request, urllib.error, os

BASE = "http://localhost:1234/v1"
OUT  = "/path/to/revision/results/v7"
SIF  = os.path.join(OUT, "cm2_network.sif")


def post(path, payload):
    req = urllib.request.Request(BASE + path,
                                 data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"},
                                 method="POST")
    with urllib.request.urlopen(req, timeout=1800) as r:
        body = r.read().decode()
    return json.loads(body) if body.strip() else {}


def get(path):
    with urllib.request.urlopen(BASE + path, timeout=600) as r:
        return json.loads(r.read().decode())


def load_network():
    res = post("/commands/network/load%20file", {"file": SIF})
    suid = res["data"]["networks"][0]
    print("loaded network SUID", suid, flush=True)
    return suid


def node_table(suid):
    rows = get("/networks/%d/tables/defaultnode/rows" % suid)
    return rows


def run(cmd, args, label):
    t0 = time.time()
    try:
        post("/commands/cluster/" + cmd, args)
    except urllib.error.HTTPError as e:
        print("FAILED %s: %s" % (label, e.read().decode()[:500]), flush=True)
        return None
    print("  %-28s %6.1f s" % (label, time.time() - t0), flush=True)
    return label


def main():
    suid = load_network()
    runs = []

    # --- MCL, inflation sweep 1.8 - 4.0 -------------------------------------
    # clusterMaker2 defaults: iterations 16, clusteringThresh 1e-15,
    # edgeCutOff 1e-15, adjustLoops true, maxResidual 1e-3.
    for infl in [1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0]:
        attr = "mcl_I%s" % str(infl).replace(".", "")
        runs.append(run("mcl", {
            "network": "SUID:%d" % suid,
            "inflation_parameter": infl,
            "adjustLoops": "true",
            "undirectedEdges": "true",
            "clusterAttribute": attr,
            "createGroups": "false",
            "showUI": "false",
            "restoreEdges": "false",
            "selectedOnly": "false",
        }, attr))

    # --- the igraph-backed clusterMaker2 algorithms -------------------------
    runs.append(run("multilevel", {"clusterAttribute": "louvain_cm2",
                                   "createGroups": "false", "showUI": "false",
                                   "isSynchronous": "true"}, "louvain_cm2"))
    runs.append(run("fastgreedy", {"clusterAttribute": "fastgreedy_cm2",
                                   "createGroups": "false", "showUI": "false",
                                   "isSynchronous": "true"}, "fastgreedy_cm2"))
    for tr in [1, 10, 100]:
        attr = "infomap_cm2_t%d" % tr
        runs.append(run("infomap", {"clusterAttribute": attr, "trials": tr,
                                    "createGroups": "false", "showUI": "false",
                                    "isSynchronous": "true"}, attr))
    runs.append(run("labelpropagation", {"clusterAttribute": "labelprop_cm2",
                                         "createGroups": "false", "showUI": "false",
                                         "isSynchronous": "true"}, "labelprop_cm2"))
    runs.append(run("glay", {"network": "SUID:%d" % suid,
                             "clusterAttribute": "glay_cm2",
                             "createGroups": "false", "showUI": "false",
                             "undirectedEdges": "true",
                             "restoreEdges": "false"}, "glay_cm2"))
    runs.append(run("leadingeigenvector", {"clusterAttribute": "leadeig_cm2",
                                           "createGroups": "false", "showUI": "false",
                                           "isSynchronous": "true"}, "leadeig_cm2"))
    for res in [0.5, 1.0, 2.0]:
        attr = "leiden_cm2_r%s" % str(res).replace(".", "")
        runs.append(run("leiden", {"clusterAttribute": attr,
                                   "resolution_parameter": res,
                                   "objective_function": "Modularity",
                                   "createGroups": "false", "showUI": "false",
                                   "isSynchronous": "true"}, attr))

    rows = node_table(suid)
    cols = sorted({k for r in rows for k in r})
    with open(os.path.join(OUT, "cm2_cytoscape_partitions.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print("wrote cm2_cytoscape_partitions.tsv with columns:", cols, flush=True)


if __name__ == "__main__":
    main()
