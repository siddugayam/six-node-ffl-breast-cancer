"""bc-GenExMiner v5.3 targeted differential expression (mode=8), GTEx+TCGA RNA-seq source
   (Gene2=G2RT), splitting method = nature of the tissue (healthy / tumour-adjacent / tumour)."""
import sys, csv, re, html, time
sys.path.insert(0, '/path/to/revision/scripts/v4')
import bcgem_lib as G

RES = '/path/to/revision/results/v4/'
GENES = ["PTEN", "SMAD4", "TGFBR2", "DICER1", "KLF4", "XIAP", "COL1A1", "COL3A1",
         "NFKB1", "RELA", "SP1", "ETS1", "VEGFA", "CCND2", "MYC", "EZH2"]

def parse(hh):
    out = {"groups": {}, "welch": "", "dtk": []}
    m = re.search(r'<table id="TableRes_Stats_NAT_1_tissu".*?</table>', hh, re.S)
    if m:
        for tr in re.findall(r'<tr[^>]*>(.*?)</tr>', m.group(0), re.S):
            c = G.cells(tr)
            if len(c) == 8 and c[0] in ("Healthy", "Tumour-adjacent", "Tumour"):
                out["groups"][c[0]] = dict(zip(["min", "q1", "median", "mean", "q3", "max", "sd"],
                                               c[1:]))
    m = re.search(r'<table id="TableRes_tissu_1".*?</table>', hh, re.S)
    if m:
        w = re.search(r'p\(Welch\)\s*([^\n]{1,25})', G.TAG.sub('', m.group(0)))
        if w:
            out["welch"] = html.unescape(w.group(1)).strip()
        for tr in re.findall(r'<tr[^>]*>(.*?)</tr>', m.group(0), re.S):
            c = G.cells(tr)
            if len(c) == 5 and c[0] in ("Healthy", "Tumour-adjacent", "Tumour"):
                out["dtk"].append((c[0], c[1].rstrip('/'), c[2], c[4]))
    return out

F = ["portal", "module", "data_source", "gene", "comparison_group", "min", "q1", "median",
     "mean", "q3", "max", "sd", "p_Welch", "DTK_comparison", "DTK_p", "direction_tumour_vs_healthy",
     "error", "retrieved_utc"]

def main():
    fh = open(RES + 'portals_bcgenexminer_diffexpr.csv', 'w', newline='')
    w = csv.DictWriter(fh, fieldnames=F); w.writeheader()
    now = lambda: time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    for g in GENES:
        s = G.session(8)
        d = {"Gene2": "G2RT", "Gene2Text": g, "Gene2Text2": g, "MP": "MPN", "Node": "NodeA",
             "ER": "ERA", "HER2": "HER2A", "PR": "PRA", "Prolif": "ProlifA", "Event": "EventAE",
             "CorrelMap": "CorrelMap_None", "Split": "SplitMed", "MSM": "MSM_All",
             "ExpressMap": "ExpressMap_NAT", "OutputFig": "OutputFigWhiskersBox",
             "Verifier": "Submit"}
        v, r = G.two_step(s, d)
        p = parse(r) if r else {"groups": {}, "welch": "", "dtk": []}
        if not p["groups"]:
            w.writerow(dict(portal="bc-GenExMiner v5.3", module="targeted differential expression (mode=8)",
                            data_source="GTEx & TCGA RNA-seq", gene=g, comparison_group="",
                            min="", q1="", median="", mean="", q3="", max="", sd="",
                            p_Welch="", DTK_comparison="", DTK_p="",
                            direction_tumour_vs_healthy="", error="RETRIEVAL_FAILED_OR_NO_DATA",
                            retrieved_utc=now()))
            print(g, "FAILED", flush=True); fh.flush(); time.sleep(3); continue
        dtk = {(a, c): (sym, pp) for a, sym, c, pp in p["dtk"]}
        sym, pp = dtk.get(("Tumour", "Healthy"), ("", ""))
        dirn = {"<<": "DOWN in tumour", ">>": "UP in tumour", "=": "no difference"}.get(sym, sym)
        for grp, st in p["groups"].items():
            w.writerow(dict(portal="bc-GenExMiner v5.3",
                            module="targeted differential expression (mode=8)",
                            data_source="GTEx & TCGA RNA-seq", gene=g, comparison_group=grp,
                            min=st["min"], q1=st["q1"], median=st["median"], mean=st["mean"],
                            q3=st["q3"], max=st["max"], sd=st["sd"], p_Welch=p["welch"],
                            DTK_comparison="Tumour vs Healthy", DTK_p=pp,
                            direction_tumour_vs_healthy=dirn, error="", retrieved_utc=now()))
        fh.flush()
        print(g, "medians:", {k: v["median"] for k, v in p["groups"].items()},
              "welch", p["welch"], "| T vs H:", sym, pp, flush=True)
        time.sleep(3)
    fh.close()
main()
