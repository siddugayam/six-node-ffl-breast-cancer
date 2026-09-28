import sys, csv, re, time, io
import requests
from pypdf import PdfReader

RES = '/path/to/revision/results/v4/'
GENES = ["PTEN", "SMAD4", "TGFBR2", "DICER1", "KLF4", "XIAP", "COL1A1", "COL3A1",
         "NFKB1", "RELA", "SP1", "ETS1", "VEGFA", "CCND2", "MYC", "EZH2"]
HOST = "http://gepia2.cancer-pku.cn/"
S = requests.Session()
S.headers.update({"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) Chrome/151.0.0.0 Safari/537.36",
                  "Referer": HOST})
now = lambda: time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())

def post(ep, d, tries=4, timeout=90):
    global S
    for i in range(tries):
        try:
            r = S.post(HOST + "assets/" + ep, data=d, timeout=timeout)
            if r.status_code == 200 and r.text.strip():
                return r
        except Exception:
            # the host silently drops idle keep-alive sockets; rebuild the session and retry
            S = requests.Session()
            S.headers.update({"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) "
                                            "Chrome/151.0.0.0 Safari/537.36", "Referer": HOST})
        time.sleep(3 * (i + 1))
    return None

# ---------- 1. expression: median TPM tumour vs normal (TCGA tumour vs TCGA+GTEx normal) ----------
def expression():
    F = ["portal", "assay", "gene", "cancer", "tumour_median_TPM", "normal_median_TPM",
         "log2FC_tumour_vs_normal", "direction", "note", "retrieved_utc"]
    fh = open(RES + "portals_gepia2_expression.csv", 'w', newline='')
    w = csv.DictWriter(fh, fieldnames=F); w.writeheader()
    import math
    for g in GENES:
        r = post("PHP2/barplot.php", {"gene": g, "log": "1", "jitter": "0.4",
                                      "pvalue": "0.01", "logfc": "1"})
        if r is None:
            w.writerow(dict(portal="GEPIA2", assay="median expression (TPM)", gene=g,
                            cancer="", tumour_median_TPM="", normal_median_TPM="",
                            log2FC_tumour_vs_normal="", direction="", note="RETRIEVAL_FAILED",
                            retrieved_utc=now())); continue
        t = re.search(r'var t=\[([^\]]*)\]', r.text)
        n = re.search(r'var n=\[([^\]]*)\]', r.text)
        ty = re.search(r'var types=\[([^\]]*)\]', r.text)
        if not (t and n and ty):
            w.writerow(dict(portal="GEPIA2", assay="median expression (TPM)", gene=g,
                            cancer="", tumour_median_TPM="", normal_median_TPM="",
                            log2FC_tumour_vs_normal="", direction="", note="UNPARSEABLE",
                            retrieved_utc=now())); continue
        tv = [float(x) for x in t.group(1).split(',')]
        nv = [float(x) for x in n.group(1).split(',')]
        tys = [x.strip().strip('"') for x in ty.group(1).split(',')]
        for c, a, b in zip(tys, tv, nv):
            lfc = math.log2((a + 1e-9) / (b + 1e-9)) if b > 0 else ''
            w.writerow(dict(portal="GEPIA2", assay="median expression (TPM)", gene=g, cancer=c,
                            tumour_median_TPM=a, normal_median_TPM=b,
                            log2FC_tumour_vs_normal=round(lfc, 4) if lfc != '' else '',
                            direction=("UP in tumour" if a > b else "DOWN in tumour"),
                            note="TCGA tumours vs TCGA+GTEx normals; medians read from the "
                                 "GEPIA2 pan-cancer bar plot",
                            retrieved_utc=now()))
        fh.flush()
        print("expr", g, "BRCA t/n:", tv[tys.index("BRCA")], nv[tys.index("BRCA")], flush=True)
        time.sleep(1)
    fh.close()

# ---------- 2. survival ----------
PDF = re.compile(r'Logrank p=([0-9.eE+\-]+)\s*HR\(high\)=([0-9.eE+\-]+)\s*p\(HR\)=([0-9.eE+\-]+)'
                 r'\s*n\(high\)=(\d+)\s*n\(low\)=(\d+)', re.S)

def survival():
    F = ["portal", "module", "gene", "dataset", "endpoint", "split", "logrank_p", "HR_high",
         "p_HR", "n_high", "n_low", "pdf_url", "error", "retrieved_utc"]
    fh = open(RES + "portals_gepia2_survival.csv", 'w', newline='')
    w = csv.DictWriter(fh, fieldnames=F); w.writeheader()
    for g in GENES:
        for meth, ename in [("os", "OS"), ("dfs", "DFS/RFS")]:
            for c1, c2, sname in [("50", "50", "median"), ("75", "25", "quartile 75/25")]:
                d = {"methodoption": meth, "dataset": "BRCA", "signature": g,
                     "highcol": "#ff0000", "lowcol": "#0000ff", "groupcutoff1": c1,
                     "groupcutoff2": c2, "axisunit": "month", "ifhr": "hr", "ifconf": "conf",
                     "signature_norm": "", "is_sub": "false", "subtype": ""}
                r = post("PHP4/survival_zf.php", d)
                row = dict(portal="GEPIA2", module="Survival Analysis (TCGA-BRCA)", gene=g,
                           dataset="BRCA", endpoint=ename, split=sname, logrank_p="", HR_high="",
                           p_HR="", n_high="", n_low="", pdf_url="", error="",
                           retrieved_utc=now())
                try:
                    j = r.json()
                    od = j.get("outdir", "")
                    if od in ("", "fail"):
                        row["error"] = "gepia returned outdir=%r" % od
                    else:
                        u = HOST + "tmp/" + od
                        row["pdf_url"] = u
                        pr = S.get(u, timeout=120)
                        txt = PdfReader(io.BytesIO(pr.content)).pages[0].extract_text()
                        txt = re.sub(r'\s+', ' ', txt)
                        m = PDF.search(txt)
                        if m:
                            (row["logrank_p"], row["HR_high"], row["p_HR"],
                             row["n_high"], row["n_low"]) = m.groups()
                        else:
                            row["error"] = "unparsed PDF text: " + txt[-160:]
                except Exception as e:
                    row["error"] = repr(e)[:200]
                w.writerow(row); fh.flush()
                print("surv", g, ename, sname, row["logrank_p"], row["HR_high"],
                      row["n_high"], row["n_low"], row["error"][:60], flush=True)
                time.sleep(1.5)
    fh.close()

expression()
survival()
