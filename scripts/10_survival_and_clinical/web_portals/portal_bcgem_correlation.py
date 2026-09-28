"""bc-GenExMiner v5.3 targeted gene-correlation (mode=4) over the 16 cross-check genes.
   Submits one request per correlation map (all patients + subtype-stratified maps) and
   downloads the portal's own CSV of pairwise Pearson r / p / n."""
import sys, re, csv, time, urllib.parse
sys.path.insert(0, '/path/to/revision/scripts/10_survival_and_clinical/web_portals')
import bcgem_lib as G

RES = '/path/to/revision/results/v4/'
GENES = ["PTEN", "SMAD4", "TGFBR2", "DICER1", "KLF4", "XIAP", "COL1A1", "COL3A1",
         "NFKB1", "RELA", "SP1", "ETS1", "VEGFA", "CCND2", "MYC", "EZH2"]
MAPS = [("CorrelMap_All", "all patients"), ("CorrelMap_ER", "by ER status"),
        ("CorrelMap_PAM50", "by PAM50 subtype"), ("CorrelMap_TNBC", "by TNBC status")]
F = ["portal", "module", "data_source", "correlation_map", "population", "gene_1", "gene_2",
     "pearson_r", "p_value", "n", "source_csv", "retrieved_utc"]

def main():
    fh = open(RES + 'portals_bcgenexminer_correlation.csv', 'w', newline='')
    w = csv.DictWriter(fh, fieldnames=F); w.writeheader()
    now = lambda: time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    for cmap, cname in MAPS:
        s = G.session(4)
        gl = ";".join(GENES)
        d = {"Gene2": "G2GS", "Gene2Text": gl, "Gene2Text2": gl, "MP": "MPN", "Node": "NodeA",
             "ER": "ERA", "HER2": "HER2A", "PR": "PRA", "Prolif": "ProlifA", "Event": "EventAE",
             "Split": "SplitMed", "ExpressMap": "ExpressMap_All", "MSM": "MSM_All",
             "CorrelMap": cmap, "OutputFig": "OutputFigHBP", "Verifier": "Submit"}
        v, r = G.two_step(s, d)
        if not r:
            print(cmap, "FAILED", flush=True); continue
        # 1) the portal's own CSV exports
        csvs = re.findall(r'href="(\./R/R_csv/[^"]+\.csv)"', r)
        got = 0
        for c in dict.fromkeys(csvs):
            url = G.B + urllib.parse.quote(c[2:], safe="/._-")
            try:
                rr = s.get(url, timeout=300)
                if rr.status_code != 200 or not rr.text.strip():
                    continue
            except Exception:
                continue
            pop = re.search(r'Gene-Correlation-Analysis_([^_]+)_', c)
            pop = pop.group(1) if pop else cname
            # portal CSV layout: ';'-separated, first col = row gene, second col = statistic name
            # ("Pearson's correlation coefficient" / "adj. p-value" / "No. patients"),
            # remaining columns = the other genes.
            lines = [l for l in rr.text.replace('\r', '').split('\n') if l.strip()]
            hdr = lines[0].split(';')
            cols = [x.strip() for x in hdr[2:]]
            block = {}
            for l in lines[1:]:
                f = [x.strip() for x in l.split(';')]
                if len(f) < 3:
                    continue
                block.setdefault(f[0], {})[f[1]] = f[2:]
            for g1, stats in block.items():
                rvals = stats.get("Pearson's correlation coefficient", [])
                pvals = stats.get('adj. p-value', [])
                nvals = stats.get('No. patients', [])
                for i, g2 in enumerate(cols):
                    if g2 == g1 or i >= len(rvals) or rvals[i] in ('NA', ''):
                        continue
                    w.writerow(dict(portal="bc-GenExMiner v5.3",
                                    module="targeted gene correlation (mode=4)",
                                    data_source="all DNA microarray data",
                                    correlation_map=cname, population=pop,
                                    gene_1=g1, gene_2=g2, pearson_r=rvals[i],
                                    p_value=pvals[i] if i < len(pvals) else '',
                                    n=nvals[i] if i < len(nvals) else '',
                                    source_csv=url, retrieved_utc=now()))
                    got += 1
        # 2) fallback: scrape the pairwise blocks straight off the page
        if got == 0:
            for m in re.finditer(r"Pearson's pairwise correlation plot for ([^:]+):\s*"
                                 r"<span class=\"geneSymbol\">([^<]+)</span>\s*<i>versus</i>\s*"
                                 r"<span class=\"geneSymbol\">([^<]+)</span>(.*?)</table>", r, re.S):
                pop, g1, g2, blk = m.groups()
                vals = re.findall(r'align="right" nowrap>([^<]+)</td>', blk)
                if len(vals) >= 3:
                    w.writerow(dict(portal="bc-GenExMiner v5.3",
                                    module="targeted gene correlation (mode=4)",
                                    data_source="all DNA microarray data", correlation_map=cname,
                                    population=pop.strip(), gene_1=g1.strip(), gene_2=g2.strip(),
                                    pearson_r=vals[0].strip(), p_value=vals[1].strip(),
                                    n=vals[2].strip().replace(' ', ' '),
                                    source_csv="scraped from result page", retrieved_utc=now()))
                    got += 1
        fh.flush()
        print(cmap, "rows", got, flush=True)
        time.sleep(4)
    fh.close()
main()
