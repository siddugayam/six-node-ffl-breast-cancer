import sys, os, csv, time, re
sys.path.insert(0, '/path/to/revision/scripts/10_survival_and_clinical/web_portals')
import ualcan_lib as U

RES = '/path/to/revision/results/v4/'
GENES = ["PTEN", "SMAD4", "TGFBR2", "DICER1", "KLF4", "XIAP", "COL1A1", "COL3A1",
         "NFKB1", "RELA", "SP1", "ETS1", "VEGFA", "CCND2", "MYC", "EZH2"]

def grouping(cats):
    j = ' | '.join(cats)
    if re.search(r'TNBC-BL1', j): return "Major subclasses (with TNBC types)"
    if re.search(r'Stage1', j): return "Individual cancer stages"
    if re.search(r'Caucasian', j): return "Patient race"
    if re.search(r'\bMale\b', j): return "Patient gender"
    if re.search(r'Yrs', j): return "Patient age"
    if re.search(r'Menopause', j): return "Menopause status"
    if re.search(r'Triple negative|TNBC', j): return "Major subclasses"
    if re.search(r'\bIDC\b', j): return "Tumor histology"
    if re.search(r'\bN0\b', j): return "Nodal metastasis status"
    if re.search(r'TP53', j): return "TP53 mutation status"
    if len(cats) == 2 and 'Normal' in j: return "Sample types"
    return "unclassified: " + j[:60]

BOXF = ["portal", "assay", "dataset", "gene", "grouping", "category", "n", "unit",
        "low_whisker", "q1", "median", "q3", "high_whisker", "url", "retrieved_utc"]
STATF = ["portal", "assay", "dataset", "gene", "grouping", "comparison", "p_value",
         "url", "retrieved_utc"]

JOBS = [
    ("expression", "TCGA-BRCA", "TCGAExResultNew2.pl?genenam={g}&ctype=BRCA",
     RES + "portals_ualcan_expression.csv"),
    ("promoter_methylation", "TCGA-BRCA", "TCGA-methyl-Result.pl?genenam={g}&ctype=BRCA",
     RES + "portals_ualcan_methylation.csv"),
    ("protein_CPTAC", "CPTAC breast", "CPTAC-Result.pl?genenam={g}&ctype=Breast",
     RES + "portals_ualcan_cptac_protein.csv"),
]
STATOUT = RES + "portals_ualcan_stats.csv"

def main():
    s = U.session()
    now = lambda: time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    sf = open(STATOUT, 'w', newline=''); sw = csv.DictWriter(sf, fieldnames=STATF); sw.writeheader()
    for assay, dataset, tmpl, outpath in JOBS:
        f = open(outpath, 'w', newline=''); w = csv.DictWriter(f, fieldnames=BOXF); w.writeheader()
        for g in GENES:
            url = U.BASE + tmpl.format(g=g)
            h = U.fetch(s, url)
            if h is None:
                print("FAIL", assay, g, flush=True)
                w.writerow(dict(portal="UALCAN", assay=assay, dataset=dataset, gene=g,
                                grouping="RETRIEVAL_FAILED", category="", n="", unit="",
                                low_whisker="", q1="", median="", q3="", high_whisker="",
                                url=url, retrieved_utc=now())); f.flush()
                continue
            if re.search(r'not (?:found|available)|no data|Invalid gene', h, re.I) and 'highcharts({' not in h:
                print("NODATA", assay, g, flush=True)
            charts = U.parse_charts(h)
            nb = 0
            for c in charts:
                if not c['boxes'] or not c['categories']:
                    continue
                gr = grouping(c['categories'])
                for cat, b in zip(c['categories'], c['boxes']):
                    nb += 1
                    w.writerow(dict(portal="UALCAN", assay=assay, dataset=dataset, gene=g,
                                    grouping=gr, category=U.label_from_cat(cat),
                                    n=U.n_from_cat(cat), unit=c['unit'],
                                    low_whisker=b['low'], q1=b['q1'], median=b['median'],
                                    q3=b['q3'], high_whisker=b['high'], url=url,
                                    retrieved_utc=now()))
                for comp, p in c['stats'].items():
                    sw.writerow(dict(portal="UALCAN", assay=assay, dataset=dataset, gene=g,
                                     grouping=gr, comparison=comp, p_value=p, url=url,
                                     retrieved_utc=now()))
            f.flush(); sf.flush()
            print(assay, g, "charts", len(charts), "boxes", nb, flush=True)
            time.sleep(1.0)
        f.close()
    sf.close()

main()
