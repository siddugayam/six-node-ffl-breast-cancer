"""Headline rows the manuscript can quote directly, each traced to its source file."""
import csv
RES = '/path/to/revision/results/v4/'
km   = list(csv.DictReader(open(RES + 'portals_kmplot_mrna_breast.csv')))
kmi  = list(csv.DictReader(open(RES + 'portals_kmplot_mirna_breast.csv')))
bc   = list(csv.DictReader(open(RES + 'portals_bcgenexminer_prognosis.csv')))
F = ["claim_group", "feature", "portal", "cohort", "endpoint", "split", "HR", "CI_95",
     "p", "n_patients", "n_events", "reading", "source_file"]
rows = []
ANCH = ["PTEN", "SMAD4", "TGFBR2", "DICER1", "KLF4"]
for g in ANCH:
    for ep in ["DMFS", "OS", "DFS"]:
        r = [x for x in bc if x['gene'] == g and x['endpoint'] == ep
             and x['population_ER'] == 'ER all' and x['population_PR'] == 'PR all'
             and x['population_HER2'] == 'HER2 all']
        if r:
            r = r[0]
            rows.append(dict(claim_group="miR-130a tumour-suppressor targets: higher mRNA = better outcome",
                feature=g, portal="bc-GenExMiner v5.3", cohort="pooled breast cancer cohorts",
                endpoint=ep, split="median", HR=r['HR'], CI_95=r['CI_95'], p=r['p_value'],
                n_patients=r['n_patients'], n_events=r['n_events'],
                reading=r['good_prognosis_direction'],
                source_file="portals_bcgenexminer_prognosis.csv"))
    for ep in ["RFS", "OS", "DMFS"]:
        r = [x for x in km if x['gene'] == g and x['endpoint'] == ep
             and x['restriction'] == 'all' and x['cutoff_method'] == 'median']
        if r:
            r = r[0]
            rows.append(dict(claim_group="miR-130a tumour-suppressor targets: higher mRNA = better outcome",
                feature="%s (%s)" % (g, r['probe']), portal="KM Plotter mRNA",
                cohort="pooled public microarray cohorts", endpoint=ep, split="median",
                HR=r['HR'], CI_95="%s - %s" % (r['CI_low'], r['CI_high']), p=r['logrank_p'],
                n_patients=r['n_total'], n_events=r['n_events'],
                reading=("high expression = good prognosis" if r['HR'] and float(r['HR']) < 1
                         else "high expression = poor prognosis"),
                source_file="portals_kmplot_mrna_breast.csv"))
for mir in ["hsa-miR-130a", "hsa-miR-130b", "hsa-miR-301a", "hsa-miR-29a"]:
    for ds in ["METABRIC", "TCGA", "GSE40267", "GSE19783"]:
        for cut in ["median", "auto_percentile"]:
            r = [x for x in kmi if x['mirna'] == mir and x['dataset'] == ds
                 and x['cutoff_method'] == cut and x['restriction'] == 'all']
            if r and r[0]['HR']:
                r = r[0]
                rows.append(dict(
                    claim_group="miR-130a survival in an independent miRNA cohort "
                                "(and its seed-family comparators)",
                    feature=mir, portal="KM Plotter breast miRNA module", cohort=ds,
                    endpoint="OS", split=cut, HR=r['HR'],
                    CI_95="%s - %s" % (r['CI_low'], r['CI_high']), p=r['logrank_p'],
                    n_patients=r['n_total'], n_events=r['n_events'],
                    reading=("higher miRNA = better OS" if float(r['HR']) < 1
                             else "higher miRNA = worse OS")
                            + ("; KM Plotter FDR " + r['FDR'] if r['FDR'] else ""),
                    source_file="portals_kmplot_mirna_breast.csv"))
with open(RES + 'portals_HEADLINE_RESULTS.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=F); w.writeheader(); w.writerows(rows)
print("wrote portals_HEADLINE_RESULTS.csv", len(rows), "rows")
