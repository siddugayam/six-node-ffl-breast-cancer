import sys, csv, time
sys.path.insert(0, '/path/to/revision/scripts/v4')
import bcgem_lib as G

OUT = '/path/to/revision/results/v4/portals_bcgenexminer_prognosis.csv'
GENES = ["PTEN", "SMAD4", "TGFBR2", "DICER1", "KLF4", "XIAP", "COL1A1", "COL3A1",
         "NFKB1", "RELA", "SP1", "ETS1", "VEGFA", "CCND2", "MYC", "EZH2"]
F = ["portal", "module", "gene", "endpoint", "context", "population_ER", "population_PR",
     "population_HER2", "p_value", "HR", "CI_95", "good_prognosis_direction", "n_patients",
     "n_events", "split", "retrieved_utc"]

def main():
    fh = open(OUT, 'w', newline=''); w = csv.DictWriter(fh, fieldnames=F); w.writeheader()
    for g in GENES:
        s = G.session(2)
        d = {"Gene2": "G2GS", "Gene2Text": g, "Gene2Text2": g, "Event": "EventAE",
             "Node": "NodeA", "ER": "ERA", "HER2": "HER2A", "PR": "PRA", "Prolif": "ProlifA",
             "CorrelMap": "CorrelMap_None", "ExpressMap": "ExpressMap_Rec", "MP": "MPO",
             "MSM": "MSM_All", "Split": "SplitMed", "SplitCustomValeur": "50",
             "Verifier": "Submit"}
        v, r = G.two_step(s, d)
        rows = G.parse_cox_tables(r) if r else []
        for x in rows:
            ep = x['table_id'].replace('TabResExh_', '')
            w.writerow(dict(portal="bc-GenExMiner v5.3", module="exhaustive prognostic (mode=2)",
                            gene=x['gene'] or g, endpoint=ep, context=x['context'],
                            population_ER=x['pop1'], population_PR=x['pop2'],
                            population_HER2=x['pop3'], p_value=x['p_value'], HR=x['HR'],
                            CI_95=x['CI'], good_prognosis_direction=x['good_prognosis'],
                            n_patients=x['n_patients'], n_events=x['n_events'],
                            split="median", retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())))
        if not rows:
            w.writerow(dict(portal="bc-GenExMiner v5.3", module="exhaustive prognostic (mode=2)",
                            gene=g, endpoint="", context="RETRIEVAL_FAILED_OR_NO_DATA",
                            population_ER="", population_PR="", population_HER2="", p_value="",
                            HR="", CI_95="", good_prognosis_direction="", n_patients="",
                            n_events="", split="median",
                            retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())))
        fh.flush()
        top = [x for x in rows if x['pop1'] == 'ER all' and x['pop2'] == 'PR all' and x['pop3'] == 'HER2 all']
        print(g, "rows", len(rows), "| overall:",
              [(x['table_id'], x['p_value'], x['HR'], x['CI'], x['n_patients']) for x in top], flush=True)
        time.sleep(3)
    fh.close()
main()
