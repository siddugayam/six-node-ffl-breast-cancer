import sys, os, csv, time, itertools
sys.path.insert(0, '/path/to/revision/scripts/v4')
import kmplot_lib as K

OUT = '/path/to/revision/results/v4/portals_kmplot_mirna_breast.csv'
MIRS = ["hsa-miR-130a", "hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "hsa-let-7b", "hsa-let-7e",
        "hsa-miR-130b", "hsa-miR-301a"]
DATASETS = ["METABRIC", "TCGA", "GSE40267", "GSE19783"]
CUTOFFS = [("auto_percentile", {"cutoff": "auto", "auto_cutoff": "percentile_unique"}),
           ("median",          {"cutoff": "cutoff_quartile", "quartile": "50"})]
# subtype restrictions, applied to the primary miRNAs
RESTR = [("all", {}),
         ("ER_pos_IHC",        {"er_status_ihc": "1"}),
         ("ER_neg_IHC",        {"er_status_ihc": "0"}),
         ("ER_pos_array",      {"er_status_array": "1"}),
         ("ER_neg_array",      {"er_status_array": "0"}),
         ("HER2_pos_IHC",      {"her2_status_ihc_fish": "1"}),
         ("HER2_neg_IHC",      {"her2_status_ihc_fish": "0"}),
         ("TNBC_stGallen",     {"molecular_subtype_stgallen": "1"}),
         ("LuminalA_stGallen", {"molecular_subtype_stgallen": "2"}),
         ("LuminalB_stGallen", {"molecular_subtype_stgallen": "3"}),
         ("HER2pos_ERneg_stGallen", {"molecular_subtype_stgallen": "4"}),
         ("LN_positive",       {"lymph_node_status": "1"}),
         ("LN_negative",       {"lymph_node_status": "0"}),
         ]

FIELDS = ["query_id", "portal", "module", "mirna", "dataset", "endpoint", "cutoff_method",
          "restriction", "HR", "CI_low", "CI_high", "logrank_p", "FDR", "n_total", "n_low",
          "n_high", "n_events", "cutoff_value", "expr_range", "median_surv_low",
          "median_surv_high", "error", "permalink", "retrieved_utc"]

def main():
    jobs = []
    qid = 0
    for mir in MIRS:
        primary = mir in ("hsa-miR-130a", "hsa-miR-29a")
        for ds in DATASETS:
            for cname, cfg in CUTOFFS:
                for rname, rcfg in RESTR:
                    if rname != "all" and not (primary and ds in ("METABRIC", "TCGA")):
                        continue
                    if rname != "all" and cname != "auto_percentile" and mir != "hsa-miR-130a":
                        continue
                    qid += 1
                    jobs.append((qid, mir, ds, cname, cfg, rname, rcfg))
    print("total jobs:", len(jobs), flush=True)
    s = K.new_session('breast_mirna')
    done = set()
    if os.path.exists(OUT):
        with open(OUT) as fh:
            for r in csv.DictReader(fh):
                done.add((r['mirna'], r['dataset'], r['cutoff_method'], r['restriction']))
    newfile = not os.path.exists(OUT)
    fh = open(OUT, 'a', newline='')
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if newfile:
        w.writeheader()
    for (qid, mir, ds, cname, cfg, rname, rcfg) in jobs:
        if (mir, ds, cname, rname) in done:
            continue
        f = dict(K.MIRNA_DEFAULTS)
        f['affyid'] = mir
        f['dataset'] = ds
        f.update(cfg)
        f.update(rcfg)
        out, hml = K.run(s, 'breast_mirna', f, tries=3)
        hr, lo, hi = K.parse_hr(out)
        n = ln = hn = ev = None
        if out.get('result_id'):
            n, ln, hn, ev = K.get_n(s, out['result_id'])
        row = dict(query_id=qid, portal="KM Plotter (kmplot.com)", module="breast cancer miRNA",
                   mirna=mir, dataset=ds, endpoint="OS", cutoff_method=cname, restriction=rname,
                   HR=hr, CI_low=lo, CI_high=hi, logrank_p=out.get('P value', ''),
                   FDR=out.get('FDR', ''), n_total=n, n_low=ln, n_high=hn, n_events=ev,
                   cutoff_value=out.get('Cutoff value used in analysis', ''),
                   expr_range=out.get('Expression range of the probe', ''),
                   median_surv_low=out.get('median_surv_low', ''),
                   median_surv_high=out.get('median_surv_high', ''),
                   error=out.get('error', ''), permalink=out.get('permalink', ''),
                   retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
        w.writerow(row); fh.flush()
        print(qid, mir, ds, cname, rname, '->', hr, lo, hi, out.get('P value', ''), 'n=', n,
              out.get('error', ''), flush=True)
        time.sleep(1.2)
    fh.close()

main()
