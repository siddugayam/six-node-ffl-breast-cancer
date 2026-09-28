import sys, os, csv, time
sys.path.insert(0, '/path/to/revision/scripts/10_survival_and_clinical/web_portals')
import kmplot_lib as K

OUT = '/path/to/revision/results/v4/portals_kmplot_mrna_breast.csv'
GENES = ["PTEN", "SMAD4", "TGFBR2", "DICER1", "KLF4", "XIAP", "COL1A1", "COL3A1",
         "NFKB1", "RELA", "SP1", "ETS1", "VEGFA", "CCND2", "MYC", "EZH2"]
ENDPOINTS = [("RFS", "0"), ("OS", "1"), ("DMFS", "2")]
CUTOFFS = [("median", {"cutoff": "cutoff_quartile", "quartile": "50"}),
           ("auto_percentile", {"cutoff": "auto", "auto_cutoff": "percentile_unique"})]
SUBSET = ["PTEN", "SMAD4", "TGFBR2", "DICER1", "KLF4", "MYC", "EZH2", "VEGFA"]
RESTR = [("all", {}),
         ("ER_pos", {"er_status_array": "1"}),
         ("ER_neg", {"er_status_array": "0"}),
         ("PAM50_basal", {"molecular_subtype_pam50": "1"}),
         ("PAM50_luminalA", {"molecular_subtype_pam50": "2"}),
         ("PAM50_luminalB", {"molecular_subtype_pam50": "3"}),
         ("PAM50_HER2", {"molecular_subtype_pam50": "4"})]

FIELDS = ["query_id", "portal", "module", "gene", "probe", "endpoint", "cutoff_method",
          "restriction", "HR", "CI_low", "CI_high", "logrank_p", "FDR", "n_total", "n_low",
          "n_high", "n_events", "cutoff_value", "expr_range", "median_surv_low",
          "median_surv_high", "error", "permalink", "retrieved_utc"]

def main():
    jobs = []
    qid = 0
    for g in GENES:
        for ename, sv in ENDPOINTS:
            for cname, cfg in CUTOFFS:
                qid += 1
                jobs.append((qid, g, ename, sv, cname, cfg, "all", {}))
    for g in SUBSET:
        for rname, rcfg in RESTR[1:]:
            qid += 1
            jobs.append((qid, g, "RFS", "0", "median",
                         {"cutoff": "cutoff_quartile", "quartile": "50"}, rname, rcfg))
    print("total jobs:", len(jobs), flush=True)
    s = K.new_session('breast')
    done = set()
    if os.path.exists(OUT):
        with open(OUT) as fh:
            for r in csv.DictReader(fh):
                done.add((r['gene'], r['endpoint'], r['cutoff_method'], r['restriction']))
    newfile = not os.path.exists(OUT)
    fh = open(OUT, 'a', newline='')
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if newfile:
        w.writeheader()
    for (qid, g, ename, sv, cname, cfg, rname, rcfg) in jobs:
        if (g, ename, cname, rname) in done:
            continue
        f = dict(K.MRNA_DEFAULTS)
        f['affyid'] = g
        f['surv'] = sv
        f.update(cfg); f.update(rcfg)
        out, hml = K.run(s, 'breast', f, tries=3)
        hr, lo, hi = K.parse_hr(out)
        n = ln = hn = ev = None
        if out.get('result_id'):
            n, ln, hn, ev = K.get_n(s, out['result_id'])
        row = dict(query_id=qid, portal="KM Plotter (kmplot.com)",
                   module="breast cancer mRNA (gene chip)", gene=g,
                   probe=K.parse_probe(hml), endpoint=ename, cutoff_method=cname,
                   restriction=rname, HR=hr, CI_low=lo, CI_high=hi,
                   logrank_p=out.get('P value', ''), FDR=out.get('FDR', ''),
                   n_total=n, n_low=ln, n_high=hn, n_events=ev,
                   cutoff_value=out.get('Cutoff value used in analysis', ''),
                   expr_range=out.get('Expression range of the probe', ''),
                   median_surv_low=out.get('median_surv_low', ''),
                   median_surv_high=out.get('median_surv_high', ''),
                   error=out.get('error', ''), permalink=out.get('permalink', ''),
                   retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
        w.writerow(row); fh.flush()
        print(qid, g, ename, cname, rname, '->', hr, lo, hi, out.get('P value', ''), 'n=', n,
              out.get('error', ''), flush=True)
        time.sleep(1.5)
    fh.close()

main()
