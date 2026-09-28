"""Re-run every kmplot row that came back without an HR, so the CSV records WHY."""
import sys, csv, time
sys.path.insert(0, '/path/to/revision/scripts/v4')
import kmplot_lib as K

RES = '/path/to/revision/results/v4/'
CUT = {"auto_percentile": {"cutoff": "auto", "auto_cutoff": "percentile_unique"},
       "median": {"cutoff": "cutoff_quartile", "quartile": "50"}}
MIR_RESTR = {"all": {}, "ER_pos_IHC": {"er_status_ihc": "1"}, "ER_neg_IHC": {"er_status_ihc": "0"},
             "ER_pos_array": {"er_status_array": "1"}, "ER_neg_array": {"er_status_array": "0"},
             "HER2_pos_IHC": {"her2_status_ihc_fish": "1"}, "HER2_neg_IHC": {"her2_status_ihc_fish": "0"},
             "TNBC_stGallen": {"molecular_subtype_stgallen": "1"},
             "LuminalA_stGallen": {"molecular_subtype_stgallen": "2"},
             "LuminalB_stGallen": {"molecular_subtype_stgallen": "3"},
             "HER2pos_ERneg_stGallen": {"molecular_subtype_stgallen": "4"},
             "LN_positive": {"lymph_node_status": "1"}, "LN_negative": {"lymph_node_status": "0"}}
GENE_RESTR = {"all": {}, "ER_pos": {"er_status_array": "1"}, "ER_neg": {"er_status_array": "0"},
              "PAM50_basal": {"molecular_subtype_pam50": "1"},
              "PAM50_luminalA": {"molecular_subtype_pam50": "2"},
              "PAM50_luminalB": {"molecular_subtype_pam50": "3"},
              "PAM50_HER2": {"molecular_subtype_pam50": "4"}}
SURV = {"RFS": "0", "OS": "1", "DMFS": "2"}

def redo(path, kind):
    rows = list(csv.DictReader(open(path)))
    fn = rows[0].keys()
    s = K.new_session('breast_mirna' if kind == 'mirna' else 'breast')
    nfix = 0
    for r in rows:
        if r['HR']:
            continue
        if kind == 'mirna':
            f = dict(K.MIRNA_DEFAULTS); f['affyid'] = r['mirna']; f['dataset'] = r['dataset']
            f.update(CUT[r['cutoff_method']]); f.update(MIR_RESTR[r['restriction']])
            cancer = 'breast_mirna'
        else:
            f = dict(K.MRNA_DEFAULTS); f['affyid'] = r['gene']; f['surv'] = SURV[r['endpoint']]
            f.update(CUT[r['cutoff_method']]); f.update(GENE_RESTR[r['restriction']])
            cancer = 'breast'
        out, h = K.run(s, cancer, f, tries=2)
        hr, lo, hi = K.parse_hr(out)
        r['HR'], r['CI_low'], r['CI_high'] = hr, lo, hi
        r['logrank_p'] = out.get('P value', '')
        r['error'] = out.get('error', '')
        if out.get('result_id'):
            n, ln, hn, ev = K.get_n(s, out['result_id'])
            r['n_total'], r['n_low'], r['n_high'], r['n_events'] = n, ln, hn, ev
        nfix += 1
        print("redo", [r.get(k) for k in ('mirna', 'gene', 'dataset', 'endpoint', 'cutoff_method',
                                          'restriction')], '->', hr, r['error'], flush=True)
        time.sleep(1.2)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(fn)); w.writeheader(); w.writerows(rows)
    print("rewrote", path, "fixed", nfix, flush=True)

redo(RES + 'portals_kmplot_mirna_breast.csv', 'mirna')
redo(RES + 'portals_kmplot_mrna_breast.csv', 'gene')
