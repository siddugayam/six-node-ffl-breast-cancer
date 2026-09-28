"""Direct replay of the exact multipart POST that kmplot.com's own form issues.
Field set captured from a real headless-Chromium submission (see portal_km_probe.py)."""
import re, html, time, json, sys
import requests

BASE = "https://kmplot.com/analysis/index.php?p=service"

MIRNA_DEFAULTS = {
    "step": "2",
    "show_plot_line_thickness": "2", "show_plot_text_size": "1.5", "show_plot_text_bold": "1",
    "show_plot_custom_main_title": "", "show_plot_pvalue_position": "topright",
    "show_plot_legend_position": "bottomleft", "show_number_at_risk": "true",
    "show_minmax": "true", "show_median": "true",
    "affyid": "", "affyid_hidden": "", "dataset": "METABRIC",
    "more_genes_param": "gene_list", "multipletesting": "false",
    "surv": "0", "quartile": "50", "cutoff": "auto", "auto_cutoff": "percentile_unique",
    "trichotomization": "false", "follow_up_threshold": "false",
    "censore_at_threshold": "true", "show_median_of_survival": "true",
    "er_status_ihc": "false", "er_status_array": "false", "pr_status_ihc": "false",
    "her2_status_ihc_fish": "false", "her2_status_array": "false",
    "molecular_subtype_stgallen": "false", "grade": "false", "lymph_node_status": "false",
    "no_treatment_most_likely_select": "1", "er_endocrine_most_likely_select": "false",
    "chemotherapy_most_likely_select": "false", "termsOfUseAccepted": "on",
}

def _text(h):
    h2 = re.sub(r'<script.*?</script>', '', h, flags=re.S | re.I)
    h2 = re.sub(r'<style.*?</style>', '', h2, flags=re.S | re.I)
    return html.unescape(re.sub('<[^>]+>', '\n', h2))

def parse_result(h):
    """Pull the labelled key/value rows out of a kmplot result page."""
    out = {}
    # the analysis summary and Results tables are <td>label</td><td>value</td>
    for m in re.finditer(r'<td[^>]*>\s*([^<]{2,60}?)\s*:?\s*</td>\s*<td[^>]*>(.*?)</td>', h, re.S):
        k = html.unescape(re.sub(r'\s+', ' ', m.group(1))).strip().rstrip(':')
        v = html.unescape(re.sub(r'\s+', ' ', re.sub('<[^>]+>', '', m.group(2)))).strip()
        if k and k not in out:
            out[k] = v
    # median survival table (two columns)
    mm = re.search(r'id="median_0"[^>]*>([^<]*)</td><td id="median_1"[^>]*>([^<]*)</td>', h)
    if mm:
        out['median_surv_low'] = mm.group(1).strip()
        out['median_surv_high'] = mm.group(2).strip()
    pl = re.search(r"index\.php\?p=view&(?:amp;)?pa_id=(\d+)&(?:amp;)?show=([A-Za-z0-9_\-]+)", h)
    if pl:
        out['permalink'] = "https://kmplot.com/analysis/index.php?p=view&pa_id=%s&show=%s" % (pl.group(1), pl.group(2))
    rid = re.search(r'result=(km_[0-9a-z_\-]+)&(?:amp;)?version=academic', h)
    if rid:
        out['result_id'] = rid.group(1)
    if 'unexpected error' in h:
        out['error'] = 'unexpected error'
    if 'The input is not valid' in h:
        out['error'] = 'input not valid'
    if 'not enough' in h.lower():
        out['error'] = 'not enough patients'
    wm = re.search(r"<i class='material-icons p-2'>warning</i>([^<]{5,200})", h)
    if wm:
        out['error'] = re.sub(r'\s+', ' ', html.unescape(wm.group(1))).strip()
    return out

def parse_hr(out):
    hr = out.get('Hazard ratio', '')
    m = re.match(r'([0-9.]+)\s*\(95% CI:\s*([0-9.]+)\s*-\s*([0-9.]+)\)', hr)
    if m:
        return m.group(1), m.group(2), m.group(3)
    m = re.match(r'([0-9.]+)', hr)
    return (m.group(1) if m else ''), '', ''

def get_n(session, result_id, retries=2):
    """Fetch the exported plot-data txt and count patients per arm."""
    if not result_id:
        return None, None, None, None
    url = ("https://kmplot.com/kmplot_commons/php/src/service/results/index.php"
           "?result=%s&version=academic&type=txt" % result_id)
    for _ in range(retries):
        try:
            r = session.get(url, timeout=180)
            if r.status_code == 200 and len(r.text) > 50:
                return parse_txt(r.text)
        except Exception:
            time.sleep(2)
    return None, None, None, None

def parse_txt(txt):
    """kmplot export: Sample \t Expression (1=high) \t Expression \t Time (months) \t Event"""
    lines = [l for l in txt.strip().split('\n') if l.strip()]
    if len(lines) < 2:
        return None, None, None, None
    lo = hi = ev = 0
    for l in lines[1:]:
        f = [x.strip().strip('"') for x in l.split('\t')]
        if len(f) < 5:
            continue
        if f[1] == '1':
            hi += 1
        elif f[1] == '0':
            lo += 1
        if f[4] == '1':
            ev += 1
    return lo + hi, lo, hi, ev

def run(session, cancer, fields, tries=3, pause=2.0):
    url = BASE + "&cancer=" + cancer
    data = dict(fields)
    files = [(k, (None, str(v))) for k, v in data.items()]
    # export_data_as_txt is submitted twice by the real form
    files.append(("export_data_as_txt", (None, "true")))
    files.append(("export_data_as_txt", (None, "true")))
    last = None
    for i in range(tries):
        try:
            r = session.post(url, files=files, timeout=300,
                             headers={"Referer": url})
            if r.status_code == 200:
                out = parse_result(r.text)
                if 'error' not in out and out.get('Hazard ratio'):
                    return out, r.text
                last = (out, r.text)
        except Exception as e:
            last = ({'error': repr(e)}, '')
        time.sleep(pause * (i + 1))
    return last if last else ({'error': 'no response'}, '')

def new_session(cancer):
    s = requests.Session()
    s.headers.update({"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                                    "(KHTML, like Gecko) Chrome/151.0.0.0 Safari/537.36"})
    s.get(BASE + "&cancer=" + cancer, timeout=120)
    return s

# ---- breast mRNA module (cancer=breast). Field set captured from a real browser submission. ----
MRNA_DEFAULTS = {
    "step": "2",
    "show_plot_line_thickness": "2", "show_plot_text_size": "1.5", "show_plot_text_bold": "1",
    "show_plot_custom_main_title": "", "show_plot_pvalue_position": "topright",
    "show_plot_legend_position": "bottomleft", "show_number_at_risk": "true",
    "show_minmax": "true", "show_median": "true",
    "redundant_publication": "true", "array_quality": "2", "propcheck": "true",
    "affyid": "", "affyid_hidden": "", "more_genes_param": "gene_list",
    "cutoff": "cutoff_quartile", "quartile": "50", "auto_cutoff": "percentile_unique",
    "trichotomization": "false", "surv": "0", "show_median_of_survival": "true",
    "follow_up_threshold": "false", "censore_at_threshold": "true",
    "probe_set_option": "jetset",
    "er_status": "false", "er_status_array": "false", "pr_status_ihc": "false",
    "her2_status_array": "false", "molecular_subtype": "false",
    "molecular_subtype_pam50": "false", "lymph_node_status": "false", "grade": "false",
    "p53_status": "false", "pietenpol_subtype": "false", "dataset": "false",
    "no_treatment_most_likely_select": "1", "er_endocrine_most_likely_select": "false",
    "chemotherapy_most_likely_select": "false",
    "b_cells": "false", "cd4_memory_t_cells": "false", "cd8_t_cells": "false",
    "macrophages": "false", "natural_killer_t_cells": "false", "regulatory_t_cells": "false",
    "termsOfUseAccepted": "on",
}
SURV_MRNA = {"0": "RFS", "1": "OS", "2": "DMFS", "3": "PPS"}

def parse_probe(h):
    m = re.search(r'The selected Affy ID is valid:\s*([^<,]+)', h)
    return m.group(1).strip() if m else ''
