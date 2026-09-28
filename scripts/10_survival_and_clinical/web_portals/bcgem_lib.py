import re, html, time
import requests, urllib3
urllib3.disable_warnings()

B = "https://bcgenex.ico.unicancer.fr/BC-GEM/"
# NOTE: this host serves an INCOMPLETE TLS chain (leaf cert is a valid Sectigo cert for
# *.ico.unicancer.fr, but the intermediate is not sent -> openssl "Verify return code: 21").
# Browsers that carry the Sectigo intermediate succeed; curl/requests do not. We therefore
# disable verification for this host only, having first confirmed the leaf certificate.
def session(mode):
    s = requests.Session(); s.verify = False
    s.headers.update({"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                                    "(KHTML, like Gecko) Chrome/151.0.0.0 Safari/537.36"})
    s.get(B + "GEM-Requete.php?mode=%d" % mode, timeout=180)
    return s

def two_step(s, payload, tries=2):
    last = ''
    for i in range(tries):
        try:
            r1 = s.post(B + "GEM-requete.php", data=payload, timeout=300)
            if 'Validation step' not in r1.text and 'Verif_Submit' not in r1.text:
                last = r1.text
                time.sleep(5); continue
            r2 = s.post(B + "GEM-requete.php", data={"Calculer": "Start analysis"}, timeout=900)
            return r1.text, r2.text
        except Exception as e:
            last = repr(e)
            time.sleep(8)
    return last, ''

TAG = re.compile(r'</?[a-zA-Z!][^>]*>')   # only real tags; a bare "< 0.0001" must survive

def cells(tr):
    return [html.unescape(re.sub(r'\s+', ' ', TAG.sub('', x))).strip()
            for x in re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', tr, re.S)]

def arrow(tr):
    m = re.search(r'fleche_(bleue|rouge|verte|orange)_(haut|bas)', tr)
    if m:
        return "high expression = good prognosis" if m.group(2) == 'haut' \
               else "low expression = good prognosis"
    return ""

def parse_cox_tables(h):
    """Every 'univariate Cox analysis' table on a bc-GenExMiner result page."""
    out = []
    for m in re.finditer(r'<table[^>]*id="(TabResExh_[A-Z]+|[^"]*)"[^>]*>(.*?)</table>', h, re.S):
        tid, tb = m.group(1), m.group(2)
        hm = re.search(r'univariate Cox analysis\s*</b>\s*<br\s*/?>\s*\(([^)]*)\)', tb, re.S)
        if not hm:
            continue
        ctx = re.sub(r'\s+', ' ', hm.group(1)).strip()
        gm = re.search(r'class="geneSymbol">([^<]*)<', tb)
        gene = gm.group(1).strip() if gm else ''
        for tr in re.findall(r'<tr[^>]*>(.*?)</tr>', tb, re.S):
            c = cells(tr)
            if len(c) < 11 or not c[0].isdigit():
                continue
            out.append(dict(table_id=tid, context=ctx, gene=gene, rank=c[0],
                            pop1=c[2], pop2=c[3], pop3=c[4], p_value=c[6], HR=c[7],
                            CI=c[8], good_prognosis=arrow(tr), n_patients=c[10].replace(' ',' '),
                            n_events=(c[11].replace(' ',' ') if len(c) > 11 else '')))
    return out

def parse_subtype_tables(h):
    """The mode=3 / mode=7 tables: rows are gene-expression signatures, not populations."""
    out = []
    for m in re.finditer(r'<table[^>]*>(.*?)</table>', h, re.S):
        tb = m.group(1)
        hm = re.search(r'univariate Cox analysis by\s*<[^>]*>\s*([^<]*)', tb, re.S)
        if 'univariate Cox analysis by' not in tb:
            continue
        head = re.sub(r'\s+', ' ', re.sub('<[^>]+>', ' ',
                      tb[:tb.find('</tr>')])).strip()
        gm = re.search(r'class="geneSymbol">([^<]*)<', tb)
        gene = gm.group(1).strip() if gm else ''
        for tr in re.findall(r'<tr[^>]*>(.*?)</tr>', tb, re.S):
            c = cells(tr)
            if len(c) < 8 or not c[0].isdigit():
                continue
            out.append(dict(context=head[:160], gene=gene, rank=c[0], signature=c[1],
                            p_value=c[3], HR=c[4], CI=c[5], good_prognosis=arrow(tr),
                            n_patients=c[7].replace(' ',' ') if len(c) > 7 else '',
                            n_events=c[8].replace(' ',' ') if len(c) > 8 else ''))
    return out
