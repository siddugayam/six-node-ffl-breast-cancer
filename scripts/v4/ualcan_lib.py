import re, html, time
import requests

BASE = "https://ualcan.path.uab.edu/cgi-bin/"

def session():
    s = requests.Session()
    s.headers.update({"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                                    "(KHTML, like Gecko) Chrome/151.0.0.0 Safari/537.36"})
    return s

def fetch(s, url, tries=3):
    for i in range(tries):
        try:
            r = s.get(url, timeout=300)
            if r.status_code == 200 and len(r.text) > 3000:
                return r.text
        except Exception:
            pass
        time.sleep(3 * (i + 1))
    return None

BOX_RE = re.compile(r'low\s*:\s*(-?[\d.]+)\s*,\s*q1\s*:\s*(-?[\d.]+)\s*,\s*median\s*:\s*(-?[\d.]+)'
                    r'\s*,\s*q3\s*:\s*(-?[\d.]+)\s*,\s*high\s*:\s*(-?[\d.]+)', re.S)

def parse_charts(h):
    """Return a list of dicts, one per Highcharts boxplot block, in page order."""
    out = []
    # each block starts at "$('#containerN...').highcharts({"
    starts = [m.start() for m in re.finditer(r"\$\('#container[^']*'\)\.highcharts\(\{", h)]
    starts.append(len(h))
    for i in range(len(starts) - 1):
        blk = h[starts[i]:starts[i + 1]]
        tm = re.search(r"title:\s*\{\s*text:\s*'([^']*)'", blk)
        title = tm.group(1) if tm else ''
        cm = re.search(r"categories:\s*\[(.*?)\]", blk, re.S)
        cats = []
        if cm:
            cats = [html.unescape(x.strip().strip("'").replace('<br>', ' ').strip())
                    for x in re.findall(r"'([^']*)'", cm.group(1))]
        ym = re.findall(r"title:\s*\{\s*text:\s*'([^']*)'", blk)
        unit = ym[2] if len(ym) > 2 else (ym[-1] if ym else '')
        boxes = BOX_RE.findall(blk)
        # stats table that follows this chart, before the next chart
        stats = {}
        for m in re.finditer(r'<td[^>]*>\s*([A-Za-z0-9_+\-]+-vs-[A-Za-z0-9_+\-]+)\s*</td>\s*'
                             r'<td[^>]*>\s*([0-9.eE+\-]*)\s*</td>', blk):
            stats[m.group(1)] = m.group(2)
        out.append(dict(title=title, categories=cats, unit=unit,
                        boxes=[dict(low=float(b[0]), q1=float(b[1]), median=float(b[2]),
                                    q3=float(b[3]), high=float(b[4])) for b in boxes],
                        stats=stats))
    return out

def n_from_cat(c):
    m = re.search(r'\(n\s*=\s*(\d+)\)', c)
    return int(m.group(1)) if m else None

def label_from_cat(c):
    return re.sub(r'\s*\(n\s*=\s*\d+\)\s*', '', c).strip()

def all_stats(h):
    """All Comparison/significance pairs anywhere on the page (fallback)."""
    d = {}
    for m in re.finditer(r'<td[^>]*>\s*([A-Za-z0-9_+\- ]+-vs-[A-Za-z0-9_+\- ]+)\s*</td>\s*'
                         r'<td[^>]*>\s*([0-9.eE+\-]*)\s*</td>', h):
        k = m.group(1).strip()
        if k not in d:
            d[k] = m.group(2).strip()
    return d
