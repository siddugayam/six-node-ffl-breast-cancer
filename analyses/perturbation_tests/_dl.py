#!/usr/bin/env python3
"""Download helper for analyses/perturbation_tests.  Every file fetched for these tests goes through fetch(), which stores it
under REV/data_cache_2026-09-27d/<subdir>/, and appends one row to DOWNLOADS.tsv:
part, source, accession, URL, date (UTC), size_bytes, md5, licence, file (relative to the cache).
A file already in the cache with a DOWNLOADS.tsv row is not fetched again.  The cap of P0 (42.9 GB) is checked
before every download against the running total of DOWNLOADS.tsv."""
import os, csv, hashlib, datetime, time, urllib.request, urllib.error, shutil

INB = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(os.path.dirname(INB))
CACHE = os.path.join(REV, 'data_cache_2026-09-27d')
LOG = os.path.join(INB, 'DOWNLOADS.tsv')
COLS = ['part', 'source', 'accession', 'url', 'date_utc', 'size_bytes', 'md5', 'licence', 'file']
CAP_BYTES = 42.9e9

def _rows():
    if not os.path.exists(LOG): return []
    return list(csv.DictReader(open(LOG), delimiter='\t'))

def total_bytes():
    return sum(int(r['size_bytes']) for r in _rows())

def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

def record(part, source, accession, url, path, licence):
    new = not os.path.exists(LOG)
    with open(LOG, 'a', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        if new: w.writerow(COLS)
        w.writerow([part, source, accession, url, datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
                    os.path.getsize(path), md5(path), licence, os.path.relpath(path, CACHE)])

def fetch(url, subdir, fname, part, source, accession, licence, tries=4, timeout=120, headers=None, max_bytes=3e9):
    """Return the local path; download only if no DOWNLOADS.tsv row exists for this file."""
    d = os.path.join(CACHE, subdir); os.makedirs(d, exist_ok=True)
    path = os.path.join(d, fname); rel = os.path.relpath(path, CACHE)
    if os.path.exists(path) and any(r['file'] == rel for r in _rows()): return path
    if total_bytes() > CAP_BYTES: raise RuntimeError('download cap reached')
    tmp = path + '.part'; err = None
    for k in range(tries):
        try:
            # NCBI answers browser-like user agents on /geo/download/ with a reCAPTCHA page; a plain client agent is served
            req = urllib.request.Request(url, headers=headers or {'User-Agent': 'curl/8.5.0'})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, 'wb') as f:
                n = 0
                while True:                                  # stop anything larger than max_bytes (default 3 GB)
                    b = r.read(1 << 20)
                    if not b: break
                    n += len(b)
                    if n > max_bytes: raise OSError(f'larger than {max_bytes / 1e9:.1f} GB; not fetched')
                    f.write(b)
            head = open(tmp, 'rb').read(16)
            if fname.endswith('.gz') and head[:2] != b'\x1f\x8b':
                raise OSError(f'not gzip data (starts {head[:15]!r})')
            os.replace(tmp, path); record(part, source, accession, url, path, licence); return path
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            err = e; time.sleep(3 * (k + 1))
    if os.path.exists(tmp): os.remove(tmp)
    raise RuntimeError(f'download failed: {url}: {err}')

def register_existing(path, part, source, accession, url, licence):
    """Record a file fetched by another tool (e.g. huggingface_hub) into the cache."""
    rel = os.path.relpath(path, CACHE)
    if not any(r['file'] == rel for r in _rows()): record(part, source, accession, url, path, licence)
