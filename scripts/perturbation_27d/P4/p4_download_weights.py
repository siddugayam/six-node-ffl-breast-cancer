#!/usr/bin/env python3
"""P4: fetch the Borzoi replicate weights 1-3 (borzoi-pytorch port of the official Calico human weights) into the
27d cache, and record each file in DOWNLOADS.tsv.  Replicate 0 was fetched on 2026-09-09 for scripts/v3/23 and is
already in the local Hugging Face cache; its revision and md5 are recorded in p4_weights.txt, not re-downloaded."""
import os, sys, glob, hashlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _dl
from huggingface_hub import snapshot_download, HfApi
HF = os.path.join(_dl.CACHE, 'huggingface')
api = HfApi(); out = []
for r in range(4):
    repo = f'johahi/borzoi-replicate-{r}'; info = api.model_info(repo); sha = info.sha
    lic = (info.card_data or {}).get('license') if info.card_data else None
    if r == 0:
        snap = glob.glob(os.path.expanduser(f'~/.cache/huggingface/hub/models--johahi--borzoi-replicate-0/snapshots/*'))
        for s in snap:
            for f in sorted(os.listdir(s)):
                p = os.path.realpath(os.path.join(s, f))
                out.append(f'replicate 0 (cached 2026-09-09): snapshot {os.path.basename(s)}; {f}; {os.path.getsize(p)} bytes; md5 {_dl.md5(p)}; hub revision now {sha}')
        continue
    d = snapshot_download(repo, revision=sha, cache_dir=HF, allow_patterns=['config.json', 'model.safetensors'])
    for f in ('config.json', 'model.safetensors'):
        p = os.path.join(d, f)
        _dl.register_existing(os.path.realpath(p), 'P4', 'Hugging Face (borzoi-pytorch port of Calico Borzoi weights)', f'{repo}@{sha}',
                              f'https://huggingface.co/{repo}/resolve/{sha}/{f}', f'{lic} (model card)')
        out.append(f'replicate {r}: revision {sha}; {f}; {os.path.getsize(p)} bytes; md5 {_dl.md5(os.path.realpath(p))}; licence {lic}')
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'p4_weights.txt'), 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
