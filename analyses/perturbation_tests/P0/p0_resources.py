#!/usr/bin/env python3
"""P0 of analyses/perturbation_tests: free disk, RAM, CPU, GPU and the download cap.  Writes p0_resources.txt (numbered lines).
Cap (fixed in the analysis plan): new downloads may use at most the smaller of 150 GB or 50 % of free space."""
import os, shutil, subprocess, datetime
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REV = os.path.dirname(os.path.dirname(INB))
CACHE = os.path.join(REV, 'data_cache')
sh = lambda c: subprocess.run(c, shell=True, capture_output=True, text=True).stdout.strip()
L = []
P = L.append
P(f'date: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}')
for name, p in (('REV', REV), ('cache', CACHE)):
    du = shutil.disk_usage(p); dev = sh(f'df --output=source "{p}" | tail -1')
    P(f'disk {name}: device {dev}; size {du.total / 1e9:.0f} GB; used {du.used / 1e9:.0f} GB; free {du.free / 1e9:.1f} GB')
free = shutil.disk_usage(CACHE).free / 1e9
P(f'download cap: min(150 GB, 50 % of {free:.1f} GB free) = {min(150, free / 2):.1f} GB')
mem = dict(l.split(':', 1) for l in open('/proc/meminfo').read().splitlines())
kb = lambda k: int(mem[k].split()[0])
P(f'RAM: total {kb("MemTotal") / 2**20:.1f} GiB; available {kb("MemAvailable") / 2**20:.1f} GiB')
P(f'CPU threads: {os.cpu_count()}; load average {os.getloadavg()[0]:.2f}')
P('GPU: ' + sh('nvidia-smi --query-gpu=name,memory.total,memory.used,utilization.gpu,power.draw,power.limit,temperature.gpu,driver_version --format=csv,noheader'))
P('CUDA (nvidia-smi): ' + sh("nvidia-smi | grep -o 'CUDA Version: [0-9.]*'"))
apps = sh('nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader')
apps = '\n'.join(', '.join([p[0], os.path.basename(p[1]), *p[2:]]) for p in (l.split(', ') for l in apps.splitlines())) if apps else ''   # program name only, no local path
P('GPU compute processes: ' + (apps or 'none'))
for pid in sh('nvidia-smi --query-compute-apps=pid --format=csv,noheader').split():
    P(f'  pid {pid}: elapsed {sh(f"ps -o etime= -p {pid}")}; command {os.path.basename(sh(f"ps -o comm= -p {pid}"))}')
P('torch: ' + sh('python3 -c "import torch; print(torch.__version__, torch.cuda.is_available())"'))
hf = os.path.expanduser('~/.cache/huggingface/hub')
for r in range(4):
    d = os.path.join(hf, f'models--johahi--borzoi-replicate-{r}')
    P(f'Borzoi replicate {r} weights cached: ' + (sh(f'du -sh "{d}" | cut -f1') if os.path.isdir(d) else 'no'))
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'p0_resources.txt'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
