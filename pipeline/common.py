"""Small, explicit IO, audit, and resource guards used by every stage."""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import threading
import time
from datetime import datetime, timezone
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]

def now():
    return datetime.now(timezone.utc).isoformat()

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def read_tsv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if not rows and not fields:
        raise ValueError(f'Refuse decorative empty table: {path}')
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='raise')
        w.writeheader()
        for row in rows:
            w.writerow({k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list)) else v
                        for k, v in row.items()})

def append_tsv(path, row):
    path = Path(path)
    rows = read_tsv(path) if path.exists() else []
    write_tsv(path, rows + [row])

def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n', encoding='utf-8')

def config():
    import yaml
    return yaml.safe_load((ROOT / 'config/design.yml').read_text())

def footprint(root=ROOT):
    total = 0
    for base, _, files in os.walk(root):
        for name in files:
            try:
                total += (Path(base) / name).stat().st_size
            except FileNotFoundError:
                pass
    return total

def budget(requested, free, used=0):
    safe = free - 1_000_000_000
    if safe <= 0:
        raise RuntimeError('Resource refusal: non-project disk reserve unavailable')
    return min(int(requested), 13_000_000_000, safe + used)

def check_resources(projected=0, pilot_size=0, reducible='candidate_count'):
    cfg = config()
    free = shutil.disk_usage(ROOT).free
    used = footprint()
    ceiling = cfg['disk_bytes']
    if (ROOT / 'config/resolved.json').exists():
        ceiling = json.loads((ROOT / 'config/resolved.json').read_text())['effective_disk_budget']
    remaining = min(ceiling - used, free - cfg['reserve_bytes'])
    adjusted = int(projected * 1.2)
    if remaining < adjusted or remaining <= 0:
        raise RuntimeError(f'Resource refusal: current_free={free}; configured_ceiling={ceiling}; '
                           f'project_used={used}; remaining_budget={remaining}; pilot_size={pilot_size}; '
                           f'projected_full={projected}; safety_adjusted={adjusted}; '
                           f'shortfall={max(0, adjusted-remaining)}; reduce={reducible}')
    return remaining

def failure(stage, command, error, resolution='unresolved', exit_status=1):
    append_tsv(ROOT / 'logs/failures.tsv', dict(timestamp=now(), stage=stage,
               command=command, error_type=type(error).__name__, exit_status=exit_status,
               message=str(error), resolution=resolution, rerun_id=''))

def fetch(url, destination, identifier='', version=''):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        check_resources(20_000_000, 1, 'download_count')
        request = Request(url, headers={'User-Agent': 'CD81ReproducibleBenchmark/0.1'})
        with urlopen(request, timeout=90) as response, destination.with_suffix(destination.suffix+'.tmp').open('wb') as f:
            shutil.copyfileobj(response, f)
        destination.with_suffix(destination.suffix+'.tmp').replace(destination)
    manifest = ROOT / 'results/download_manifest.tsv'
    rows = read_tsv(manifest) if manifest.exists() else []
    record = dict(resource=str(destination.relative_to(ROOT)), identifier=identifier,
                  source_url=url, download_timestamp=now(), size_bytes=destination.stat().st_size,
                  sha256=sha256(destination), version=version, notes='Original bytes retained; cache timestamp retained on rerun')
    if not any(r['resource'] == record['resource'] for r in rows):
        write_tsv(manifest, rows + [record])
    return destination

def run_command(command, cwd=ROOT, timeout=3600):
    completed = subprocess.run([str(x) for x in command], cwd=cwd, text=True,
                               capture_output=True, timeout=timeout)
    if completed.returncode:
        raise RuntimeError(f'{command}: exit={completed.returncode}\n{completed.stdout}\n{completed.stderr}')
    return completed.stdout

class MeasuredStage:
    """Sample process-tree RSS and project bytes, also use Windows peak WS counters."""
    def __init__(self, name, command):
        self.name, self.command = name, command
        self.done = threading.Event()
        self.peak_rss = 0
        self.peak_disk = 0

    def sample(self):
        import psutil
        try:
            parent = psutil.Process()
            procs = [parent] + parent.children(recursive=True)
            rss = 0
            for process in procs:
                try:
                    info = process.memory_info()
                    rss += max(info.rss, getattr(info, 'peak_wset', 0))
                except psutil.Error:
                    pass
            self.peak_rss = max(self.peak_rss, rss)
            self.peak_disk = max(self.peak_disk, footprint())
        except psutil.Error:
            pass

    def __enter__(self):
        check_resources()
        self.start = now()
        self.t0 = time.monotonic()
        self.before = footprint()
        self.sample()
        def monitor():
            while not self.done.wait(1):
                self.sample()
                if self.peak_rss > config()['ram_bytes']:
                    # Stop children first; leave a useful refusal rather than OOM.
                    import psutil
                    for child in psutil.Process().children(recursive=True):
                        try:
                            child.terminate()
                        except psutil.Error:
                            pass
        self.thread = threading.Thread(target=monitor, daemon=True)
        self.thread.start()
        return self

    def __exit__(self, kind, error, tb):
        self.done.set()
        self.thread.join()
        self.sample()
        record = dict(stage=self.name, command=self.command, start_time=self.start, end_time=now(),
                      elapsed_seconds=round(time.monotonic()-self.t0, 4), exit_status=1 if error else 0,
                      peak_rss_bytes=self.peak_rss, disk_before_bytes=self.before,
                      disk_peak_bytes=self.peak_disk, disk_after_bytes=footprint(),
                      threads=config()['threads'], host_class=platform.platform(),
                      notes='1-second process-tree RSS sampling; Windows peak working-set counters; 1-second logical project bytes')
        append_tsv(ROOT / 'logs/resource_usage.tsv', record)
        write_tsv(ROOT / 'results/resource_summary.tsv', read_tsv(ROOT / 'logs/resource_usage.tsv'))
        if error:
            failure(self.name, self.command, error)
