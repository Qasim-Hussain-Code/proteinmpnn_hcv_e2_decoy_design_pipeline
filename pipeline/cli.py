from __future__ import annotations
import argparse
import importlib
import json
import os
from pathlib import Path
import shutil
import sys
from .common import ROOT, MeasuredStage, budget, config, now, write_json

STAGES = ['sources', 'structures', 'humanize', 'ground_truth', 'sequences',
          'diversity', 'states', 'benchmark', 'pilot', 'design', 'discovery',
          'freeze', 'heldout', 'sensitivity', 'figures', 'report', 'verify', 'smoke']

def main():
    parser = argparse.ArgumentParser(description='CPU fixed-backbone CD81 benchmark')
    parser.add_argument('stage', choices=['configure', 'all'] + STAGES)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--ram', type=float, default=14)
    parser.add_argument('--disk', type=float, default=13)
    parser.add_argument('--seed', type=int, default=20261004)
    parser.add_argument('--yes', action='store_true')
    args = parser.parse_args()
    os.chdir(ROOT)
    if args.stage == 'configure':
        free = shutil.disk_usage(ROOT).free
        resolved = dict(threads=args.threads, ram_bytes=int(args.ram*1e9), seed=args.seed,
                        current_free_bytes=free, effective_disk_budget=budget(args.disk*1e9, free),
                        reserve_bytes=1_000_000_000, timestamp=now(), python=sys.executable)
        if args.threads < 1 or args.ram <= 0 or args.disk <= 0:
            raise ValueError('Threads, RAM and disk must be positive')
        write_json(ROOT / 'config/resolved.json', resolved)
        from .common import write_tsv
        write_tsv(ROOT/'results/configuration.tsv',[dict(parameter=k,value=v) for k,v in resolved.items()])
        (ROOT / 'project.conf').write_text(''.join(f'{k.upper()}={json.dumps(v)}\n' for k,v in resolved.items()
                                                 if k not in {'timestamp', 'python'})+
                                          f'PYTHON="{sys.executable.replace(chr(92), chr(47))}"\n', encoding='utf-8')
        import yaml
        cfg = config()
        cfg.update(threads=args.threads, ram_bytes=resolved['ram_bytes'], master_seed=args.seed)
        (ROOT / 'config/design.yml').write_text(yaml.safe_dump(cfg, sort_keys=False), encoding='utf-8')
        print(json.dumps(resolved, indent=2))
        return
    stages = STAGES[:-1] if args.stage == 'all' else [args.stage]
    for stage in stages:
        module_name = {'sources':'sources','structures':'structures','humanize':'structures',
                       'ground_truth':'literature','sequences':'diversity','diversity':'diversity',
                       'states':'scoring','benchmark':'scoring','pilot':'design','design':'design',
                       'discovery':'selection','freeze':'selection','heldout':'selection',
                       'sensitivity':'selection','figures':'reporting','report':'reporting',
                       'verify':'verification','smoke':'verification'}[stage]
        print(f'START {stage}', flush=True)
        with MeasuredStage(stage, ' '.join(sys.argv)):
            if stage in {'structures','humanize','ground_truth','sequences','diversity','states','benchmark'} and (ROOT/'results/candidate_freeze.tsv').exists():
                from .selection import require_freeze
                require_freeze()
                print('Frozen scientific inputs retained; use a fresh checkout for a new run.',flush=True)
            else:
                getattr(importlib.import_module('pipeline.'+module_name), stage)()
        print(f'DONE {stage}', flush=True)

if __name__ == '__main__':
    main()
