"""Copy a clean scientific template without overwriting an existing analysis."""
import argparse
from pathlib import Path
import shutil
import subprocess
import yaml
from .common import ROOT,config,run_command

def main():
    parser=argparse.ArgumentParser(description='Create a new analysis directory containing code and predeclared settings, without frozen results, environment, weights or raw downloads')
    parser.add_argument('--destination',required=True)
    args=parser.parse_args(); target=Path(args.destination).resolve()
    if target.exists():
        raise RuntimeError('Destination already exists; refuse overwriting an analysis')
    target.mkdir(parents=True)
    for name in ['pipeline','scripts','tests']:
        shutil.copytree(ROOT/name,target/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for name in ['.gitignore','LICENSE','CITATION.cff','pyproject.toml','requirements.lock.txt','run_all.sh']:
        shutil.copyfile(ROOT/name,target/name)
    (target/'config').mkdir()
    shutil.copyfile(ROOT/'config/upstream_lock.json',target/'config/upstream_lock.json')
    cfg=config(); cfg['full_count_per_cell']=None; cfg['scoring_workers']=1
    (target/'config/design.yml').write_text(yaml.safe_dump(cfg,sort_keys=False))
    (target/'docs').mkdir(); shutil.copyfile(ROOT/'docs/methods_notes.md',target/'docs/methods_notes.md')
    (target/'data').mkdir(); shutil.copyfile(ROOT/'data/README.md',target/'data/README.md')
    name=run_command(['git','config','user.name']).strip(); email=run_command(['git','config','user.email']).strip()
    run_command(['git','init'],cwd=target)
    run_command(['git','add','.gitignore'],cwd=target)
    run_command(['git','-c','user.name='+name,'-c','user.email='+email,'commit','-m','add_gitignore'],cwd=target)
    print(f'Created clean analysis template: {target}')
    print('Create its own Python 3.14 environment, install bootstrap smoke dependencies, configure resources, install locked scientific dependencies, then run the pipeline. This starts a new live-query analysis, not a replay of the frozen dataset.')

if __name__=='__main__':
    main()
