"""Reproduce the canonical installation with telemetry and wheel hashes."""
import json
import os
from pathlib import Path
import shutil
import sys
import zipfile
from .common import ROOT, MeasuredStage, append_tsv, check_resources, now, run_command, sha256, write_tsv

def main():
    work=ROOT/'.cache/install_validation'
    if (ROOT/'results/installation_report.json').exists() and (ROOT/'results/wheel_manifest.tsv').exists() and (ROOT/'results/installation_validation_complete.json').exists():
        return
    wheels=work/'wheels'; environment=work/'venv'
    check_resources(2_000_000_000,1,'optional_installation_validation')
    wheels.mkdir(parents=True,exist_ok=True)
    temporary=work/'tmp'; temporary.mkdir(exist_ok=True)
    for key in ['TEMP','TMP','TMPDIR','PIP_CACHE_DIR']:
        os.environ[key]=str(temporary)
    with MeasuredStage('installation_validation','Download locked binary wheels; fresh project-local venv; offline install; CPU check'):
        run_command([sys.executable,'-m','pip','download','--no-cache-dir','--only-binary=:all:','--dest',wheels,
                     '-r',ROOT/'requirements.lock.txt','--extra-index-url','https://download.pytorch.org/whl/cpu'],timeout=1800)
        run_command([sys.executable,'-m','venv',environment])
        python=environment/('Scripts/python.exe' if sys.platform=='win32' else 'bin/python')
        run_command([python,'-m','pip','install','--no-index','--find-links',wheels,'-r',ROOT/'requirements.lock.txt',
                     '--report',ROOT/'results/installation_report.json'],timeout=1800)
        run_command([python,'-c','import torch; assert torch.version.cuda is None; print(torch.__version__)'])
        report=json.loads((ROOT/'results/installation_report.json').read_text())
        rows=[]
        for wheel in sorted(wheels.glob('*.whl')):
            with zipfile.ZipFile(wheel) as z:
                metadata=z.read(next(n for n in z.namelist() if n.endswith('.dist-info/METADATA'))).decode('utf-8')
                licensefiles=[n for n in z.namelist() if '/licenses/' in n or n.endswith('/LICENSE') or n.endswith('/LICENSE.txt')]
            rows.append(dict(wheel=wheel.name,size_bytes=wheel.stat().st_size,sha256=sha256(wheel),
                             installed_distribution_metadata_license=[line for line in metadata.splitlines() if line.startswith(('License:','License-Expression:'))],
                             license_files=licensefiles,cpu_torch_verified=True,date_checked=now()))
            append_tsv(ROOT/'results/download_manifest.tsv',dict(resource=str(wheel.relative_to(ROOT)),identifier=wheel.name,
                       source_url='https://download.pytorch.org/whl/cpu' if wheel.name.startswith('torch-') else 'https://pypi.org/simple/',
                       download_timestamp=now(),size_bytes=wheel.stat().st_size,sha256=sha256(wheel),version=wheel.name.split('-')[1],
                       notes='Canonical locked binary-wheel installation reproduced with telemetry; exact file hash retained; disposable wheel removed after validation'))
        write_tsv(ROOT/'results/wheel_manifest.tsv',rows)
        from .wheel_sources import wheel_sources
        wheel_sources()
        from .common import write_json
        write_json(ROOT/'results/installation_validation_complete.json',dict(project_local_temp=True,cpu_only=True,wheel_count=len(rows),timestamp=now()))
    # All disposable targets are absolute and proven project-local before deletion.
    resolved=work.resolve()
    if not resolved.is_relative_to(ROOT.resolve()) or resolved==ROOT.resolve():
        raise RuntimeError('Refuse cleanup outside intended project cache')
    shutil.rmtree(resolved)

if __name__=='__main__':
    main()
