"""Primary source snapshots; no third-party code is vendored in Git."""
import importlib.metadata
import json
import platform
import shutil
import sys
import zipfile
from .common import ROOT, fetch, now, run_command, sha256, write_json, write_tsv

def sources():
    existing=ROOT/'results/software_manifest.tsv'
    if existing.exists() and (ROOT/'results/candidate_freeze.tsv').exists():
        from .common import read_tsv
        for row in read_tsv(existing):
            if row['weight_hash'] and sha256(ROOT/row['model_or_weight_file'])!=row['weight_hash']:
                raise RuntimeError('Pinned executable/model hash changed after freeze')
        return
    rows = []
    for name, purpose in [('numpy','numerical arrays'), ('biopython','mmCIF parsing, sequence alignment and SASA'),
                          ('scipy','distance and statistics'),('matplotlib','figures'),('psutil','resource measurement'),
                          ('pytest','tests'),('pyyaml','configuration'),('torch','CPU inverse folding')]:
        dist = importlib.metadata.distribution(name)
        rows.append(dict(tool=name,purpose=purpose,version=dist.version,source_url=f'https://pypi.org/project/{name}/',
                         git_commit='',release_or_tag=dist.version,license=dist.metadata.get('License-Expression') or dist.metadata.get('License','See installed license'),
                         license_url=f'https://pypi.org/project/{name}/',model_or_weight_file='',weight_hash='',date_checked=now(),
                         local_or_remote='local',redistributable='review upstream terms',notes='Installed metadata preserved in requirements.lock.txt'))
    rows.append(dict(tool='Python',purpose='runtime',version=platform.python_version(),source_url='https://www.python.org/',
                     git_commit='',release_or_tag=platform.python_version(),license='PSF',license_url='https://docs.python.org/3/license.html',
                     model_or_weight_file='',weight_hash='',date_checked=now(),local_or_remote='pre-existing',redistributable='yes',notes='Host installation'))
    for name,command,url,terms in [
        ('Git',['git','--version'],'https://git-scm.com/about.html','GPL-2.0'),
        ('Bash',['C:/Program Files/Git/bin/bash.exe' if sys.platform=='win32' else 'bash','--version'],'https://www.gnu.org/software/bash/','GPL-3.0-or-later'),
        ('GCC',['g++','--version'],'https://gcc.gnu.org/','GPL-3.0-or-later; runtime library exception 3.1')]:
        version=run_command(command).splitlines()[0]
        rows.append(dict(tool=name,purpose='Local execution and source compilation',version=version,source_url=url,
                         git_commit='',release_or_tag=version,license=terms,license_url=url,model_or_weight_file='',weight_hash='',
                         date_checked=now(),local_or_remote='pre-existing',redistributable='not redistributed',notes='Host command version recorded; no system tool downloaded'))
    shellcheck=ROOT/'vendor/shellcheck/shellcheck.exe'
    if sys.platform=='win32' and not shellcheck.exists():
        archive=fetch('https://github.com/koalaman/shellcheck/releases/download/v0.11.0/shellcheck-v0.11.0.zip',
                      ROOT/'data/raw/shellcheck-v0.11.0.zip','shellcheck','v0.11.0')
        with zipfile.ZipFile(archive) as bundle:
            executable=next(n for n in bundle.namelist() if n.endswith('shellcheck.exe'))
            shellcheck.parent.mkdir(parents=True,exist_ok=True)
            shellcheck.write_bytes(bundle.read(executable))
    binary=str(shellcheck) if shellcheck.exists() else shutil.which('shellcheck')
    if binary:
        rows.append(dict(tool='ShellCheck',purpose='Bash static analysis',version=run_command([binary,'--version']).split('version: ')[-1].splitlines()[0],
                         source_url='https://github.com/koalaman/shellcheck',git_commit='',release_or_tag='v0.11.0' if shellcheck.exists() else 'host installation',
                         license='GPL-3.0-or-later',license_url='https://github.com/koalaman/shellcheck/blob/v0.11.0/LICENSE',
                         model_or_weight_file=str(shellcheck.relative_to(ROOT)) if shellcheck.exists() else '',weight_hash=sha256(shellcheck) if shellcheck.exists() else '',
                         date_checked=now(),local_or_remote='local ignored',redistributable='not redistributed',notes='Official release or host installation; static analysis only'))
    for name in ['ProteinMPNN','EvoEF2']:
        folder = ROOT / 'vendor' / name
        commit = run_command(['git','rev-parse','HEAD'],cwd=folder).strip()
        source_url = 'https://github.com/'+('dauparas/ProteinMPNN' if name == 'ProteinMPNN' else 'tommyhuangthu/EvoEF2')
        version = (folder/'VERSION').read_text().strip() if (folder/'VERSION').exists() else commit
        weights = list(folder.glob('*model_weights/v_48_020.pt')) if name == 'ProteinMPNN' else [folder/('EvoEF2_local.exe' if (folder/'EvoEF2_local.exe').exists() else 'EvoEF2.exe')]
        for weight in weights:
            rows.append(dict(tool=name,purpose='fixed-backbone inverse folding' if name=='ProteinMPNN' else 'sidechain modeling and physical scoring',
                             version=version,source_url=source_url,git_commit=commit,release_or_tag='HEAD verified on retrieval',license='MIT' if name=='ProteinMPNN' else 'MIT LICENSE conflicts with academic-use README',
                             license_url=source_url+'/blob/'+commit+'/LICENSE',model_or_weight_file=str(weight.relative_to(ROOT)),weight_hash=sha256(weight),date_checked=now(),
                             local_or_remote='local, ignored',redistributable='not redistributed',
                             notes='Weights bundled in official MIT repository; no separate weight license found; do not infer separate grant' if name=='ProteinMPNN' else 'Local academic computation; static compilation: g++ -O3 -static -o EvoEF2_local.exe src/*.cpp; binary/code excluded due to terms ambiguity'))
        for path in [folder/'LICENSE',folder/'README.md']+weights:
            from .common import append_tsv
            append_tsv(ROOT/'results/download_manifest.tsv',dict(resource=str(path.relative_to(ROOT)),identifier=name,
                       source_url=source_url+'/blob/'+commit+'/'+str(path.relative_to(folder)).replace('\\','/'), download_timestamp=now(),
                       size_bytes=path.stat().st_size,sha256=sha256(path),version=commit,notes='Retrieved by Git; selective files audited'))
    import torch
    if torch.version.cuda is not None or '+cpu' not in torch.__version__:
        raise RuntimeError('Refuse CUDA or unverified CPU PyTorch build')
    write_tsv(ROOT/'results/software_manifest.tsv',rows)
    (ROOT/'requirements.lock.txt').write_text(run_command([__import__('sys').executable,'-m','pip','freeze']),encoding='utf-8')
    write_json(ROOT/'config/upstream_lock.json',{r['tool']:r['git_commit'] for r in rows if r['git_commit']})
    resources=[('PDB','https://www.wwpdb.org/about/usage','CC0/public archive usage; original files retrieved, not bulk redistributed'),
               ('GenBank','https://www.ncbi.nlm.nih.gov/home/about/policies/','Public database; submitter rights may apply; raw records not redistributed'),
               ('PMC articles','https://pmc.ncbi.nlm.nih.gov/about/copyright/','Article-specific copyright; original prose and figures not redistributed'),
               ('Original code','LICENSE','MIT; author from existing Git identity')]
    write_tsv(ROOT/'results/license_audit.tsv',[dict(resource=n,license_source=u,license=l,date_checked=now(),redistribution_status='derived factual tables only' if n!='Original code' else 'MIT') for n,u,l in resources])
    tree=[]
    for name in ['ProteinMPNN','EvoEF2']:
        folder=ROOT/'vendor'/name
        commit=run_command(['git','rev-parse','HEAD'],cwd=folder).strip()
        for path in sorted(folder.rglob('*')):
            if path.is_file() and '.git' not in path.parts:
                tree.append(dict(repository=name,git_commit=commit,relative_path=str(path.relative_to(folder)),size_bytes=path.stat().st_size,sha256=sha256(path),
                                 origin='locally compiled' if path.name=='EvoEF2_local.exe' else 'upstream Git checkout'))
    write_tsv(ROOT/'results/upstream_file_manifest.tsv',tree)
