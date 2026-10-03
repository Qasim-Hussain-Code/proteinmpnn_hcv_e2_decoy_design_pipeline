"""Primary source snapshots; no third-party code is vendored in Git."""
import importlib.metadata
import json
import platform
from .common import ROOT, fetch, now, run_command, sha256, write_json, write_tsv

def sources():
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
    for name in ['ProteinMPNN','EvoEF2']:
        folder = ROOT / 'vendor' / name
        commit = run_command(['git','rev-parse','HEAD'],cwd=folder).strip()
        source_url = 'https://github.com/'+('dauparas/ProteinMPNN' if name == 'ProteinMPNN' else 'tommyhuangthu/EvoEF2')
        version = (folder/'VERSION').read_text().strip() if (folder/'VERSION').exists() else commit
        weights = list(folder.glob('*model_weights/v_48_020.pt')) if name == 'ProteinMPNN' else [folder/'EvoEF2.exe']
        for weight in weights:
            rows.append(dict(tool=name,purpose='fixed-backbone inverse folding' if name=='ProteinMPNN' else 'sidechain modeling and physical scoring',
                             version=version,source_url=source_url,git_commit=commit,release_or_tag='HEAD verified on retrieval',license='MIT' if name=='ProteinMPNN' else 'MIT LICENSE conflicts with academic-use README',
                             license_url=source_url+'/blob/'+commit+'/LICENSE',model_or_weight_file=str(weight.relative_to(ROOT)),weight_hash=sha256(weight),date_checked=now(),
                             local_or_remote='local, ignored',redistributable='not redistributed',
                             notes='Weights bundled in official MIT repository; no separate weight license found; do not infer separate grant' if name=='ProteinMPNN' else 'Local academic computation; binary/code excluded from repository due to ambiguity'))
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
