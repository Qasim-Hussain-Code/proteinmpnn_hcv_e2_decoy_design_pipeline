"""Resolve exact official wheel URLs and verify them against measured install hashes."""
import html
import json
import re
from urllib.parse import unquote,urljoin
from .common import ROOT,fetch,read_tsv,write_tsv

def wheel_sources():
    rows=[]
    for wheel in read_tsv(ROOT/'results/wheel_manifest.tsv'):
        filename=wheel['wheel']; name,version=filename.split('-')[:2]
        if name=='torch':
            source='https://download.pytorch.org/whl/cpu/torch/'
            page=fetch(source,ROOT/'data/raw/torch_cpu_wheel_index.html','torch_CPU_index')
            links=re.findall(r'href="([^"]+)"',page.read_text(encoding='utf-8'))
            matches=[html.unescape(urljoin(source,link)) for link in links if unquote(link.split('#')[0].split('/')[-1])==filename]
            match=next((url for url in matches if '#sha256='+wheel['sha256'] in url),None)
            if not match:
                raise RuntimeError(f'Official CPU wheel URL/hash not resolved: {filename}')
            url=match.split('#')[0]
        else:
            source=f'https://pypi.org/pypi/{name}/{version}/json'
            snapshot=fetch(source,ROOT/f'data/raw/wheel_metadata/{name}-{version}.json',name,version)
            data=json.loads(snapshot.read_text(encoding='utf-8'))
            match=next((r for r in data['urls'] if r['filename']==filename and r['digests']['sha256']==wheel['sha256']),None)
            if not match:
                raise RuntimeError(f'Official wheel filename/hash mismatch: {filename}')
            url=match['url']
        rows.append(dict(wheel=filename,version=version,official_download_url=url,sha256=wheel['sha256'],
                         verified_against='Official filename and published SHA256 equal the measured installation wheel',metadata_source=source))
    write_tsv(ROOT/'results/wheel_source_urls.tsv',rows)
    manifest=read_tsv(ROOT/'results/download_manifest.tsv'); lookup={r['wheel']:r for r in rows}
    for row in manifest:
        if row['identifier'] in lookup and '.whl' in row['resource']:
            row['source_url']=lookup[row['identifier']]['official_download_url']
            row['notes']+='; exact official URL verified against published SHA256'
    write_tsv(ROOT/'results/download_manifest.tsv',manifest)

if __name__=='__main__':
    wheel_sources()
