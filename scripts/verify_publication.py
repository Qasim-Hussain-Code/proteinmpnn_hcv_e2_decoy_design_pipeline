"""Audit the retained publication tree without requiring ignored model caches."""
from pathlib import Path
import argparse
from collections import Counter
import hashlib
import json
import re
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from pipeline.common import read_tsv,write_json,sha256,now
from pipeline.selection import require_freeze
from pipeline.followup import snapshot
from pipeline.scoring import predicted_direction
from pipeline.verification import language_errors,species_errors

def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT,text=True).splitlines()
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only',action='store_true',help='Print the audit without rewriting the retained report')
    args=parser.parse_args()
    checks=[]
    def check(name,value,detail=None):checks.append(dict(check=name,passed=bool(value),detail=detail))
    require_freeze();snapshot();check('original_scientific_freeze_and_archive',True)
    first=git('rev-list','--max-parents=0','HEAD')[0]
    check('gitignore_first_commit',git('ls-tree','--name-only',first)==['.gitignore'])
    subjects=git('log','--format=%s','HEAD');check('two_or_three_word_commit_messages',all(re.fullmatch(r'[a-z]+(?:_[a-z]+){1,2}',s) for s in subjects),subjects)
    check('no_ai_coauthor',not any(re.search(r'Co-authored-by:.*(?:codex|openai|assistant)',s,re.I) for s in git('log','--format=%B','HEAD')))
    files=git('ls-files');forbidden=('.venv/','vendor/','data/raw/','data/work/','.cache/','followup/states/','followup/logs/')
    check('no_bulk_or_regenerable_caches',not any(p.startswith(forbidden) or (p.startswith('followup/models/') and (p.endswith('_fixed.pdb') or p.endswith('/humanized.pdb'))) for p in files))
    largest=sorted([(p,(ROOT/p).stat().st_size) for p in files],key=lambda x:-x[1])[:5]
    check('tracked_file_limit_50mb',all(n<=50_000_000 for p,n in largest),largest)
    objects=git('rev-list','--objects','HEAD');query='\n'.join(line.split(' ',1)[0] for line in objects)+'\n'
    info=subprocess.check_output(['git','cat-file','--batch-check=%(objectname) %(objecttype) %(objectsize)'],input=query.encode(),cwd=ROOT).decode().splitlines()
    blobs=[(line.split()[0],int(line.split()[2])) for line in info if ' blob ' in line]
    check('all_history_blob_limit_50mb',all(size<=50_000_000 for oid,size in blobs))
    # High-confidence credential signatures only. Never print matched contents.
    patterns=[re.compile(rb'gh[opusr]_[A-Za-z0-9_]{20,}'),re.compile(rb'github_pat_[A-Za-z0-9_]{30,}'),re.compile(rb'AKIA[A-Z0-9]{16}'),re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----')]
    found=[]
    for oid,size in blobs:
        if size>2_000_000:continue
        data=subprocess.check_output(['git','cat-file','blob',oid],cwd=ROOT)
        if any(p.search(data) for p in patterns):found.append(oid)
    check('history_credential_scan',not found,found)
    for row in read_tsv(ROOT/'followup/artifact_manifest.tsv'):check('artifact:'+row['path'],sha256(ROOT/row['path'])==row['sha256'])
    for row in read_tsv(ROOT/'followup/code_manifest.tsv'):check('code:'+row['path'],sha256(ROOT/row['path'])==row['sha256'])
    check('species_provenance',not species_errors(read_tsv(ROOT/'results/structure_provenance.tsv')))
    for path in [ROOT/'README.md',ROOT/'followup/README.md',*(ROOT/'docs').glob('*.md'),ROOT/'publication/README.md']:
        check('language:'+path.relative_to(ROOT).as_posix(),not language_errors(path.read_text(encoding='utf-8')))
    for r in read_tsv(ROOT/'followup/benchmark.tsv'):check('direction:'+r['model']+':'+r['mutation']+':'+r['mode'],predicted_direction(float(r['delta']),.5)==r['predicted'])
    dev=json.loads((ROOT/'followup/development_sequences.json').read_text());reserve=json.loads((ROOT/'followup/quarantine/future_reserve.json').read_text());panel=read_tsv(ROOT/'followup/balanced_development_panel.tsv');represented={r['haplotype'] for r in panel}
    check('expanded_counts',len(dev)==526 and len(reserve)==6 and len(panel)==306)
    observed=sum(r['haplotype'] in represented for r in dev)/len(dev);reported=json.loads((ROOT/'followup/balanced_panel_summary.json').read_text())
    check('coverage_recomputed',abs(observed-reported['coverage'])<1e-12 and observed>=.9)
    check('coverage_by_genotype',all(sum(r['genotype']==g and r['haplotype'] in represented for r in dev)/sum(r['genotype']==g for r in dev)>=.9 for g in {r['genotype'] for r in dev}))
    check('reserve_components_disjoint',not ({r['component'] for r in reserve}&{r['component'] for r in dev}))
    truth=read_tsv(ROOT/'followup/ground_truth.tsv');check('corrected_assay_counts',sum(r['included']=='true' for r in truth)==34 and sum(r['included']=='true' and r['followup_role']!='external_source_challenge' for r in truth)==28)
    check('soluble_double_mutant_attribution',not any(r['included']=='true' and r['source_id']=='higginbottom2000' and 'GST' in r['cd81_construct'] and r['mutation'] in {'F186L','E188K'} for r in truth))
    check('negative_benchmark_visible','Scoring remains unreliable' in (ROOT/'README.md').read_text())
    check('original_and_data_licenses',(ROOT/'LICENSE').read_text().startswith('MIT License') and len(read_tsv(ROOT/'publication/data_licenses.tsv'))==5)
    report=dict(timestamp=now(),passed=all(r['passed'] for r in checks),check_count=len(checks),failed=[r for r in checks if not r['passed']],largest_tracked_files=largest,history_blobs_scanned=len(blobs),scientific_interpretation='Reproducible computational results; failed biological benchmark',resource_note='Historical unmeasured follow-up peaks are not claimed.')
    if not args.check_only:
        write_json(ROOT/'publication/publication_checks.json',report)
    print(json.dumps(report,indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
