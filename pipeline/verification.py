"""Offline smoke test and final checks that fail on scientific provenance errors."""
from __future__ import annotations
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from .common import ROOT, config, footprint, read_tsv, run_command, sha256, write_json, write_tsv

def smoke():
    from .diversity import deterministic_split,entropy
    from .selection import freeze_ids,heldout_read
    from .scoring import mutation_syntax,predicted_direction
    from .structures import read_structure,interface
    fixture=ROOT/'tests/fixtures/tiny.pdb'; structure=read_structure(fixture)
    assert interface(structure[0]['E'],structure[0]['R'])['atom_contacts']==1
    data=[dict(accession=str(i),genotype='1',subtype='1a',haplotype=h,sequence_hash=h) for i,h in enumerate(['AAAA','CCCC','DDDD','EEEE','FFFF','GGGG'])]
    split=deterministic_split(data,20261004)
    assert {r['split'] for r in split}=={'discovery','held_out'}
    assert entropy('AACC')==1
    assert mutation_syntax('T','R',163,'A','T')=='TR163A'
    assert predicted_direction(2)=='reduced'
    metrics=[dict(candidate_id='synthetic_1',reference_score=1,worst_delta_wt=3,median_delta_wt=2),dict(candidate_id='synthetic_2',reference_score=2,worst_delta_wt=0,median_delta_wt=0)]
    selected,front=freeze_ids(metrics,{c:dict(hydrophobic_sasa_fraction=.4) for c in ['synthetic_1','synthetic_2']},{},1)
    assert selected=={'synthetic_1':['single_state'],'synthetic_2':['escape_aware']}
    (ROOT/'data/work').mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=ROOT/'data/work') as directory:
        root=Path(directory); secret=root/'heldout.json'; secret.write_text('[]')
        try:
            heldout_read(secret,root)
        except RuntimeError:
            pass
        else:
            raise AssertionError('Held-out invisibility failed')
    write_json(ROOT/'results/smoke_test.json',dict(passed=True,fixture='synthetic geometry and analytical score fixture',
               internet=False,gpu=False,external_models=False,scientific_outputs=False,seed=20261004))
    print('Offline synthetic end-to-end fixture passed; no scientific score or sequence result inferred.')

def language_errors(text):
    problems=[]
    if '\u2014' in text:
        problems.append('em dash')
    banned=['it is worth noting','it is important to note',"in today's rapidly evolving",'plays a crucial role','serves as a testament','paving the way for','in conclusion','delve','seamless','underscores','showcases']
    for phrase in banned:
        if re.search(r'\b'+re.escape(phrase)+r'\b',text,re.I):
            problems.append(phrase)
    for line in text.splitlines():
        for term in ['binding energy','binding affinity','free energy','high affinity','binder','neutralizer']:
            if re.search(r'\b'+term+r'\b',line,re.I):
                # Negative statements and explicitly experimental/source-based
                # language are permitted; affirmative score-as-affinity is not.
                if not re.search(r'not |does not |cannot |experimental|rigorous|no |primary|paper|published|not estimate',line,re.I):
                    problems.append('unsupported score language: '+line)
    if re.search('[\U0001F300-\U0001FAFF]',text):
        problems.append('emoji')
    return problems

def species_errors(rows):
    errors=[]
    for row in rows:
        if 'CD81' in row['description'].upper():
            expected={'7MWX':'Saguinus oedipus','3X0E':'Homo sapiens'}.get(row['pdb_id'])
            if expected and row['organism']!=expected:
                errors.append(f'{row["pdb_id"]} receptor species incorrect')
    return errors

def verify():
    from .selection import require_freeze
    findings=[]
    def check(name,passed,details=''):
        findings.append(dict(check=name,passed=passed,details=details))
    check('species',not species_errors(read_tsv(ROOT/'results/structure_provenance.tsv')))
    for path in [ROOT/'README.md',*list((ROOT/'docs').glob('*.md'))]:
        errors=language_errors(path.read_text(encoding='utf-8')); check('language_'+path.name,not errors,';'.join(errors))
    require_freeze(); check('freeze_integrity',True)
    freeze=json.loads((ROOT/'results/candidate_freeze_manifest.json').read_text())
    from .structures import sequence,read_structure
    frozen_rows=read_tsv(ROOT/'results/candidate_freeze.tsv')
    check('frozen_model_sequences',all(sequence(read_structure(ROOT/f'data/processed/frozen_candidates/{r["candidate_id"]}.pdb')[0]['R'])==r['sequence'] for r in frozen_rows))
    if (ROOT/'results/cleanup_summary.json').exists():
        cleanup=json.loads((ROOT/'results/cleanup_summary.json').read_text())
        manifest=ROOT/'results/disposable_model_manifest.tsv'
        check('cleanup_manifest_integrity',sha256(manifest)==cleanup['manifest_sha256'])
        removed=read_tsv(manifest)
        check('retained_model_hashes',all(sha256(ROOT/r['retained_model'])==r['sha256'] for r in removed if r['retained_model']))
        check('disposable_models_removed',all(not (ROOT/r['path']).exists() for r in removed))
    successful_heldout=[r for r in read_tsv(ROOT/'logs/resource_usage.tsv') if r['stage']=='heldout' and r['exit_status']=='0']
    check('freeze_precedes_successful_heldout',bool(successful_heldout) and all(r['start_time']>freeze['selection_timestamp'] for r in successful_heldout))
    for row in read_tsv(ROOT/'results/leakage_audit.tsv'):
        check(row['check'],row['passed'].lower()=='true',row['value'])
    for pdb in ['7MWX','3X0E']:
        source=ROOT/f'data/raw/{pdb}.cif'
        hashes={r['sha256'] for r in read_tsv(ROOT/'results/structure_provenance.tsv') if r['pdb_id']==pdb}
        check('source_unchanged_'+pdb,sha256(source) in hashes)
    states=read_tsv(ROOT/'results/e2_state_manifest.tsv')
    for row in states:
        check('state_hash_'+row['haplotype_id']+'_'+row['split'],sha256(ROOT/row['path'])==row['structure_hash'])
    import torch
    check('cpu_torch',torch.version.cuda is None,torch.__version__)
    resource=read_tsv(ROOT/'logs/resource_usage.tsv'); peak=max(int(r['disk_peak_bytes']) for r in resource)
    check('sampled_disk_hard_ceiling',peak<=13_000_000_000,str(peak))
    check('rss_ceiling',max(int(r['peak_rss_bytes']) for r in resource)<=config()['ram_bytes'])
    import shutil
    check('disk_reserve',shutil.disk_usage(ROOT).free>=config()['reserve_bytes'])
    check('figures',all((ROOT/r['file']).exists() and sha256(ROOT/r['file'])==r['sha256'] for r in read_tsv(ROOT/'results/figure_manifest.tsv')))
    check('traceability_sources',all((ROOT/r['source_file']).exists() for r in read_tsv(ROOT/'results/readme_traceability.tsv')))
    shellcheck=shutil.which('shellcheck') or (str(ROOT/'vendor/shellcheck/shellcheck.exe') if (ROOT/'vendor/shellcheck/shellcheck.exe').exists() else None)
    check('shellcheck_available',shellcheck is not None,'unavailable: Bash syntax checked; shellcheck not falsely marked passed')
    for path in [ROOT/'run_all.sh',*list((ROOT/'scripts').glob('*.sh'))]:
        result=subprocess.run(['C:/Program Files/Git/bin/bash.exe' if sys.platform=='win32' else 'bash','-n',str(path)],capture_output=True,text=True)
        check('bash_syntax_'+path.name,result.returncode==0,result.stderr)
        if shellcheck:
            result=subprocess.run([shellcheck,'-x',str(path)],cwd=ROOT,capture_output=True,text=True)
            check('shellcheck_'+path.name,result.returncode==0,result.stdout)
    files=run_command(['git','ls-files']).splitlines()
    largest=sorted([(f,(ROOT/f).stat().st_size) for f in files if (ROOT/f).exists()],key=lambda p:-p[1])[:5]
    write_tsv(ROOT/'results/largest_tracked_files.tsv',[dict(path=p,size_bytes=n) for p,n in largest])
    check('tracked_size_limit',all((ROOT/f).stat().st_size<=50_000_000 for f in files if (ROOT/f).exists()))
    check('ignored_bulk_files',not any(f.startswith(('.venv/','vendor/','data/raw/','data/work/','.cache/')) for f in files))
    generated=read_tsv(ROOT/'results/generated_sequences.tsv.gz')
    from .design import constraints
    policy=read_tsv(ROOT/'config/fixed_positions.tsv'); scaffold=json.loads((ROOT/'data/processed/humanization.json').read_text())
    check('independent_sequence_constraints',all(constraints(r['sequence'],scaffold['human_sequence'],scaffold['human_positions'],
        [int(p['sequence_index'])-1 for p in policy if p['design_arm']==r['design_arm'] and p['designable']=='false'],config()['disulfide_positions'])==r['filter_status'] for r in generated))
    from .common import validate_schema
    validate_schema(generated,['candidate_id','sequence','seed','effective_seed','ProteinMPNN_model','filter_status'])
    validate_schema(read_tsv(ROOT/'results/candidate_freeze.tsv'),['candidate_id','sequence','config_hash','software_manifest_hash','selection_timestamp'])
    validate_schema(read_tsv(ROOT/'results/score_dictionary.tsv'),['score_name','tool','version','mathematical_or_algorithmic_definition','direction','units_or_unitless','what_it_can_support','what_it_cannot_support'])
    check('required_schemas',True)
    check('effective_seeds',all(r['seed']==r['effective_seed'] and int(r['seed'])!=0 for r in generated))
    cells={(r['design_arm'],r['ProteinMPNN_model'],r['seed'],r['temperature']) for r in generated}
    check('matched_generation_budget',len(cells)==40 and all(sum((r['design_arm'],r['ProteinMPNN_model'],r['seed'],r['temperature'])==cell for r in generated)==config()['full_count_per_cell'] for cell in cells))
    for name in ['discovery_scores.tsv.gz','heldout_scores.tsv']:
        scores=read_tsv(ROOT/'results'/name)
        expected_ids={r['haplotype_id'] for r in scores}
        wt=[r for r in scores if r['candidate_id']=='wild_type' and r['score_status']=='passed']
        check('wild_type_comparator_'+name,{r['haplotype_id'] for r in wt}==expected_ids and len(wt)==len(expected_ids))
    from .structures import read_structure
    from Bio.SeqUtils import seq1
    base=read_structure(ROOT/'data/processed/humanized_repaired.pdb')
    state_mutations_valid=True
    for row in states:
        model=read_structure(ROOT/row['path'])
        for mutation in row['mutations_from_reference'].split(',') if row['mutations_from_reference'] else []:
            wt,pos,new=mutation[0],int(mutation[1:-1]),mutation[-1]
            state_mutations_valid &= seq1(base[0]['E'][(' ',pos,' ')].resname)==wt and seq1(model[0]['E'][(' ',pos,' ')].resname)==new
    check('state_mutation_reference_and_identity',bool(state_mutations_valid))
    software=read_tsv(ROOT/'results/software_manifest.tsv')
    check('software_model_hashes',all(not r['weight_hash'] or sha256(ROOT/r['model_or_weight_file'])==r['weight_hash'] for r in software))
    check('downloaded_source_hashes',all(not (ROOT/r['resource']).exists() or sha256(ROOT/r['resource'])==r['sha256'] for r in read_tsv(ROOT/'results/download_manifest.tsv')))
    check('fixed_positions',all(r['filter_status'] not in {'fixed_position_violation','disulfide_violation'} for r in generated))
    check('benchmark_warning_propagates',all('failed' in r['warning'].lower() for r in read_tsv(ROOT/'results/final_computational_priorities.tsv')))
    required=['structure_provenance','humanization_metrics','software_manifest','download_manifest','cd81_mutation_benchmark','hcv_sequence_manifest',
              'e2_position_conservation','e2_interface_residues','e2_haplotypes','leakage_audit','e2_state_manifest','design_filter_funnel','candidate_freeze',
              'heldout_scores','heldout_metrics','seed_variance','sensitivity_analysis','resource_summary','score_dictionary','readme_traceability','scientific_discomfort']
    check('required_results',all((ROOT/f'results/{name}.tsv').exists() and read_tsv(ROOT/f'results/{name}.tsv') for name in required))
    check('compressed_results',all((ROOT/f'results/{name}.tsv.gz').exists() for name in ['generated_sequences','discovery_scores']))
    write_tsv(ROOT/'results/verification.tsv',findings)
    failures=[r for r in findings if not r['passed'] and r['check']!='shellcheck_available']
    if failures:
        raise RuntimeError(f'Verification failed: {failures}')
    print(json.dumps(dict(passed=True,shellcheck_available=bool(shellcheck),largest_tracked_files=largest),indent=2))
