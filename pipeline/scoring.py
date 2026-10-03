"""Identical preparation and EvoEF2 scoring for WT, controls, and candidates."""
from __future__ import annotations
import copy
import json
import os
from pathlib import Path
import re
import shutil
import time
from collections import Counter
from scipy.stats import binomtest
from Bio.SeqUtils import seq1
from .common import ROOT, config, failure, now, read_tsv, run_command, sha256, write_json, write_tsv
from .structures import interface, read_structure, save_pdb, sequence

def evo_binary():
    folder=ROOT/'vendor/EvoEF2'
    for name in ['EvoEF2_local.exe','EvoEF2','EvoEF2.exe']:
        if (folder/name).exists():
            return folder/name
    raise RuntimeError('EvoEF2 executable missing; run scripts/02_install.sh')

def evo(command,pdb,work=None,extra=()):
    work=Path(work or ROOT/'data/work/evo'/str(os.getpid()))
    work.mkdir(parents=True,exist_ok=True)
    # Upstream's PDB path buffer is short. Pass a local filename, never a long
    # Windows absolute path; preserve the input bytes in a disposable workspace.
    local=work/Path(pdb).name
    if Path(pdb).resolve()!=local.resolve():
        shutil.copyfile(pdb,local)
    cmd=[evo_binary(),f'--command={command}',f'--pdb={local.name}',*extra]
    output=run_command(cmd,cwd=work)
    log=ROOT/'logs/evo'/f'{command}_{Path(pdb).stem}_{os.getpid()}.txt'
    log.parent.mkdir(parents=True,exist_ok=True); log.write_text(output,encoding='utf-8')
    if 'invalid' in output.lower() and command.startswith('Compute'):
        raise RuntimeError(f'Refuse scoring incomplete sidechains: {log}')
    return output

def repair(pdb,destination):
    destination=Path(destination)
    if destination.exists():
        model=read_structure(destination)
        if not list(model.get_atoms()) or any('CA' not in r for chain in model[0] for r in chain if r.id[0]==' '):
            raise RuntimeError(f'Incomplete repaired model cache: {destination}')
        return destination
    work=ROOT/'data/work/repair'/str(os.getpid()); work.mkdir(parents=True,exist_ok=True)
    evo('RepairStructure',pdb,work)
    output=work/(Path(pdb).stem+'_Repair.pdb')
    if not output.exists():
        raise RuntimeError(f'RepairStructure did not create {output}')
    destination.parent.mkdir(parents=True,exist_ok=True); shutil.move(str(output),destination)
    return destination

def mutation_syntax(wild_type,chain,position,new,reference):
    if reference != wild_type:
        raise ValueError(f'Reference residue mismatch: expected {wild_type}, got {reference} at {chain}{position}')
    if new not in 'ACDEFGHIKLMNPQRSTVWY':
        raise ValueError('Noncanonical mutation')
    return f'{wild_type}{chain}{position}{new}'

def mutate(pdb,mutations,destination):
    destination=Path(destination)
    if destination.exists():
        check=read_structure(destination)
        for chain,pos,wt,new in mutations:
            if seq1(check[0][chain][(' ',int(pos),' ')].resname)!=new:
                raise RuntimeError(f'Cached mutation identity mismatch: {destination}')
        return destination
    if not mutations:
        shutil.copyfile(pdb,destination); return destination
    structure=read_structure(pdb)
    syntax=[]
    for chain,pos,wt,new in mutations:
        reference=seq1(structure[0][chain][(' ',int(pos),' ')].resname)
        syntax.append(mutation_syntax(wt,chain,pos,new,reference))
    work=ROOT/'data/work/mutate'/str(os.getpid()); work.mkdir(parents=True,exist_ok=True)
    mutation_file=work/'individual_list.txt'; mutation_file.write_text(','.join(syntax)+';\n')
    evo('BuildMutant',pdb,work,extra=[f'--mutant_file={mutation_file.name}'])
    output=work/(Path(pdb).stem+'_Model_0001.pdb')
    if not output.exists():
        output=work/(Path(pdb).stem+'_Model_1.pdb')
    if not output.exists():
        raise RuntimeError('BuildMutant output missing')
    repair(output,destination)
    output.unlink()
    check=read_structure(destination)
    for chain,pos,wt,new in mutations:
        if seq1(check[0][chain][(' ',int(pos),' ')].resname)!=new:
            raise RuntimeError('BuildMutant silently failed mutation identity check')
    return destination

def energy(pdb,kind='ComputeBinding'):
    output=evo(kind,pdb,extra=['--split=E,R'] if kind=='ComputeBinding' else [])
    values=re.findall(r'^Total\s*=\s*([-+\d.eE]+)',output,re.M)
    if len(values)!=1:
        raise RuntimeError(f'Expected exactly one total from {kind}: found {values}')
    return float(values[0])

def combine(receptor_pdb,e2_pdb,destination):
    # Preserve the exact already-written coordinate records. Entity deepcopy
    # follows parent links and needlessly copies a whole complex each state.
    e2lines=[line for line in Path(e2_pdb).read_text().splitlines(True) if line.startswith('ATOM') and line[21]=='E']
    rlines=[line for line in Path(receptor_pdb).read_text().splitlines(True) if line.startswith('ATOM') and line[21]=='R']
    Path(destination).write_text(''.join(e2lines)+'TER\n'+''.join(rlines)+'TER\nEND\n')
    return destination

def states():
    base=repair(ROOT/'data/processed/humanized_unrepaired.pdb',ROOT/'data/processed/humanized_repaired.pdb')
    audit=interface(read_structure(base)[0]['E'],read_structure(base)[0]['R'])
    write_tsv(ROOT/'results/repair_audit.tsv',[dict(input='humanized_unrepaired.pdb',output='humanized_repaired.pdb',
              clashes_after=audit['clashes'],contacts_after=audit['atom_contacts'],sha256=sha256(base),
              num_of_runs=1,method='EvoEF2 verified source default; no backbone minimization')])
    build_states('discovery')

def build_states(split):
    if split=='held_out':
        from .selection import require_freeze
        require_freeze()
    base=ROOT/'data/processed/humanized_repaired.pdb'
    haplotypes=read_tsv(ROOT/f'data/processed/{split}_panel.tsv')
    mapping={int(r['h77_position']):int(r['structure_position']) for r in read_tsv(ROOT/'results/e2_numbering.tsv') if r['h77_position']}
    structure=read_structure(base)
    rows=[]
    folder=ROOT/f'data/processed/states/{split}'; folder.mkdir(parents=True,exist_ok=True)
    reference=folder/'reference.pdb'; shutil.copyfile(base,reference)
    rows.append(dict(haplotype_id='reference',split=split,genotype='experimental_template',subtype_distribution='not_population_sample',number_of_sequences=0,
                     mutations_from_reference='',modeled_positions='',unmodeled_positions='',repair_status='repaired',clash_count_before='',
                     clash_count_after=interface(structure[0]['E'],structure[0]['R'])['clashes'],structure_hash=sha256(reference),path=str(reference.relative_to(ROOT))))
    for row in haplotypes:
        mutations=[]
        for key,new in json.loads(row['residues']).items():
            hpos=int(key)
            if hpos not in mapping:
                raise RuntimeError('State contains unmapped interface residue')
            pos=mapping[hpos]; wt=seq1(structure[0]['E'][(' ',pos,' ')].resname)
            if wt!=new:
                mutations.append(('E',pos,wt,new))
        target=folder/(row['haplotype_id']+'.pdb')
        mutate(base,mutations,target)
        model=read_structure(target); audit=interface(model[0]['E'],model[0]['R'])
        rows.append(dict(haplotype_id=row['haplotype_id'],split=split,genotype=row['genotypes'],subtype_distribution=row['subtypes'],number_of_sequences=row['count'],
                         mutations_from_reference=','.join(f'{wt}{pos}{new}' for chain,pos,wt,new in mutations),modeled_positions=','.join(str(pos) for chain,pos,wt,new in mutations),
                         unmodeled_positions='',repair_status='BuildMutant plus RepairStructure verified source defaults',clash_count_before='not_defined_for_missing_sidechains',
                         clash_count_after=audit['clashes'],structure_hash=sha256(target),path=str(target.relative_to(ROOT))))
    manifest=ROOT/'results/e2_state_manifest.tsv'
    prior=[r for r in read_tsv(manifest) if r['split']!=split] if manifest.exists() else []
    write_tsv(manifest,prior+rows)

def predicted_direction(delta,tolerance=0.5):
    return 'enhanced' if delta < -tolerance else 'reduced' if delta>tolerance else 'approximately_neutral'

def benchmark():
    base=repair(ROOT/'data/processed/humanized_unrepaired.pdb',ROOT/'data/processed/humanized_repaired.pdb')
    wt_score=energy(base)
    rows=[]; cache={}
    for truth in read_tsv(ROOT/'config/cd81_mutation_ground_truth.tsv'):
        if truth['included']!='true':
            continue
        mut=truth['mutation']; match=re.fullmatch(r'([A-Z])(\d+)([A-Z])',mut)
        if not match:
            continue
        wt,pos,new=match.groups(); pos=int(pos)
        if mut not in cache:
            out=ROOT/f'data/processed/controls/{mut}.pdb'; out.parent.mkdir(parents=True,exist_ok=True)
            mutate(base,[('R',pos,wt,new)],out)
            cache[mut]=(energy(out),interface(read_structure(out)[0]['E'],read_structure(out)[0]['R'])['clashes'])
        score,clashes=cache[mut]; delta=score-wt_score; pred=predicted_direction(delta,config()['score_neutral_tolerance'])
        exp=truth['normalized_direction']
        rows.append(dict(source_id=truth['source_id'],mutation=mut,wild_type_score=wt_score,mutant_score=score,delta_score=delta,predicted_direction=pred,
                         experimental_direction=exp,concordant=str(pred==exp).lower() if exp in {'enhanced','reduced','approximately_neutral'} else 'not_classifiable',
                         assay_context=truth['assay_type'],construct=truth['cd81_construct'],strain=truth['hcv_e2_or_virus_strain'],
                         structural_notes=f'Hybrid model; repaired identically; clashes={clashes}',doi=truth['doi_or_pmid']))
    write_tsv(ROOT/'results/cd81_mutation_benchmark.tsv',rows)
    binding=[r for r in rows if 'entry' not in r['assay_context'] and r['experimental_direction'] in {'enhanced','reduced','approximately_neutral'}]
    unique={(r['mutation'],r['experimental_direction']):r for r in binding}
    directional=[r for r in unique.values() if r['experimental_direction']!='approximately_neutral']
    correct=sum(r['concordant']=='true' for r in directional)
    positive=any(r['mutation']=='T163A' and r['concordant']=='true' for r in directional)
    passed=bool(directional) and correct/len(directional)>0.5 and positive
    counts=Counter(r['experimental_direction'] for r in directional)
    ci=binomtest(correct,len(directional)).proportion_ci(method='exact') if directional else None
    write_json(ROOT/'results/benchmark_summary.json',dict(assay_rows=len(rows),unique_mutations=len(cache),
               classifiable_binding_controls=len(unique),directional_binding_controls=len(directional),concordant_directional=correct,
               positive_control_recovered=positive,benchmark_passed=passed,
               concordance_ci95=[ci.low,ci.high] if ci else [],
               concordance_interval_method='Clopper-Pearson exact; distinct mutation-direction controls; assay dependence remains',
               experimental_direction_counts=dict(counts),
               majority_direction_baseline_correct=max(counts.values(),default=0),majority_direction_baseline_total=len(directional),
               warning='Model outputs only; categorical validation does not establish quantitative affinity' if passed else 'Interface scoring failed directional/positive-control gate; no improved-binding interpretation allowed',
               quantitative_correlation='Not calculated: qualitative observations; incompatible assays not pooled'))
