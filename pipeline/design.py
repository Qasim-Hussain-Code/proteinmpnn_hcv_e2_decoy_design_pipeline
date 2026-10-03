"""Official upstream ProteinMPNN invocation and independent constraint checks."""
from __future__ import annotations
from collections import Counter
import copy
import json
import os
from pathlib import Path
import re
import sys
import time
import numpy as np
from Bio import SeqIO
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from .common import ROOT, check_resources, config, footprint, now, read_tsv, run_command, sha256, write_json, write_tsv
from .structures import interface, read_structure, save_pdb, sequence

def fixed_policy():
    model=read_structure(ROOT/'data/processed/humanized_unrepaired.pdb')
    receptor=model[0]['R']; e2=model[0]['E']
    distance=interface(receptor,e2)['min_distances']; rows=[]
    wt=sequence(receptor)
    for arm in ['A','B']:
        for i,residue in enumerate(receptor):
            position=residue.id[1]; d=distance[position]
            essential=position in config()['disulfide_positions'] or position in config()['indispensable_positions']
            # All positions with complete backbone can be redesigned, subject to
            # receptor-surface preservation in A and cautious reliable contacts in B.
            reliable=all(atom in residue for atom in ['N','CA','C','O'])
            fixed=essential or not reliable or (arm=='A' and d<=config()['interface_cutoff'])
            rows.append(dict(position=position,sequence_index=i+1,wild_type=wt[i],designable=str(not fixed).lower(),design_arm=arm,
                             reason='disulfide' if position in config()['disulfide_positions'] else 'experimental_recognition_constraint' if essential else
                                    'incomplete_backbone' if not reliable else 'Arm A preserves geometric interface' if fixed else 'predeclared sequence redesign',
                             source='Primary ground truth plus experimental coordinates',structural_distance_to_e2=d,
                             experimental_evidence='Conservative preservation across assay disagreement' if position in config()['indispensable_positions'] else '',
                             notes='Policy computed without natural-sequence conservation or held-out scores'))
    write_tsv(ROOT/'config/fixed_positions.tsv',rows)
    return rows

def constraints(sequence_new,wt,positions,fixed,required_cysteines,allow_extra=False):
    if len(sequence_new)!=len(wt):
        return 'invalid_length'
    if any(sequence_new[i]!=wt[i] for i in fixed):
        return 'fixed_position_violation'
    if any(sequence_new[positions.index(p)]!='C' for p in required_cysteines):
        return 'disulfide_violation'
    if not allow_extra and sequence_new.count('C')!=wt.count('C'):
        return 'unexpected_cysteine'
    if any(a not in 'ACDEFGHIKLMNPQRSTVWY' for a in sequence_new):
        return 'noncanonical_sequence'
    return 'passed'

def generate(arm,model,seed,temp,count,tag):
    folder=ROOT/f'data/work/mpnn/{tag}'
    fasta=folder/'seqs/humanized_unrepaired.fa'
    policy=read_tsv(ROOT/'config/fixed_positions.tsv')
    fixed=[int(r['sequence_index']) for r in policy if r['design_arm']==arm and r['designable']=='false']
    fixedfile=ROOT/f'data/processed/mpnn_fixed_{arm}.jsonl'
    fixedfile.write_text(json.dumps({'humanized_unrepaired':{'R':fixed,'E':[]}})+'\n')
    weights=ROOT/f'vendor/ProteinMPNN/{"vanilla" if model=="standard" else "soluble"}_model_weights'
    command=[sys.executable,ROOT/'vendor/ProteinMPNN/protein_mpnn_run.py','--pdb_path',ROOT/'data/processed/humanized_unrepaired.pdb',
             '--pdb_path_chains','R','--path_to_model_weights',weights,'--out_folder',folder,'--num_seq_per_target',count,
             '--sampling_temp',temp,'--seed',seed,'--batch_size',1,'--fixed_positions_jsonl',fixedfile,
             '--omit_AAs','CX','--model_name','v_48_020','--save_score',1]
    if model=='soluble':
        command.append('--use_soluble_model')
    if not fasta.exists():
        os.environ['CUDA_VISIBLE_DEVICES']=''
        os.environ['OMP_NUM_THREADS']=str(config()['threads'])
        os.environ['MKL_NUM_THREADS']=str(config()['threads'])
        command=[arg.as_posix() if isinstance(arg,Path) else arg for arg in command]
        output=run_command(command,timeout=7200)
        log=ROOT/f'logs/mpnn/{tag}.txt'; log.parent.mkdir(parents=True,exist_ok=True); log.write_text(output)
    records=list(SeqIO.parse(fasta,'fasta'))
    if not records or f'seed={seed}' not in records[0].description:
        raise RuntimeError('Effective ProteinMPNN seed not verified in official output header')
    if len(records)-1!=count:
        raise RuntimeError('Official MPNN sample count mismatch')
    wt=json.loads((ROOT/'data/processed/humanization.json').read_text())['human_sequence']
    positions=json.loads((ROOT/'data/processed/humanization.json').read_text())['human_positions']
    rows=[]
    scores=np.load(folder/'scores/humanized_unrepaired.npz')['score']
    # upstream output includes designed chains only for this one-chain design.
    for i,record in enumerate(records[1:]):
        seq=str(record.seq).split('/')[0]
        status=constraints(seq,wt,positions,[x-1 for x in fixed],config()['disulfide_positions'])
        if status in {'fixed_position_violation','disulfide_violation'}:
            raise RuntimeError(f'Critical MPNN enforcement failure: {status}')
        mutations=[f'{a}{p}{b}' for a,b,p in zip(wt,seq,positions) if a!=b]
        pa=ProteinAnalysis(seq)
        longest=max((len(s) for s in re.split('0',''.join('1' if a==b else '0' for a,b in zip(wt,seq)))),default=0)
        rows.append(dict(candidate_id='',sequence=seq,length=len(seq),design_arm=arm,ProteinMPNN_model=model,
                         seed=seed,effective_seed=seed,temperature=temp,generation_index=i+1,ProteinMPNN_score=float(scores[i]),
                         filter_status=status,sequence_identity_percent=100*sum(a==b for a,b in zip(wt,seq))/len(wt),mutation_count=len(mutations),
                         mutation_list=','.join(mutations),longest_unchanged_stretch=longest,net_charge_ph7=pa.charge_at_pH(7),pI=pa.isoelectric_point(),
                         hydrophobic_fraction=sum(a in 'AVILMFWY' for a in seq)/len(seq),sequence_sha256=__import__('hashlib').sha256(seq.encode()).hexdigest(),
                         upstream_fasta_sha256=sha256(fasta),generation_tag=tag))
    return rows

def pilot():
    if (ROOT/'results/mpnn_pilot.json').exists():
        return
    fixed_policy(); before=footprint(); t0=time.monotonic()
    rows=generate('A','standard',config()['generation_seeds'][0],0.1,config()['pilot_count'],'pilot')
    elapsed=time.monotonic()-t0; growth=footprint()-before
    write_tsv(ROOT/'results/pilot_sequences.tsv.gz',rows)
    # Budget is chosen from throughput before any generated candidate score.
    # One engineering cap: at most 20 minutes projected MPNN sampling.
    per_cell=max(1,min(25,int(1200/max(elapsed,0.1)*len(rows)/40)))
    projection=growth*per_cell*40/max(1,len(rows))+50_000_000
    check_resources(projection,len(rows),'full_count_per_cell')
    import yaml
    cfg=config(); cfg['full_count_per_cell']=per_cell
    (ROOT/'config/design.yml').write_text(yaml.safe_dump(cfg,sort_keys=False))
    write_json(ROOT/'results/mpnn_pilot.json',dict(count=len(rows),unique=len({r['sequence'] for r in rows}),elapsed_seconds=elapsed,
               disk_growth_bytes=growth,output_bytes_per_sequence=growth/max(1,len(rows)),chosen_count_per_cell=per_cell,
               projected_generation_count=per_cell*40,projected_elapsed_seconds=elapsed*per_cell*40/len(rows),
               budget_rule='Throughput-limited 20 minute MPNN generation engineering cap, 25 per model/arm/seed/temperature cell maximum. Chosen before scores. Pilot excluded from candidate pool.',
               model_loading_included=True,threads=config()['threads']))

def design():
    if (ROOT/'results/generated_sequences.tsv.gz').exists():
        return
    if config()['full_count_per_cell'] is None:
        raise RuntimeError('Pilot must precede full budget selection')
    pilotdata=json.loads((ROOT/'results/mpnn_pilot.json').read_text())
    check_resources(pilotdata['output_bytes_per_sequence']*pilotdata['projected_generation_count']+50_000_000,
                    pilotdata['count'],'full_count_per_cell')
    rows=[]
    for arm in ['A','B']:
        for model in ['standard','soluble']:
            for seed in config()['generation_seeds']:
                for temp in config()['temperatures']:
                    tag=f'{arm}_{model}_{seed}_{temp}'
                    rows.extend(generate(arm,model,seed,temp,config()['full_count_per_cell'],tag))
                    print(f'Generated {len(rows)} sequences',flush=True)
    seen=set()
    for i,row in enumerate(rows):
        row['candidate_id']=f'decoy_{i+1:04d}'
        row['duplicate']=str(row['sequence'] in seen).lower(); seen.add(row['sequence'])
    write_tsv(ROOT/'results/generated_sequences.tsv.gz',rows)
    counts=Counter(r['filter_status'] for r in rows)
    write_tsv(ROOT/'results/design_filter_funnel.tsv',[dict(transition='generated',count=len(rows)),
              dict(transition='valid_length',count=sum(r['filter_status']!='invalid_length' for r in rows)),
              dict(transition='fixed_positions_and_disulfides',count=sum(r['filter_status'] not in {'invalid_length','fixed_position_violation','disulfide_violation'} for r in rows)),
              dict(transition='canonical_and_no_extra_cysteine',count=counts['passed']),
              dict(transition='unique_passed',count=len({r['sequence'] for r in rows if r['filter_status']=='passed'}))])
    seeds=[]
    for arm in ['A','B']:
        for model in ['standard','soluble']:
            for seed in config()['generation_seeds']:
                cell=[r for r in rows if r['design_arm']==arm and r['ProteinMPNN_model']==model and r['seed']==seed]
                seeds.append(dict(design_arm=arm,model=model,seed=seed,count=len(cell),unique=len({r['sequence'] for r in cell}),
                                  mean_nll=np.mean([r['ProteinMPNN_score'] for r in cell]),mean_identity=np.mean([r['sequence_identity_percent'] for r in cell]),
                                  mean_hydrophobic_fraction=np.mean([r['hydrophobic_fraction'] for r in cell])))
    write_tsv(ROOT/'results/seed_variance.tsv',seeds)
