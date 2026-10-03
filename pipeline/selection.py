"""Discovery-explicit scoring; immutable candidate freeze; gated test evaluation."""
from __future__ import annotations
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time
from functools import lru_cache
import numpy as np
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1
from .common import ROOT, check_resources, config, failure, now, read_tsv, sha256, write_json, write_tsv
from .scoring import build_states, combine, energy, mutate
from .structures import interface, read_structure, save_pdb, sequence

def require_freeze(root=ROOT):
    table=Path(root)/'results/candidate_freeze.tsv'
    manifest=Path(root)/'results/candidate_freeze_manifest.json'
    if not table.exists() or not manifest.exists():
        raise RuntimeError('Leakage refusal: held-out data inaccessible before candidate freeze')
    data=json.loads(manifest.read_text())
    if sha256(table)!=data['candidate_table_sha256']:
        raise RuntimeError('Frozen candidate table changed: invalidate run explicitly')
    for name,key in [('config/design.yml','config_hash'),('results/software_manifest.tsv','software_manifest_hash')]:
        if sha256(Path(root)/name)!=data[key]:
            raise RuntimeError(f'Frozen input changed: {name}')
    for name,expected in data.get('method_file_hashes',{}).items():
        if sha256(Path(root)/name)!=expected:
            raise RuntimeError(f'Frozen scientific implementation changed: {name}')
    for name,expected in data.get('input_file_hashes',{}).items():
        if sha256(Path(root)/name)!=expected:
            raise RuntimeError(f'Frozen selection input changed: {name}')
    return data

def heldout_read(path,root=ROOT):
    require_freeze(root)
    return json.loads(Path(path).read_text())

def install_visibility_guard(root=ROOT):
    """Audit actual Python file opens, not just a README or filename convention."""
    def audit(event,args):
        if event!='open' or not isinstance(args[0],(str,bytes)):
            return
        path=str(args[0]).replace('\\','/').lower()
        mode=args[1]
        if any(term in path for term in ['quarantine/','heldout_sequences','held_out_panel','states/held_out/']):
            if mode is None or 'r' in str(mode):
                require_freeze(root)
    sys.addaudithook(audit)

def pareto_indices(array):
    a=np.asarray(array,float)
    return [i for i in range(len(a)) if not any(np.all(a[j]<=a[i]) and np.any(a[j]<a[i]) for j in range(len(a)) if j!=i)]

def summarize(scores,wt):
    scores=np.asarray(scores,float); wt=np.asarray(wt,float); delta=scores-wt
    return dict(median_score=float(np.median(scores)),worst_score=float(scores.max()),best_score=float(scores.min()),
                iqr=float(np.quantile(scores,.75)-np.quantile(scores,.25)),std=float(scores.std()),
                fraction_better_wt=float(np.mean(delta<0)),fraction_within_wt_tolerance=float(np.mean(delta<=config()['wt_tolerance'])),
                worst_delta_wt=float(delta.max()),median_delta_wt=float(np.median(delta)))

def score_one_task(task):
    candidate,states,split=task
    if split=='discovery':
        install_visibility_guard()
    os.environ['OMP_NUM_THREADS']='1'; os.environ['MKL_NUM_THREADS']='1'
    return score_pool([candidate],states,split,_serial=True)

@lru_cache(maxsize=64)
def cached_state_structure(path):
    return read_structure(path)

def score_pool(candidates,discovery_state_files,split='discovery',_serial=False):
    """Caller must pass explicit state records. No test files read in selection."""
    if split!='discovery':
        require_freeze()
    if split=='discovery' and any(r['split']!='discovery' for r in discovery_state_files):
        raise RuntimeError('Held-out state passed to discovery selection')
    if len(candidates)>2 and not _serial:
        from concurrent.futures import ProcessPoolExecutor
        output=[]; bios=[]
        checkpoint=ROOT/f'data/work/{split}_score_checkpoint.tsv.gz'
        with ProcessPoolExecutor(max_workers=config().get('scoring_workers',1)) as executor:
            tasks=[(candidate,discovery_state_files,split) for candidate in candidates]
            for n,(rows,bio) in enumerate(executor.map(score_one_task,tasks)):
                output.extend(rows); bios.extend(bio)
                if n%10==0:
                    print(f'Scored {n+1}/{len(candidates)} {split} candidates',flush=True)
                    write_tsv(checkpoint,output)
        return output,bios
    base=ROOT/'data/processed/humanized_repaired.pdb'
    wt=read_structure(base)[0]['R']; positions=[r.id[1] for r in wt]; wild=sequence(wt)
    work=ROOT/'data/work/scoring'/str(os.getpid()); work.mkdir(parents=True,exist_ok=True)
    output=[]; bios=[]
    receptor_models=ROOT/'data/work/receptor_models'; receptor_models.mkdir(parents=True,exist_ok=True)
    for n,candidate in enumerate(candidates):
        candidate_id=candidate['candidate_id']; target=receptor_models/(candidate_id+'.pdb')
        mutations=[('R',p,a,b) for p,a,b in zip(positions,wild,candidate['sequence']) if a!=b]
        # Identical reference complex packing for WT and each generated sequence;
        # viral state E2 sidechains retained independently. No state-specific
        # receptor repacking, explicitly a fixed-conformation approximation.
        try:
            if candidate_id=='wild_type':
                shutil.copyfile(base,target)
            else:
                mutate(base,mutations,target)
            structural=read_structure(target)
            monomer=structural.copy(); monomer[0].detach_child('E')
            mono=work/(candidate_id+'_monomer.pdb'); save_pdb(monomer,mono)
            stability=energy(mono,'ComputeStability')
            ShrakeRupley(n_points=100,probe_radius=1.4).compute(monomer,level='R')
            residues=list(monomer[0]['R']); total=sum(r.sasa for r in residues)
            hydrophobic=sum(r.sasa for r in residues if seq1(r.resname) in 'AVILMFWY')
            patch=max((sum(r.sasa for r in residues if seq1(r.resname) in 'AVILMFWY' and np.linalg.norm(r['CA'].coord-center['CA'].coord)<8)
                       for center in residues),default=0)
            bios.append(dict(candidate_id=candidate_id,stability_score=stability,sasa_angstrom2=total,
                             hydrophobic_sasa_fraction=hydrophobic/max(total,1),hydrophobic_patch_proxy_angstrom2=patch,
                             preparation='Reference-complex BuildMutant + RepairStructure; E2 conformation then recombined per state'))
            mono.unlink()
            for state in discovery_state_files:
                combined=work/'score_complex.pdb'; combine(target,ROOT/state['path'],combined)
                audit=interface(cached_state_structure(str(ROOT/state['path']))[0]['E'],structural[0]['R'])
                score=energy(combined)
                output.append(dict(candidate_id=candidate_id,haplotype_id=state['haplotype_id'],split=split,interaction_score=score,
                                   genotype=state['genotype'],number_of_sequences=state['number_of_sequences'],clash_count=audit['clashes'],
                                   contact_count=audit['atom_contacts'],score_status='passed',failure=''))
        except Exception as exc:
            failure(f'{split}_score',candidate_id,exc)
            for state in discovery_state_files:
                output.append(dict(candidate_id=candidate_id,haplotype_id=state['haplotype_id'],split=split,interaction_score='',
                                   genotype=state['genotype'],number_of_sequences=state['number_of_sequences'],clash_count='',contact_count='',score_status='failed',failure=str(exc)))
        if n%10==0 and not _serial:
            print(f'Scored {n+1}/{len(candidates)} {split} candidates',flush=True)
    return output,bios

def wt_candidate():
    return dict(candidate_id='wild_type',sequence=sequence(read_structure(ROOT/'data/processed/humanized_repaired.pdb')[0]['R']))

def discovery():
    install_visibility_guard()
    output=ROOT/'results/discovery_scores.tsv.gz'
    if output.exists():
        cached=read_tsv(output)
        pool=[r for r in read_tsv(ROOT/'results/generated_sequences.tsv.gz') if r['filter_status']=='passed' and r['duplicate']=='false']
        states=[r for r in read_tsv(ROOT/'results/e2_state_manifest.tsv') if r['split']=='discovery']
        keys={(r['candidate_id'],r['haplotype_id']) for r in cached}
        expected={(c['candidate_id'],s['haplotype_id']) for c in [wt_candidate()]+pool for s in states}
        if len(cached)!=len(expected) or keys!=expected or not (ROOT/'results/candidate_biophysics.tsv').exists() or not (ROOT/'results/model_comparison.tsv').exists():
            raise RuntimeError('Incomplete discovery artifacts; refuse silent cache skip')
        return
    pool=[r for r in read_tsv(ROOT/'results/generated_sequences.tsv.gz') if r['filter_status']=='passed' and r['duplicate']=='false']
    states=[r for r in read_tsv(ROOT/'results/e2_state_manifest.tsv') if r['split']=='discovery']
    from .common import MeasuredStage,footprint
    before=footprint(); t0=time.monotonic()
    with MeasuredStage('scoring_pilot','WT + first unique candidate on reference state'):
        pilot,bio=score_pool([wt_candidate(),pool[0]],states[:1])
    pilot_elapsed=time.monotonic()-t0; pilot_growth=max(1,footprint()-before)
    if any(r['score_status']=='failed' for r in pilot):
        raise RuntimeError('Scoring pilot failed; refuse full discovery operation')
    projection=pilot_growth*(len(pool)+1)/2+len(pool)*len(states)*2000
    check_resources(projection,2,'candidate_count')
    write_json(ROOT/'results/scoring_projection.json',dict(pilot_candidates=2,pilot_states=1,full_candidates=len(pool)+1,
               full_states=len(states),pilot_elapsed_seconds=pilot_elapsed,pilot_growth_bytes=pilot_growth,
               projected_disk_bytes=projection,safety_margin=1.2,
               note='Receptor-only reusable models, scored complexes disposable; no state-specific receptor repacking'))
    rows,bios=score_pool([wt_candidate()]+pool,states)
    write_tsv(output,rows); write_tsv(ROOT/'results/candidate_biophysics.tsv',bios)
    bmap={r['candidate_id']:r for r in bios}; model_summary=[]
    for arm in ['A','B']:
        for model in ['standard','soluble']:
            for seed in config()['generation_seeds']:
                cell=[r for r in pool if r['design_arm']==arm and r['ProteinMPNN_model']==model and int(r['seed'])==seed and r['candidate_id'] in bmap]
                if cell:
                    model_summary.append(dict(design_arm=arm,model=model,seed=seed,count=len(cell),
                         mean_hydrophobic_sasa_fraction=np.mean([float(bmap[r['candidate_id']]['hydrophobic_sasa_fraction']) for r in cell]),
                         mean_hydrophobic_patch_proxy=np.mean([float(bmap[r['candidate_id']]['hydrophobic_patch_proxy_angstrom2']) for r in cell]),
                         mean_stability_score=np.mean([float(bmap[r['candidate_id']]['stability_score']) for r in cell]),
                         mean_charge=np.mean([float(r['net_charge_ph7']) for r in cell]),mean_pI=np.mean([float(r['pI']) for r in cell])))
    write_tsv(ROOT/'results/model_comparison.tsv',model_summary)

def freeze_ids(metrics,bio,candidate_rows,count=5):
    eligible=[r for r in metrics if r['candidate_id'] in bio]
    front=pareto_indices([[float(r['worst_delta_wt']),float(bio[r['candidate_id']]['hydrophobic_sasa_fraction'])] for r in eligible])
    pareto={eligible[i]['candidate_id'] for i in front}
    single=sorted(eligible,key=lambda r:(float(r['reference_score']),r['candidate_id']))[:count]
    escape=sorted(eligible,key=lambda r:(float(r['worst_delta_wt']),float(r['median_delta_wt']),r['candidate_id']))[:count]
    selected=defaultdict(list)
    for label,items in [('single_state',single),('escape_aware',escape)]:
        for r in items:
            selected[r['candidate_id']].append(label)
    return selected,pareto

def freeze():
    install_visibility_guard()
    if (ROOT/'results/candidate_freeze.tsv').exists():
        require_freeze(); return
    rows=read_tsv(ROOT/'results/discovery_scores.tsv.gz')
    wt={r['haplotype_id']:float(r['interaction_score']) for r in rows if r['candidate_id']=='wild_type' and r['score_status']=='passed'}
    grouped=defaultdict(list)
    for r in rows:
        if r['candidate_id']!='wild_type':
            grouped[r['candidate_id']].append(r)
    metrics=[]
    for cid,items in sorted(grouped.items()):
        if any(r['score_status']!='passed' for r in items) or len(items)!=len(wt):
            continue
        natural=[r for r in items if r['haplotype_id']!='reference']
        m=summarize([float(r['interaction_score']) for r in natural],[wt[r['haplotype_id']] for r in natural])
        m.update(candidate_id=cid,reference_score=next(float(r['interaction_score']) for r in items if r['haplotype_id']=='reference'),
                 max_clashes=max(int(r['clash_count']) for r in items),selection_status='scoreable',
                 genotype_stratified={g:summarize([float(r['interaction_score']) for r in natural if g in r['genotype'].split(',')],
                       [wt[r['haplotype_id']] for r in natural if g in r['genotype'].split(',')]) for g in sorted({g for r in natural for g in r['genotype'].split(',')})})
        metrics.append(m)
    write_tsv(ROOT/'results/discovery_metrics.tsv',metrics)
    bio={r['candidate_id']:r for r in read_tsv(ROOT/'results/candidate_biophysics.tsv')}
    candidates={r['candidate_id']:r for r in read_tsv(ROOT/'results/generated_sequences.tsv.gz')}
    wt_clashes={r['haplotype_id']:int(r['clash_count']) for r in rows if r['candidate_id']=='wild_type' and r['score_status']=='passed'}
    eligible_ids={cid for cid,items in grouped.items() if all(r['score_status']=='passed' and int(r['clash_count'])<=wt_clashes[r['haplotype_id']]+config()['max_excess_clashes'] for r in items)}
    eligible=[m for m in metrics if m['candidate_id'] in eligible_ids]
    selected,pareto=freeze_ids(eligible,bio,candidates,config()['freeze_per_strategy'])
    frozen=[]; timestamp=now(); cfg_hash=sha256(ROOT/'config/design.yml'); sw_hash=sha256(ROOT/'results/software_manifest.tsv')
    for m in metrics:
        cid=m['candidate_id']
        if cid in selected:
            frozen.append(dict(**candidates[cid],selection_rule=','.join(selected[cid]),discovery_metrics=m,biophysics=bio[cid],pareto_eligible=cid in pareto,
                               selection_timestamp=timestamp,config_hash=cfg_hash,software_manifest_hash=sw_hash))
    write_tsv(ROOT/'results/candidate_freeze.tsv',frozen)
    write_json(ROOT/'results/candidate_freeze_manifest.json',dict(candidate_table_sha256=sha256(ROOT/'results/candidate_freeze.tsv'),
               candidate_ids=sorted(selected),selection_timestamp=timestamp,config_hash=cfg_hash,software_manifest_hash=sw_hash,
               discovery_inputs={'scores':sha256(ROOT/'results/discovery_scores.tsv.gz'),'states':sha256(ROOT/'data/processed/discovery_panel.tsv')},
               method_file_hashes={f'pipeline/{name}.py':sha256(ROOT/f'pipeline/{name}.py') for name in ['common','structures','diversity','scoring','design','selection']},
               input_file_hashes={name:sha256(ROOT/name) for name in ['results/discovery_scores.tsv.gz','data/processed/discovery_panel.tsv',
                   'results/generated_sequences.tsv.gz','results/candidate_biophysics.tsv','config/fixed_positions.tsv','config/split_manifest.tsv',
                   'config/cd81_mutation_ground_truth.tsv','results/benchmark_summary.json','data/processed/humanized_repaired.pdb']},
               rule='Same unique candidate pool: reference interaction-score ascending; escape maximum paired WT delta ascending, then median paired delta, then neutral identifier. Pareto annotated separately; not an extra selection gate.',
               heldout_opened=False))
    funnel=read_tsv(ROOT/'results/design_filter_funnel.tsv')
    funnel.extend([dict(transition='scoreable',count=len(metrics)),dict(transition='no_severe_excess_clash',count=len(eligible)),dict(transition='pareto_eligible',count=len(pareto)),
                   dict(transition='discovery_selected_and_frozen',count=len(frozen))])
    write_tsv(ROOT/'results/design_filter_funnel.tsv',funnel)

def heldout():
    require_freeze()
    if (ROOT/'results/heldout_metrics.tsv').exists():
        expected=['heldout_scores.tsv','heldout_summary.json','final_computational_priorities.tsv']
        if not all((ROOT/'results'/name).exists() for name in expected):
            raise RuntimeError('Incomplete held-out artifacts; refuse silent cache skip')
        completed=read_tsv(ROOT/'results/heldout_metrics.tsv')
        if {r['candidate_id'] for r in completed}!={r['candidate_id'] for r in read_tsv(ROOT/'results/candidate_freeze.tsv')}:
            raise RuntimeError('Incomplete held-out candidate evaluation')
        return
    from .diversity import panel
    rows=heldout_read(ROOT/'data/processed/quarantine/heldout_sequences.json')
    if not rows:
        raise RuntimeError('No independent held-out sequences after cluster split')
    positions=json.loads((ROOT/'results/diversity_summary.json').read_text())['interface_positions']
    hpanel=panel(rows,positions,config()['coverage_target'],config()['max_heldout_states'])
    write_tsv(ROOT/'data/processed/held_out_panel.tsv',hpanel)
    build_states('held_out')
    states=[r for r in read_tsv(ROOT/'results/e2_state_manifest.tsv') if r['split']=='held_out' and r['haplotype_id']!='reference']
    frozen=read_tsv(ROOT/'results/candidate_freeze.tsv')
    scores,bios=score_pool([wt_candidate()]+frozen,states,'held_out')
    wt={r['haplotype_id']:float(r['interaction_score']) for r in scores if r['candidate_id']=='wild_type' and r['score_status']=='passed'}
    for row in scores:
        row['paired_delta_wt']=float(row['interaction_score'])-wt[row['haplotype_id']] if row['score_status']=='passed' else ''
    write_tsv(ROOT/'results/heldout_scores.tsv',scores)
    metrics=[]
    for candidate in frozen:
        items=[r for r in scores if r['candidate_id']==candidate['candidate_id']]
        if any(r['score_status']=='failed' for r in items):
            continue
        m=summarize([float(r['interaction_score']) for r in items],[wt[r['haplotype_id']] for r in items])
        genotype={g:summarize([float(r['interaction_score']) for r in items if g in r['genotype'].split(',')],
                             [wt[r['haplotype_id']] for r in items if g in r['genotype'].split(',')]) for g in sorted({g for r in items for g in r['genotype'].split(',')})}
        metrics.append(dict(candidate_id=candidate['candidate_id'],selection_rule=candidate['selection_rule'],**m,genotype_stratified=genotype,
                            failure_count=0,number_of_states=len(items)))
    write_tsv(ROOT/'results/heldout_metrics.tsv',metrics)
    # Equal-state comparison; state bootstrap, preserve shared candidates as paired.
    single=[c['candidate_id'] for c in frozen if 'single_state' in c['selection_rule']]
    escape=[c['candidate_id'] for c in frozen if 'escape_aware' in c['selection_rule']]
    stateids=sorted(wt)
    lookup={(r['candidate_id'],r['haplotype_id']):float(r['paired_delta_wt']) for r in scores if r['score_status']=='passed'}
    matrix_single=np.array([[lookup[c,s] for s in stateids] for c in single]); matrix_escape=np.array([[lookup[c,s] for s in stateids] for c in escape])
    effect=float(np.mean(matrix_escape.max(axis=1))-np.mean(matrix_single.max(axis=1)))
    rng=np.random.default_rng(config()['bootstrap_seed']); boots=[]
    for i in range(config()['bootstrap_replicates']):
        idx=rng.integers(0,len(stateids),len(stateids))
        boots.append(float(np.mean(matrix_escape[:,idx].max(axis=1))-np.mean(matrix_single[:,idx].max(axis=1))))
    from scipy.stats import spearmanr
    dm={r['candidate_id']:r for r in read_tsv(ROOT/'results/discovery_metrics.tsv')}
    rank=spearmanr([float(dm[r['candidate_id']]['worst_delta_wt']) for r in metrics],[float(r['worst_delta_wt']) for r in metrics])
    write_json(ROOT/'results/heldout_summary.json',dict(states=len(stateids),sequences=len(rows),frozen=len(frozen),single_state_candidates=len(single),escape_aware_candidates=len(escape),
               shared_candidates=len(set(single)&set(escape)),worst_delta_effect_escape_minus_single=effect,ci95=np.quantile(boots,[.025,.975]).tolist(),
               bootstrap_unit='held-out structural haplotype state, equally weighted; correlated components remain a limitation',bootstrap_replicates=config()['bootstrap_replicates'],bootstrap_seed=config()['bootstrap_seed'],
               bootstrap_method='paired percentile state bootstrap of mean candidate maximum paired WT delta',rank_spearman=float(rank.statistic),
               no_experimental_binding_measurements=True,benchmark_warning=json.loads((ROOT/'results/benchmark_summary.json').read_text())['warning']))
    # A computational priority list is an annotation on the immutable frozen set.
    write_tsv(ROOT/'results/final_computational_priorities.tsv',[dict(**c,heldout_metrics=next(r for r in metrics if r['candidate_id']==c['candidate_id']),
               warning=json.loads((ROOT/'results/benchmark_summary.json').read_text())['warning']) for c in frozen])
    funnel=read_tsv(ROOT/'results/design_filter_funnel.tsv'); funnel.append(dict(transition='heldout_evaluated',count=len(metrics)))
    write_tsv(ROOT/'results/design_filter_funnel.tsv',funnel)
    final=ROOT/'data/processed/frozen_candidates'; final.mkdir(parents=True,exist_ok=True)
    for c in frozen:
        shutil.copyfile(ROOT/f'data/work/receptor_models/{c["candidate_id"]}.pdb',final/(c['candidate_id']+'.pdb'))

def sensitivity():
    require_freeze()
    from .diversity import panel
    discovery=json.loads((ROOT/'data/processed/discovery_sequences.json').read_text())
    positions=json.loads((ROOT/'results/diversity_summary.json').read_text())['interface_positions']
    rows=[]
    for target in config()['coverage_sensitivity']:
        selected=panel(discovery,positions,target,config()['max_discovery_states'])
        rows.append(dict(parameter='discovery_coverage_target',value=target,number_states=len(selected),achieved_coverage=selected[0]['achieved_coverage'] if selected else 0,
                         interpretation='Coverage selection only; candidates remain frozen; no additional structural score inferred'))
    for cutoff in config()['interface_sensitivity']:
        count=sum(float(r['cutoff_angstrom'])==cutoff for r in read_tsv(ROOT/'results/e2_interface_residues.tsv'))
        rows.append(dict(parameter='geometric_interface_cutoff',value=cutoff,number_states=count,achieved_coverage='',interpretation='Residue count sensitivity; no retroactive redesign'))
    for tolerance in [0.0,0.5,1.0]:
        b=read_tsv(ROOT/'results/cd81_mutation_benchmark.tsv')
        from .scoring import predicted_direction
        distinct={r['mutation']:r for r in b if 'entry' not in r['assay_context'] and r['experimental_direction'] in {'enhanced','reduced'}}
        concordance=sum(predicted_direction(float(r['delta_score']),tolerance)==r['experimental_direction'] for r in distinct.values())
        rows.append(dict(parameter='benchmark_neutral_tolerance',value=tolerance,number_states=concordance,achieved_coverage=len(distinct),interpretation='Categorical benchmark sensitivity; no threshold tuning'))
    write_tsv(ROOT/'results/sensitivity_analysis.tsv',rows)
