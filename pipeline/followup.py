"""Isolated, reproducible follow-up. Never edits the frozen experiment."""
from __future__ import annotations
import argparse
import copy
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from urllib.request import Request, urlopen
import numpy as np
from Bio.PDB import MMCIFParser, Superimposer
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.SeqUtils import seq1, seq3
from scipy.spatial.distance import cdist
from .common import ROOT, now, read_tsv, sha256, write_json, write_tsv, check_resources
from .structures import entity_metadata, interface, map_alignment, read_structure, save_pdb, standard_residues, sequence, heavy
from .scoring import evo_binary, combine, predicted_direction

OUT = ROOT / 'followup'
RAW = ROOT / 'data/raw/followup'
WORK = ROOT / 'data/work/followup'

def protocol():
    return json.loads((OUT/'protocol.json').read_text())

def snapshot():
    """Hash all original tracked science and outputs, including the old README."""
    path=OUT/'original_snapshot.tsv'
    if path.exists():
        rows=read_tsv(path)
        bad=[r['path'] for r in rows if not (ROOT/r['path']).exists() or sha256(ROOT/r['path'])!=r['sha256']]
        if bad:
            raise RuntimeError(f'Original archive changed: {bad}')
        return rows
    files=subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
    rows=[dict(path=p,sha256=sha256(ROOT/p)) for p in files if not p.startswith(('followup/','pipeline/followup','tests/test_followup')) and (ROOT/p).is_file()]
    write_tsv(path,rows)
    return rows

def lock():
    snapshot()
    manifests=[(ROOT/'results/download_manifest.tsv','resource'),(OUT/'downloads.tsv','path')]
    for file,key in manifests:
        if file.exists():
            for row in read_tsv(file):
                if (file.parent==OUT or row[key].replace('\\','/') in {'data/raw/7MWX.cif','data/raw/3X0E.cif','data/raw/P60033.fasta','data/raw/H77.gb'}) and sha256(ROOT/row[key])!=row['sha256']:
                    raise RuntimeError(f'Source snapshot changed: {row[key]}')
    hashes={p.name:sha256(p) for p in [OUT/'protocol.json',Path(__file__)]}
    file=OUT/'protocol_lock.json'
    if file.exists():
        recorded=json.loads(file.read_text())
        if recorded['protocol_sha256']!=hashes['protocol.json']:
            raise RuntimeError('Follow-up protocol changed after lock')
    else:
        write_json(file,dict(timestamp=now(),protocol_sha256=hashes['protocol.json'],initial_code_sha256=hashes['followup.py'],labels_consulted=True,interpretation='Exploratory diagnostic analysis, not a blinded validation'))

def fetch(url,destination,identifier):
    destination=Path(destination); destination.parent.mkdir(parents=True,exist_ok=True)
    manifest=OUT/'downloads.tsv'
    rows=read_tsv(manifest) if manifest.exists() else []
    name=str(destination.relative_to(ROOT))
    prior=next((r for r in rows if r['path']==name),None)
    if destination.exists() and prior and sha256(destination)!=prior['sha256']:
        raise RuntimeError(f'Source hash mismatch {name}')
    if not destination.exists():
        check_resources(25_000_000)
        temporary=destination.with_suffix(destination.suffix+'.tmp')
        for attempt in range(4):
            try:
                request=Request(url,headers={'User-Agent':'CD81Followup/0.1 (reproducible academic analysis)'})
                with urlopen(request,timeout=45) as response, temporary.open('wb') as stream:
                    shutil.copyfileobj(response,stream)
                if not temporary.stat().st_size:
                    raise RuntimeError('Empty download')
                temporary.replace(destination); break
            except Exception:
                temporary.unlink(missing_ok=True)
                if attempt==3: raise
                time.sleep(2*(attempt+1))
    if not prior:
        rows.append(dict(path=name,identifier=identifier,source_url=url,retrieved_utc=now(),bytes=destination.stat().st_size,sha256=sha256(destination)))
        write_tsv(manifest,rows)
    return destination

def prepare():
    lock()
    for identifier in ['1G8Q','5TCX']:
        fetch(f'https://files.rcsb.org/download/{identifier}.cif',RAW/f'{identifier}.cif',identifier)
    fetch('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=15670777&retmode=xml',RAW/'PMID15670777.xml','PMID15670777')
    curate()

def curate():
    rows=[]; audit=[]
    for i,original in enumerate(read_tsv(ROOT/'config/cd81_mutation_ground_truth.tsv'),1):
        row=dict(original,original_row_id=i,followup_role='original_evidence')
        if row['source_id']=='higginbottom2000' and row['cd81_construct']=='GST human CD81 LEL' and row['mutation'] in {'F186L','E188K'}:
            row['included']='false'; row['exclusion_reason']='Soluble construct was F186L+E188K, not the reported single mutation'
            row['followup_role']='curation_correction_excluded'
            audit.append(dict(original_row_id=i,mutation=row['mutation'],action='exclude misassigned soluble single mutant',evidence='Methods, Cloning and expression of CD81 LEL fusion proteins; Results, Fig 4',source=row['primary_source']))
        if row['source_id']=='higginbottom2000' and row['mutation']=='D196E':
            row['notes']+='; abstract reports no effect whereas Results describe weaker/intermediate binding; direction follows assay-specific Results and remains uncertain'
        rows.append(row)
    template=next(r for r in rows if r['source_id']=='higginbottom2000' and r['mutation']=='F186L' and 'GST' in r['cd81_construct'])
    rows.append(dict(template,mutation='F186L+E188K',wild_type_residue='F,E',position_full_length_cd81='186,188',position_structure='186,188',included='true',exclusion_reason='',normalized_direction='reduced',reported_direction='reduced',measurement='double mutant has no E2 binding signal',original_row_id='',followup_role='corrected_double_mutant',notes='Single cell-surface F186L and E188K observations retained separately; no double-mutant effect attributed to either single substitution',evidence_anchor='Methods soluble GST constructs; Results and Fig 4'))
    for mutation,direction in [('K124T','binding_retained_unquantified'),('V146E','binding_retained_unquantified'),('F150S','reduced'),('T166I','reduced'),('C157S','reduced'),('C190R','reduced')]:
        rows.append(dict(source_id='drummer2005',citation='Drummer, Wilson and Poumbourios, BBRC 328:251-257 (2005)',doi_or_pmid='10.1016/j.bbrc.2004.12.160; PMID15670777',cd81_construct='recombinant human CD81 LEL; exact boundaries not established from abstract',assay_type='recombinant LEL E2 binding; abstract-level evidence',hcv_e2_or_virus_strain='not specified in retrieved abstract',mutation=mutation,wild_type_residue=mutation[0],position_full_length_cd81=mutation[1:-1],position_structure=mutation[1:-1],measurement='retained binding' if direction.startswith('binding_retained') else 'loss of E2 binding',measurement_units='qualitative',reported_direction=direction,normalized_direction=direction,replicate_information='not specified in abstract',notes='All six mutations abolished recombinant LEL dimerization; folding/dimerization effects are not captured by fixed-backbone interface scoring; retained does not establish WT-equivalent affinity',included='true',exclusion_reason='',primary_source='https://pubmed.ncbi.nlm.nih.gov/15670777/',evidence_anchor='PubMed abstract; soluble and full-length effects are not pooled',original_row_id='',followup_role='external_source_challenge'))
    write_tsv(OUT/'ground_truth.tsv',rows)
    write_tsv(OUT/'curation_audit.tsv',audit)

def invoke(command,pdb,tag,extra=()):
    folder=WORK/'evo'/str(os.getpid()); folder.mkdir(parents=True,exist_ok=True)
    target=folder/'m.pdb'; shutil.copyfile(pdb,target)
    # Resolve junctions before launching: upstream constructs library paths in
    # short fixed buffers, so a long alias can corrupt its own data filenames.
    binary=evo_binary().resolve()
    args=[str(binary),f'--command={command}','--pdb=m.pdb',*extra]
    result=subprocess.run(args,cwd=folder,capture_output=True,text=True,timeout=240)
    log=OUT/'logs'/f'{tag}_{command}.txt'; log.parent.mkdir(parents=True,exist_ok=True)
    log.write_text(result.stdout+'\nSTDERR\n'+result.stderr,encoding='utf-8')
    if result.returncode or not result.stdout.strip():
        raise RuntimeError(f'EvoEF2 failed: {tag}; {log}')
    if command.startswith('Compute') and 'invalid' in result.stdout.lower():
        raise RuntimeError(f'Incomplete model refused: {log}')
    write_json(log.with_suffix('.json'),dict(command=args,input_sha256=sha256(pdb),binary_sha256=sha256(binary),timestamp=now(),exit_status=result.returncode))
    return result.stdout,folder

def repair(pdb,destination,tag):
    if not Path(destination).exists():
        (WORK/'evo'/str(os.getpid())/'m_Repair.pdb').unlink(missing_ok=True)
        _,folder=invoke('RepairStructure',pdb,tag)
        generated=folder/'m_Repair.pdb'
        if not generated.exists(): raise RuntimeError('RepairStructure output missing')
        shutil.copyfile(generated,destination)
    return Path(destination)

def mutations(text):
    parsed=[]
    for item in text.split('+'):
        match=re.fullmatch(r'([A-Z])(\d+)([A-Z])',item)
        if not match: raise ValueError(f'Invalid mutation {text}')
        a,pos,b=match.groups(); parsed.append(('R',int(pos),a,b))
    return parsed

def mutate(pdb,changes,destination,tag):
    model=read_structure(pdb)
    for chain,pos,wt,new in changes:
        if seq1(model[0][chain][(' ',pos,' ')].resname)!=wt:
            raise ValueError(f'Wrong WT identity {chain}:{pos}')
    if not destination.exists():
        folder=WORK/'evo'/str(os.getpid()); folder.mkdir(parents=True,exist_ok=True)
        (folder/'individual_list.txt').write_text(','.join(f'{wt}{chain}{pos}{new}' for chain,pos,wt,new in changes)+';\n')
        (folder/'m_Model_0001.pdb').unlink(missing_ok=True)
        _,folder=invoke('BuildMutant',pdb,tag,['--mutant_file=individual_list.txt'])
        built=folder/'m_Model_0001.pdb'
        if not built.exists(): raise RuntimeError('BuildMutant output missing')
        repair(built,destination,tag+'_repair')
    model=read_structure(destination)
    for chain,pos,wt,new in changes:
        if seq1(model[0][chain][(' ',pos,' ')].resname)!=new:
            raise RuntimeError(f'Mutation identity failure {tag}')
    return destination

def energy(pdb,tag):
    log=OUT/'logs'/f'{tag}_ComputeBinding.txt'
    if log.exists():
        output=log.read_text()
        provenance=log.with_suffix('.json')
        if provenance.exists() and json.loads(provenance.read_text())['input_sha256']!=sha256(pdb):
            raise RuntimeError(f'Score cache input mismatch {tag}')
    else:
        output,_=invoke('ComputeBinding',pdb,tag,['--split=E,R'])
    if 'invalid' in output.lower():raise RuntimeError(f'Incomplete model score refused {tag}')
    if len(re.findall(r'^\s*Total\s*=',output,re.M))!=1:raise RuntimeError(f'Ambiguous total {tag}')
    terms={k:float(v) for k,v in re.findall(r'^\s*(\w+)\s*=\s*([-+\d.eE]+)',output,re.M)}
    if 'Total' not in terms: raise RuntimeError(f'Score parse failed {tag}')
    return terms

def human_reference():
    return ''.join((ROOT/'data/raw/P60033.fasta').read_text().splitlines()[1:])

def numbered(chain,require_human=False):
    residues=standard_residues(chain)
    mapping,_=map_alignment(human_reference(),sequence(chain))
    out={mapping[i]:r for i,r in enumerate(residues) if i in mapping and 113<=mapping[i]<=201}
    if require_human and any(seq1(r.resname)!=human_reference()[p-1] for p,r in out.items()):
        raise RuntimeError('Experimental human sequence mismatch; do not silently humanize it')
    return out

def make_models():
    lock(); check_resources(protocol()['maximum_new_model_bytes'])
    complex_model=MMCIFParser(QUIET=True,auth_chains=False).get_structure('bound',str(ROOT/'data/raw/7MWX.cif'))
    rows=[]; source_rows=[]; positions=[]
    for spec in protocol()['model_specs']:
        ident=spec['id']; folder=OUT/'models'/ident; folder.mkdir(parents=True,exist_ok=True)
        native=folder/'native.pdb'; prepared=folder/'wt.pdb'
        if spec.get('kind')=='archive':
            shutil.copyfile(ROOT/'data/processed/humanized_unrepaired.pdb',native)
            shutil.copyfile(ROOT/'data/processed/humanized_repaired.pdb',prepared)
            rmsd=float(json.loads((ROOT/'data/processed/humanization.json').read_text())['rmsd']); nfit=81; label='experimental human backbone, archival all-CA rigid fit'
        else:
            echain=complex_model[0][spec.get('e2','A')]; target=numbered(complex_model[0][spec.get('cd81','E')])
            model=Structure(ident); m=Model(0); model.add(m); e=Chain('E'); r=Chain('R'); m.add(e);m.add(r)
            for residue in standard_residues(echain): e.add(residue.copy())
            if spec.get('kind')=='bound_template':
                for p,residue in target.items():
                    new=residue.copy(); new.id=(' ',p,' '); new.resname=seq3(human_reference()[p-1]).upper()
                    # Restore the source identities below, then explicitly mutate their sidechains.
                    r.add(new)
                rmsd=0.;nfit=0;label='tamarin bound backbone with human sequence; diagnostic, not experimental human coordinates'
                for p,residue in target.items():
                    wt=seq1(residue.resname); new=human_reference()[p-1]
                    if wt!=new: r[(' ',p,' ')].resname=residue.resname
                save_pdb(model,native)
                changes=[('R',p,seq1(residue.resname),human_reference()[p-1]) for p,residue in target.items() if seq1(residue.resname)!=human_reference()[p-1]]
                mutated=folder/'humanized.pdb'; mutate(native,changes,mutated,ident+'_humanize')
                shutil.copyfile(mutated,prepared)
                source_rows.append(dict(model=ident,pdb_id='7MWX',label_chain=spec['cd81'],organism='Saguinus oedipus',modified_positions=[p for c,p,a,b in changes],notes=label))
            else:
                pdb=spec['pdb']; path=(ROOT/'data/raw/3X0E.cif') if pdb=='3X0E' else RAW/f'{pdb}.cif'
                meta=entity_metadata(MMCIF2Dict(str(path)))
                entity=next((v for v in meta.values() if spec['chain'] in v['chains']),None)
                if not entity or entity['organism']!='Homo sapiens' or 'CD81' not in entity['description'].upper():
                    raise RuntimeError(f'Species/entity verification failure {spec}')
                human=MMCIFParser(QUIET=True,auth_chains=False).get_structure(pdb,str(path))[0][spec['chain']]
                source=numbered(human,True); common=sorted(set(source)&set(target))
                if spec['fit']=='core': common=[p for p in common if any(a<=p<=b for a,b in protocol()['core_fit_human_positions'])]
                sup=Superimposer();sup.set_atoms([target[p]['CA'] for p in common],[source[p]['CA'] for p in common]);rmsd=float(sup.rms);nfit=len(common)
                for p,residue in source.items():
                    new=residue.copy();new.id=(' ',p,' ');sup.apply(list(new.get_atoms()));r.add(new)
                label='experimental human backbone, '+spec['fit']+' rigid CA fit'
                save_pdb(model,native); repair(native,prepared,ident+'_wt')
                source_rows.append(dict(model=ident,pdb_id=pdb,label_chain=spec['chain'],organism=entity['organism'],modified_positions=[],notes=label))
        structure=read_structure(prepared); pre=interface(read_structure(native)[0]['E'],read_structure(native)[0]['R']); post=interface(structure[0]['E'],structure[0]['R'])
        r=structure[0]['R']; distances={str(a)+'-'+str(b):float(np.linalg.norm(r[(' ',a,' ')]['SG'].coord-r[(' ',b,' ')]['SG'].coord)) for a,b in [(156,190),(157,175)]}
        row=dict(model=ident,description=label,ca_rmsd=rmsd,fit_atoms=nfit,observed_receptor_residues=len(standard_residues(r)),missing_positions=[p for p in range(113,202) if (' ',p,' ') not in r],contacts_before=pre['atom_contacts'],contacts_after=post['atom_contacts'],clashes_before=pre['clashes'],clashes_after=post['clashes'],disulfide_sg_distances=distances,wt_sha256=sha256(prepared),path=str(prepared.relative_to(ROOT)))
        rows.append(row)
        for p in [124,146,150,157,162,163,166,171,181,182,184,186,188,190,196]:
            if (' ',p,' ') in r:
                aa=heavy([r[(' ',p,' ')]]); bb=heavy(structure[0]['E']); ds=cdist([a.coord for a in aa],[b.coord for b in bb]);i,j=np.unravel_index(ds.argmin(),ds.shape)
                positions.append(dict(model=ident,cd81_position=p,min_distance=float(ds[i,j]),cd81_atom=aa[i].name,e2_structure_position=bb[j].parent.id[1],e2_atom=bb[j].name))
        print('prepared',ident,'clashes',post['clashes'],flush=True)
        write_tsv(OUT/'model_audit.tsv',rows)
    write_tsv(OUT/'structure_sources.tsv',source_rows)
    write_tsv(OUT/'control_contacts.tsv',positions)

def benchmark():
    lock(); rows=read_tsv(OUT/'ground_truth.tsv'); scores=[]
    unique=sorted({r['mutation'] for r in rows if r['included']=='true'})
    for audit in read_tsv(OUT/'model_audit.tsv'):
        ident=audit['model'];base=ROOT/audit['path'];folder=base.parent
        wt=energy(base,ident+'_wt');scores.append(dict(model=ident,mutation='WT',mode='repacked',**wt))
        for text in unique:
            changes=mutations(text)
            if any((' ',p,' ') not in read_structure(base)[0][c] for c,p,a,b in changes):
                continue
            key=text.replace('+','_');target=folder/(key+'.pdb');mutate(base,changes,target,ident+'_'+key)
            total=energy(target,ident+'_'+key)
            scores.append(dict(model=ident,mutation=text,mode='repacked',**total))
            fixed=folder/(key+'_fixed.pdb');combine(target,base,fixed)
            terms=energy(fixed,ident+'_'+key+'_fixed')
            scores.append(dict(model=ident,mutation=text,mode='fixed_e2',**terms))
        write_tsv(OUT/'energy_terms.tsv',scores)
        print('scored',ident,flush=True)
    summarize_benchmark(scores,rows)

def summarize_benchmark(scores,truth):
    lookup={(r['model'],r['mutation'],r['mode']):r for r in scores}; outputs=[];summaries=[]
    labels={(r['mutation'],r['normalized_direction']) for r in truth if r['included']=='true' and 'entry' not in r['assay_type'] and r['normalized_direction'] in {'enhanced','reduced'} and '+' not in r['mutation'] and r['followup_role']=='original_evidence'}
    extra={(r['mutation'],r['normalized_direction']) for r in truth if r['followup_role']=='external_source_challenge' and r['normalized_direction']=='reduced'}
    for model in sorted({r['model'] for r in scores}):
        wt=float(lookup[model,'WT','repacked']['Total'])
        for mode in ['repacked','fixed_e2']:
            for row in truth:
                key=(model,row['mutation'],mode)
                if row['included']!='true' or key not in lookup:continue
                score=lookup[key];delta=float(score['Total'])-wt
                outputs.append(dict(model=model,mutation=row['mutation'],mode=mode,source_id=row['source_id'],assay=row['assay_type'],role=row['followup_role'],expected=row['normalized_direction'],delta=delta,predicted=predicted_direction(delta,.5),vdwrep_delta=float(score.get('interD_vdwrep',0))-float(lookup[model,'WT','repacked'].get('interD_vdwrep',0))))
            for tolerance in protocol()['diagnostic_tolerance_grid']:
                correct=sum(predicted_direction(float(lookup[model,mut,mode]['Total'])-wt,tolerance)==label for mut,label in labels if (model,mut,mode) in lookup)
                n=sum((model,mut,mode) in lookup for mut,label in labels)
                external=sum(predicted_direction(float(lookup[model,mut,mode]['Total'])-wt,tolerance)==label for mut,label in extra if (model,mut,mode) in lookup)
                summaries.append(dict(model=model,mode=mode,tolerance=tolerance,original_single_correct=correct,original_single_total=n,always_reduced_baseline=sum(label=='reduced' for mut,label in labels if (model,mut,mode) in lookup),external_reduced_correct=external,external_reduced_total=sum((model,mut,mode) in lookup for mut,label in extra),wt_total=wt,wt_vdwrep=float(lookup[model,'WT','repacked'].get('interD_vdwrep',0)),interpretation='diagnostic only; no model or threshold selected by these labels'))
    write_tsv(OUT/'benchmark.tsv',outputs);write_tsv(OUT/'benchmark_summary.tsv',summaries)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('stage',choices=['prepare','models','benchmark','sequences','coverage','states','controls','report','verify']);args=parser.parse_args()
    start=time.monotonic()
    if args.stage in ['sequences','coverage','states']:
        from . import followup_diversity
        getattr(followup_diversity,args.stage)()
    elif args.stage in ['controls','report','verify']:
        from . import followup_report
        getattr(followup_report,args.stage)()
    else:
        {'prepare':prepare,'models':make_models,'benchmark':benchmark}[args.stage]()
    snapshot();print('complete',args.stage,'seconds',round(time.monotonic()-start,2),flush=True)

if __name__=='__main__':main()
