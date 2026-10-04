"""Expanded convenience cohort, leakage-aware future reserve, and state panel."""
from __future__ import annotations
from collections import Counter, defaultdict
import hashlib
import json
import random
import time
from urllib.parse import urlencode
from Bio import SeqIO
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
from .common import ROOT, now, read_tsv, write_json, write_tsv, sha256, check_resources
from .structures import interface, map_alignment, standard_residues, sequence, read_structure
from .diversity import EUTIL, e2_translation, metadata, ambiguous_reason, panel
from .followup import OUT, RAW, protocol, fetch, lock, mutate

def sequences():
    lock(); queries=[];records={};manifest=[];proteins=[]
    # Retain all archived retrievals, including their documented exclusions.
    for path in sorted((ROOT/'data/raw').glob('hcv_*.gb')):
        stratum=path.stem.removeprefix('hcv_')
        for record in SeqIO.parse(path,'genbank'):
            records.setdefault(record.id,(record,stratum,'archived_retrieval'))
    for stratum in protocol()['new_sequence_queries']:
        term=f'"{stratum}"[All Fields]' if stratum not in {'7','8'} else f'("genotype {stratum}"[All Fields] OR "{stratum}a"[All Fields])'
        query=f'(txid11103[Organism:exp] OR "Orthohepacivirus hominis"[Organism]) AND 8000:11000[SLEN] AND {term} NOT synthetic[Title] NOT chimeric[Title]'
        url=EUTIL+'esearch.fcgi?'+urlencode(dict(db='nuccore',term=query,retmode='json',retmax=protocol()['maximum_records_per_query']))
        source=fetch(url,RAW/f'esearch_{stratum}.json',stratum)
        result=json.loads(source.read_text())['esearchresult'];ids=result['idlist']
        queries.append(dict(stratum=stratum,query=query,accepted_query=result.get('querytranslation',''),warnings=result.get('errorlist',{}),available_count=result['count'],requested=len(ids),uids=','.join(ids),sampling=protocol()['sampling'],search_sha256=sha256(source)))
        if ids:
            time.sleep(.4)
            url=EUTIL+'efetch.fcgi?'+urlencode(dict(db='nuccore',id=','.join(ids),rettype='gb',retmode='text'))
            batch=fetch(url,RAW/f'hcv_{stratum}.gb',','.join(ids)); received=list(SeqIO.parse(batch,'genbank'))
            if len(received)!=len(ids):raise RuntimeError(f'Incomplete retrieval for {stratum}: {len(received)}/{len(ids)}')
            for record in received: records.setdefault(record.id,(record,stratum,'expanded_retrieval'))
        write_tsv(OUT/'queries.tsv',queries)
        print('retrieved stratum',stratum,'records',len(ids),flush=True)
        time.sleep(.4)
    for accession,(record,stratum,cohort) in sorted(records.items()):
        protein,reason=e2_translation(record)
        if protein: reason=ambiguous_reason(protein,.01)
        title=record.description.lower()
        journal=' '.join(getattr(r,'journal','') for r in record.annotations.get('references',[])).lower()
        if 'patent' in journal: reason='patent_record'
        elif any(t in title for t in ['chimeric','synthetic construct','infectious clone','cell culture adapted']):reason='laboratory_construct_or_adaptation'
        meta=metadata(record,stratum)
        manifest.append(dict(accession=accession,cohort=cohort,description=record.description,**meta,protein_length=len(protein or ''),status='excluded' if reason else 'pending_mapping',reason=reason,record_sequence_sha256=hashlib.sha256(str(record.seq).encode()).hexdigest()))
        if protein and not reason:proteins.append(dict(accession=accession,protein=protein,cohort=cohort,**meta))
    write_tsv(OUT/'sequence_manifest.tsv',manifest)
    write_json(RAW/'polyproteins.json',proteins)

def expanded_positions():
    old=sorted(int(k) for k in json.loads(read_tsv(ROOT/'data/processed/discovery_panel.tsv')[0]['residues']))
    positions=set(old);rows=[dict(h77_position=p,source='archival humanized geometric interface',e2_chain='A',structure_position='',modeled_in_archival='') for p in old]
    complex_model=MMCIFParser(QUIET=True,auth_chains=False).get_structure('bound',str(ROOT/'data/raw/7MWX.cif'))
    h77=json.loads((ROOT/'data/processed/h77_reference.json').read_text())['protein']
    for eid,rid in [('A','E'),('B','H')]:
        e=complex_model[0][eid];mapping,_=map_alignment(h77[383:746],sequence(e),384)
        pmap={r.id[1]:mapping[i] for i,r in enumerate(standard_residues(e)) if i in mapping}
        contacts=interface(e,complex_model[0][rid])
        for pos,dist in contacts['min_distances'].items():
            if dist<=5.:
                if pos not in pmap:raise RuntimeError('Unmapped bound interface residue')
                p=pmap[pos];positions.add(p);rows.append(dict(h77_position=p,source='experimental tamarin geometric interface',e2_chain=eid,structure_position=pos,modeled_in_archival=''))
    for row in read_tsv(ROOT/'results/e2_literature_interface_annotations.tsv'):
        p=int(row['h77_position']);positions.add(p);rows.append(dict(h77_position=p,source='curated experimental E2 binding annotation',e2_chain='',structure_position='',modeled_in_archival=''))
    mapping={int(r['h77_position']):int(r['structure_position']) for r in read_tsv(ROOT/'results/e2_numbering.tsv') if r['h77_position']}
    for row in rows:row['modeled_in_archival']=row['h77_position'] in mapping
    write_tsv(OUT/'interface_positions.tsv',rows)
    return old,sorted(positions),mapping

def future_split(rows,old_accessions,old_haplotypes,old_hashes,old_indices,seed):
    """Connected projected/expanded neighbors; all exposed clusters stay in development."""
    haplotypes=sorted({r['haplotype'] for r in rows}); parent={h:h for h in haplotypes}
    def find(h):
        while parent[h]!=h:parent[h]=parent[parent[h]];h=parent[h]
        return h
    projected={h:''.join(h[i] for i in old_indices) for h in haplotypes}
    for i,a in enumerate(haplotypes):
        for b in haplotypes[i+1:]:
            if sum(x!=y for x,y in zip(a,b))<=1 or sum(x!=y for x,y in zip(projected[a],projected[b]))<=1:
                aa,bb=find(a),find(b)
                if aa!=bb:parent[max(aa,bb)]=min(aa,bb)
    groups=defaultdict(list)
    for row in rows:groups[find(row['haplotype'])].append(row)
    exposed={}
    for g,items in groups.items():
        exposed[g]=any(r['accession'] in old_accessions or r['sequence_hash'] in old_hashes or any(sum(x!=y for x,y in zip(projected[r['haplotype']],h))<=1 for h in old_haplotypes) for r in items)
    eligible=sorted(g for g in groups if not exposed[g]);rng=random.Random(seed);rng.shuffle(eligible)
    n=round(.2*len(eligible))
    if len(eligible)>=2:n=max(1,n)
    reserved=set(eligible[:n])
    result=[]
    for g,items in sorted(groups.items()):
        for r in items:result.append(dict(r,component=hashlib.sha256(g.encode()).hexdigest()[:16],split='future_reserve' if g in reserved else 'development',exposed_component=exposed[g]))
    return result

def curve(rows,positions):
    chosen=panel(rows,positions,target=1.,maximum=len(rows))
    represented=0;out=[]
    for i,row in enumerate(chosen,1):
        represented+=int(row['count']);out.append(dict(states=i,accessions_represented=represented,total_accessions=len(rows),coverage=represented/max(1,len(rows)),haplotype_id=row['haplotype_id']))
    return chosen,out

def coverage():
    lock();old_positions,positions,mapping=expanded_positions()
    h77=json.loads((ROOT/'data/processed/h77_reference.json').read_text())['protein']; ref=h77[383:746]
    mapped=[];reasons={}
    for record in json.loads((RAW/'polyproteins.json').read_text()):
        fragment=record['protein'][350:800];m,_=map_alignment(ref,fragment,384);aligned={p:fragment[i] for i,p in m.items()}
        missing=[p for p in positions if p not in aligned];identity=sum(aligned.get(p)==h77[p-1] for p in range(384,747))/363
        reason='expanded_interface_missing' if missing else 'severe_alignment_failure' if identity<.35 else ''
        if not reason:reason=ambiguous_reason(''.join(aligned[p] for p in positions),.01)
        if not reason:
            indices={v:k for k,v in m.items()}
            if any(p+1 in positions and indices[p+1]-indices[p]!=1 for p in positions):reason='unmodelable_interface_insertion'
        reasons[record['accession']]=reason
        if not reason:
            seq=''.join(aligned.get(p,'-') for p in range(384,747));hashed=hashlib.sha256(seq.encode()).hexdigest()[:16]
            mapped.append(dict(accession=record['accession'],cohort=record['cohort'],genotype=record['genotype'],subtype=record['subtype'],e2_sequence=seq,sequence_hash=hashed,haplotype=''.join(aligned[p] for p in positions)))
    old_rows=json.loads((ROOT/'data/processed/discovery_sequences.json').read_text())+json.loads((ROOT/'data/processed/quarantine/heldout_sequences.json').read_text())
    indices=[positions.index(p) for p in old_positions]
    split=future_split(mapped,{r['accession'] for r in read_tsv(ROOT/'results/hcv_sequence_manifest.tsv')},{r['haplotype'] for r in old_rows},{r['sequence_hash'] for r in old_rows},indices,protocol()['reserve_seed'])
    if not split:raise RuntimeError('No eligible expanded sequences')
    dev=[r for r in split if r['split']=='development'];reserve=[r for r in split if r['split']=='future_reserve']
    write_json(OUT/'development_sequences.json',dev)
    # Sequence reserve is not consumed by modeling or selection in this follow-up.
    quarantine=OUT/'quarantine';quarantine.mkdir(exist_ok=True);write_json(quarantine/'future_reserve.json',reserve)
    write_tsv(OUT/'split_manifest.tsv',[{k:v for k,v in r.items() if k not in {'haplotype','e2_sequence'}} for r in split])
    manifest=read_tsv(OUT/'sequence_manifest.tsv')
    for r in manifest:
        if r['accession'] in reasons:r['reason']=reasons[r['accession']];r['status']='excluded' if r['reason'] else 'included'
    write_tsv(OUT/'sequence_manifest.tsv',manifest)
    selected,progress=curve(dev,positions);write_tsv(OUT/'coverage_curve.tsv',progress)
    target=next((r['states'] for r in progress if r['coverage']>=.9),len(selected))
    # Preserve panel semantics and avoid truncating at the old cap.
    chosen=panel(dev,positions,.9,target);write_tsv(OUT/'development_panel.tsv',chosen)
    prior=json.loads((ROOT/'data/processed/discovery_sequences.json').read_text());_,old_curve=curve(prior,old_positions);write_tsv(OUT/'archival_coverage_curve.tsv',old_curve)
    counts=[]
    for split_name,items in [('retrieved',manifest),('included',split),('development',dev),('future_reserve',reserve)]:
        for genotype,n in sorted(Counter(r['genotype'] for r in items).items()):
            counts.append(dict(cohort=split_name,genotype=genotype,count=n))
    write_tsv(OUT/'genotype_counts.tsv',counts)
    represented={r['haplotype'] for r in chosen};by_genotype=[]
    for genotype in sorted({r['genotype'] for r in dev}):
        items=[r for r in dev if r['genotype']==genotype];n=sum(r['haplotype'] in represented for r in items)
        by_genotype.append(dict(genotype=genotype,development_accessions=len(items),represented=n,coverage=n/len(items)))
    write_tsv(OUT/'genotype_coverage.tsv',by_genotype)
    old24=next((r['coverage'] for r in old_curve if r['states']==24),1)
    old90=next(r['states'] for r in old_curve if r['coverage']>=.9)
    summary=dict(retrieved=len(manifest),included=len(split),excluded=len(manifest)-len(split),development=len(dev),future_reserve=len(reserve),unique_e2=len({r['sequence_hash'] for r in split}),expanded_interface_positions=positions,original_interface_positions=old_positions,unmodeled_interface_positions=[p for p in positions if p not in mapping],development_haplotypes=len({r['haplotype'] for r in dev}),states_for_90_percent=target,achieved_coverage=float(chosen[0]['achieved_coverage']),archival_24_state_coverage=old24,archival_states_for_90_percent=old90,unknown_genotype=sum(r['genotype']=='unknown' for r in split),genotypes_missing=sorted(set(map(str,range(1,9)))-{r['genotype'] for r in split}),reserve_note='Novel components only; archived test is spent. A reserve from the same retrieval process is not an independent experimental validation.',coverage_note='Exact interface haplotype coverage of the bounded eligible accession cohort, not world population or therapeutic protection')
    write_json(OUT/'coverage_summary.json',summary)
    print(json.dumps(summary,sort_keys=True),flush=True)

def states():
    lock();check_resources(protocol()['maximum_new_model_bytes'])
    _,positions,mapping=expanded_positions();base=ROOT/'data/processed/humanized_repaired.pdb';s=read_structure(base);rows=[]
    folder=OUT/'states';folder.mkdir(exist_ok=True)
    start=time.monotonic();panel_rows=read_tsv(OUT/'development_panel.tsv')
    modeled=[]
    for row in panel_rows:
        residues=json.loads(row['residues']);changes=[];unmodeled=[]
        for key,new in residues.items():
            p=int(key)
            if p not in mapping:unmodeled.append(p);continue
            pos=mapping[p];wt=seq1(s[0]['E'][(' ',pos,' ')].resname)
            if wt!=new:changes.append(('E',pos,wt,new))
        target=folder/(row['haplotype_id']+'.pdb')
        if changes:mutate(base,changes,target,'state_'+row['haplotype_id'])
        elif not target.exists():__import__('shutil').copyfile(base,target)
        audit=interface(read_structure(target)[0]['E'],read_structure(target)[0]['R'])
        rows.append(dict(haplotype_id=row['haplotype_id'],count=row['count'],genotypes=row['genotypes'],mutation_count=len(changes),mutations=[f'{wt}{pos}{new}' for c,pos,wt,new in changes],unmodeled_h77_positions=unmodeled,clashes=audit['clashes'],sha256=sha256(target),path=str(target.relative_to(ROOT))))
        write_tsv(OUT/'state_manifest.tsv',rows)
        if len(rows)==3:
            elapsed=time.monotonic()-start;projected=int(sum(target.stat().st_size for target in folder.glob('*.pdb'))/3*len(panel_rows));check_resources(projected,3,'expanded_state_count')
            write_json(OUT/'state_pilot.json',dict(states=3,seconds=elapsed,projected_seconds=elapsed/3*len(panel_rows),projected_bytes=projected,planned_states=len(panel_rows),safety_factor=1.2))
        if len(rows)%10==0:print('built states',len(rows),'/',len(panel_rows),flush=True)
    write_json(OUT/'state_summary.json',dict(states=len(rows),total_seconds=time.monotonic()-start,bytes=sum((ROOT/r['path']).stat().st_size for r in rows),model='Archival experimental human CD81 backbone; fixed E2 backbone',claims='Coordinates prepared only. No candidate re-ranking or affinity claim; glycans, backbone adaptation and viral indels not modeled.'))

def balanced_panel():
    """Extend the original global panel to reach the same target within every recorded genotype."""
    dev=json.loads((OUT/'development_sequences.json').read_text());positions=json.loads((OUT/'coverage_summary.json').read_text())['expanded_interface_positions']
    base=read_tsv(OUT/'development_panel.tsv');chosen={r['haplotype'] for r in base};grouped=defaultdict(list)
    for row in dev:grouped[row['haplotype']].append(row)
    for genotype in sorted({r['genotype'] for r in dev}):
        items=[r for r in dev if r['genotype']==genotype]
        while sum(r['haplotype'] in chosen for r in items)/len(items)<.9:
            eligible=set(r['haplotype'] for r in items)-chosen
            h=min(eligible,key=lambda h:(-sum(r['genotype']==genotype for r in grouped[h]),-len(grouped[h]),h));chosen.add(h)
    coverage=sum(r['haplotype'] in chosen for r in dev)/len(dev)
    output=[]
    for h in sorted(chosen):
        items=grouped[h];output.append(dict(haplotype_id='hap_'+hashlib.sha256(h.encode()).hexdigest()[:12],haplotype=h,count=len(items),genotypes=','.join(sorted({r['genotype'] for r in items})),subtypes=Counter(r['subtype'] for r in items),residues=dict(zip(map(str,positions),h)),coverage_target=.9,achieved_coverage=coverage))
    write_tsv(OUT/'balanced_development_panel.tsv',output)
    write_tsv(OUT/'balanced_genotype_coverage.tsv',[dict(genotype=g,accessions=sum(r['genotype']==g for r in dev),represented=sum(r['genotype']==g and r['haplotype'] in chosen for r in dev),coverage=sum(r['genotype']==g and r['haplotype'] in chosen for r in dev)/sum(r['genotype']==g for r in dev)) for g in sorted({r['genotype'] for r in dev})])
    # The initial panel is retained. Build only its supplement, in the same model.
    original=OUT/'development_panel.tsv';temporary=OUT/'balanced_only_panel.tsv'
    write_tsv(temporary,[r for r in output if r['haplotype'] not in {r['haplotype'] for r in base}])
    build_supplement(read_tsv(temporary))
    write_json(OUT/'balanced_panel_summary.json',dict(states=len(output),additional_states=len(output)-len(base),coverage=coverage,minimum_genotype_coverage=min(float(r['coverage']) for r in read_tsv(OUT/'balanced_genotype_coverage.tsv')),note='Secondary coverage extension after observing genotype imbalance; no experimental or candidate scores used; genotype 8 structural exclusion persists'))

def build_supplement(panel_rows):
    _,positions,mapping=expanded_positions();base=ROOT/'data/processed/humanized_repaired.pdb';s=read_structure(base);rows=read_tsv(OUT/'state_manifest.tsv');folder=OUT/'states'
    for row in panel_rows:
        changes=[];unmodeled=[]
        for key,new in json.loads(row['residues']).items():
            p=int(key)
            if p not in mapping:unmodeled.append(p);continue
            pos=mapping[p];wt=seq1(s[0]['E'][(' ',pos,' ')].resname)
            if wt!=new:changes.append(('E',pos,wt,new))
        target=folder/(row['haplotype_id']+'.pdb')
        if changes:mutate(base,changes,target,'state_'+row['haplotype_id'])
        elif not target.exists():__import__('shutil').copyfile(base,target)
        if row['haplotype_id'] not in {r['haplotype_id'] for r in rows}:
            rows.append(dict(haplotype_id=row['haplotype_id'],count=row['count'],genotypes=row['genotypes'],mutation_count=len(changes),mutations=[f'{wt}{pos}{new}' for c,pos,wt,new in changes],unmodeled_h77_positions=unmodeled,clashes=interface(read_structure(target)[0]['E'],read_structure(target)[0]['R'])['clashes'],sha256=sha256(target),path=str(target.relative_to(ROOT))))
            write_tsv(OUT/'state_manifest.tsv',rows)
    print('balanced states',len(rows),flush=True)
