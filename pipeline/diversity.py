"""Public sequence retrieval, H77 mapping, cluster split, discovery-only panel."""
from __future__ import annotations
from collections import Counter, defaultdict
import hashlib
import json
import math
import random
import re
import time
from urllib.parse import urlencode
from Bio import SeqIO
import numpy as np
from .common import ROOT, config, failure, fetch, now, read_tsv, sha256, write_json, write_tsv
from .structures import map_alignment, read_structure, sequence, standard_residues

EUTIL='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
STRATA=['1a','1b','2a','2b','3a','4a','5a','6a']

def e2_translation(record):
    cds=[f for f in record.features if f.type=='CDS' and 'translation' in f.qualifiers]
    proteins=[f.qualifiers['translation'][0] for f in cds]
    if not proteins:
        return None,'no_annotated_CDS_translation'
    protein=max(proteins,key=len)
    if len(protein)<700:
        return None,'insufficient_polyprotein_coverage'
    if '*' in protein:
        return None,'stop_codon'
    # Do not slice at fixed absolute indices in non-H77 genotypes.
    return protein,''

def ambiguous_reason(seq,max_fraction=0.01):
    if '*' in seq:
        return 'stop_codon'
    if sum(a not in 'ACDEFGHIKLMNPQRSTVWY-' for a in seq)/max(1,len(seq))>max_fraction:
        return 'excessive_ambiguous_amino_acids'
    return ''

def entropy(values):
    count=Counter(a for a in values if a in 'ACDEFGHIKLMNPQRSTVWY')
    n=sum(count.values())
    return -sum((v/n)*math.log2(v/n) for v in count.values()) if n else 0.0

def metadata(record,stratum):
    source=next((f for f in record.features if f.type=='source'),None)
    q=source.qualifiers if source else {}
    texts=' '.join([record.description]+[v for vs in q.values() for v in vs])
    match=re.search(r'(?:genotype|subtype)\s*[:=]?\s*([1-8][a-z]?)\b',texts,re.I)
    subtype=match.group(1).lower() if match else 'unknown'
    # Query enrichment is retained explicitly and never mislabeled as a verified genotype.
    return dict(genotype=subtype[0] if subtype!='unknown' else 'unknown',subtype=subtype,
                query_stratum=stratum,genotype_evidence='record source/definition' if match else 'unclassified; query term alone insufficient',
                geography=';'.join(q.get('geo_loc_name',q.get('country',[]))),collection_date=';'.join(q.get('collection_date',[])))

def sequences():
    rows=[]; dataset=[]; seen=set(); proteins=[]
    reference=fetch(EUTIL+'efetch.fcgi?'+urlencode(dict(db='nuccore',id='AF009606.1',rettype='gb',retmode='text')),
                    ROOT/'data/raw/H77.gb','AF009606.1')
    record=SeqIO.read(reference,'genbank'); hprotein,reason=e2_translation(record)
    if reason:
        raise RuntimeError('H77 reference CDS unavailable')
    write_json(ROOT/'data/processed/h77_reference.json',dict(accession=record.id,protein=hprotein,e2_start=384,e2_end=746,
               numbering='H77 polyprotein amino acids, one-based',sha256=sha256(reference)))
    for stratum in STRATA:
        query=f'(txid11103[Organism:exp] OR "Orthohepacivirus hominis"[Organism]) AND 8000:11000[SLEN] AND "{stratum}"[All Fields] NOT patent[filter] NOT synthetic[Title] NOT chimeric[Title]'
        params=dict(db='nuccore',term=query,retmode='json',retmax=40)
        search=fetch(EUTIL+'esearch.fcgi?'+urlencode(params),ROOT/f'data/raw/esearch_{stratum}.json',stratum)
        data=json.loads(search.read_text())['esearchresult']; ids=data['idlist']
        dataset.append(dict(dataset='GenBank nuccore',stratum=stratum,exact_query=query,accepted_query=data.get('querytranslation',''),source_query_warnings=data.get('errorlist',{}),download_date=now(),records_returned=data['count'],
                            records_requested=len(ids),sampling='first 40 NCBI default-order records per query; convenience sample',numbering='H77 AF009606.1',ids=','.join(ids)))
        if not ids:
            continue
        time.sleep(0.4)
        batch=fetch(EUTIL+'efetch.fcgi?'+urlencode(dict(db='nuccore',id=','.join(ids),rettype='gb',retmode='text')),
                    ROOT/f'data/raw/hcv_{stratum}.gb',','.join(ids))
        for record in SeqIO.parse(batch,'genbank'):
            if record.id in seen:
                continue
            seen.add(record.id)
            protein,reason=e2_translation(record)
            if protein:
                reason=ambiguous_reason(protein,config()['max_ambiguous_fraction'])
            title=record.description.lower()
            if any(term in title for term in ['chimeric','synthetic construct','infectious clone','cell culture adapted']):
                reason='laboratory_construct_or_adaptation'
            meta=metadata(record,stratum)
            row=dict(accession=record.id,sequence_length=len(record.seq),protein_length=len(protein or ''),
                     **meta,inclusion_status='excluded' if reason else 'pending_mapping',exclusion_reason=reason,
                     accession_multiplicity=1,download_date=now(),description=record.description)
            rows.append(row)
            if protein and not reason:
                proteins.append(dict(accession=record.id,protein=protein,**meta))
    write_tsv(ROOT/'config/datasets.tsv',dataset)
    write_tsv(ROOT/'results/hcv_sequence_manifest.tsv',rows)
    write_json(ROOT/'data/processed/hcv_polyproteins.json',proteins)

def deterministic_split(rows,seed,heldout_fraction=0.2,near_hamming=1):
    """Connected components prevent exact and Hamming-neighbor haplotype leakage."""
    haplotypes=sorted({r['haplotype'] for r in rows})
    parent={h:h for h in haplotypes}
    def find(h):
        while parent[h]!=h:
            parent[h]=parent[parent[h]]; h=parent[h]
        return h
    for i,a in enumerate(haplotypes):
        for b in haplotypes[i+1:]:
            if len(a)==len(b) and sum(x!=y for x,y in zip(a,b))<=near_hamming:
                ra,rb=find(a),find(b)
                if ra!=rb:
                    parent[max(ra,rb)]=min(ra,rb)
    groups=defaultdict(list)
    for row in rows:
        groups[find(row['haplotype'])].append(row)
    strata=defaultdict(list)
    for group,items in groups.items():
        modal=Counter((r['genotype'],r['subtype']) for r in items).most_common(1)[0][0]
        strata[modal].append(group)
    rng=random.Random(seed); assignments={}; fallback=[]
    for stratum,clusters in sorted(strata.items()):
        clusters=sorted(clusters); rng.shuffle(clusters)
        if len(clusters)>=5:
            n=max(1,round(len(clusters)*heldout_fraction))
            for i,g in enumerate(clusters):
                assignments[g]='held_out' if i<n else 'discovery'
        else:
            fallback.extend(clusters)
    fallback=sorted(fallback); rng.shuffle(fallback)
    n=round(len(fallback)*heldout_fraction)
    if len(fallback)>=2:
        n=max(1,n)
    for i,g in enumerate(fallback):
        assignments[g]='held_out' if i<n else 'discovery'
    return [dict(**row,cluster_id=hashlib.sha256(find(row['haplotype']).encode()).hexdigest()[:16],split=assignments[find(row['haplotype'])]) for row in rows]

def panel(rows,positions,target=0.9,maximum=24):
    grouped=defaultdict(list)
    for row in rows:
        grouped[row['haplotype']].append(row)
    remaining=set(grouped); chosen=[]; total=len(rows); represented=0
    genotypes=sorted({r['genotype'] for r in rows})
    # Round-robin genotypes; within each choose most frequent remaining haplotype.
    while remaining and len(chosen)<maximum and represented/max(total,1)<target:
        progressed=False
        for genotype in genotypes:
            eligible=[h for h in remaining if any(r['genotype']==genotype for r in grouped[h])]
            if not eligible:
                continue
            h=sorted(eligible,key=lambda h:(-sum(r['genotype']==genotype for r in grouped[h]),-len(grouped[h]),h))[0]
            remaining.remove(h); chosen.append(h); represented+=len(grouped[h]); progressed=True
            if len(chosen)>=maximum or represented/max(total,1)>=target:
                break
        if not progressed:
            break
    output=[]
    for h in chosen:
        items=grouped[h]
        output.append(dict(haplotype_id='hap_'+hashlib.sha256(h.encode()).hexdigest()[:12],haplotype=h,count=len(items),
                           genotypes=','.join(sorted({r['genotype'] for r in items})),subtypes=Counter(r['subtype'] for r in items),
                           residues=dict(zip(map(str,positions),h)),coverage_target=target,achieved_coverage=represented/max(total,1)))
    return output

def diversity():
    h77=json.loads((ROOT/'data/processed/h77_reference.json').read_text())['protein']
    ref=h77[383:746]
    structure=read_structure(ROOT/'data/processed/humanized_unrepaired.pdb')
    e2=standard_residues(structure[0]['E']); mapping,score=map_alignment(ref,sequence(structure[0]['E']),384)
    numbering=[]
    for i,r in enumerate(e2):
        pos=mapping.get(i)
        numbering.append(dict(structure_position=r.id[1],structure_residue=__import__('Bio.SeqUtils',fromlist=['seq1']).seq1(r.resname),
                              h77_position=pos or '',h77_reference_residue=h77[pos-1] if pos else '',mapping_method='Biopython global affine-gap sequence alignment',
                              reference_accession='AF009606.1',mapping_status='mapped' if pos else 'unmapped'))
    write_tsv(ROOT/'results/e2_numbering.tsv',numbering)
    positions=[]; interface_rows=read_tsv(ROOT/'results/e2_interface_residues.tsv')
    pmap={r['structure_position']:r['h77_position'] for r in numbering}
    for row in interface_rows:
        row['h77_position']=pmap[int(row['e2_structure_position'])]
        if not row['h77_position']:
            raise RuntimeError('Interface residue lacks reliable H77 mapping')
        if float(row['cutoff_angstrom'])==config()['interface_cutoff']:
            positions.append(int(row['h77_position']))
    positions=sorted(set(positions)); write_tsv(ROOT/'results/e2_interface_residues.tsv',interface_rows)
    mapped=[]; reasons={}; unique={}
    for record in json.loads((ROOT/'data/processed/hcv_polyproteins.json').read_text()):
        protein=record['protein']
        # Long polyproteins: compare E2 neighborhood only, with generous flanks.
        fragment=protein[350:800]
        m,_=map_alignment(ref,fragment,384)
        aligned={p:fragment[i] for i,p in m.items()}
        missing=[p for p in positions if p not in aligned]
        identity=sum(aligned.get(p)==h77[p-1] for p in range(384,747))/363
        reason='interface_coverage_missing' if missing else 'severe_alignment_failure' if identity<0.35 else ''
        if not reason:
            reason=ambiguous_reason(''.join(aligned.get(p,'-') for p in positions),config()['max_ambiguous_fraction'])
        # Structural side-chain substitution cannot represent interface indels.
        if not reason:
            iface_indices=sorted(i for i,p in m.items() if p in positions)
            for p in positions:
                if p+1 in positions:
                    indices={v:k for k,v in m.items()}
                    if indices[p+1]-indices[p]!=1:
                        reason='unmodelable_interface_insertion'; break
        reasons[record['accession']]=reason
        if not reason:
            e2seq=''.join(aligned.get(p,'-') for p in range(384,747))
            key=hashlib.sha256(e2seq.encode()).hexdigest()[:16]
            unique[key]=e2seq
            mapped.append(dict(accession=record['accession'],genotype=record['genotype'],subtype=record['subtype'],
                               e2_sequence=e2seq,sequence_hash=key,haplotype=''.join(aligned[p] for p in positions)))
    split=deterministic_split(mapped,config()['master_seed'],config()['heldout_fraction'],config()['near_haplotype_hamming'])
    if not split:
        raise RuntimeError('No eligible natural sequences; preserve exclusions')
    # This is the quarantine boundary. Downstream selection receives discovery only.
    write_tsv(ROOT/'config/split_manifest.tsv',[{k:v for k,v in r.items() if k not in {'e2_sequence','haplotype'}} for r in split])
    discovery=[r for r in split if r['split']=='discovery']; heldout=[r for r in split if r['split']=='held_out']
    write_json(ROOT/'data/processed/discovery_sequences.json',discovery)
    quarantine=ROOT/'data/processed/quarantine'; quarantine.mkdir(parents=True,exist_ok=True)
    write_json(quarantine/'heldout_sequences.json',heldout)
    manifest=read_tsv(ROOT/'results/hcv_sequence_manifest.tsv'); multiplicity=Counter(r['sequence_hash'] for r in mapped)
    lookup={r['accession']:r for r in mapped}
    for row in manifest:
        if row['accession'] in reasons:
            row['exclusion_reason']=reasons[row['accession']]; row['inclusion_status']='excluded' if row['exclusion_reason'] else 'included'
            if row['accession'] in lookup:
                row['accession_multiplicity']=multiplicity[lookup[row['accession']]['sequence_hash']]
    write_tsv(ROOT/'results/hcv_sequence_manifest.tsv',manifest)
    # Explicit no-leakage audit; metadata split construction may inspect input sequences,
    # but no held-out score or conservation informs design or state selection.
    dh={r['haplotype'] for r in discovery}; hh={r['haplotype'] for r in heldout}
    overlap=len(dh&hh); near=sum(sum(a!=b for a,b in zip(x,y))<=config()['near_haplotype_hamming'] for x in dh for y in hh)
    write_tsv(ROOT/'results/leakage_audit.tsv',[dict(check='interface_haplotype_overlap',value=overlap,passed=overlap==0),
              dict(check='near_interface_haplotype_overlap',value=near,passed=near==0),
              dict(check='duplicate_E2_sequence_overlap',value=len({r['sequence_hash'] for r in discovery}&{r['sequence_hash'] for r in heldout}),passed=not ({r['sequence_hash'] for r in discovery}&{r['sequence_hash'] for r in heldout}))])
    if overlap or near:
        raise RuntimeError('Leakage audit failed')
    dpanel=panel(discovery,positions,config()['coverage_target'],config()['max_discovery_states'])
    write_tsv(ROOT/'data/processed/discovery_panel.tsv',dpanel)
    write_tsv(ROOT/'results/e2_haplotypes.tsv',dpanel)
    conservation=[]
    # Discovery-only conservation. Held-out conservation can be added after freeze.
    for p in range(384,747):
        values=[r['e2_sequence'][p-384] for r in discovery]; count=Counter(values); canonical=sum(v for a,v in count.items() if a in 'ACDEFGHIKLMNPQRSTVWY')
        conservation.append(dict(h77_position=p,number_observed=canonical,amino_acid_frequencies={a:v/canonical for a,v in count.items() if a in 'ACDEFGHIKLMNPQRSTVWY'} if canonical else {},
               shannon_entropy=entropy(values),consensus_residue=count.most_common(1)[0][0] if count else '',reference_residue=h77[p-1],
               gap_fraction=count['-']/max(1,len(values)),ambiguous_fraction=sum(v for a,v in count.items() if a not in 'ACDEFGHIKLMNPQRSTVWY-')/max(1,len(values)),
               genotype_specific_frequencies={g:dict(Counter(r['e2_sequence'][p-384] for r in discovery if r['genotype']==g)) for g in sorted({r['genotype'] for r in discovery})},
               split='discovery',is_interface=p in positions))
    write_tsv(ROOT/'results/e2_position_conservation.tsv',conservation)
    write_json(ROOT/'results/diversity_summary.json',dict(retrieved=len(manifest),included=len(mapped),excluded=sum(r['inclusion_status']=='excluded' for r in manifest),
               unique_e2_sequences=len(unique),discovery=len(discovery),held_out=len(heldout),discovery_haplotypes=len(dh),heldout_haplotypes=len(hh),
               interface_positions=positions,discovery_states=len(dpanel),clusters=len({r['cluster_id'] for r in split}),
               split_rule='Connected components at interface Hamming distance <=1, stratify modal genotype/subtype if >=5 clusters; small strata pooled deterministic fallback',
               caveat='Not population prevalence. First 40 per genotype-enriched query. Unknown metadata remain unknown.',
               versions=dict(biopython=__import__('Bio').__version__),alignment_command='python -m pipeline.cli diversity'))
