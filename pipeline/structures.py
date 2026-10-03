"""Entity-based species verification, residue numbering, and experimental hybrid."""
from __future__ import annotations
from collections import defaultdict
import copy
import json
import numpy as np
from Bio import Align
from Bio.PDB import MMCIFParser, PDBIO, Select, Superimposer
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from scipy.spatial.distance import cdist
from .common import ROOT, config, fetch, now, sha256, write_json, write_tsv

def column(cif, key):
    val = cif.get(key, [])
    return [val] if isinstance(val, str) else val

def entity_metadata(cif):
    organisms = {}
    for category, id_key, org_key in [('_entity_src_gen','entity_id','pdbx_gene_src_scientific_name'),
                                     ('_entity_src_nat','entity_id','pdbx_organism_scientific')]:
        for eid, org in zip(column(cif,category+'.'+id_key), column(cif,category+'.'+org_key)):
            organisms[eid] = org
    entities = {}
    for eid, description in zip(column(cif,'_entity.id'),column(cif,'_entity.pdbx_description')):
        entities[eid] = dict(description=description, organism=organisms.get(eid,'unknown'),chains=[])
    for chain,eid in zip(column(cif,'_struct_asym.id'),column(cif,'_struct_asym.entity_id')):
        entities[eid]['chains'].append(chain)
    return entities

def standard_residues(chain):
    return [r for r in chain if r.id[0]==' ' and seq1(r.resname,custom_map={'UNK':'X'})!='X']

def sequence(chain):
    return ''.join(seq1(r.resname) for r in standard_residues(chain))

def map_alignment(reference, query, start=1):
    aligner = Align.PairwiseAligner(mode='global',match_score=2,mismatch_score=-1,
                                   open_gap_score=-8,extend_gap_score=-0.5)
    alignment = aligner.align(reference,query)[0]
    mapping = {}
    for (a,b),(c,d) in zip(*alignment.aligned):
        for i,j in zip(range(a,b),range(c,d)):
            mapping[j] = i+start
    return mapping, alignment.score

def heavy(chain):
    return [a for r in chain for a in r.get_atoms() if a.element not in {'H','D'}]

def interface(chain_a,chain_b,cutoff=5.0):
    aa,bb = heavy(chain_a),heavy(chain_b)
    distances = cdist([a.coord for a in aa],[a.coord for a in bb])
    hits = np.argwhere(distances<=cutoff)
    pairs = {(aa[i].parent.id[1],bb[j].parent.id[1]) for i,j in hits}
    mins = {r.id[1]:float(np.min([distances[i].min() for i,a in enumerate(aa) if a.parent is r])) for r in chain_a}
    return dict(atom_contacts=int(len(hits)),residue_contacts=len(pairs),
                clashes=int(np.sum(distances<config()['clash_distance'])),min_distances=mins,
                pairs=sorted(pairs))

def read_structure(pdb):
    from Bio.PDB import PDBParser
    return PDBParser(QUIET=True).get_structure('model',str(pdb))

def save_pdb(structure,path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path))

def structures():
    rows=[]
    for identifier in ['7MWX','3X0E']:
        file=fetch(f'https://files.rcsb.org/download/{identifier}.cif',ROOT/f'data/raw/{identifier}.cif',identifier)
        cif=MMCIF2Dict(str(file))
        meta=entity_metadata(cif)
        structure=MMCIFParser(QUIET=True,auth_chains=False).get_structure(identifier,str(file))
        observed={}
        for chain in structure[0]:
            observed[chain.id]=sequence(chain)
        authors=defaultdict(set)
        for label,author in zip(column(cif,'_atom_site.label_asym_id'),column(cif,'_atom_site.auth_asym_id')):
            authors[label].add(author)
        revisions=list(zip(column(cif,'_pdbx_audit_revision_history.major_revision'),column(cif,'_pdbx_audit_revision_history.minor_revision')))
        version='.'.join(revisions[-1]) if revisions else 'unknown'
        method=';'.join(column(cif,'_exptl.method'))
        resolution=';'.join(column(cif,'_refine.ls_d_res_high'))
        for eid,entity in meta.items():
            if 'CD81' in entity['description'].upper():
                expected='Saguinus oedipus' if identifier=='7MWX' else 'Homo sapiens'
                if entity['organism']!=expected:
                    raise RuntimeError(f'Critical species failure: {identifier}: {entity}')
            rows.append(dict(pdb_id=identifier,entity_id=eid,description=entity['description'],organism=entity['organism'],
                             label_chains=','.join(entity['chains']),author_chains=','.join(sorted({a for c in entity['chains'] for a in authors[c]})),
                             version=version,method=method,resolution_angstrom=resolution,download_date=now(),sha256=sha256(file),
                             biological_assemblies=column(cif,'_pdbx_struct_assembly.id'),assembly_chains=column(cif,'_pdbx_struct_assembly_gen.asym_id_list'),
                             sequence={c:observed.get(c,'') for c in entity['chains']},
                             missing_residues={k:v for k,v in cif.items() if k.startswith('_pdbx_unobs_or_zero_occ_residues.')},
                             alternate_conformations=sorted(set(column(cif,'_atom_site.label_alt_id')))))
        write_json(ROOT/f'data/processed/{identifier}_metadata.json',dict(entities=meta,author_chains={k:sorted(v) for k,v in authors.items()},
                    assemblies={k:v for k,v in cif.items() if k.startswith('_pdbx_struct_assembly')},
                    citation={k:v for k,v in cif.items() if k.startswith('_citation.')}))
    write_tsv(ROOT/'results/structure_provenance.tsv',rows)

def humanize():
    complex_path=ROOT/'data/raw/7MWX.cif'
    human_path=ROOT/'data/raw/3X0E.cif'
    cmeta=json.loads((ROOT/'data/processed/7MWX_metadata.json').read_text())
    hmeta=json.loads((ROOT/'data/processed/3X0E_metadata.json').read_text())
    cs=MMCIFParser(QUIET=True,auth_chains=False).get_structure('complex',str(complex_path))
    hs=MMCIFParser(QUIET=True,auth_chains=False).get_structure('human',str(human_path))
    cd81ids=[c for e in cmeta['entities'].values() if e['organism']=='Saguinus oedipus' and 'CD81' in e['description'].upper() for c in e['chains']]
    e2ids=[c for e in cmeta['entities'].values() if 'eE2' in e['description'] or ('precursor' in e['description'] and 'virus' in e['organism'].lower()) for c in e['chains']]
    if not e2ids:
        e2ids=[c for e in cmeta['entities'].values() if e['description'].endswith('eE2') for c in e['chains']]
    hids=[c for e in hmeta['entities'].values() if e['organism']=='Homo sapiens' and 'CD81' in e['description'].upper() for c in e['chains']]
    # Resolve receptor/E2 pair by geometry, then lexical tie break; no chain assumptions.
    contacts=[(interface(cs[0][e],cs[0][c])['atom_contacts'],e,c) for e in e2ids for c in cd81ids]
    _,eid,cid=sorted(contacts,key=lambda x:(-x[0],x[1],x[2]))[0]
    hc=hs[0][sorted(hids)[0]]
    tc=cs[0][cid]
    hr,tr=standard_residues(hc),standard_residues(tc)
    mapping,_=map_alignment(sequence(tc),sequence(hc))
    matched=[(tr[j-1],hr[i]) for i,j in mapping.items() if 'CA' in tr[j-1] and 'CA' in hr[i]]
    sup=Superimposer()
    sup.set_atoms([a['CA'] for a,b in matched],[b['CA'] for a,b in matched])
    identity=sum(a.resname==b.resname for a,b in matched)/len(matched)
    # Human full-length numbering independently checked against UniProt P60033.
    reference=fetch('https://rest.uniprot.org/uniprotkb/P60033.fasta',ROOT/'data/raw/P60033.fasta','P60033')
    human_ref=''.join(reference.read_text().splitlines()[1:])
    hnumber,_=map_alignment(human_ref,sequence(hc))
    for i,r in enumerate(hr):
        if i not in hnumber or human_ref[hnumber[i]-1]!=seq1(r.resname):
            raise RuntimeError(f'Human CD81 reference mismatch at {r.id}: abort rather than silently humanizing sequence')
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain
    hybrid=Structure('humanized_hybrid'); model=Model(0); hybrid.add(model)
    echain=Chain('E'); rchain=Chain('R'); model.add(echain); model.add(rchain)
    for r in standard_residues(cs[0][eid]):
        echain.add(copy.deepcopy(r))
    for i,r in enumerate(hr):
        r=copy.deepcopy(r); r.id=(' ',hnumber[i],' ')
        sup.apply(list(r.get_atoms())); rchain.add(r)
    output=ROOT/'data/processed/humanized_unrepaired.pdb'; save_pdb(hybrid,output)
    before=interface(cs[0][eid],tc); after=interface(echain,rchain)
    excluded=[r.id[1] for r in tr if r not in [a for a,b in matched]]
    write_tsv(ROOT/'results/humanization_metrics.tsv',[dict(alignment_method='Global sequence alignment then least-squares C-alpha Kabsch; all mapped observed CA atoms',
               atoms_used=len(matched),matched_residues=[(a.id[1],hnumber[hr.index(b)]) for a,b in matched],sequence_identity=identity,
               ca_rmsd_angstrom=float(sup.rms),excluded_tamarin_residues=excluded,
               excluded_human_alignment_positions=[hnumber[i] for i,r in enumerate(hr) if r not in [b for a,b in matched]],
               missing_human_positions=[p for p in range(min(hnumber.values()),max(hnumber.values())+1) if p not in hnumber.values()],
               e2_label_chain=eid,tamarin_label_chain=cid,human_label_chain=hc.id,
               e2_author_chain=cmeta['author_chains'][eid],tamarin_author_chain=cmeta['author_chains'][cid],
               modeled_chain_mapping={'E':eid,'R':hc.id},contacts_before=before['atom_contacts'],contacts_after=after['atom_contacts'],
               residue_contacts_before=before['residue_contacts'],residue_contacts_after=after['residue_contacts'],
               initial_steric_clashes=after['clashes'],sidechain_repair_required=after['clashes']>0,
               excluded_reason='Unobserved or unmatched residues; no coordinates invented',structure_sha256=sha256(output))])
    write_tsv(ROOT/'results/human_residue_mapping.tsv',[dict(sequence_index=i+1,full_length_position=hnumber[i],structure_original_position=r.id[1],wild_type=seq1(r.resname)) for i,r in enumerate(hr)])
    hetero=[]
    for chain in cs[0]:
        for r in chain:
            if r.id[0]!=' ':
                hetero.append(dict(label_chain=chain.id,author_chains=cmeta['author_chains'][chain.id],residue=r.resname,
                                  position=r.id[1],atoms=len(list(r.get_atoms())),reason='Protein-only EvoEF2/ProteinMPNN model',
                                  near_selected_e2=min((float(np.linalg.norm(a.coord-b.coord)) for a in r.get_atoms() for b in heavy(cs[0][eid])),default=None)))
    write_tsv(ROOT/'results/removed_heteroatoms.tsv',hetero)
    positions=[]
    for cutoff in config()['interface_sensitivity']+[config()['interface_cutoff']]:
        audit=interface(echain,rchain,cutoff)
        for pos,distance in audit['min_distances'].items():
            if distance<=cutoff:
                residue=echain[(' ',pos,' ')]
                positions.append(dict(e2_structure_position=pos,reference_residue=seq1(residue.resname),cutoff_angstrom=cutoff,
                                      min_distance_angstrom=distance,definition='Any protein heavy atom within cutoff of humanized CD81',h77_position='pending_sequence_mapping'))
    write_tsv(ROOT/'results/e2_interface_residues.tsv',positions)
    write_json(ROOT/'data/processed/humanization.json',dict(e2_chain=eid,tamarin_chain=cid,human_chain=hc.id,
              rmsd=float(sup.rms),human_sequence=sequence(rchain),human_positions=[r.id[1] for r in rchain],e2_sequence=sequence(echain)))
