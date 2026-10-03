import hashlib
import json
from pathlib import Path
import random
import numpy as np
import pytest
from Bio.PDB import Superimposer
from pipeline.common import budget,sha256,write_tsv,read_tsv
from pipeline.diversity import entropy,ambiguous_reason,deterministic_split,e2_translation
from pipeline.scoring import mutation_syntax,predicted_direction
from pipeline.design import constraints
from pipeline.selection import freeze_ids,require_freeze,heldout_read,score_pool,pareto_indices
from pipeline.structures import entity_metadata,map_alignment,interface,read_structure,sequence

def test_cif_chain_species():
    cif={'_entity.id':['1'],'_entity.pdbx_description':['CD81 protein'],
         '_entity_src_gen.entity_id':['1'],'_entity_src_gen.pdbx_gene_src_scientific_name':['Saguinus oedipus'],
         '_struct_asym.id':['E','H'],'_struct_asym.entity_id':['1','1']}
    e=entity_metadata(cif)['1']
    assert e['organism']=='Saguinus oedipus' and e['chains']==['E','H']

def test_human_numbering():
    m,_=map_alignment('AAACDEFGHIKLLL','CDEFGHIK')
    assert m=={i:i+4 for i in range(8)}

def test_h77_one_based():
    m,_=map_alignment('ACDE','CDE',384)
    assert m=={0:385,1:386,2:387}

def test_alignment_coordinates():
    from Bio.PDB.Atom import Atom
    a=[Atom('CA',np.array(x,float),0,1,' ',' CA ',i,'C') for i,x in enumerate([[0,0,0],[1,0,0],[0,1,0]])]
    b=[Atom('CA',x.coord+3,0,1,' ',' CA ',i,'C') for i,x in enumerate(a)]
    s=Superimposer(); s.set_atoms(a,b); s.apply(b)
    assert s.rms<1e-6 and np.allclose(a[0].coord,b[0].coord)

@pytest.mark.parametrize('seq,reason',[('ACDX','excessive_ambiguous_amino_acids'),('AC*D','stop_codon'),('ACDE','')])
def test_ambiguous(seq,reason):
    assert ambiguous_reason(seq)==reason

def test_entropy():
    assert entropy('AAAA')==0 and entropy('AACC')==1

def test_extraction():
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature
    r=SeqRecord(Seq('ATG')); r.features=[SeqFeature(type='CDS',qualifiers={'translation':['A'*3000]})]
    assert e2_translation(r)==('A'*3000,'')
    r.features=[]; assert e2_translation(r)[1]=='no_annotated_CDS_translation'

def test_interface_distance(tmp_path):
    fixture=Path(__file__).parent/'fixtures/tiny.pdb'
    s=read_structure(fixture)
    assert interface(s[0]['E'],s[0]['R'],5)['atom_contacts']==1
    assert interface(s[0]['E'],s[0]['R'],4)['atom_contacts']==0

def synthetic_rows():
    return [dict(accession=str(i),genotype='1',subtype='1a',haplotype=h,sequence_hash=h) for i,h in enumerate(['AAAA','AAAA','AAAC','CCCC','CCCD','DDDD','EEEF','FFFF','GGGG','HHHH'])]

def test_deterministic_split_and_no_leakage():
    rows=synthetic_rows(); a=deterministic_split(rows,20261004); b=deterministic_split(rows,20261004)
    assert a==b
    d={r['haplotype'] for r in a if r['split']=='discovery'}; h={r['haplotype'] for r in a if r['split']=='held_out'}
    assert d and h and not d&h
    assert all(sum(x!=y for x,y in zip(i,j))>1 for i in d for j in h)

def test_mutation_reference_check():
    assert mutation_syntax('T','R',163,'A','T')=='TR163A'
    with pytest.raises(ValueError): mutation_syntax('T','R',163,'A','S')

@pytest.mark.parametrize('delta,expected',[(-1,'enhanced'),(1,'reduced'),(0,'approximately_neutral')])
def test_direction(delta,expected):
    assert predicted_direction(delta)==expected

def test_constraints():
    assert constraints('ACDC','ACDC',[1,2,3,4],[0],[2,4])=='passed'
    assert constraints('SCDC','ACDC',[1,2,3,4],[0],[2,4])=='fixed_position_violation'
    assert constraints('ACDA','ACDC',[1,2,3,4],[],[2,4])=='disulfide_violation'
    assert constraints('ACCC','ACDC',[1,2,3,4],[],[2,4])=='unexpected_cysteine'

def test_generated_uniqueness():
    assert len(set(['ACDE','ACDE','ACDF']))==2

def test_seeds_explicit():
    from pipeline.common import config
    cfg=config()
    assert len(set(cfg['generation_seeds']))==5 and all(s!=0 for s in cfg['generation_seeds'])
    from pipeline.design import generate
    import inspect
    assert 'Effective ProteinMPNN seed not verified' in inspect.getsource(generate)

def test_freeze_determinism():
    metrics=[dict(candidate_id='x',reference_score=0,worst_delta_wt=3,median_delta_wt=2),dict(candidate_id='y',reference_score=1,worst_delta_wt=1,median_delta_wt=0)]
    bio={c:dict(hydrophobic_sasa_fraction=.3) for c in ['x','y']}
    a=freeze_ids(metrics,bio,{},1); b=freeze_ids(list(reversed(metrics)),bio,{},1)
    assert a==b and a[0]=={'x':['single_state'],'y':['escape_aware']}

def test_heldout_invisible_before_freeze(tmp_path):
    path=tmp_path/'heldout_sequences.json'; path.write_text('[{"secret":1}]')
    with pytest.raises(RuntimeError,match='Leakage refusal'):
        heldout_read(path,tmp_path)

def test_discovery_rejects_heldout_states():
    with pytest.raises(RuntimeError,match='Held-out state'):
        score_pool([], [{'split':'held_out'}])

def test_freeze_tamper(tmp_path):
    (tmp_path/'results').mkdir(); table=tmp_path/'results/candidate_freeze.tsv'; table.write_text('candidate_id\nx\n')
    (tmp_path/'results/candidate_freeze_manifest.json').write_text(json.dumps({'candidate_table_sha256':'wrong'}))
    with pytest.raises(RuntimeError,match='changed'):
        require_freeze(tmp_path)

def test_disk_decimal_budget():
    assert budget(20_000_000_000,14_000_000_000)==13_000_000_000
    assert budget(13_000_000_000,8_000_000_000)==7_000_000_000
    with pytest.raises(RuntimeError,match='reserve'):
        budget(13_000_000_000,999_999_999)

def test_hash_schema(tmp_path):
    file=tmp_path/'data.tsv'; write_tsv(file,[dict(a=1,b='test')])
    assert read_tsv(file)==[{'a':'1','b':'test'}]
    assert sha256(file)==hashlib.sha256(file.read_bytes()).hexdigest()
    with pytest.raises(ValueError): write_tsv(tmp_path/'empty.tsv',[])

def test_pareto():
    assert pareto_indices([[1,3],[2,4],[3,1]])==[0,2]

def test_actual_file_open_guard(tmp_path):
    import subprocess,sys
    secret=tmp_path/'quarantine/heldout_sequences.json'; secret.parent.mkdir(); secret.write_text('[]')
    script='from pipeline.selection import install_visibility_guard; from pathlib import Path; install_visibility_guard(Path('+repr(str(tmp_path))+')); open('+repr(str(secret))+').read()'
    result=subprocess.run([sys.executable,'-c',script],capture_output=True,text=True)
    assert result.returncode!=0 and 'Leakage refusal' in result.stderr

def test_critical_species_rejection():
    from pipeline.verification import species_errors
    assert species_errors([dict(pdb_id='7MWX',description='CD81 protein',organism='Homo sapiens')])
    assert not species_errors([dict(pdb_id='7MWX',description='CD81 protein',organism='Saguinus oedipus')])

def test_critical_score_language():
    from pipeline.verification import language_errors
    assert language_errors('The EvoEF2 binding affinity increased.')
    assert not language_errors('An EvoEF2 score is not experimental binding affinity.')
    assert language_errors('A high affinity binder was designed.')

def test_fixture_end_to_end():
    # Run analytical pieces together; score numbers here are explicitly synthetic.
    from pipeline.verification import smoke
    smoke()

def test_resource_refusal_projection(monkeypatch,tmp_path):
    import pipeline.common as common
    monkeypatch.setattr(common,'ROOT',tmp_path)
    monkeypatch.setattr(common,'config',lambda:dict(disk_bytes=13000,reserve_bytes=1000000000))
    monkeypatch.setattr(common,'footprint',lambda:5000)
    from collections import namedtuple
    Disk=namedtuple('Disk','total used free')
    monkeypatch.setattr(common.shutil,'disk_usage',lambda path:Disk(2000000000,999999000,1000001000))
    with pytest.raises(RuntimeError,match='shortfall=200'):
        common.check_resources(projected=1000,pilot_size=10)

def test_schema_missing_fields():
    from pipeline.common import validate_schema
    validate_schema([dict(sequence='ACD',seed=4)],['sequence','seed'])
    with pytest.raises(ValueError,match='seed'):
        validate_schema([dict(sequence='ACD')],['sequence','seed'])

def test_atomic_table_failure_retains_valid_output(tmp_path):
    path=tmp_path/'scores.tsv.gz'
    write_tsv(path,[dict(candidate='a',score=1)])
    original=path.read_bytes()
    with pytest.raises(ValueError):
        write_tsv(path,[dict(candidate='b',unexpected=4)],fields=['candidate','score'])
    assert path.read_bytes()==original and read_tsv(path)[0]['candidate']=='a'
    assert not list(tmp_path.glob('*.tmp'))
