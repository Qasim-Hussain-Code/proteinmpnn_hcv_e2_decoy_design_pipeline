"""Scientific boundaries and leakage prevention in the isolated follow-up."""
from pathlib import Path
import pytest
from pipeline.common import ROOT, read_tsv
from pipeline.followup import mutations, snapshot
from pipeline.followup_diversity import future_split, curve

def test_double_mutant_keeps_two_changes():
    assert mutations('F186L+E188K')==[('R',186,'F','L'),('R',188,'E','K')]
    with pytest.raises(ValueError):mutations('F186L;E188K')

def test_exposed_neighbor_bridge_cannot_enter_future_reserve():
    rows=[dict(accession='old',haplotype='AAAA',sequence_hash='1'),dict(accession='bridge',haplotype='AAAC',sequence_hash='2'),dict(accession='novel',haplotype='AACC',sequence_hash='3'),dict(accession='other1',haplotype='CCCC',sequence_hash='4'),dict(accession='other2',haplotype='GGGG',sequence_hash='5')]
    result=future_split(rows,{'old'},{'AAAA'},{'1'},[0,1,2,3],17)
    assert all(r['split']=='development' for r in result if r['accession'] in {'old','bridge','novel'})
    assert sum(r['split']=='future_reserve' for r in result)==1

def test_projected_interface_neighbors_stay_together():
    rows=[dict(accession='a',haplotype='AACCCC',sequence_hash='1'),dict(accession='b',haplotype='AAGGGG',sequence_hash='2'),dict(accession='c',haplotype='TTAAAA',sequence_hash='3')]
    result=future_split(rows,{'a'},set(),set(),[0,1],42)
    assert all(r['split']=='development' for r in result if r['accession'] in {'a','b'})

def test_coverage_does_not_stop_at_old_cap():
    rows=[dict(haplotype=f'{i:03d}',genotype='1',subtype='1a') for i in range(40)]
    selected,progress=curve(rows,[1,2,3])
    assert len(selected)==40 and progress[-1]['coverage']==1
    assert next(r['states'] for r in progress if r['coverage']>=.9)==36

def test_corrected_soluble_controls_are_not_misattributed():
    rows=read_tsv(ROOT/'followup/ground_truth.tsv')
    soluble=[r for r in rows if r['source_id']=='higginbottom2000' and 'GST' in r['cd81_construct'] and r['included']=='true']
    assert not any(r['mutation'] in {'F186L','E188K'} for r in soluble)
    assert any(r['mutation']=='F186L+E188K' and r['normalized_direction']=='reduced' for r in soluble)
    assert all(r['normalized_direction']=='binding_retained_unquantified' for r in rows if r['source_id']=='drummer2005' and r['mutation'] in {'K124T','V146E'})

def test_original_archive_unchanged():
    assert len(snapshot())>100

def test_evo_launch_resolves_long_junction_alias(monkeypatch,tmp_path):
    from types import SimpleNamespace
    import pipeline.followup as followup
    binary=tmp_path/'evo.exe';binary.write_bytes(b'fixture executable')
    class LongAlias:
        def resolve(self):return binary
        def __str__(self):return 'long_alias/'*30+'evo.exe'
    captured=[]
    def run(args,**kwargs):
        captured.append(args)
        return SimpleNamespace(returncode=0,stdout='Total = 0.00\n',stderr='')
    monkeypatch.setattr(followup,'evo_binary',lambda:LongAlias())
    monkeypatch.setattr(followup,'OUT',tmp_path/'outputs')
    monkeypatch.setattr(followup,'WORK',tmp_path/'work')
    monkeypatch.setattr(followup.subprocess,'run',run)
    followup.invoke('ComputeBinding',ROOT/'tests/fixtures/tiny.pdb','fixture')
    assert captured[0][0]==str(binary)
